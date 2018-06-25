from ..abstract import EIWrapper, AtomicOperation
from ..rnaseq.mikado import MikadoStats
from ..rnaseq.mikado.workflow import Mikado
from ..repeats.workflow import RepeatMasking
import os


class FlnWrapper(EIWrapper):

    def __init__(self, mikado: Mikado, masker: RepeatMasking):

        super().__init__()

        if mikado.stats:
            # Start extracting the data ...
            sequence_extractor = MikadoSequenceExtractor(mikado.stats)
            self.add_edge(mikado, sequence_extractor)
            indexer = MikadoSequenceIndexer(sequence_extractor)
            self.add_edge(sequence_extractor, indexer)
            splitter = SplitMikadoFasta(indexer)
            self.add_edge(indexer, splitter)
            fln_chunks = [FLN(splitter, chunk) for chunk in range(1, splitter.chunks + 1)]
            self.add_edges_from([(splitter, chunk) for chunk in fln_chunks])
            concat = ConcatenateFLN(fln_chunks)
            self.add_edges_from([(chunk, concat) for chunk in fln_chunks])
            fln_filter = FilterFLN(concat)
            self.add_edge(concat, fln_filter)
            self.categories = dict((cat, ExtractFLNClass(fln_filter, mikado.stats, cat))
                                   for cat in ("Gold", "Silver", "Bronze"))
            self.add_edges_from([(mikado.stats, self.categories[cat]) for cat in self.categories])
            self.add_edges_from([(fln_filter, self.categories[cat]) for cat in self.categories])
            self.add_edges_from([(masker, self.categories[cat]) for cat in self.categories])
            # TODO: devise a final flag
            self.add_final_flag()


class MikadoSequenceExtractor(AtomicOperation):

    def __init__(self, stats: MikadoStats):

        super().__init__()
        self.input = {"loci": stats.input["link"]}
        self.input.update(stats.output)
        self.input["genome"] = self.genome
        self.configuration = stats.configuration
        self.outdir = stats.outdir
        self.output = {"transcripts": os.path.join(self.outdir, "mikado.loci.transcripts.fasta")}

    @property
    def loader(self):
        return ["gffread"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        genome = self.genome
        outdir = self.outdir

        cmd = "{load} gffread -g {genome} -w {output[transcripts]} -C {input[loci]}".format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def rulename(self):
        return "extract_mikado_fasta"


class MikadoSequenceIndexer(AtomicOperation):

    def __init__(self, extractor: MikadoSequenceExtractor):
        super().__init__()
        self.input = extractor.output
        self.configuration = extractor.configuration
        self.outdir = extractor.outdir
        self.output = {"fai": self.input["transcripts"] + ".fai"}
        self.log = os.path.join(self.outdir, "mikado.faidx.log")

    @property
    def rulename(self):
        return "faidx_mikado_transcripts"

    @property
    def loader(self):
        return ["samtools"]

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} samtools faidx {input[transcripts]} > {log} 2>&1".format(**locals())

        return cmd


class SplitMikadoFasta(AtomicOperation):

    def __init__(self, indexer: MikadoSequenceIndexer):
        super().__init__()
        self.input = indexer.output
        self.input.update(indexer.input)
        self.configuration = indexer.configuration
        self.log = os.path.join(os.path.dirname(self.outdir), "split.log")
        self.output["split_flag"] = os.path.join(self.outdir, "split.done")
        self.output["chunks"] = ["{outprefix}_{chunk}.fasta".format(
            outprefix=self.outprefix, chunk=str(_).zfill(3)) for _ in range(1, self.chunks + 1)]

    @property
    def chunks(self):
        raise NotImplementedError()

    @property
    def rulename(self):
        return "split_mikado_loci_fasta"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):

        load = self.load
        log = self.log
        input, output = self.input, self.output
        chunks = self.chunks
        outprefix = self.outprefix

        cmd = "{load} split_fasta.py -m {chunks} {input[fa]} {outprefix} && touch {output[split_flag]}".format(
            **locals()
        )

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "abinitio", "1-FLN", "Chunks")

    @property
    def outprefix(self):
        return os.path.join(self.outdir, "chunk")


class FLN(AtomicOperation):

    def __init__(self, mikado_split: SplitMikadoFasta, chunk):

        super().__init__()
        self.__chunk = None
        self.__split = mikado_split
        self.configuration = mikado_split.configuration
        self.chunk = chunk
        self.input = {"chunk": self.fasta_chunk, "split_flag": self.__split.output["split_flag"]}
        self.output = {"ptseq": os.path.join(self.outdir, "pt_seqs"),
                       "flag": os.path.join(self.outdir, "fln.done")}
        self.touch = True
        self.log = os.path.join(self.outdir, "fln.log")

    @property
    def chunk(self):
        return self.__chunk

    @chunk.setter
    def chunk(self, chunk):

        if not isinstance(chunk, int):
            raise ValueError(chunk)
        if ("{outprefix}_{chunk}.fasta".format(outprefix=self.__split.outprefix, chunk=str(chunk).zfill(3))
                not in self.input["chunks"]):
            raise ValueError("{} not found among inputs!")
        self.__chunk = chunk

    @property
    def fasta_chunk(self):
        return os.path.join(self.__split.outdir, "{outprefix}_{chunk}.fasta".format(
            outprefix=self.__split.outprefix, chunk=str(self.chunk).zfill(3)))

    @property
    def loader(self):
        return ["full_lengther_next"]

    @property
    def rulename(self):
        return "fln_chunk_{}".format(self.chunk)

    @property
    def cmd(self):

        taxon = self.taxon
        dbs = self.fln_dbs
        outdir = self.outdir
        log = os.path.basename(self.log)
        fasta_link_src = os.path.relpath(os.path.abspath(self.input["fasta_chunk"]),
                                         start=os.path.abspath(self.outdir))

        cmd = "{load} mkdir -p {outdir} && cd {outdir} && ln -s {fasta_link_src} chunk.fasta && "
        cmd += "full_lengther_next -w {threads} -g {taxon} -a {dbs} -f chunk.fasta 2>&1 > {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "abinitio", "1-FLN", "FLN", "Chunk_{}".format(
            str(self.chunk).zfill(3)))

    @property
    def taxon(self):
        raise NotImplementedError()

    @property
    def fln_dbs(self):
        raise NotImplementedError()


class ConcatenateFLN(AtomicOperation):

    def __init__(self, runs: [FLN]):

        assert runs
        super().__init__()
        self.configuration = runs[0].configuration
        self.input["flns"] = [run.output["ptseq"] for run in runs]
        self.output["fln_table"] = os.path.join(self.outdir, "fln_table.txt")

    @property
    def threads(self):
        return 1

    @property
    def rulename(self):
        return "concatenate_flns"

    @property
    def loader(self):
        return []

    @property
    def cmd(self):

        input, output = self.input, self.output
        input = " ".join(self.input["flns"])
        outdir = self.outdir
        cmd = """mkdir -p {outdir} && cat {input} | """
        cmd += """awk 'NR==1 || $1!="Query_id"' > {output[fln_table]}""".format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "abinitio", "1-FLN", "output")


class FilterFLN(AtomicOperation):

    def __init__(self, concat: ConcatenateFLN):
        super().__init__()

        self.configuration = concat.configuration
        self.input = concat.output
        self.outdir = concat.outdir
        self.output = {"table": self.outprefix + ".table.txt",
                       "list": self.outprefix + ".list.txt"}
        self.log = os.path.join(self.outdir, "filter_fln.log")

    @property
    def outprefix(self):

        return os.path.join(self.outdir, "fln_filtered")

    @property
    def threads(self):
        return 1

    @property
    def rulename(self):
        return "filter_fln"

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        input, log = self.input, self.log

        cmd = "{load} "
        cmd += "mkdir -p {outdir} && filter_fln.py {input[table]} {input[loci]} {outprefix} > {log} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd


class ExtractFLNClass(AtomicOperation):

    def __init__(self, fln_filter: FilterFLN, stats: MikadoStats, category):

        super().__init__()
        self.configuration = fln_filter.configuration
        self.input = fln_filter.output
        self.input["loci"] = stats.input["loci"]
        self.outdir = fln_filter.outdir
        self.__category = None
        self.category = category
        self.output = {"list": os.path.join(self.outdir, "{}_models.txt".format(self.category)),
                       "gff3": os.path.join(self.outdir, "{}_models.gff3".format(self.category)),
                       }

    @property
    def category(self):
        return self.__category

    @category.setter
    def category(self, category):
        if category not in ("Gold", "Silver", "Bronze"):
            raise ValueError("Invalid category! Expected one of {}, received: {}".format(
                ", ".join(("Gold", "Silver", "Bronze")), category
            ))
        self.__category = category

    @property
    def rulename(self):
        return "extract_fln_{}".format(self.category.lower())

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        category = self.category

        cmd = "{load} "
        cmd += """awk '{{if ($4=="{category}") {{OFS="\\t"; print $1,$2}} }} {input[list]} > {output[list]} && """
        cmd += """ mikado util grep {output[list]} {input[loci]} {output[gff3]}"""
        cmd = cmd.format(**locals())

        return cmd

    @property
    def threads(self):
        return 1

