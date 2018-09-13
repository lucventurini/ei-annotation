#!/usr/bin/env python3

from ..abstract import EIWrapper, AtomicOperation
from ..rnaseq.mikado.pick import MikadoStats
from ..rnaseq.mikado import Mikado
from ..rnaseq.mikado.abstract import MikadoOp
from ..preparation import SanitizeProteinBlastDB, BlastxIndex, DiamondIndex
# from ..repeats.__init__ import RepeatMasking
import os
import abc
import functools


class FlnWrapper(EIWrapper):

    __final_rulename__ = "fln_done"

    def __init__(self, mikado: Mikado):  #, masker: RepeatMasking):

        super().__init__()

        if mikado.stats:
            # Start extracting the data ...
            self.configuration = mikado.configuration
            if self.dbs:
                self.sanitizer = SanitizeProteinBlastDB(self.configuration, dbs=self.dbs)
                self.blast_index = BlastxIndex(self.sanitizer)
                self.add_edge(self.sanitizer, self.blast_index)
            else:
                self.blast_index = None
            sequence_extractor = MikadoSequenceExtractor(mikado.stats)
            if mikado.stats.is_long is True:
                diamond_index = DiamondIndex(sequence_extractor, key="proteins", rule_suffix="fln_long")
            else:
                diamond_index = DiamondIndex(sequence_extractor, key="proteins", rule_suffix="fln")
            self.self_blast = SelfDiamondP(diamond_index, is_long=mikado.stats.is_long)

            self.add_edge(mikado, sequence_extractor)
            self.add_edge(sequence_extractor, diamond_index)
            self.add_edge(diamond_index, self.self_blast)
            indexer = MikadoSequenceIndexer(sequence_extractor)
            self.add_edge(sequence_extractor, indexer)
            splitter = SplitMikadoFasta(indexer)
            self.add_edge(indexer, splitter)
            fln_chunks = [FLN(splitter, chunk, self.blast_index) for chunk in range(1, splitter.chunks + 1)]

            self.add_edges_from([(splitter, chunk) for chunk in fln_chunks])
            if self.blast_index:
                self.add_edges_from([(self.blast_index, chunk) for chunk in fln_chunks])

            concat = ConcatenateFLN(fln_chunks)
            self.add_edges_from([(chunk, concat) for chunk in fln_chunks])
            convert_mikado = ConvertMikadoToBed12(mikado.stats)
            self.add_edge(mikado, convert_mikado)
            self.fln_filter = FilterFLN(concat, convert_mikado, self.self_blast)
            self.add_edge(concat, self.fln_filter)
            self.add_edge(convert_mikado, self.fln_filter)
            self.add_edge(self.self_blast, self.fln_filter)
            self.training = FlnCategoryStats(self.fln_filter, "Training", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.training)
            self.gold = FlnCategoryStats(self.fln_filter, "Gold", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.gold)
            self.silver = FlnCategoryStats(self.fln_filter, "Silver", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.silver)
            self.bronze = FlnCategoryStats(self.fln_filter, "Bronze", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.bronze)

            if mikado.stats.is_long:
                self.__final_rulename__ += "_long"

            self.add_final_flag()
        else:
            pass

    @property
    def flag_name(self):
        return os.path.join(self.fln_filter.outdir, "fln.done")

    @property
    def dbs(self):
        return self.configuration["mikado_homology"]["prot_dbs"]


class FLNOp(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, ancestor, is_long=None):

        super().__init__()
        if is_long is None:
            self.__is_long = ancestor.is_long
        else:
            self.__is_long = is_long
        self.configuration = ancestor.configuration
        self.input.update(ancestor.output)

    @property
    @abc.abstractmethod
    def _rulename(self):
        pass

    @property
    def is_long(self):
        return self.__is_long

    @property
    def rulename(self):
        return "{rulename}{long}".format(rulename=self._rulename,
                                         long="_long" if self.is_long else "")

    @property
    def flndir(self):
        return os.path.join("abinitio", "1-FLN{long}".format(long="-long-reads" if self.is_long else ""))

    @property
    def toolname(self):
        return "full_lengther_next"


class SelfDiamondP(FLNOp):

    def __init__(self, index: DiamondIndex, is_long):

        super().__init__(index, is_long=is_long)
        self.input = index.output
        self.input["query"] = index.input["db"]
        self.outdir = os.path.dirname(self.input["query"])
        self.output["blast_txt"] = os.path.join(
            self.outdir,
            os.path.splitext(os.path.basename(self.input["query"]))[0] + ".dmnd.txt")
        self.__loader = index.loader
        self.log = os.path.join(self.outdir, "diamond_proteins.log")

    @property
    def fmt(self):
        line = "6 qseqid sseqid pident qstart qend sstart send "
        line += "qlen slen length nident mismatch positive gapopen gaps evalue bitscore"
        return line

    @property
    def loader(self):
        return self.__loader

    @property
    def cmd(self):
        load, input, output = self.load, self.input, self.output
        threads = self.threads
        log = self.log
        fmt = self.fmt
        cmd = "{load} "
        cmd += "diamond blastp --threads {threads} --outfmt {fmt} --compress 0 "
        cmd += " --out {output[blast_txt]} --db {input[db]} --query {input[query]} --sensitive "
        evalue = 1
        cmd += " --evalue {evalue} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def _rulename(self):
        return "self_diamond_blastp"


class MikadoSequenceExtractor(FLNOp):

    def __init__(self, stats: MikadoStats):

        super().__init__(stats)
        self.input["loci"] = stats.input["link"]
        self.input["genome"] = self.genome
        self.outdir = stats.outdir
        self.output = {"transcripts": os.path.join(self.outdir, "mikado.loci.transcripts.fasta"),
                       "proteins": os.path.join(self.outdir, "mikado.loci.proteins.fasta"),
                       "cds": os.path.join(self.outdir, "mikado.loci.cds.fasta")}

    @property
    def loader(self):
        return ["gffread"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        genome = self.genome
        outdir = self.outdir
        cmd = "{load} gffread -g {genome} -w {output[transcripts]}  "
        cmd += " -x {output[cds]} -y {output[proteins]} "
        cmd += " -C {input[loci]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "extract_mikado_fasta"

    @property
    def is_small(self):
        return True


class MikadoSequenceIndexer(FLNOp):

    def __init__(self, extractor: MikadoSequenceExtractor):
        super().__init__(extractor)
        self.outdir = extractor.outdir
        self.output = {"fai": self.input["transcripts"] + ".fai"}
        self.log = os.path.join(self.outdir, "mikado.faidx.log")

    @property
    def _rulename(self):
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

    @property
    def is_small(self):
        return True


class SplitMikadoFasta(FLNOp):

    def __init__(self, indexer: MikadoSequenceIndexer):
        super().__init__(indexer)
        self.input.update(indexer.input)
        self.log = os.path.join(os.path.dirname(self.outdir), "split.log")
        self.output["split_flag"] = os.path.join(self.outdir, "split.done")
        self.output["chunks"] = ["{outprefix}_{chunk}.fasta".format(
            outprefix=self.outprefix, chunk=str(_).zfill(3)) for _ in range(1, self.chunks + 1)]

    @property
    def chunks(self):
        return self.configuration.get("full_lengther_next", dict()).get("chunks", 1)

    @property
    def _rulename(self):
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

        cmd = "{load} split_fasta.py -m {chunks} {input[transcripts]} {outprefix} && touch {output[split_flag]}".format(
            **locals()
        )

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        # TODO: Change
        return os.path.join(self.configuration["outdir"], self.flndir, "Chunks")

    @property
    def outprefix(self):
        return os.path.join(self.outdir, "chunk")

    @property
    def is_small(self):
        return True


@functools.lru_cache(typed=True, maxsize=4)
def get_fln_version(loader):

    return True


class FLN(FLNOp):

    def __init__(self,
                 mikado_split: SplitMikadoFasta,
                 chunk,
                 blastx_index: BlastxIndex=None):

        super().__init__(mikado_split)
        self.__chunk = None
        self.__split = mikado_split
        self.chunk = chunk
        # self.input["fasta_chunk"] = self.fasta_chunk
        self.output = {"dbannotated": os.path.join(self.outdir, "fln_results", "pt_seqs"),
                       "flag": os.path.join(self.outdir, "fln.done")}
        self.touch = False
        self.log = os.path.join(self.outdir, "fln.log")
        if blastx_index is not None:
            self.input["fln_user"] = blastx_index.output["db"]

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
    def user_db(self):
        if self.input.get("fln_user", None):
            prefix = os.path.splitext(self.input.get("fln_user", None))[0]
            loc = os.path.relpath(prefix, start=self.outdir)
            return " -u {loc} ".format(**locals())
        else:
            return " "

    @property
    def fasta_chunk(self):
        return os.path.join(self.__split.outdir, "{outprefix}_{chunk}.fasta".format(
            outprefix=os.path.basename(self.__split.outprefix), chunk=str(self.chunk).zfill(3)))

    @property
    def loader(self):
        return ["full_lengther_next"]

    @property
    def _rulename(self):
        return "fln_chunk_{}".format(self.chunk)

    @property
    def cmd(self):

        taxon = self.taxon
        dbs = self.fln_dbs
        outdir = self.outdir
        log = os.path.basename(self.log)
        fasta_link_src = os.path.relpath(os.path.abspath(self.fasta_chunk),
                                         start=os.path.abspath(self.outdir))
        load = self.load
        threads = self.threads
        output = self.output
        db_verify = os.path.relpath(self.output["dbannotated"],
                                    start=self.outdir)
        flag = os.path.basename(self.output["flag"])

        dbs = "-a {dbs}".format(**locals()) if get_fln_version(self.load) else ""
        user_db = self.user_db

        cmd = "{load} mkdir -p {outdir} && cd {outdir} &&"
        cmd += " if [[ ! -e chunk.fasta ]]; then ln -s {fasta_link_src} chunk.fasta; fi && "
        cmd += "full_lengther_next {user_db} -w {threads} --taxon_group={taxon} {dbs} -f chunk.fasta 2>&1 > {log}"
        cmd += " && [[ -s {db_verify} ]] && touch {flag}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], self.flndir,
                            "FLN", "Chunk_{}".format(
            str(self.chunk).zfill(3)))

    @property
    def taxon(self):
        assert self.configuration["programs"]["full_lengther_next"]["taxon"] in ["fungi",
                                                                 "human", "invertebrates",
                                                                 "mammals", "plants", "rodents",
                                                                 "vertebrates"]

        return self.configuration["programs"]["full_lengther_next"]["taxon"]

    @property
    def fln_dbs(self):
        if "dbs" not in self.configuration["programs"]["full_lengther_next"]:
            return "stnp"
        else:
            assert all(item in ("s", "t", "n", "p", "c") for item
                       in self.configuration["programs"]["full_lengther_next"]["dbs"])
            return "".join(self.configuration["programs"]["full_lengther_next"]["dbs"])


class ConcatenateFLN(FLNOp):

    def __init__(self, runs: [FLN]):

        assert runs
        super().__init__(runs[0])
        self.configuration = runs[0].configuration
        self.input = {"dbs": [run.output["dbannotated"] for run in runs]}
        self.output["table"] = os.path.join(self.outdir, "fln_table.txt")

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "concatenate_flns"

    @property
    def loader(self):
        return []

    @property
    def cmd(self):

        input, output = self.input, self.output
        input = " ".join(self.input["dbs"])
        outdir = self.outdir
        cmd = """mkdir -p {outdir} && cat {input} | """
        cmd += """awk 'NR==1 || $1!="Query_id"' > {output[table]}"""
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], self.flndir, "output")

    @property
    def is_small(self):
        return True


class ConvertMikadoToBed12(MikadoOp):

    def __init__(self, mikado: MikadoStats):

        super().__init__(is_long=mikado.is_long)
        self.configuration = mikado.configuration
        self.input = mikado.output
        self.input.update(mikado.input)
        self.output = {"bed12": os.path.splitext(mikado.input["link"])[0] + ".bed12"}
        self.log = os.path.join(self.mikado_dir, "logs", "convert_to_bed12.log")

    @property
    def _rulename(self):
        return "convert_mikado_to_bed12"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        log = self.log
        cmd = "{load} mikado util convert -of bed12 {input[loci]} {output[bed12]} 2> {log} > {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def is_small(self):
        return True


class FilterFLN(FLNOp):

    def __init__(self, concat: ConcatenateFLN,
                 to_bed12: ConvertMikadoToBed12,
                 self_blast: SelfDiamondP):
        super().__init__(concat)
        self.input.update(to_bed12.output)
        self.input.update(self_blast.output)
        # self.input.update(to_bed12.input)
        self.outdir = concat.outdir
        self.output = {"table": self.outprefix + ".table.txt",
                       # "list": self.outprefix + ".list.txt",
                       }
        for category in ['Training', 'Gold', 'Silver', 'Bronze']:
            self.output[category] = self.outprefix + ".{category}.gff3".format(**locals())

        self.log = os.path.join(self.outdir, "filter_fln.log")

    @property
    def outprefix(self):
        # TODO: change
        return os.path.join(self.outdir, "output")

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "filter_fln"

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        input, log = self.input, self.log
        outprefix = self.outprefix
        mikado_loci = os.path.splitext(self.input["bed12"])[0]

        cmd = "{load} "
        cmd += "mkdir -p {outdir} && filter_fln.py {input[table]} {mikado_loci} {input[blast_txt]}"
        cmd += " {outprefix} > {log} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def is_small(self):
        return False


class FlnCategoryStats(AtomicOperation):

    def __init__(self, fln: FilterFLN, category, is_long=False):

        super().__init__()
        self.is_long = is_long
        self.input = fln.output
        # self.input.update(index.input)
        self.category = category
        self.outdir = fln.outdir
        self.configuration = fln.configuration
        self.output = {"stats": os.path.splitext(self.input[self.category])[0] + ".stats"}
        self.message = "Calculating statistics for FLN filtering, {self.category} category".format(**locals())
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return "stats_fln_{category}{long}".format(category=self.category, long="" if not self.is_long else "_long")

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["link"])

    @property
    def cmd(self):

        load = self.load
        input, output, log = self.input, self.output, self.log
        inp_file = self.input[self.category]
        cmd = "{load}  if [ -s {inp_file} ]; then "
        cmd += " mikado util stats {inp_file} {output[stats]} > {log} 2>&1; else "
        cmd += " touch {output[stats]}; fi"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def is_small(self):
        return True
