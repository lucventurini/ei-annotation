from ..abstract import AtomicOperation  # , EIWrapper
from .abstract import ChunkProteins, ProteinChunkAligner, FilterAlignments, _get_value, ProteinWrapper
from ..repeats import RepeatMasking
from ..preparation import FaidxProtein, SanitizeProteinBlastDB
from ..rnaseq.alignments.portcullis import PortcullisWrapper
import os


class MKVTreeIndex(AtomicOperation):

    def __init__(self, masked: RepeatMasking):

        super().__init__()
        self.configuration = masked.configuration
        self.masked = masked
        self.input["genome"] = self.masked_genome
        self.output["genome"] = os.path.join(self.outdir, os.path.basename(self.masked_genome))
        self.output["flag"] = os.path.join(self.outdir, "vtree_index.done")
        self.touch = True
        self.log = os.path.join(self.outdir, "mkvtree.log")

    @property
    def rulename(self):
        return "mkvtree_index"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["vmatch", "gth"]

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", "vmatch")

    @property
    def cmd(self):

        outdir = self.outdir
        load = self.load
        masked = os.path.basename(self.masked_genome)
        link_dna_dest = "{}.dna".format(os.path.basename(self.masked_genome))
        link_src = os.path.relpath(self.masked_genome, start=self.outdir)
        log = os.path.relpath(self.log, start=self.outdir)
        cmd = "{load} mkdir -p {outdir} && cd {outdir} && "
        cmd += "ln -sf {link_src} {masked} && ln -sf {link_src} {link_dna_dest} && "
        cmd += "mkvtree -v -dna -allout -pl -db {link_dna_dest} > {log} 2> {log}"
        cmd = cmd.format(**locals())
        return cmd


class GTH(ProteinChunkAligner):

    ___toolname__ = "gth"

    def __init__(self, chunks: ChunkProteins, chunk, vtree: MKVTreeIndex):
        super().__init__(chunks, chunk, vtree.masked)
        self.input.update(vtree.output)
        self.output["gff3"] = os.path.join(self.outdir, "{tool}.{dbname}_{chunk}.gff3").format(
            chunk=self.chunk, tool=self.toolname, dbname=self.dbname)
        self.log = os.path.join(self.logdir, "{tool}.{dbname}_{chunk}.log".format(
            dbname=self.dbname, tool=self.toolname, chunk=self.chunk))

    @property
    def loader(self):
        return ["gth", "ei-annotation"]

    @property
    def threads(self):
        # GTH relies on parallelising the input data, rather than on multithreading
        return 1

    @property
    def toolname(self):
        return self.___toolname__

    @property
    def species(self):
        if "species" in self.configuration["programs"][self.toolname]:
            return " -species {} ".format(self.configuration["programs"][self.toolname]["species"])
        else:
            return " "

    @property
    def coverage(self):
        if self._coverage_value:
            return " -gcmincoverage {} ".format(self._coverage_value)
        else:
            return " "

    @property
    def extra(self):
        return self.configuration["programs"].get(self.toolname, {}).get("extra", " ")

    @property
    def gcintron(self):
        if _get_value(self.configuration, self.dbname, "max_intron_middle"):
            return "-gcmaxgapwidth {}".format(_get_value(self.configuration, self.dbname, "max_intron_middle"))
        else:
            return " "

    @property
    def source(self):
        return "{dbname}.{cov}cov{hamming}id".format(dbname=self.dbname,
                                                     cov=int(round(self.coverage)),
                                                     hamming=int(round(self.hamming)))

    @property
    def hamming(self):
        val = _get_value(self.configuration, self.dbname, "hamming")
        if val and isinstance(val, int):
            return " -prhdist {val}".format(**locals())
        else:
            return ""

    @property
    def complete(self):
        val = _get_value(self.configuration, self.dbname, "prhdist")
        if val:
            return " -startcodon -finalstopcodon"
        else:
            return " "

    @property
    def seedlength(self):
        val = _get_value(self.configuration, self.dbname, "seedlength")
        if val and isinstance(val, int):
            return " -prseedlength {val} ".format(**locals())
        else:
            return ""

    @property
    def cmd(self):

        # gth -intermediate  -introncutout -first 10 -species rice -gff3out -gcmincoverage 50 -paralogs
        # -o eitest/proteins/alignments/chunk_003.gff3 -force
        # -genomic eitest/inputs/reference/genome.fa -protein eitest/proteins/chunks/chunk_003.stops.fasta

        load = self.load
        gcintron = self.gcintron  # Max intron length
        cmd = "{load} gth -intermediate  -introncutout {gcintron} "
        coverage = self.coverage
        hamming = self.hamming
        complete = self.complete
        species = self.species
        extra = self.extra
        seedlength = self.seedlength
        cmd += " {species} -gff3out {coverage} {hamming} {seedlength} {complete} {extra} -paralogs "
        input, output = self.input, self.output
        logdir, log = os.path.dirname(self.log), self.log
        cmd += " -genomic {input[genome]} -protein {input[fasta]} 2> {log} "
        cmd += """ | gth_correct.py > {output[gff3]} 2> {log}"""
        cmd = cmd.format(**locals())
        return cmd


class CollapseGTH(AtomicOperation):

    def __init__(self, runs: [GTH]):
        super().__init__()
        self.configuration = runs[0].configuration
        assert len(set([_.dbname for _ in runs])) == 1
        self.dbname = runs[0].dbname
        self.input["chunks"] = [run.output["gff3"] for run in runs]
        self.output["gff3"] = os.path.join(self.outdir, "{tool}.{dbname}.gff3").format(tool="gth", dbname=self.dbname)
        self.log = os.path.join(self.logdir, "collapse_{dbname}.log".format(dbname=self.dbname))

    @property
    def loader(self):
        return ["genometools"]

    @property
    def rulename(self):
        return "collapse_gth_{}".format(self.dbname)

    @property
    def cmd(self):
        inputs = " ".join(self.input["chunks"])
        output = self.output
        log = self.log
        cmd = "gt gff3 -sort -tidy -addids no -retainids -o {output[gff3]}"
        cmd += " -force {inputs} 2> {log} > {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "output")

    @property
    def logdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "logs")

    @property
    def threads(self):
        return 1

    @property
    def is_small(self):
        return True


class GTHProteinWrapper(ProteinWrapper):

    def __init__(self,
                 masker: RepeatMasking, portcullis: PortcullisWrapper):

        super().__init__(masker, portcullis)

    def execute_protein(self, db):
        vtree = MKVTreeIndex(self.masker)
        sanitised = SanitizeProteinBlastDB(self.configuration, db)
        faidx = FaidxProtein(sanitised)
        self.add_edge(sanitised, faidx)
        chunk_proteins = ChunkProteins(sanitised)
        self.add_edge(faidx, chunk_proteins)
        self.add_edge(self.masker, vtree)
        chunks = [GTH(chunk_proteins, chunk, vtree) for chunk in range(1, chunk_proteins.chunks + 1)]
        self.add_edges_from([(chunk_proteins, chunk) for chunk in chunks])
        self.add_edges_from([(vtree, chunk) for chunk in chunks])
        collapsed = CollapseGTH(chunks)
        self.add_edges_from([(chunk, collapsed) for chunk in chunks])
        filterer = FilterAlignments(collapse=collapsed,
                                    portcullis=self.portcullis,
                                    masker=self.masker)
        self.add_edge(self.masker, filterer)
        self.add_edge(self.portcullis, filterer)
        self.add_edge(collapsed, filterer)
        self._add_to_proteins(filterer)
