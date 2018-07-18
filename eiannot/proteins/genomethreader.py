from .chunking import ChunkProteins
from ..abstract import AtomicOperation, EIWrapper
from .abstract import ProteinChunkAligner
from ..repeats import RepeatMasking
from ..preparation import FaidxProtein, SanitizeProteinBlastDB
from ..rnaseq.alignments.portcullis import PortcullisWrapper
import os


class GTH(ProteinChunkAligner):

    ___toolname__ = "gth"

    def __init__(self, chunks: ChunkProteins, chunk, masked: RepeatMasking):
        super().__init__(chunks, chunk, masked)
        self.output["gff3"] = os.path.join(self.outdir, "{chunk}.{tool}.gff3").format(chunk=self.chunk,
                                                                                     tool=self.toolname)
        self.log = os.path.join(self.logdir, "{tool}.{chunk}.log".format(tool=self.toolname,
                                                                         chunk=self.chunk))

    @property
    def loader(self):
        return ["gth"]

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
        if "coverage" in self.configuration["homology"]:
            return " -gcmincoverage {} ".format(self.configuration["homology"]["coverage"])
        else:
            return " "

    @property
    def extra(self):
        return self.configuration["programs"].get(self.toolname, {}).get("extra", " ")

    @property
    def gcintron(self):
        return "-gcmaxgapwidth {}".format(self.max_intron)

    @property
    def cmd(self):

        # gth -intermediate  -introncutout -first 10 -species rice -gff3out -gcmincoverage 50 -paralogs
        # -o eitest/proteins/alignments/chunk_003.gff3 -force
        # -genomic eitest/inputs/reference/genome.fa -protein eitest/proteins/chunks/chunk_003.stops.fasta

        load = self.load
        gcintron = self.gcintron
        coverage = self.coverage
        cmd = "{load} gth -intermediate  -introncutout -first 10 {gcintron} "
        species = self.species
        extra = self.extra
        cmd += " {species} -gff3out {coverage} {extra} -paralogs "
        input, output = self.input, self.output
        cmd += " -genomic {input[genome]} -protein {input[fasta]} "
        cmd += """ | awk -F "\\t" '{{OFS="\\t"; if ($3=="gene") {{$3="match"}} else {{if ($3=="exon") {{$3="match_part"}} }}; print $0}}' """
        cmd += " > {output[gff3]}"

        cmd = cmd.format(**locals())
        return cmd


class CollapseGTH(AtomicOperation):

    def __init__(self, runs: [GTH]):
        super().__init__()
        self.configuration = runs[0].configuration
        self.input["chunks"] = [run.output["gff3"] for run in runs]
        self.output["gff3"] = os.path.join(self.outdir, "{tool}.gff3").format(tool="gth")

    @property
    def loader(self):
        return ["genometools"]

    @property
    def rulename(self):
        return "collapse_gth"

    @property
    def cmd(self):
        inputs = " ".join(self.input["chunks"])
        output = self.output
        cmd = "gt gff3 -sort -tidy -addids no -retainids -o {output[gff3]}"
        cmd += " -force {inputs}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "output")


class FilterGTH(AtomicOperation):

    __rulename__ = "filter_exonerate_alignments"

    def __init__(self, collapse: CollapseGTH, portcullis: PortcullisWrapper, masker: RepeatMasking):

        super().__init__()
        self.configuration = collapse.configuration
        self.input = collapse.output
        self.input["junctions"] = portcullis.junctions
        self.input["fai"] = masker.fai
        self.input["genome"] = self.masked_genome
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "filter_gth.log")
        self.output["gff3"] = os.path.join(self.outdir, "gth.filtered.gff3")

    @property
    def rulename(self):
        return self.__rulename__

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]  # TODO: this will need to be updated with the proper details

    @property
    def min_intron(self):
        return self.configuration["homology"]["min_intron"]

    @property
    def max_intron_middle(self):
        return self.configuration["homology"]["max_intron_middle"]

    @property
    def max_intron_ends(self):
        return self.configuration["homology"]["max_intron_ends"]

    @property
    def identity(self):
        return self.configuration["homology"]["identity"]

    @property
    def coverage(self):
        return self.configuration["homology"]["coverage"]

    @property
    def cmd(self):

        load = self.load
        mini = self.min_intron
        maxe, maxm = self.max_intron_ends, self.max_intron_middle
        genome = self.masked_genome
        outdir = self.outdir
        logdir = os.path.dirname(self.log)
        min_coverage, min_identity = self.coverage, self.identity
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && "
        cmd += " filter_exonerate.py -minI {mini} -maxE {maxe} -maxM {maxm} -j {input[junctions]} -g {genome} "
        cmd += " -minid {min_identity} -mincov {min_coverage} {input[gff3]} {output[gff3]} 2> {log} > {log}"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "output")


class GTHProteinWrapper(EIWrapper):

    __final_rulename__ = "gth_done"

    def __init__(self,
                 masker: RepeatMasking, portcullis: PortcullisWrapper):

        super().__init__()
        self.configuration = masker.configuration
        sanitised = SanitizeProteinBlastDB(self.configuration)

        if sanitised.protein_dbs and self.execute and self.use_gth:
            faidx = FaidxProtein(sanitised)
            self.add_edge(sanitised, faidx)
            chunk_proteins = ChunkProteins(sanitised)
            self.add_edge(faidx, chunk_proteins)
            chunks = [GTH(chunk_proteins, chunk, masker) for chunk in range(1, chunk_proteins.chunks + 1)]
            self.add_edges_from([(chunk_proteins, chunk) for chunk in chunks])
            self.add_edges_from([(masker, chunk) for chunk in chunks])
            collapsed = CollapseGTH(chunks)
            self.add_edges_from([(chunk, collapsed) for chunk in chunks])
            filterer = FilterGTH(collapse=collapsed, portcullis=portcullis,
                                 masker=masker)
            self.add_edge(masker, filterer)
            self.add_edge(portcullis, filterer)
            self.add_edge(collapsed, filterer)
            self.add_final_flag()
            assert self.exit

    @property
    def flag_name(self):
        return os.path.join(self.exit.outdir, "proteins.done")

    @property
    def execute(self):
        return self.configuration["homology"].get("execute", True)

    @property
    def use_gth(self):
        return not self.configuration["homology"].get("use_exonerate", False)
