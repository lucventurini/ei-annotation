from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats, FilterGF
import itertools


class ScallopWrapper(ShortAssemblerWrapper):

    __toolname__ = "scallop"
    __tag__ = "RPKM"  # We are not going to use it here

    def __init__(self, aln_wrapper):

        super().__init__(aln_wrapper)

        if len(self.runs) > 0 and len(self.bams) > 0:
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                scallop = Scallop(bam, run, create_link=False)
                scallop_filter = FilterGF(scallop,
                                          tag=self.tag,
                                          monoexonic_threshold=self.mono_threshold,
                                          multiexonic_threshold=self.multi_threshold)
                self.add_edge(scallop, scallop_filter)
                stat = AsmStats(scallop_filter)
                self.add_edge(aln_wrapper, scallop)
                self.add_edge(scallop_filter, stat)
                self.add_to_gf(stat)
                continue


class Scallop(ShortAssembler):

    __toolname__ = "scallop"

    def __init__(self, bam, run, create_link):

        super().__init__(bam, run, create_link=create_link)

    @property
    def loader(self):
        return ["scallop"]

    @property
    def threads(self):
        return 1

    @property
    def strand(self):
        if self.sample.strandedness in ("f", "fr-secondstrand"):
            return "--library_type second"
        elif self.sample.strandedness in ("r", "fr-firststrand"):
            return "--library_type first"
        elif self.sample.strandedness == "fr-unstranded" or self.sample.strandedness is None:
            return "--library_type unstranded"
        else:
            raise ValueError(self.sample.strandedness)

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        extra = self.extra
        threads = self.threads
        outdir = self.gfdir
        link_src = self.link_src
        input = self.input
        output = self.output
        cmd += "mkdir -p {outdir} && "
        strand = self.strand
        log = self.log
        cmd += "scallop -i {input[bam]} -o {output[gf]} {strand} {extra} > {log} 2>&1"
        if self._create_link is True:
            cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def input_reads(self):
        return ""
