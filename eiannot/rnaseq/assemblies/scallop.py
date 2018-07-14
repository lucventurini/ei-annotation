from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats
from ...abstract import AtomicOperation, EIWrapper, ShortSample
import os
import itertools


class ScallopWrapper(ShortAssemblerWrapper):

    __toolname__ = "scallop"

    def __init__(self, aln_wrapper):

        super().__init__(aln_wrapper)

        if len(self.runs) > 0 and len(self.bams) > 0:
            scallops = []
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                scallop = Scallop(bam, run)
                stat = AsmStats(scallop)
                self.add_edge(aln_wrapper, scallop)
                self.add_edge(scallop, stat)
                self.add_to_gf(stat)
                scallops.append(stat)
                continue


class Scallop(ShortAssembler):

    __toolname__ = "scallop"

    def __init__(self, bam, run):

        super().__init__(bam, run)

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
        cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def input_reads(self):
        return ""
