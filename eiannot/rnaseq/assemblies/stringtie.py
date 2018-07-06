from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats
from ...abstract import AtomicOperation, EIWrapper, ShortSample
import os
import itertools
import re


class StringtieWrapper(ShortAssemblerWrapper):

    __toolname__ = "stringtie"

    def __init__(self, aln_wrapper):
        super().__init__(aln_wrapper)

        if len(self.bams) > 0 and len(self.runs) > 0:
            stringties = []
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                stringtie = Stringtie(bam, run)
                stringties.append(stringtie)
                stat = AsmStats(stringtie)
                self.add_edge(stringtie, stat)
                self.add_to_gf(stat)
            # flag = StringtieFlag(stringties, self.outdir)
            # self.add_node(flag)
            # self.add_edges_from([(stringtie, flag) for stringtie in stringties])
            return

    @property
    def toolname(self):
        return "stringtie"


class Stringtie(ShortAssembler):

    __toolname__ = "stringtie"

    def __init__(self, bam, run):

        super().__init__(bam, run)

    @property
    def loader(self):
        return ["stringtie"]

    @property
    def cmd(self):
        load = self.load
        cmd = " {load} "
        extra = self.extra
        # Default options for Stringtie. We do not want to overcome those provided by
        if not re.search("\-f ([0-9]|\.)*", extra):
            extra += " -f 0.05 "
        if not re.search("\-m [0-9]*", extra):
            extra += " -m 200 "

        threads = self.threads
        outdir = self.gfdir
        link_src = self.link_src
        input = self.input
        output = self.output
        log = self.log
        gtf = os.path.join(self.gfdir, "transcripts.gtf")
        if self.ref_transcriptome is not None:
            trans = " -G {ref_transcriptome} ".format(ref_transcriptome=self.ref_transcriptome)
        else:
            trans = " "
        alrun = self.alrun
        cmd += " stringtie {input[bam]} -l Stringtie_{alrun} {extra} "
        cmd += " {trans} -o {output[gf]} -p {threads} > {log} 2>&1 && "
        cmd += "ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def strand(self):
        return ''

    @property
    def input_reads(self):
        return ""
