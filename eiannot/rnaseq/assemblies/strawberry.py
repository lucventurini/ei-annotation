from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats, FilterGF
import os
import itertools


class StrawberryWrapper(ShortAssemblerWrapper):

    __toolname__ = "strawberry"
    __tag__ = "TPM"

    def __init__(self, aln_wrapper):
        super().__init__(aln_wrapper)

        if len(self.bams) > 0 and len(self.runs) > 0:
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                berry = Strawberry(bam, run, create_link=False)
                straw_filter = FilterGF(berry,
                                        tag=self.tag,
                                        monoexonic_threshold=self.mono_threshold,
                                        multiexonic_threshold=self.multi_threshold)
                stat = AsmStats(straw_filter)
                self.add_edge(berry, straw_filter)
                self.add_edge(straw_filter, stat)
                self.add_to_gf(stat)
            return


class Strawberry(ShortAssembler):

    __toolname__ = "strawberry"

    def __init__(self, bam, run, create_link=False):

        super().__init__(bam, run, create_link=create_link)
        self.output["gf"] = os.path.join(self.gfdir, "assembled_transcripts.gtf")

    @property
    def loader(self):
        return ["strawberry"]

    @property
    def cmd(self):
        load = self.load
        cmd = " {load} mkdir -p {gfdir} && "
        extra = self.extra
        # # Default options for Stringtie. We do not want to overcome those provided by
        # if not re.search("\-f ([0-9]|\.)*", extra):
        #     extra += " -f 0.05 "
        # if not re.search("\-m [0-9]*", extra):
        #     extra += " -m 200 "

        threads = self.threads
        outdir = self.gfdir
        link_src = self.link_src
        input = self.input
        output = self.output
        log = self.log
        gtf = os.path.join(self.gfdir, "transcripts.gtf")
        gfdir = self.gfdir
        min_intron, max_intron = self.min_intron, self.max_intron
        if self.ref_transcriptome is not None:
            trans = " -g {ref_transcriptome} ".format(ref_transcriptome=self.ref_transcriptome)
        else:
            trans = " "
        alrun = self.alrun
        strand = self.strand
        cmd += " strawberry {extra} {trans} -o {gfdir} -j {min_intron} -J {max_intron} {strand} "
        cmd += " -p {threads} {input[bam]} > {log} 2>&1 "
        if self._create_link is True:
            cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def strand(self):
        if not self.sample.stranded:
            return ''
        elif self.sample.strandedness == "fr-firststrand":
            return "--rf"
        elif self.sample.strandedness == "fr-secondstrand":
            return "--fr"

    @property
    def input_reads(self):
        return ""
