from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats
import os
import itertools


class CufflinksWrapper(ShortAssemblerWrapper):

    def __init__(self, align_wrapper):
        super().__init__(align_wrapper)

        if len(self.bams) > 0 and len(self.runs) > 0:
            cuffs = []
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                cufflinks = Cufflinks(bam, run)
                cuffs.append(cufflinks)
                stat = AsmStats(cufflinks)
                self.add_edge(cufflinks, stat)
                self.add_to_gf(stat)
                # self.add_to_gf(cuffs)


    @property
    def toolname(self):
        return "cufflinks"


class Cufflinks(ShortAssembler):

    def __init__(self, bam, run):

        super().__init__(bam, run)
        # Cufflinks has a set name for the output
        self.output["gf"] = os.path.join(self.gfdir, "transcripts.gtf")

    @property
    def toolname(self):
        return "cufflinks"

    @property
    def loader(self):
        return ["cufflinks"]

    @property
    def cmd(self):
        load = self.load
        strand = self.strand
        cmd = "{load}"
        extra = self.extra
        threads = self.threads
        outdir = self.gfdir
        link_src = self.link_src
        input = self.input
        output = self.output
        log = self.log
        gtf = os.path.join(self.gfdir, "transcripts.gtf")
        min_intron = self.min_intron
        max_intron = self.max_intron
        if self.ref_transcriptome is not None:
            trans = "--GTF-guide={ref_transcriptome}".format(ref_transcriptome=self.ref_transcriptome)
        else:
            trans = ""
        cmd += "cufflinks --output-dir={outdir} --num-threads={threads} {trans} {strand}"
        cmd += " --min-intron-length={min_intron} --max-intron-length={max_intron}"
        cmd += " --no-update-check {extra} {input[bam]} > {log} 2>&1 && "
        cmd += "ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def strand(self):
        if self.sample.strandedness:
            return '--library-type={}'.format(self.sample)
        else:
            return ""

    @property
    def input_reads(self):
        raise NotImplementedError()
