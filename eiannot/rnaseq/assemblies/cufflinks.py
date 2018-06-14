from .abstract import ShortAssembler, ShortAssemblerWrapper
from ...abstract import AtomicOperation, EIWrapper, ShortSample
import os
import itertools


class CufflinksWrapper(ShortAssemblerWrapper):

    def __init__(self, align_wrapper):
        super().__init__(align_wrapper)

        if len(self.bams) > 0 and len(self.runs) > 0:
            cuffs = []
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                cufflinks = Cufflinks(bam, run, self.configuration, self.outdir)
                cuffs.append(cufflinks)
                self.add_to_gf(cuffs)
            flag = CufflinksFlag(cuffs, self.outdir)
            self.add_edges_from([(cuff, flag) for cuff in cuffs])
            return

    @property
    def toolname(self):
        return "cufflinks"


class CufflinksFlag(AtomicOperation):

    def __init__(self, cufflinks, outdir):
        super().__init__()
        self.input = {"gtfs": [cuff.output["link"] for cuff in cufflinks]}
        self.touch = True
        self.output = {"flag": os.path.join(outdir, "cufflinks.done")}


class Cufflinks(ShortAssembler):

    def __init__(self, bam, run, configuration, outdir):

        super().__init__(bam, run, configuration, outdir)
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
