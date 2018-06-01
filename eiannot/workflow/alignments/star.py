from .. import AtomicOperation, EIWrapper, ShortSample, LongSample
from . import IndexBuilder, ShortAligner
import os
import itertools


class StarIndex(IndexBuilder):

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)

        self.output = {"index": os.path.join(self._outdir, "SAindex")}
        self.touch = False
        self.message = "Indexing genome with star"

    @property
    def toolname(self):
        return "star"

    @property
    def loader(self):
        return ["star"]

    def cmd(self):
        load = self.load
        threads = self.__configuration["threads"]
        input = self.input
        log = self.log
        align_dir = os.path.abspath(os.path.dirname(self.output["index"]))
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        if "ref_transcriptome" in self.input:
            trans = "--sjdbGTFfile {input[ref_transcriptome]}"
        else:
            trans = ""
        extra = self.__configuration["programs"]["star"]["index"]
        cmd = "{load}"
        cmd += "cd {align_dir} && "
        cmd += "STAR --runThreadN {threads} --runMode GenomeGenerate --genomeDir . {trans}"
        cmd += "--genomeFastaFiles {input[genome]} {extra} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class StarAligner(ShortAligner):

    def __init__(self, index, sample, run, configuration, outdir, ref_transcriptome=None):

        super(StarAligner, self).__init__(configuration=configuration,
                                          sample=sample, run=run,
                                          outdir=outdir, index=index, ref_transcriptome=ref_transcriptome)

    @property
    def compression_option(self):
        if self.sample.suffix == ".gz":
            return "--readFilesCommand zcat"
        elif self.sample.suffix in (".bz", ".bz2"):
            return "--readFilesCommand bzcat"
        else:
            return ""

    @property
    def input_reads(self):
        read1 = os.path.abspath(self.input["read1"])
        read2 = os.path.abspath(self.input["read2"])

        snippet = "--readFilesIn {read1}"
        if read2:
            snippet += " {read2}"
        snippet = snippet.format(**locals())
        return snippet

    @property
    def cmd(self):

        load = self.load
        index = os.path.dirname(self.input["index"])
        outdir = self.bamdir
        cmd = "{load}"
        cmd += "cd {outdir} &&"
        cmd += "STAR --runThreadN {threads} --runMode alignReads --genomeDir {index}"
        rfc = self.compression_option
        infiles = self.input_reads
        cmd += "{rfc} --runMode alignReads {infiles}"
        cmd += "--outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM XS NM MD"
        cmd += "--outSAMstrandField intronMotif"
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += "--alignIntronMin {min_intron} --alignIntronMax {max_intron} --alignMatesGapMax {max_intron}"
        if self.ref_transcriptome:
            ref_transcriptome = self.ref_transcriptome
            cmd += "--sjdbGTFfile {ref_transcriptome}"
        extra = self.extra
        bamdir = self.bamdir
        log = self.log
        cwd = os.getcwd()  # TODO: this is probably *not* what we need
        cmd += "--outFileNamePrefix {bamdir} {extra} > {log} 2>&1"
        cmd += " && cd {cwd} && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["star"]

    @property
    def toolname(self):
        return "star"
