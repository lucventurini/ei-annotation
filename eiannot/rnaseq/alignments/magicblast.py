from .abstract import LongAligner, LongWrapper, IndexBuilder
from ...abstract import LongSample
from ...preparation import PrepareWrapper
import os
import itertools
import re


__doc__ = """This module implements long-read alignment for MagicBlast, a new implementation by NCBI which is
splice-aware."""


class MagicBlastWrapper(LongWrapper):

    __toolname__ = "magicblast"

    def __init__(self, prepare_flag: PrepareWrapper):
        super().__init__(prepare_flag)

        if len(self.runs) > 0 and len(self.samples) > 0:
            # First, we have to collect the runs that we have to perform.
            indexer = self.indexer(prepare_flag.configuration)
            self.add_edge(prepare_flag.exit, indexer)
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                mb_run = MagicBlast(indexer=indexer, sample=sample, run=run)
                self.add_edge(indexer, mb_run)
                mb_convert = MagicBlastConvert(aligner=mb_run)
                self.add_edge(mb_run, mb_convert)
                self.add_to_gfs(mb_convert)

    @property
    def _do_stats(self):
        return False

    @property
    def prebuilt(self):
        return False

    @property
    def indexer(self):
        """Contrary to other aligners, MiniMap2 does not require a precomputed index"""
        return MagicBlastIndexer


class MagicBlastIndexer(IndexBuilder):

    __toolname__ = "magicblast"

    def __init__(self, configuration):

        super().__init__(configuration)
        self.input["genome"] = self.genome
        self.output = {"flag": os.path.join(self.outdir, "magicblast_index.done")}

    @property
    def message(self):
        return "Building the MagicBlast index for {} in {}".format(self.input["genome"],
                                                                   self.outdir)

    @property
    def cmd(self):
        input, output = self.input, self.output
        load = self.load
        index = self.index
        log = self.log
        cmd = "{load} makeblastdb -dbtype nucl -parse_seqids -in {input[genome]} -out {index} 2> {log} > {log} "
        cmd += "&& touch {output[flag]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["magicblast"]

    @property
    def threads(self):
        return 1

    @property
    def index(self):
        return os.path.join(self.outdir, os.path.basename(self.genome))


class MagicBlast(LongAligner):

    __toolname__ = "magicblast"

    def __init__(self, indexer: MagicBlastIndexer, sample: LongSample, run):
        super().__init__(indexer=indexer, sample=sample, run=run)
        self.output = {"bam": self.bam}
        self.log = os.path.join(os.path.dirname(self.bam), "magicblast.log")

    @property
    def loader(self):
        return ["magicblast", "samtools"]

    @property
    def cmd(self):

        load = self.load
        extra = self.extra
        input, output = self.input, self.output
        outdir = self.outdir
        log = self.log
        max_intron, strand_option = self.max_intron, self.strand_option
        threads = self.threads
        genome = self.genome
        cmd = "{load} mkdir -p {outdir} && "
        infmt = self.infmt
        index = self.indexer.index
        threads = self.threads
        cmd += "magicblast -query {input[read1]} {infmt} -db {index} -outfmt sam -num_threads {threads} "
        identity = self.identity
        log = self.log
        extra = self.extra
        max_intron = self.max_intron
        strand = self.strand_option
        cmd += "{identity} {extra} -max_intron_length {max_intron} -splice T -reftype genome {strand}"
        cmd += " | samtools view -bS - | "
        cmd += "samtools sort -@ {threads} --reference {genome} -T {output[bam]}.sort -o {output[bam]} -"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def strand_option(self):
        if self.sample.stranded:
            return "-fr"
        else:
            return ""

    @property
    def infmt(self):
        if self.sample.fileformat in ("asn1", "asn1b", "fastc", "fastq", "fasta"):
            return "-infmt {}".format(self.sample.fileformat)
        else:
            return "-infmt fasta"

    @property
    def identity(self):
        if "identity" in self.configuration["programs"][self.toolname]:
            return "-perc_identity {}".format(self.configuration["programs"][self.toolname]["identity"])
        else:
            return ""

    @property
    def rulename(self):
        return "{toolname}_{sample.label}_{run}".format(sample=self.sample, run=self.run, toolname=self.toolname)

    @property
    def suffix(self):
        return ".bam"

    @property
    def bam(self):
        return os.path.join(self.outdir, self.toolname, "{label}-{run}", "magicblast.bam").format(
            label=self.sample.label,
            run=self.run
        )


class MagicBlastConvert(LongAligner):

    __toolname__ = "magicblast"

    def __init__(self, aligner: MagicBlast):

        super().__init__(indexer=aligner.indexer, sample=aligner.sample, run=aligner.run)
        self.aligner = aligner
        self.input = self.aligner.output
        self.output = {"gf": self.gtf,
                       "link": self.link,}

    @property
    def rulename(self):
        return self.aligner.rulename + "_convert"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):

        load = self.load
        input, output = self.input, self.output
        outdir = self.outdir
        cmd = "{load} mkdir -p {outdir} && bam2gtf.py {input[bam]} {output[gf]}"
        # Now link
        link_dir = os.path.dirname(self.link)
        link_src = os.path.relpath(self.output["gf"], start=os.path.dirname(self.output["link"]))
        link_dest = os.path.basename(self.link)
        cmd += " && mkdir -p {link_dir} && cd {link_dir} && ln -sf {link_src} {link_dest}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def suffix(self):
        return ".gtf"

    @property
    def is_small(self):
        return True

    @property
    def gtf(self):
        return os.path.join(self.outdir, self.toolname, "{label}-{run}",
                            "magicblast.gtf").format(label=self.sample.label, run=self.run)
