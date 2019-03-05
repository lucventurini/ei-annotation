from .abstract import LongSample, LongAligner, LongWrapper, IndexBuilder
from ...preparation import PrepareWrapper
import os
import itertools
import re


__doc__ = """Minimap2 is a strange aligner, in the sense that it does not use - like HISAT2, STAR, etc. - a single index file
which can then be queried using different parameters. Rather, it tries to index everything in memory for every run,
and requires to re-index if one of the key parameters (-k, -w or -H) is changed by any one of the runs.
This requires, for the sake of consistency, to dicth pre-indexing altogether."""


_nanopore_types = ("ont-direct",)


class MiniMap2Wrapper(LongWrapper):

    __toolname__ = "minimap2"

    def __init__(self, prepare_flag: PrepareWrapper):
        super().__init__(prepare_flag)

        if len(self.runs) > 0 and len(self.samples) > 0:
            # First, we have to collect the runs that we have to perform.
            indexer = self.indexer(prepare_flag.configuration, prepare_flag)
            self.add_edge(prepare_flag.exit, indexer)
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                mini_run = MiniMap2(indexer=indexer, sample=sample, run=run)
                self.add_edge(indexer, mini_run)
                mini_convert = Minimap2Convert(aligner=mini_run)
                self.add_edge(mini_run, mini_convert)
                self.add_to_gfs(mini_convert)

    @property
    def _do_stats(self):
        return False

    @property
    def prebuilt(self):
        return False

    @property
    def indexer(self):
        """Contrary to other aligners, MiniMap2 does not require a precomputed index"""
        return MiniMap2SpliceIndexer


class MiniMap2SpliceIndexer(IndexBuilder):

    __toolname__ = "minimap2"

    def __init__(self, configuration, prepare_flag: PrepareWrapper):

        """This is a mock index."""

        self.configuration = configuration
        super().__init__(configuration)
        self.input["genome"] = self.genome
        self.input["fai"] = prepare_flag.fai.output["fai"]
        self.output["index"] = os.path.join(self.outdir, "genome.fa")
        self.log = os.path.join(self.outdir, "index.log")

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", "minimap2", "splice")

    @property
    def loader(self):
        return []

    @property
    def cmd(self):
        threads = self.threads
        load = self.load
        log = self.log
        input, output = self.input, self.output
        outdir = self.outdir
        genome = os.path.abspath(self.input["genome"])
        cmd = "mkdir -p {outdir} && ln -rs {genome} {output[index]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def index(self):
        return self.output["index"]

    @property
    def threads(self):
        return 1

    @property
    def local(self):
        return True

    @property
    def is_small(self):
        return True


class MiniMap2(LongAligner):

    __toolname__ = "minimap2"

    def __init__(self, indexer: MiniMap2SpliceIndexer, sample: LongSample, run):
        super().__init__(indexer=indexer, sample=sample, run=run)
        self.input["index"] = self.indexer.output["index"]
        self.output = {"bam": self.bam}
        self.log = os.path.join(os.path.dirname(self.bam), "minimap.log")

    @property
    def loader(self):
        return ["minimap2", "ei-annotation", "samtools"]

    @property
    def strand_option(self):
        if self.sample.stranded:
            return "-u f"
        else:
            return "-u b"

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
        noncanonical_cost = self.noncanonical_cost
        cmd = "{load} mkdir -p {outdir} && minimap2 -x splice -c --cs=long {extra}"
        cmd += " -G {max_intron} {strand_option} -t {threads} "
        cmd += " -I $(determine_genome_size.py -G {genome}) -a "
        cmd += " {noncanonical_cost} {genome} {input[read1]} 2> {log} | "  # -C 5 :> cost for non-canonical splicing site
        cmd += " samtools view -bS - | "
        cmd += "samtools sort -@ {threads} --reference {genome} -T {output[bam]}.sort -o {output[bam]} -"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def rulename(self):
        return "minimap2_{sample.label}_{run}".format(sample=self.sample, run=self.run)

    @property
    def noncanonical_cost(self):
        if self.sample.type in ("cdna", "est"):
            return " -C 5 "
        else:
            return " "

    @property
    def suffix(self):
        return ".bam"

    @property
    def bam(self):
        return os.path.join(self.outdir, "minimap2", "{label}-{run}", "minimap2.bam").format(
            label=self.sample.label,
            run=self.run
        )


class Minimap2Convert(LongAligner):

    __toolname__ = "minimap2"

    def __init__(self, aligner: MiniMap2):

        super().__init__(indexer=aligner.indexer, sample=aligner.sample, run=aligner.run)
        self.aligner = aligner
        self.input = self.aligner.output
        self.output = {"gf": self.bed12,
                       "link": self.link,
                       "paf": os.path.join(
                           os.path.dirname(self.bed12),
                           re.sub("\.bed12", "", os.path.basename(self.bed12)) + ".paf.gz")
                       }

    @property
    def rulename(self):
        return self.aligner.rulename + "_convert"

    @property
    def loader(self):
        return ["minimap2", "ei-annotation", "mikado", "samtools"]

    @property
    def cmd(self):

        load = self.load
        input, output = self.input, self.output
        cmd = "{load} k8 $(which paftools.js) sam2paf <(samtools view -h {input[bam]}) | gzip -c - > {output[paf]}"
        cmd += " && k8 $(which paftools.js) splice2bed -m <(samtools view -h {input[bam]}) | "
        # Needed to correct for the fact that minimap2 BED12
        cmd += " correct_bed12_mappings.py > {output[gf]}"
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
        return ".bed12"

    @property
    def bed12(self):
        return os.path.join(self.outdir, "minimap2", "{label}-{run}", "minimap2.bed12").format(
            label=self.sample.label,
            run=self.run
        )

    @property
    def is_small(self):
        return False
