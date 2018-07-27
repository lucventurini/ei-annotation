from .abstract import LongAligner, LongWrapper, IndexBuilder, IndexLinker
from ...preparation import PrepareWrapper
import os
import itertools
import re


class MiniMap2Wrapper(LongWrapper):

    __toolname__ = "minimap2"

    def __init__(self, prepare_flag: PrepareWrapper):
        super().__init__(prepare_flag)

        if len(self.runs) > 0 and len(self.samples) > 0:
            indexer = self.indexer(prepare_flag.configuration, prepare_flag)
            self.add_edge(prepare_flag.exit, indexer)
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                mini_run = MiniMap2(indexer=indexer, sample=sample, run=run)
                self.add_edge(indexer, mini_run)
                self.add_to_gfs(mini_run)

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

        self.configuration = configuration
        super().__init__(configuration, self.outdir)
        self.input["genome"] = self.genome
        self.input["fai"] = prepare_flag.fai.output["fai"]
        self.output["index"] = os.path.join(self.outdir, "genome.idx")
        self.log = os.path.join(self.outdir, "index.log")

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", "minimap2", "splice")

    @property
    def loader(self):
        return ["minimap2", "ei-annotation"]

    @property
    def cmd(self):
        threads = self.threads
        load = self.load
        log = self.log
        input, output = self.input, self.output
        outdir = self.outdir
        cmd = "{load} mkdir -p {outdir} && minimap2 -x splice -I $(determine_genome_size.py -G {input[genome]}) "
        cmd += " -t {threads} -d {output[index]} {input[genome]} > {log} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def index(self):
        return self.output["index"]


class MiniMap2(LongAligner):

    __toolname__ = "minimap2"

    def __init__(self, indexer, sample, run):
        super().__init__(indexer=indexer, sample=sample, run=run)
        # self.input["index"] = self.indexer.output["index"]
        self.output = {"link": self.link,
                       "gf": self.bed12,
                       "paf": os.path.join(
                           os.path.dirname(self.bed12),
                           re.sub("\.bed12", "", os.path.basename(self.bed12)) + ".paf.gz")}
        self.log = os.path.join(os.path.dirname(self.bed12), "minimap.log")

    @property
    def loader(self):
        return ["minimap2", "ei-annotation"]

    @property
    def strand_option(self):
        if self.sample.stranded:
            return " -u f"
        else:
            return "-u b"

    @property
    def cmd(self):

        load = self.load
        extra = self.extra
        input, output = self.input, self.output
        type_args = self.type_args
        outdir = self.outdir
        log = self.log
        paf = re.sub("\.gz$", "", self.output["paf"])
        max_intron, strand_option = self.max_intron, self.strand_option
        threads = self.threads
        genome = self.genome

        cmd = "{load} mkdir -p {outdir} && minimap2 -x splice -c --cs=long {extra} {type_args}"
        cmd += " -G {max_intron} {strand_option} -t {threads} "
        cmd += " -I $(determine_genome_size.py -G {genome}) "
        cmd += " -C 5 {input[index]} {input[read1]} 2> {log} > {paf} "  # -C 5 :> cost for non-canonical splicing site
        cmd += " && k8 $(which paftools.js) splice2bed -m {paf} | "
        # Needed to correct for the fact that minimap2 BED12
        cmd += " correct_bed12_mappings.py > {output[gf]} && gzip {paf}"
        # Now link
        link_dir = os.path.dirname(self.link)
        link_src = os.path.relpath(self.output["gf"], start=os.path.dirname(self.output["link"]))
        link_dest = os.path.basename(self.link)
        cmd += " && mkdir -p {link_dir} && cd {link_dir} && ln -sf {link_src} {link_dest}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return ".bed12"

    @property
    def type_args(self):
        if self.sample.type in ("ont-direct",):
            return "-k 14"
        else:
            return ""

    @property
    def rulename(self):
        return "minimap2_{sample.label}_{run}".format(sample=self.sample, run=self.run)

    @property
    def bed12(self):
        return os.path.join(self.outdir, "minimap2", "{label}-{run}", "minimap2.bed12").format(
            label=self.sample.label,
            run=self.run
        )
