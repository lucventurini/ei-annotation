from .abstract import LongAligner, LongWrapper
import os
import itertools


class MiniMap2Wrapper(LongWrapper):

    __toolname__ = "minimap2"

    def __init__(self, prepare_flag):
        super().__init__(prepare_flag)

        if len(self.runs) > 0 and len(self.samples) > 0:
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                mini_run = MiniMap2(indexer=prepare_flag, sample=sample, run=run)
                self.add_to_gfs(mini_run)
                self.add_edge(prepare_flag, mini_run)

    @property
    def prebuilt(self):
        return False

    @property
    def indexer(self):
        """Contrary to other aligners, MiniMap2 does not require a precomputed index"""
        return None


class MiniMap2(LongAligner):

    __toolname__ = "minimap2"

    def __init__(self, indexer, sample, run):
        super().__init__(indexer=indexer, sample=sample, run=run)
        self.input["genome"] = self.genome
        self.output = {"link": self.link,
                       "gf": self.bed12}
        self.log = os.path.join(os.path.dirname(self.bed12), "minimap.log")

    @property
    def loader(self):
        return ["minimap2"]

    @property
    def cmd(self):

        load = self.load
        extra = self.extra
        input, output = self.input, self.output
        type_args = self.type_args
        outdir = self.outdir
        log = self.log

        cmd = "{load} mkdir -p {outdir} && minimap2 -x splice -c --cs=long {extra} {type_args}"
        cmd += " -C 5 {input[genome]} {input[read1]} 2> {log} |"  # -C 5 :> cost for non-canonical splicing site
        cmd += "k8 $(which paftools.js) splice2bed -m - > {output[gf]} "
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
