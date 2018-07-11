from .abstract import LongAligner, LongWrapper
import os


class MiniMap2Wrapper(LongWrapper):

    __toolname__ = "minimap2"

    def __init__(self, prepare_flag):
        super().__init__(prepare_flag)



    @property
    def prebuilt(self):
        return False

    @property
    def indexer(self):
        """Contrary to other aligners, MiniMap2 does not require a precomputed index"""
        return None


class MiniMap2(LongAligner):

    __toolname__ = "minimap2"

    def __init__(self, sample, run):
        super().__init__( sample=sample, run=run)
        self.output = {"link": self.link,
                       "gf": os.path.join(self.outdir, "minimap." + self.suffix)}

    @property
    def loader(self):
        return ["minimap2"]

    @property
    def cmd(self):

        load = self.load
        extra = self.extra
        input, output = self.input, self.output
        type_args = self.type_args

        cmd = "{load} minimap2 -x splice -c --cs=long {extra} {type_args}"
        cmd += " -C 5 {input[genome]} {input[read1]} |"  # -C 5 :> cost for non-canonical splicing site
        cmd += "k8 $(which paftools.js) splice2bed -m - > {output[bed]} "
        # Now link
        link_src = os.path.relpath(self.output["gf"], start=os.path.dirname(self.output["link"]))
        link_dest = os.path.basename(self.link)
        cmd += " && mkdir -p {outdir} && cd {outdir} && ln -sf {link_src} {link_dest}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "bed12"

    @property
    def type_args(self):
        if self.sample.type in ("ont-direct",):
            return "-k 14"





