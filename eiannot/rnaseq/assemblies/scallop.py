from .abstract import ShortAssembler, ShortAssemblerWrapper
from ...abstract import AtomicOperation, EIWrapper, ShortSample
import os
import itertools


class ScallopWrapper(ShortAssemblerWrapper):

    def __init__(self, aln_wrapper):

        super().__init__(aln_wrapper)

        if len(self.runs) > 0 and len(bams) > 0:
            scallops = []
            for bam, run in itertools.product(bams, range(len(self.runs))):
                scallop = Scallop(bam, run, configuration, self.outdir)
                scallops.append(scallop)
                self.add_to_gf(scallop)
                continue
            flag = ScallopFlag(scallops, self.outdir)
            self.add_edges_from([(scallop, flag) for scallop in scallops])


class ScallopFlag(AtomicOperation):

    def __init__(self, scallops, outdir):
        super().__init__()
        self.input = {"gtfs": [scallop.output["link"] for scallop in scallops]}
        self.touch = True
        self.output = {"flag": os.path.join(outdir, "scallop.done")}


class Scallop(ShortAssembler):

    def __init__(self, bam, run, configuration, outdir):

        super().__init__(bam, run, configuration, outdir)

    @property
    def toolname(self):
        return "scallop"

    @property
    def loader(self):
        return ["scallop"]

    @property
    def threads(self):
        return 1

    @property
    def strand(self):
        if self.sample.strandedness in ("f", "fr-secondstrand"):
            return "--library_type second"
        elif self.sample.strandedness in ("r", "fr-secondstrand"):
            return "--library_type first"
        elif self.sample.strandedness == "fr-unstranded" or self.sample.strandedness is None:
            return "--library_type unstranded"
        else:
            raise ValueError(self.sample.strandedness)

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        extra = self.extra
        threads = self.threads
        outdir = self.gfdir
        link_src = self.link_src
        input = self.input
        output = self.output
        cmd += "mkdir -p {outdir} && "
        strand = self.strand
        log = self.log
        cmd += "scallop -i {input[bam]} -o {output[gf]} {strand} {extra} > log 2>&1"
        cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"
