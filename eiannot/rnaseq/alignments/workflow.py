# from .bam import BamStats, BamIndex, BamSort
from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper, GmapLongReads
from .tophat2 import TopHat2Wrapper
from .star import StarWrapper, StarLongWrapper
from .abstract import ShortAlnFlag
from ...abstract import EIWrapper, AtomicOperation
import os


class ShortAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarWrapper,
                "hisat": HisatWrapper,
                "tophat2": TopHat2Wrapper,
                "gsnap": GsnapWrapper}

    def __init__(self, prepare_wrapper: EIWrapper):

        super().__init__()
        self.__bams = []
        self.configuration = prepare_wrapper.configuration
        self.__prepare_flag = prepare_wrapper
        instances = []
        flags = []
        for wrapper in self.wrappers.values():
            instance = wrapper(self.configuration, self.__prepare_flag)
            instance.finalise()
            if len(instance.bams) > 0:
                instances.append(instance)
                flags.append(instance.exit)
                self.merge([instance])
                self.__bams.extend(instance.bams)

        final_flag = AlnFlag(flags, outdir=self.outdir, prepare_wrapper=prepare_wrapper)
        self.add_node(final_flag)  # We have to add the node, otherwise the workflow will be empty
        self.add_edges_from([(flag, final_flag) for flag in flags])
        assert self.exit == final_flag, self.nodes

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments"))

    @property
    def samples(self):
        return self.configuration["short_reads"]

    @property
    def bams(self):
        return self.__bams.copy()


class LongAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarLongWrapper,
                # "gmap": G  # TODO: finish GMAP here
                }

    def __init__(self, prepare_wrapper):

        super().__init__()
        stats = []
        self.__gf_rules = []
        self.configuration = prepare_wrapper.configuration
        prepare_flag = prepare_wrapper.output

        for wrapper in self.wrappers.values():
            instance = wrapper(self.configuration, prepare_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            self.__gf_rules.extend(instance.gfs)
            for gf in instance.gfs:
                stat = LongAlignerStats(gf)
                self.add_edge(gf, stat)
                stats.append(stat)

        flag = LongAlignersFlag(stats, outdir=self.outdir)
        self.add_node(flag)
        self.add_edges_from([(stat, flag) for stat in stats])
        pass

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments"))

    @property
    def gfs(self):
        return self.__gf_rules


class LongAlignerStats(AtomicOperation):

    def __init__(self, aligner):
        super().__init__()
        self.input = aligner.output
        self.outdir = aligner.outdir
        self.output = {"stats": os.path.splitext(self.input["gf"])[0] + ".stats"}
        self.message = "Calculating statistics for: {input[gf]}".format(input=self.input)
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["gf"])

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[gf]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class LongAlignersFlag(AtomicOperation):

    def __init__(self, stats_runs: [LongAlignerStats], outdir=None):

        super().__init__()
        self.touch = True
        if stats_runs:
            outdir = os.path.dirname(os.path.dirname(stats_runs[0].output["stats"]))
        else:
            assert outdir is not None
        self.output["flag"] = os.path.join(outdir, "long_reads.done")

    @property
    def rulename(self):
        return "aln_all"

    @property
    def loader(self):
        return []


class AlnFlag(AtomicOperation):

    def __init__(self, flags: [ShortAlnFlag], prepare_wrapper=None, outdir=None):

        super().__init__()
        self.touch = True
        if flags:
            self.input = {"flags": [flag.output["flag"] for flag in flags]}
            outdir = os.path.dirname(os.path.dirname(flags[0].output["flag"]))
        else:
            assert outdir is not None
            assert prepare_wrapper is not None
            self.input = prepare_wrapper.output
        self.output["flag"] = os.path.join(outdir, "short_reads.done")

    @property
    def rulename(self):
        return "aln_all"

    @property
    def loader(self):
        return []
