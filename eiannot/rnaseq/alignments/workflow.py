# from .bam import BamStats, BamIndex, BamSort
from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper, GmapLongWrapper
from .tophat2 import TopHat2Wrapper
from .star import StarWrapper, StarLongWrapper
from .abstract import AlnFlag
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
        for tool, wrapper in self.wrappers.items():
            instance = wrapper(self.configuration, self.__prepare_flag)
            instance.finalise()
            # print(tool, instance.bams)
            if len(instance.bams) > 0:
                instances.append(instance)
                flags.append(instance.exit)
                self.merge([instance])
                self.__bams.extend(instance.bams)

        final_flag = ShortAlignersFlag(flags, outdir=self.outdir, prepare_wrapper=prepare_wrapper)
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

    wrappers = {
        "star": StarLongWrapper,
        "gmap": GmapLongWrapper,
    }

    def __init__(self, prepare_wrapper):

        super().__init__()
        self.__gf_rules = []
        self.configuration = prepare_wrapper.configuration
        self.__prepare_flag = prepare_wrapper
        flags = []

        for wrapper in self.wrappers.values():
            instance = wrapper(self.__prepare_flag)
            instance.finalise()
            assert instance.exit
            self.merge([instance])
            self.__gf_rules.extend(instance.gfs)
            flags.append(instance.exit)

        # print(*[flag.rulename for flag in flags])
        flag = LongAlignersFlag(flags, outdir=self.outdir)
        self.add_node(flag)
        self.add_edges_from([(f_flag, flag) for f_flag in flags])
        # print(flag.input)
        self.add_flag_to_inputs(prepare_wrapper, "prep_flag", "fai")

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "4-long-alignments"))

    @property
    def gfs(self):
        return self.__gf_rules


class LongAlignersFlag(AtomicOperation):

    def __init__(self, stats_runs: [AlnFlag], outdir=None):

        super().__init__()
        self.touch = True
        if stats_runs:
            outdir = os.path.dirname(stats_runs[0].output["flag"])
            self.input["flags"] = [stat.output["flag"] for stat in stats_runs]
        else:
            assert outdir is not None
        self.output["flag"] = os.path.join(outdir, "long_reads.done")

    @property
    def rulename(self):
        return "long_aln_all"

    @property
    def loader(self):
        return []


class ShortAlignersFlag(AtomicOperation):

    def __init__(self, flags: [AlnFlag], prepare_wrapper=None, outdir=None):

        super().__init__()
        self.touch = True
        if flags:
            self.input = {"flags": [flag.output["flag"] for flag in flags]}
            outdir = os.path.dirname(flags[0].output["flag"])
        else:
            assert outdir is not None
            assert prepare_wrapper is not None
            self.input = prepare_wrapper.output
        self.output["flag"] = os.path.join(outdir, "short_reads.done")

    @property
    def rulename(self):
        return "short_aln_all"

    @property
    def loader(self):
        return []
