# from .bam import BamStats, BamIndex, BamSort
from ...abstract import FinalFlag
from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper, GmapLongWrapper
from .tophat2 import TopHat2Wrapper
from .star import StarWrapper, StarLongWrapper
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

        self.add_final_flag(os.path.join(self.outdir, "short_reads.done"), "short_aln_all")

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
        self.add_final_flag(os.path.join(self.outdir, "long_reads.done"), "long_aln_all")

        # print(flag.input)
        self.add_flag_to_inputs(prepare_wrapper, "prep_flag", "fai")

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "4-long-alignments"))

    @property
    def gfs(self):
        return self.__gf_rules



class ShortAlignersFlag(AtomicOperation):

    def __init__(self, flags: [FinalFlag], prepare_wrapper=None, outdir=None):

        super().__init__()
        self.touch = True
        if flags:
            self.input = {"flags": [flag.output["flag"] for flag in flags]}
            outdir = os.path.dirname(flags[0].output["flag"])
            self.configuration = flags[0].configuration
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
