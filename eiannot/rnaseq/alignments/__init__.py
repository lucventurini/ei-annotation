from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper, GmapLongWrapper
from .tophat2 import TopHat2Wrapper
from .star import StarWrapper, StarLongWrapper
from ...abstract import EIWrapper
from ...preparation import PrepareWrapper
from .minimap2 import MiniMap2Wrapper
import os


class ShortAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarWrapper,
                "hisat": HisatWrapper,
                "tophat2": TopHat2Wrapper,
                "gsnap": GsnapWrapper}

    def __init__(self, prepare_wrapper: PrepareWrapper):

        super().__init__()
        self.__bams = []
        self.configuration = prepare_wrapper.configuration
        self.__prepare_flag = prepare_wrapper

        instances = []
        flags = []
        for tool, wrapper in self.wrappers.items():
            instance = wrapper(self.configuration, self.prepare_wrapper)
            instance.finalise()

            # print(tool, instance.bams)
            if len(instance.bams) > 0:
                instances.append(instance)
                flags.append(instance.exit)
                self.merge([instance])
                self.__bams.extend(instance.bams)

        self.add_flag_to_inputs(prepare_wrapper, "prepare_flag", "flag")
        self.add_final_flag()

    __final_rulename__ = "short_aln_all"

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "short_reads.done")

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments"))

    @property
    def samples(self):
        return self.configuration["short_reads"]

    @property
    def bams(self):
        return self.__bams.copy()

    @property
    def prepare_wrapper(self):
        return self.__prepare_flag


class LongAlignmentsWrapper(EIWrapper):

    wrappers = {
        "star": StarLongWrapper,
        "gmap": GmapLongWrapper,
        "minimap2": MiniMap2Wrapper
    }

    def __init__(self, prepare_wrapper):

        super().__init__()
        self.__gf_rules = []
        self.configuration = prepare_wrapper.configuration
        self.__prepare_flag = prepare_wrapper
        flags = []

        for wrapper in self.wrappers.values():
            instance = wrapper(self.prepare_wrapper)
            instance.finalise()
            assert instance.exit
            self.merge([instance])
            self.__gf_rules.extend(instance.gfs)
            flags.append(instance.exit)

        # print(*[flag.rulename for flag in flags])
        self.add_final_flag()

        # print(flag.input)
        self.add_flag_to_inputs(prepare_wrapper, "prepare_flag", "flag")

    __final_rulename__ = "long_aln_all"

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "long_reads.done")

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "4-long-alignments"))

    @property
    def gfs(self):
        return self.__gf_rules

    @property
    def prepare_wrapper(self):
        return self.__prepare_flag