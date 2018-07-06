from .abstract import ShortAssemblerWrapper
from .class2 import Class2Wrapper
from .stringtie import StringtieWrapper
from .cufflinks import CufflinksWrapper
from .trinity_gg import TrinityGGWrapper
from .scallop import ScallopWrapper
from ...abstract import EIWrapper, EIWorfkflow
from ..alignments.workflow import ShortAlignmentsWrapper
from .strawberry import StrawberryWrapper
import os


class AssemblyWrapper(EIWrapper):

    wrappers = {"cufflinks": CufflinksWrapper,
                "scallop": ScallopWrapper,
                "trinity": TrinityGGWrapper,
                "class2": Class2Wrapper,
                "stringtie": StringtieWrapper,
                "strawberry": StrawberryWrapper}

    def __init__(self, short_alignments: ShortAlignmentsWrapper):

        super().__init__(configuration=short_alignments.configuration)

        self.__gfs = []
        for tool, wrapper in self.wrappers.items():
            instance = wrapper(short_alignments)
            if instance.gfs:
                self.add_edge(short_alignments, instance)
                instance.add_flag_to_inputs(short_alignments, "aln_flag", "flag")
                self.__add_to_gfs(instance)

        # print(self.gfs)
        self.add_final_flag(os.path.join(self.outdir, "all.done"), rulename="asm_all")

    @property
    def gfs(self):
        return self.__gfs

    def __add_to_gfs(self, wrapper: ShortAssemblerWrapper):
        self.__gfs.extend(wrapper.gfs)

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))
