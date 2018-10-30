from .abstract import ShortAssemblerWrapper
from .class2 import Class2Wrapper
from .stringtie import StringtieWrapper
from .cufflinks import CufflinksWrapper
from .trinity_gg import TrinityGGWrapper
from .scallop import ScallopWrapper
from ...abstract import EIWrapper, EIWorfkflow
from ..alignments.__init__ import ShortAlignmentsWrapper
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
        self.__short_alignments = set()
        for tool, wrapper in self.wrappers.items():
            instance = wrapper(short_alignments)
            if instance.gfs:
                self.add_edge(short_alignments, instance)
                instance.add_flag_to_inputs(short_alignments, "aln_flag", "flag")
                self.__add_to_bams(instance)
                self.__add_to_gfs(instance)

        self.__short_alignments = list(self.__short_alignments)
        self.add_final_flag()

    __final_rulename__ = "asm_all"

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "all.done")

    @property
    def gfs(self):
        return self.__gfs

    @property
    def short_alignments(self):
        return self.__short_alignments

    def __add_to_gfs(self, wrapper: ShortAssemblerWrapper):
        self.__gfs.extend(wrapper.gfs)

    def __add_to_bams(self, wrapper: ShortAssemblerWrapper):
        self.__short_alignments.update(set(wrapper.bams))

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))
