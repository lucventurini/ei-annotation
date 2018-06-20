from .abstract import AsmStats, AsmFlag, ShortAssemblerWrapper
from .class2 import Class2Wrapper
from .stringtie import StringtieWrapper
from .cufflinks import CufflinksWrapper
from .trinity_gg import TrinityGGWrapper
from .scallop import ScallopWrapper
from ...abstract import EIWrapper, EIWorfkflow
from ..alignments.workflow import ShortAlignmentsWrapper
import os


class AssemblyWrapper(EIWrapper):

    wrappers = {"cufflinks": CufflinksWrapper,
                "scallop": ScallopWrapper,
                "trinity": TrinityGGWrapper,
                "class2": Class2Wrapper,
                "stringtie": StringtieWrapper}

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
        final_flag = AsmFlag(self.gfs, outdir=self.outdir)
        self.add_node(final_flag)
        self.add_edges_from([(gf, final_flag) for gf in self.gfs])
        # print(self.exit.input, self.output["flag"])

    @property
    def gfs(self):
        return self.__gfs

    def __add_to_gfs(self, wrapper: ShortAssemblerWrapper):
        self.__gfs.extend(wrapper.gfs)

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))
