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

        stats = []
        self.__gfs = []
        for wrapper in self.wrappers.values():
            instance = wrapper(short_alignments)
            instance.add_flag_to_inputs()
            self.merge([instance])
            self.__add_to_gfs(instance)
            stats.extend([AsmStats(rule) for rule in instance.gfs])
            self.add_edges_from([(instance.exit, stat) for stat in stats])

        final_flag = AsmFlag(stats, outdir=self.outdir)
        self.add_node(final_flag)
        self.add_edges_from([(stat, final_flag) for stat in stats])

    @property
    def gfs(self):
        return self.__gfs

    def __add_to_gfs(self, wrapper: ShortAssemblerWrapper):
        self.__gfs.extend(wrapper.gfs)

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))
