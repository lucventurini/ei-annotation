from .abstract import AsmStats, AsmFlag
from .class2 import Class2Wrapper
from .stringtie import StringtieWrapper
from .cufflinks import CufflinksWrapper
from .trinity_gg import TrinityGGWrapper
from .scallop import ScallopWrapper
from ... import EIWrapper, EIWorfkflow


class AssemblyWrapper(EIWrapper):

    wrappers = {"cufflinks": CufflinksWrapper,
                "scallop": ScallopWrapper,
                "trinity": TrinityGGWrapper,
                "class2": Class2Wrapper,
                "stringtie": StringtieWrapper}

    def __init__(self, configuration, bams, aln_flag):

        super().__init__()
        stats = []
        for wrapper in self.wrappers.values():
            instance = wrapper(configuration, bams, aln_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            stats.extend([AsmStats(rule) for rule in instance.gfs])
            self.add_edges_from([(instance.exit, stat) for stat in stats])

        final_flag = AsmFlag(stats)
        self.add_edges_from([(stat, final_flag) for stat in stats])
