from .bam import BamStats, BamIndex, BamSort, AlnFlag
from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper
from .tophat2 import TopHat2Wrapper
from .portcullis import PortcullisWrapper
from .star import StarWrapper
from ... import EIWrapper, EIWorfkflow


class ShortAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarWrapper,
                "hisat": HisatWrapper,
                "tophat2": TopHat2Wrapper,
                "gsnap": GsnapWrapper}

    def __init__(self, configuration, bams, aln_flag):

        super().__init__()
        stats = []
        for wrapper in self.wrappers.values():
            instance = wrapper(configuration, bams, aln_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            for bam in instance.bams:
                # TODO: here add the entrance from the "preparation" stage
                # self.add_edges_from()
                sorter = BamSort(bam)
                self.add_edge(bam, sorter)
                indexer = BamIndex(sorter)
                self.add_edge(sorter, indexer)
                stater = BamStats(indexer)
                self.add_edge(indexer, stater)
                stats.append(stater)

        final_flag = AlnFlag(stats)
        self.add_edges_from([(stat, final_flag) for stat in stats])
