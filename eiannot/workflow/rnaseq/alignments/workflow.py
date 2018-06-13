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
        bams = []
        for wrapper in self.wrappers.values():
            instance = wrapper(configuration, bams, aln_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            bams.extend(instance.bams)
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
        if self.execute_portcullis is True:
            self.portcullis = PortcullisWrapper(self.configuration, bams)
            self.merge(self.portcullis)
            flag = self.portcullis.flag
            self.add_edges_from([stat, flag] for stat in stats)

    @property
    def execute_portcullis(self):

        # TODO: probably we should push this somewhere else
        return self.configuration["portcullis"]["execute"]

    @property
    def portcullis_junctions(self):
        return self.portcullis.output["bed"]


class LongAlignmentsWrapper(EIWrapper):



    def __init__(self):
        pass
