from ...abstract import EIWrapper
from ..alignments.__init__ import LongAlignmentsWrapper
from ..alignments.portcullis import PortcullisWrapper
from ..assemblies import AssemblyWrapper
from .prepare import MikadoConfig, MikadoPrepare
from .serialise import MikadoSerialise
from .homology import MikadoHomologyWrapper
from .pick import MikadoPick, IndexMikado, MikadoStats
from .orfs import Prodigal, TransDecoderWrapper, GTCDS  # TransdecoderLongOrf, TransdecoderPred, TransdecoderOrfBlast
import os
import networkx as nx


# TODO: we have to execute Mikado probably TWICE, once on only long reads, once on long and short reads.


class Mikado(EIWrapper):

    __final_rulename__ = "mikado_done"
    __final_long_rulename__ = "mikado_done_long"

    def __init__(self,
                 assemblies: AssemblyWrapper,
                 long_alignments: LongAlignmentsWrapper,
                 portcullis: PortcullisWrapper,
                 only_long=False):

        super().__init__()
        self.__indexer, self.__picker, self.__stats = None, None, None
        self.assemblies = assemblies
        self.long_alignments = long_alignments
        self.short_alignments = assemblies.short_alignments
        self.portcullis = portcullis

        self.configuration = assemblies.configuration
        if only_long is True:
            execute = (len(self.long_alignments.gfs) > 0)
        else:
            execute = (len(self.assemblies.gfs) > 0 or len(self.long_alignments.gfs) > 0)

        self.configurer = MikadoConfig(portcullis_wrapper=self.portcullis,
                                       assemblies=self.assemblies,
                                       long_aln_wrapper=long_alignments,
                                       is_long=only_long)
        if execute:

            self.add_edges_from([step, self.configurer] for step in
                                [self.assemblies, self.long_alignments, self.portcullis])

            self.preparer = MikadoPrepare(self.configurer)
            self.add_edge(self.configurer, self.preparer)
            if self.orf_caller == "Prodigal":
                self.orfs = Prodigal(self.preparer)
                self.add_edge(self.preparer, self.orfs)
            elif self.orf_caller == "GTCDS":
                self.orfs = GTCDS(self.preparer)
                self.add_edge(self.preparer, self.orfs)
            elif self.orf_caller == "Transdecoder":
                self.trans_wrapper = TransDecoderWrapper(self.preparer)
                self.add_edge(self.preparer, self.trans_wrapper)
                self.orfs = self.trans_wrapper.orfs
            else:
                self.orfs = None
            self.homologies = MikadoHomologyWrapper(self.preparer)
            if self.homologies.execute is False:
                self.homologies = None
            else:
                self.add_edge(self.preparer, self.homologies)
            self.serialiser = MikadoSerialise(prepare=self.preparer,
                                              homology=self.homologies,
                                              orfs=self.orfs,
                                              portcullis=self.portcullis)
            self.add_edges_from([_, self.serialiser] for _ in [self.preparer,
                                                               self.homologies, self.orfs] if _ is not None)
            self.picker = MikadoPick(self.serialiser)
            self.add_edge(self.serialiser, self.picker)
            self.indexer = IndexMikado(self.picker)
            self.add_edge(self.picker, self.indexer)
            self.__stats = MikadoStats(self.indexer)
            self.add_edge(self.indexer, self.__stats)

        if only_long:
            self.__final_rulename__ = self.__final_long_rulename__

        self.add_final_flag()

        try:
            _ = self.exit.rulename
        except ValueError:
            nodes = [_.rulename for _ in self if not self.adj[_]]
            raise ValueError([(_, [n.rulename for n in nx.ancestors(self.graph, _)]) for _ in nodes])

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "mikado.done")

    @property
    def outdir(self):
        return self.configurer.mikado_dir

    @property
    def gfs(self):
        return self.assemblies.gfs

    @property
    def junctions(self):
        return self.portcullis.junctions

    @property
    def orf_caller(self):
        if self.configuration.get("orfs", dict()).get("execute", True):
            if self.configuration.get("mikado", dict()).get("use_prodigal", True):
                return "Prodigal"
            elif self.configuration.get("mikado", dict()).get("use_gt", False):
                return "GTCDS"
            else:
                return "Transdecoder"
        else:
            return None

    @property
    def indexer(self):
        return self.__indexer

    @indexer.setter
    def indexer(self, indexer):

        if not isinstance(indexer, IndexMikado):
            raise TypeError(type(indexer))
        self.__indexer = indexer

    @property
    def picker(self):
        return self.__picker

    @picker.setter
    def picker(self, picker):

        if not isinstance(picker, MikadoPick):
            raise TypeError(type(picker))
        self.__picker = picker

    @property
    def stats(self):
        return self.__stats

    @stats.setter
    def stats(self, stats):

        if not isinstance(stats, MikadoPick):
            raise TypeError(type(stats))
        self.__stats = stats
