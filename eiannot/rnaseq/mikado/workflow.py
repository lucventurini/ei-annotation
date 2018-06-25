from . import EIWrapper, AssemblyWrapper, LongAlignmentsWrapper, PortcullisWrapper
from . import MikadoConfig, MikadoPrepare, MikadoHomologyWrapper, MikadoPick
from .orfs import Prodigal, TransdecoderLongOrf, TransdecoderPred
from . import FaidxGenome, MikadoSerialise
from . import IndexMikado, MikadoStats


class Mikado(EIWrapper):
    def __init__(self,
                 assemblies: AssemblyWrapper,
                 long_alignments: LongAlignmentsWrapper,
                 portcullis: PortcullisWrapper):

        super().__init__()
        self.__indexer, self.__picker, self.__stats = None, None, None
        self.configuration = assemblies.configuration
        self.assemblies = assemblies
        self.long_alignments = long_alignments
        self.portcullis = portcullis

        self.configurer = MikadoConfig(portcullis_wrapper=self.portcullis,
                                       assemblies=self.assemblies,
                                       long_aln_wrapper=long_alignments)
        self.add_edges_from([step, self.configurer] for step in
                            [self.assemblies, self.long_alignments, self.portcullis])

        self.preparer = MikadoPrepare(self.configurer)
        self.add_edge(self.configurer, self.preparer)
        if self.orf_caller == "Prodigal":
            self.orfs = Prodigal(self.preparer)
            self.add_edge(self.preparer, self.orfs)
        elif self.orf_caller == "Transdecoder":
            longorf = TransdecoderLongOrf(self.preparer)
            self.add_edge(self.preparer, longorf)
            self.orfs = TransdecoderPred(self.preparer, longorf)
            self.add_edge(longorf, self.orfs)
        else:
            self.orfs = None
        self.homologies = MikadoHomologyWrapper(self.preparer)
        if self.homologies.execute is False:
            self.homologies = None
        else:
            self.add_edge(self.preparer, self.homologies)

        self.faidx_genome = FaidxGenome(None, self.configuration)
        self.serialiser = MikadoSerialise(prepare=self.preparer,
                                          homology=self.homologies,
                                          orfs=self.orfs,
                                          faidx=self.faidx_genome,
                                          portcullis=self.portcullis)
        self.add_edges_from([_, self.serialiser] for _ in [self.faidx_genome, self.preparer,
                                                           self.homologies, self.orfs] if _ is not None)

        self.picker = MikadoPick(self.serialiser)
        self.add_edge(self.serialiser, self.picker)
        self.indexer = IndexMikado(self.picker)
        self.add_edge(self.picker, self.indexer)
        stats = MikadoStats(self.indexer)
        self.add_edge(self.indexer, stats)
        try:
            _ = self.exit.rulename
        except ValueError:
            nodes = [_.rulename for _ in self if not self.adj[_]]
            raise ValueError([(_, [n.rulename for n in nx.ancestors(self.graph, _)]) for _ in nodes])

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado")

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