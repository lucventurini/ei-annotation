from ..abstract import AtomicOperation
import abc
import os
from ..repeats import RepeatMasking
from .chunking import ChunkProteins


class ProteinChunkAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self,
                 chunks: ChunkProteins,
                 chunk,
                 masked: RepeatMasking):
        super().__init__()
        self._chunks = chunks
        self.configuration = chunks.configuration
        self.__chunk = None
        self.chunk = chunk
        self.input["flag"] = chunks.output["flag"]
        self.input["genome"] = self.masked_genome
        self.input["repeat_flag"] = masked.exit.output["flag"]
        assert self.input["fasta"] in chunks.output["chunks"]

    @property
    def rulename(self):
        return "{toolname}_protein_{chunk}".format(toolname=self.toolname,
                                                   chunk=str(self.chunk).zfill(3))

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def chunk(self):
        return self.__chunk

    @chunk.setter
    def chunk(self, chunk):
        assert isinstance(chunk, int)
        fasta = os.path.join(self._chunks.outdir, "chunk_{}.fasta".format(str(chunk).zfill(3)))
        assert fasta in self._chunks.output["chunks"], (fasta, self._chunks.output["chunks"])
        self.input["fasta"] = fasta
        self.__chunk = chunk

    @property
    @abc.abstractmethod
    def cmd(self):
        pass

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")

    @property
    def logdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "logs")