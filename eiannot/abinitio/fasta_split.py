from ..abstract import AtomicOperation
from ..preparation import FaidxGenome
from .abstract import augustus_root_dir
import os


class FastaMSplitter(AtomicOperation):

    """This step of the pipeline has the role to create a database of chunks for running Augustus."""

    def __init__(self, faidx: FaidxGenome):

        super().__init__()
        self.input["genome"] = self.genome
        self.input["fai"] = faidx.output["fai"]
        self.output["chunk_db"] = os.path.join(self._root_dir,
                                               augustus_root_dir,
                                               "chunks.db")

    @property
    def num_chunks(self):
        # TODO put this into the configuration
        return self.configuration.get("abinitio", dict()).get("chunks", 1000)

    @property
    def chunks(self):
        return range(self.num_chunks)

    @property
    def minsize(self):
        # TODO put this into the configuration
        return self.configuration.get("abinitio", dict()).get("minsize", 5 * 10**6)

    @property
    def minoverlap(self):
        # TODO put this into the configuration
        return self.configuration.get("abinitio", dict()).get("minsize", 5 * 10**5)

    @property
    def loader(self):
        return ["ei-annotation"]

    @property
    def rulename(self):
        return "genome_chunker_augustus"

    @property
    def cmd(self):
        """Command to split out the genome into manageable chunks. *THE* problem
        that we have to solve is that we cannot predict beforehand the number of chunks,
        and that the splitting script from Augustus actually tries to keep chromosomes together.
        The command we end up using should - of course - try to have overlapping sequences as much as humanly possible.
        """

        load = self.load
        chunks = self.num_chunks
        cmd = "{load} split_genome_fasta.py -n {chunks} "
        minsize, minoverlap = self.minsize, self.minoverlap
        cmd += "-ms {minsize} -minO {minoverlap} {input[genome]} {output[chunk_db]}"
        cmd = cmd.format(**locals())
        return cmd
