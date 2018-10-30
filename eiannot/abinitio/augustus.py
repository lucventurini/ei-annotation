from ..abstract import AtomicOperation, EIWrapper
from .train_augustus import TrainAugustusWrapper
from ..rnaseq.mikado import Mikado
from ..repeats import RepeatMasking
from ..preparation import FaidxGenome
from .converters import ConvertToHints
from ..proteins import ProteinWrapper
import pkg_resources
import os
import sys


class AugustusWrapper(EIWrapper):

    def __init__(self,
                 mikado: Mikado,
                 mikado_long: Mikado,
                 rmasker: RepeatMasking,
                 proteins: ProteinWrapper,
                 run: int):

        super().__init__(configuration=mikado.configuration)

        converter = ConvertToHints(mikado=mikado,
                                   portcullis=mikado.portcullis,
                                   alignments=mikado.long_alignments,
                                   repeats=rmasker,
                                   mikado_long=mikado_long,
                                   run=run,
                                   proteins=proteins)

        pass


# We have to split out the FASTA first


class FastaMSplitter(AtomicOperation):

    def __init__(self, faidx: FaidxGenome):

        super().__init__()
        self.input["genome"] = self.genome
        self.input["fai"] = faidx.output["fai"]
        self.output["chunk_db"] = pass

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
        cmd = "{load} split_genome_fasta.py -n {chunks} "
        minsize, minoverlap = self.minsize, self.minoverlap
        cmd += "-ms {minsize} -minO {minoverlap} {input[genome]} {output[chunk_db]}"
        cmd = cmd.format(**locals())
        return cmd


class RunAugustus(AtomicOperation):

    def __init__(self,
                 fastam: FastaMSplitter,
                 trained: TrainAugustusWrapper,
                 chunk: int,
                 aug_run: int):

        super().__init__()
        self.input.update(fastam.output)


        self.__chunk = chunk
        self.__run = aug_run

        pass

    @property
    def run(self):
        return self.__run

    @property
    def chunk(self):
        return self.__chunk

    @property
    def rulename(self):
        return "run_augustus_run{}_{}".format(self.run, self.chunk)

    @property
    def species(self):

        if self.configuration.get("training", {}).get("train", True):
            return self.input["trained"]
        else:
            model = self.configuration.get("training", {}).get("model_species", None)
            assert model is not None
            return model

    @property
    def loader(self):
        return ["ei-annotation", "augustus", "mikado"]

    @property
    def extrinsic(self):
        if self.configuration["abinitio"]["extrinsic"] is None:
            return pkg_resources.resource_filename("eiannot.configuration", "extrinsic.ei_augustus_generic.cfg")
        else:
            return self.configuration["abinitio"]["extrinsic"]

    def cmd(self):

        cmd = "{load}"
        load = self.load
        species = self.species
        extrinsic = self.extrinsic
        input, output = self.input
        log = self.log
        chunk = self.chunk
        if input["hints"] is not None:
            hints = "--hintsfile={input[hints]}"
        else:
            hints = ""
        cmd += "execute_augustus.py {input[genome]} {input[chunk_db]} {chunk} {output[gtf]} "
        # Execute augustus will take as last argument the augustus command string
        cmd += "\"augustus --species {species}  --UTR=on --extrinsicCfgFile={extrinsic} "
        cmd += "--stopCodonExcludedFromCDS=true --genemodel=partial "
        cmd += " --alternatives-from-evidence=true {hints} --noInFrameStop=true "
        cmd += " --allow_hinted_splicesites=atac --errfile={log} \""
        cmd = cmd.format(**locals())
        return cmd
