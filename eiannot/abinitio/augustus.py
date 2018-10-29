from ..abstract import AtomicOperation, EIWrapper
from .train_augustus import TrainAugustusWrapper
from ..rnaseq.mikado import Mikado
from ..repeats import RepeatMasking
from ..proteins import ProteinWrapper
import os
import sys


class AugustusWrapper(EIWrapper):

    def __init__(self,
                 mikado: Mikado,
                 mikado_long: Mikado,
                 rmasker: RepeatMasking):

        super().__init__(configuration=mikado.configuration)

        pass


# We have to split out the FASTA first


class FastaMSplitter(AtomicOperation):

    def __init__(self):

        super().__init__()
        self.input["genome"] = self.genome

    @property
    def loader(self):
        return ["augustus"]

    @property
    def cmd(self):
        """Command to split out the genome into manageable chunks. *THE* problem
        that we have to solve is that we cannot predict beforehand the number of chunks,
        and that the splitting script from Augustus actually tries to keep chromosomes together.
        The command we end up using should - of course - try to have overlapping sequences as much as humanly possible.
        """
        pass


class Augustus(AtomicOperation):

    def __init__(self):

        pass

    @property
    def species(self):

        if self.configuration.get("training", {}).get("train", True):
            return self.input["trained"]
        else:
            model = self.configuration.get("training", {}).get("model_species", None)
            assert model is not None
            return model

    def cmd(self):

        cmd = "{load}"
        load = self.load
        species = self.species
        extrinsic = self.extrinsic
        input, output = self.input
        log = self.log
        cmd += "augustus --species {species}  --UTR=on --extrinsicCfgFile={extrinsic} "
        cmd += "--stopCodonExcludedFromCDS=true --genemodel=partial "
        cmd += " --alternatives-from-evidence=true --hintsfile={input[hints]} --noInFrameStop=true "
        cmd += " --allow_hinted_splicesites=atac --outfile={output[gtf]} --errfile={log} "
        cmd += "{input[genome]}"
        return cmd
