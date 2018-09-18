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
