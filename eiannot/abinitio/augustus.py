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


class PrepareProteins(AtomicOperation):

    def __init__(self, proteins: ProteinWrapper, db: str):

        super.__init__()
        self.configuration=proteins.configuration
        self.__db = db

    @property
    def rulename(self):
        return "convert_to_hints_protein_{db}".format(db=self.db)

    @property
    def db(self):
        return self.__db

