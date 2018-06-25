from ..abstract import AtomicOperation, EIWrapper
from ..rnaseq.mikado.workflow import Mikado
from ..rnaseq.alignments.portcullis import PortcullisWrapper
from ..repeats.workflow import RepeatMasking
from ..proteins.workflow import ExonerateProteinWrapper


class ConvertToHints(EIWrapper):

    def __init__(self):
        pass


class ConvertMikado(AtomicOperation):

    pass

class ConvertRepeats(AtomicOperation):

    pass

class ConvertProteins(AtomicOperation):
    pass

class ConvertJunctions(AtomicOperation):

    pass


class GetCoverage(AtomicOperation):
    pass