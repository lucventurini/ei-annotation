from .preparation.prepare import parse_samplesheet
from .rnaseq import *


__name__ = "eiannot"
__version__ = "0.1"


def load_pre_cmd(*args):
    """
    Used to prefix a shell command that utilises some external software with another command used to load that software
    """
    cc = "set +u && "
    for arg in args:
        cc += "{} && ".format(arg)
    return cc
