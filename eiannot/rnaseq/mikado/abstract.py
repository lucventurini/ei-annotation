from ...abstract import AtomicOperation
import os
import abc

__modes__ = ("permissive", "stringent", "nosplit", "split", "lenient")


class MikadoOp(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, is_long=False):

        super().__init__()
        self.__is_long = is_long

    @property
    def is_long(self):
        assert isinstance(self.__is_long, bool)
        return self.__is_long

    @property
    def mikado_dir(self):
        if self.is_long is False:
            return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado")
        else:
            return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado-long-reads")

    @property
    @abc.abstractmethod
    def _rulename(self):
        pass

    @property
    def rulename(self):
        if self.is_long:
            long = "_long"
        else:
            long = ""
        rule = self._rulename

        return "{rule}{long}".format(**locals())
