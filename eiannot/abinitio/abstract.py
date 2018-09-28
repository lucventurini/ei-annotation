from ..abstract import AtomicOperation, EIWrapper
import os
import abc


augustus_root_dir = outdir = os.path.join("abinitio")


class AugustusMethod(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, configuration):

        super().__init__()
        self.configuration = configuration

    @property
    def _augustus_dir(self):
        return os.path.join(self._root_dir, augustus_root_dir)

    @property
    @abc.abstractmethod
    def __subfolder(self):
        pass
