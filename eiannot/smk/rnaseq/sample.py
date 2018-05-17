import os
import itertools


class Sample:

    """This class defines a sample, together with its reads."""

    # TODO: fill this
    __orientations = ["f", "r", "ff", "fr", "rf"]

    def __init__(self,
                 name,
                 read1,
                 read2,
                 strandness=None,
                 is_long=False,
                 orientation="fr"):

        self.__name = name

        if not read1 or not os.path.exists(read1):
            raise OSError("At least the read 1 file must be present on the filesystem")
        self.__read1 = read1
        if (read2) and not os.path.exists(read2):
            raise OSError("Read 2 file not found")
        self.__read2 = read2
        self.__is_long = False
        if not isinstance(is_long, bool):
            raise TypeError("The is_long parameter must be boolean")
        self.__is_long = is_long
        if isinstance(orientation, bytes):
            orientation = orientation.decode()
        elif not isinstance(orientation, str):
            raise TypeError("Orientation must be a string")
        orientation = orientation.lower()
        if orientation not in self.__orientations:
            raise ValueError("Invalid orientation: {}".format(orientation))
        if strandness not in ["firststrand", "secondstrand", "unstranded", None]:
            raise ValueError("Invalid strandness: {}")
        self.__strandness = strandness

    @property
    def name(self):
        return self.__name

    @property
    def read1(self):
        return self.__read1

    @property
    def read2(self):
        return self.__read2

    @property
    def paired(self):
        return (self.read2)

    @property
    def is_long(self):
        return (self.__is_long)

    @property
    def strandness(self):
        return self.__strandness

    @property
    def orientation(self):
        return self.orientation
