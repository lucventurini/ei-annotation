import os


class Sample:

    """This class defines a sample, together with its reads."""

    # TODO: fill this
    __strandedness = {}

    def __init__(self, read1, read2, strand_specific=False, is_long=False):

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
        if strand_specific not in self.__strandedness:
            raise KeyError("Invalid strandedness specified")
        self.__strand_specific = strand_specific

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
