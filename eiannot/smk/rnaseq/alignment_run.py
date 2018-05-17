from .sample import Sample
from frozendict import frozendict
import os


class AlignmentRuns:

    __aln_abrv = {"tophat": ("tph", "bam"),
                  "star": ("sta", "bam"),
                  "gsnap": ("gsp", "bam"),
                  "hisat": ("hst", "bam"),
                  "gmap": ("gmp", "gff")}

    def __init__(self, align_methods, outdir, samples: dict[Sample]):

        self.__methods = align_methods
        self.__outdir = outdir
        self.__align_runs = dict.fromkeys(self.__methods)
        self.__all_runs = dict()

        # The idea here is to create a data structure which will link a filename to the sample
        # it needs to be created.

        for aln in self.__methods:
            self.__align_runs[aln] = dict()
            for sample in samples:  # Here samples will be a dictionary of Sample objects
                for index, setting in enumerate(self.__methods[aln]):
                    # "Key" here is the output file
                    key = "{aln}-{sample}-{index}.{suff}".format(aln=aln,
                                                                 suff=self.__aln_abrv.get(aln, ["", ".bam"])[1],
                                                                 sample=sample, index=index)
                    key = os.path.join(self.outdir, key)
                    self.__align_runs[aln][key] = dict()
                    self.__align_runs[aln][key]["setting"] = setting
                    self.__align_runs[aln][key]["sample"] = sample
                    self.__align_runs[aln][key]["aln_string"] = "{aln}-{sample}-{index}".format(
                        aln=aln,
                        sample=sample, index=index)
                    self.__align_runs[aln][key]["abrv_string"] = "{abrv}-{sample}-{index}".format(
                        abrv=self.__aln_abrv.get(aln, [aln])[0],
                        sample=sample, index=index)
                    self.__all_runs[key] = (aln, key)

        # Freeze the dictionary
        self.__align_runs = frozendict(self.__align_runs)

    @property
    def methods(self):

        return self.__methods.keys()

    def __iter__(self):
        return iter(self.__all_runs.keys())

    def __getitem__(self, key):
        if key in self.__all_runs:
            return self.__align_runs[self.__all_runs[key][0]][key]
        elif key in self.__align_runs:
            return self.__all_runs[key]
        else:
            raise KeyError

    @property
    def outdir(self):
        return self.__outdir
