from .sample import Sample
from .alignment_run import AlignmentRuns
import os


class AsmRuns:

    __asm_abrv = {"cufflinks": ("cuf", "gtf"),
                  "stringtie": ("stn", "gtf"),
                  "class": ("cls", "gtf"),
                  "trinity": ("trn", "gff"),
                  "scallop": ("scl", "gtf")
                  }

    __gtf_methods = set()

    def __init__(self, methods, asm_dir, aln_runs: AlignmentRuns):

        self.__methods = methods
        self.__all_runs = dict()
        self.__asm_runs = dict.fromkeys(methods)
        self.__outdir = asm_dir

        # This is a dictionary which has the name of the method to be used, and the settings
        for method in self.__methods:
            abrv, suffix = self.__asm_abrv.get(method, (method, "gtf"))
            self.__asm_runs[method] = dict()
            for aln in aln_runs:  # This will iterate over the names
                aln_string = aln_runs[aln]["aln_string"]
                aln_abrv_string = aln_runs[aln]["abrv_string"]
                sample = aln_runs[aln]["sample"]
                for index, setting in enumerate(self.methods[method]):
                    key = "{asm}-{index}-{aln_string}.{suffix}".format(asm=method,
                                                                       index=index,
                                                                       aln_string=aln_string,
                                                                       suffix=suffix)
                    key = os.path.join(self.outdir, key)
                    self.__asm_runs[method][key] = dict()
                    self.__asm_runs[method][key]["setting"] = setting
                    self.__asm_runs[method][key]["sample"] = sample
                    self.__asm_runs[method][key]["abrv"] = "{abrv}-{index}-{aln_abrv_string}".format(
                        abrv=abrv, index=index, aln_abrv_string=aln_abrv_string)
                    self.__all_runs[key] = (method, key)

    @property
    def methods(self):
        return self.__methods.keys()

    @property
    def output_array(self):
        """This returns a comma separated list of tags"""
        return ",".join(self.keys())

    @property
    def label_array(self):
        """This returns a comma separated list of tags for Mikado"""
        arr = []
        for key in self:
            arr.append(self[key]["abrv"])
        return ",".join(arr)

    @property
    def outdir(self):
        return self.__outdir

    def __iter__(self):
        return iter(self.__all_runs.keys())

    def __getitem__(self, key):
        if key in self.__all_runs:
            return self.__asm_runs[self.__all_runs[key][0]][key]
        elif key in self.__asm_runs:
            return self.__all_runs[key]
        else:
            raise KeyError

    def keys(self):
        return self.__all_runs.keys()