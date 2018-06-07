import abc
from eiannot.workflow import ShortSample, LongSample, AtomicOperation, Sample
import os
import re


class ShortAssembler(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, bam, run, configuration, outdir, ref_transcriptome=None):

        super(ShortAssembler, self).__init__()
        self.__configuration, self.__sample, self.__run = None, None, None
        self.configuration = configuration
        self.sample = self.get_sample()
        self.run = run
        self.__ref_transcriptome = None
        self.ref_transcriptome = ref_transcriptome
        self._outdir = outdir
        self.input["bam"] = bam
        self.log = os.path.join(self._outdir, "{toolname}-{sample}-{run}-{alrun}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run,
            alrun=self.alrun
        ))
        self.output["gf"] = os.path.join(
          self._outdir, "output", "{toolname}-{sample}-{run}-{alrun}.{suffix}".format(
                toolname=self.toolname, sample=self.sample, run=self.run,
                suffix=self.suffix, alrun=self.alrun
            )
        )
        self.message = "Using {toolname} to assemble (run: {run}): {input[bam]}".format(
            input=self.input, run=run, toolname=self.toolname)

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    def get_sample(self):
        bam = os.path.basename(self.input["bam"])
        groups = re.search("(.*)-([^-]*)-([^-]*)", os.path.splitext(bam)[0]).groups()
        if groups is None:
            raise ValueError(bam)
        return "{groups[1]}".format(**locals())

    @property
    def alrun(self):
        """This method will infer the alignment run from the name of the BAM file."""
        bam = os.path.basename(self.input["bam"])
        groups = re.search("(.*)-([^-]*)-([^-]*)", os.path.splitext(bam)[0]).groups()
        if groups is None:
            raise ValueError(bam)
        return "{groups[0]}-{groups[2]}".format(**locals())

    @property
    def gfdir(self):
        return os.path.join(self._outdir, self.toolname,
                            "{sample}-{run}".format(sample=self.sample.label, run=self.run))

    @property
    def link(self):
        return os.path.join(self._outdir, "output", "{toolname}-{sample}-{run}-".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run)
        )

    @property
    def sample(self):
        return self.__sample

    @sample.setter
    def sample(self, sample):
        if not isinstance(sample, ShortSample):
            raise TypeError
        self.__sample = sample

    @property
    def run(self):
        return self.__run

    @run.setter
    def run(self, run):
        if not isinstance(run, int):
            raise TypeError
        runs = self.__configuration["programs"][self.toolname]["runs"]
        if run not in range(len(runs)):
            raise ValueError
        self.__run = run

    @property
    def extra(self):
        return self.__configuration["programs"][self.toolname]["runs"][self.run]

    @property
    def ref_transcriptome(self):
        return self.__ref_transcriptome

    @ref_transcriptome.setter
    def ref_transcriptome(self, ref):
        if not os.path.exists(os.path.abspath(ref)):
            raise ValueError
        self.__ref_transcriptome = ref
        self.input["ref_transcriptome"] = ref

    @property
    def index(self):
        return self.input["index"]

    @property
    @abc.abstractmethod
    def strand(self):
        pass

    @property
    @abc.abstractmethod
    def input_reads(self):
        pass

    @property
    def min_intron(self):
        return max(self.configuration["reference"]["min_intron"], 20)

    @property
    def max_intron(self):
        return self.configuration["reference"]["max_intron"]

    @property
    @abc.abstractmethod
    def suffix(self):
        pass

    @property
    def link_src(self):
        return os.path.join("..", self.toolname,
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run),
                            os.path.basename(self.output["bam"]))

    @property
    def rulename(self):
        # We have to give a unique name to each alignment
        return "{toolname}-{sample}-{run}".format(sample=self.sample.label,
                                                  run=self.run, toolname=self.toolname)