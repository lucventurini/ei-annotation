import abc
from eiannot.workflow import ShortSample, LongSample, AtomicOperation, Sample, EIWrapper
import os
from .bam import BamStats


class ShortAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, index, sample, run, configuration, outdir, ref_transcriptome=None):

        super(ShortAligner, self).__init__()
        self.__configuration, self.__sample, self.__run = None, None, None
        self.input = {"index": index}
        self.configuration = configuration
        self.sample = sample
        self.run = run
        self.input["index"] = index
        self.__ref_transcriptome = None
        self.ref_transcriptome = ref_transcriptome
        self._outdir = outdir
        self.input["index"] = index
        self.input["read1"] = self.sample.read1
        if self.sample.read2 is not None:
            self.input["read2"] = self.sample.read2
        self.log = os.path.join(self._outdir, "{toolname}-{sample}-{run}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run
        ))

        self.message = "Aligning input with {toolname} (sample {sample} - run {run})".format(
            toolname=self.toolname,
            sample=self.sample.label,
            run=self.run
        )

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    @property
    def bamdir(self):
        return os.path.join(self._outdir, self.toolname,
                            "{sample}-{run}".format(sample=self.sample.label, run=self.run))

    @property
    def link(self):
        return os.path.join(self._outdir, "output", "{toolname}-{sample}-{run}.bam".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run)
        )

    @property
    def sample(self):
        return self.__sample

    @sample.setter
    def sample(self, sample):
        if isinstance(sample, ShortSample):
            if sample.label not in self.configuration["short_reads"]:
                raise KeyError("Sample {sample.label} not found in the configuration!".format(**locals()))
            self.__sample = sample
        elif isinstance(sample, (str, bytes)):
            if sample not in self.configuration["short_reads"]:
                raise KeyError("Sample {sample} not found in the configuration!".format(**locals()))
            self.__sample = self.configuration["short_reads"][sample]

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
    def link_src(self):
        return os.path.join("..", self.toolname,
                            "{sample}-{run}".format(sample=self.sample.label, run=self.run),
                            os.path.basename(self.output["bam"]))

    @property
    def rulename(self):
        # We have to give a unique name to each alignment
        return "{toolname}-{sample}-{run}".format(sample=self.sample.label,
                                                  run=self.run, toolname=self.toolname)

    @property
    def sample(self):
        return self.__sample


class LongAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, index, sample, run, configuration, ref_transcriptome=None):

        super().__init__()
        self.__configuration, self.__sample, self.__run = None, None, None
        self.input = {"index": index}
        self.configuration = configuration
        self.__sample = None
        self.sample = sample
        self.run = run
        self.input["index"] = index
        self.__ref_transcriptome = None
        self.ref_transcriptome = ref_transcriptome
        self.input["index"] = index
        self.input["read1"] = self.sample.read1
        if self.sample.read2 is not None:
            self.input["read2"] = self.sample.read2
        self.log = os.path.join(self._outdir, "{toolname}-{sample}-{run}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run
        ))

        self.message = "Aligning input with {toolname} (sample {sample} - run {run})".format(
            toolname=self.toolname,
            sample=self.sample.label,
            run=self.run
        )

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    @property
    @abc.abstractmethod
    def suffix(self):
        pass

    @property
    def sample(self):
        return self.__sample

    @sample.setter
    def sample(self, sample):
        if isinstance(sample, LongSample):
            if sample.label not in self.configuration["long_reads"]:
                raise KeyError("Sample {sample.label} not found in the configuration!".format(**locals()))
            self.__sample = sample
        elif isinstance(sample, (str, bytes)):
            if sample not in self.configuration["long_reads"]:
                raise KeyError("Sample {sample} not found in the configuration!".format(**locals()))
            self.__sample = self.configuration["long_reads"][sample]


class IndexBuilder(AtomicOperation, metaclass=abc.ABCMeta):

    """Abstract class for building the index of a genome, given a tool"""

    def __init__(self, configuration, outdir):
        super(IndexBuilder, self).__init__()
        self.input = {"genome": self.genome}
        if configuration.get("reference", dict()).get("transcriptome", ""):
            self.input["ref_transcriptome"] = os.path.abspath(configuration["reference"]["transcriptome"])
        self._outdir = os.path.join(outdir, "rnaseq", "2-alignments", "index", self.toolname)
        self.log = os.path.join(outdir, "rnaseq", "2-alignments", "index", "log", "{}.log".format(self.toolname))
        self.__configuration = configuration


    @property
    def extra(self):
        return self.__configuration["programs"][self.toolname]["index"]

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    @property
    def rulename(self):
        return "{toolname}_index".format(toolname=self.toolname)


class ShortWrapper(EIWrapper, metaclass=abc.ABCMeta):

    def __init__(self, configuration, prepare_flag):
        self.__prepare_flag = prepare_flag
        super().__init__()
        self.__bam_rules = set()
        self.configuration = configuration

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    def add_to_bams(self, rule):
        if not isinstance(rule, ShortAligner):
            raise TypeError
        if "link" not in rule.output:
            raise KeyError("Link not found for rule {}".format(rule.rulename))
        self.__bam_rules.add(rule.rulename)

    @property
    def bams(self):
        """This property will the output files (ie BAMs) of the wrapper."""
        return self.__bam_rules

    @property
    def runs(self):
        return self.configuration["programs"][self.toolname]["runs"]

    def add_flag_to_inputs(self):
        for rule in self:
            rule.input["aln_flag"] = self.__prepare_flag

    @property
    @abc.abstractmethod
    def indexer(self):
        assert isinstance(self.__indexer, IndexBuilder)
        return self.__indexer(self.configuration, self.outdir)

    @property
    def samples(self):
        return self.configuration["short_reads"]

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["out_dir"], "rnaseq", "1-alignments"))
