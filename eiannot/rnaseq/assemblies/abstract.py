import abc
from ...abstract import ShortSample, LongSample, AtomicOperation, Sample, EIWrapper
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
        if bam is not None:
            self.input["bam"] = bam
        self.log = os.path.join(self._outdir, "{toolname}-{sample}-{run}-{alrun}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run,
            alrun=self.alrun
        ))
        self.output["link"] = os.path.join(
          self._outdir, "output", "{toolname}-{sample}-{run}-{alrun}.{suffix}".format(
                toolname=self.toolname, sample=self.sample, run=self.run,
                suffix=self.suffix, alrun=self.alrun
            )
        )
        self.output["gf"] = os.path.join(self.gfdir, os.path.basename(self.output["link"]))

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
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run,
                                                            alrun=self.alrun))

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
    @abc.abstractmethod
    def suffix(self):
        pass

    @property
    def link_src(self):
        return os.path.join("..", self.toolname,
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run, alrun=self.alrun),
                            os.path.basename(self.output["bam"]))

    @property
    def rulename(self):
        # We have to give a unique name to each alignment
        return "{toolname}-{run}-{alrun}".format(alrun=self.alrun,
                                                 run=self.run, toolname=self.toolname)


class ShortAssemblerWrapper(EIWrapper, metaclass=abc.ABCMeta):

    def __init__(self, aln_wrapper):
        super().__init__(configuration=aln_wrapper.configuration)
        self.__gf_rules = set()
        self.aln_flag = aln_wrapper.output
        self.bams = aln_wrapper.bams

    @property
    @abc.abstractmethod
    def toolname(self):
        pass

    def add_to_gf(self, rule):
        if not isinstance(rule, ShortAssembler):
            raise TypeError
        if "link" not in rule.output:
            raise KeyError("Link not found for rule {}".format(rule.rulename))
        self.__gf_rules.add(rule.rulename)

    @property
    def gfs(self):
        """This property will the output files (ie GTF/GFFs) of the wrapper."""
        return self.__gf_rules

    @property
    def runs(self):
        return self.configuration["programs"][self.toolname]["runs"]

    def add_flag_to_inputs(self):
        for rule in self:
            rule.input["aln_flag"] = self.aln_flag

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["out_dir"], "rnaseq", "2-assemblies"))


class AsmStats(AtomicOperation):

    def __init__(self, asm_run: ShortAssembler):

        super().__init__()
        self.input["gf"] = asm_run.output["link"]
        self.output["stats"] = self.input["gf"] + ".stats"
        self.message = "Computing assembly stats for: {input[gf]}"
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return "asm_gf_stats_{}".format(os.path.basename(os.path.splitext(self.input["gf"])[0]))

    @property
    def cmd(self):
        load, output, input = self.load, self.output, self.input
        cmd = "{load} mikado util stats {input[gf]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["mikado"]

    @property
    def threads(self):
        return 1


class AsmFlag(AtomicOperation):

    def __init__(self, stats_runs: [AsmStats]):

        super().__init__()
        self.touch = True
        outdir = os.path.dirname(os.path.dirname(stats_runs[0].output["stats"]))
        self.output["flag"] = os.path.join(outdir, "all.done")

    @property
    def rulename(self):
        return "asm_all"

    @property
    def loader(self):
        return []
