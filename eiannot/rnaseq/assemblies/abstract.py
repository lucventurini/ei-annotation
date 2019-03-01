import abc
from ...abstract import ShortSample, LongSample, AtomicOperation, Sample, EIWrapper
from ..alignments.bam import BamStats
import os
import re
import networkx as nx


class ShortAssembler(AtomicOperation, metaclass=abc.ABCMeta):

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError("Wrapper {} does not have a defined toolname!".format(cls.__name__))
        cls.__final_rulename__ = "{toolname}_flag".format(toolname=cls.__toolname__)
        super().__init_subclass__()

    def __init__(self, bam: BamStats, run, create_link=True, ref_transcriptome=None):

        super(ShortAssembler, self).__init__()
        self.__sample, self.__run = None, None
        assert create_link in (True, False)  # This must be boolean
        self._create_link = create_link
        self.configuration = bam.configuration
        self.sample = bam.sample
        self.run = run
        self.alrun = bam.align_run
        self.__ref_transcriptome = None
        self.ref_transcriptome = ref_transcriptome
        if bam is not None:
            self.input["bam"] = bam.input["bam"]
        self.log = os.path.join(self.outdir, "{toolname}-{sample}-{run}-{alrun}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run,
            alrun=self.alrun
        ))
        if self._create_link is True:
            self.output["link"] = self.link
        self.output["gf"] = os.path.join(self.gfdir, os.path.basename(self.link))

        self.message = "Using {toolname} to assemble (run: {run}): {input[bam]}".format(
            input=self.input, run=run, toolname=self.toolname)

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def gfdir(self):
        return os.path.join(self.outdir, self.toolname,
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run,
                                                            alrun=self.alrun))

    @property
    def link(self):
        return os.path.join(self.outdir, "output", "{toolname}-{sample}-{run}-{alrun}.{suffix}".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run,
            suffix=self.suffix, alrun=self.alrun
        ))

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
        runs = self.configuration["programs"][self.toolname]["runs"]
        if run not in range(len(runs)):
            raise ValueError
        self.__run = run

    @property
    def extra(self):
        return self.configuration["programs"][self.toolname]["runs"][self.run]

    @property
    def ref_transcriptome(self):
        return self.__ref_transcriptome

    @ref_transcriptome.setter
    def ref_transcriptome(self, ref):
        if ref is None:
            return
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
    @abc.abstractmethod
    def suffix(self):
        pass

    @property
    def link_src(self):
        return os.path.join("..", self.toolname,
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run, alrun=self.alrun),
                            os.path.basename(self.output["gf"]))

    @property
    def rulename(self):
        # We have to give a unique name to each alignment
        return "{toolname}-{run}-{alrun}".format(alrun=self.alrun,
                                                 run=self.run, toolname=self.toolname)

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))

    @property
    def label(self):
        alrun = re.sub("--", "-", re.sub(self.sample.label, '', self.alrun))
        alrun = re.sub(".sorted", "", alrun)
        return "{sample}-{toolname}-{run}-{alrun}".format(
            sample=self.sample.label,
            run=self.run,
            toolname=self.toolname,
            alrun=alrun)


class ShortAssemblerWrapper(EIWrapper, metaclass=abc.ABCMeta):

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError("Wrapper {} does not have a defined toolname!".format(cls.__name__))
        if not hasattr(cls, "__tag__"):
            raise NotImplementedError("Wrapper {} does not have a defined abundance tag!".format(cls.__name__))

        cls.__final_rulename__ = "{toolname}_flag".format(toolname=cls.__toolname__)
        super().__init_subclass__()

    def __init__(self, aln_wrapper):
        super().__init__(configuration=aln_wrapper.configuration)
        self.__gf_rules = set()
        self.aln_flag = aln_wrapper.exit
        self.__bams = []
        self.bams = aln_wrapper.bams
        self.configuration = aln_wrapper.configuration

    @property
    def toolname(self):
        return self.__toolname__

    def add_to_gf(self, rule):
        if not isinstance(rule, AsmStats):
            raise TypeError
        if "gf" not in rule.input:
            raise KeyError("Link not found for rule {}".format(rule.rulename))
        self.__gf_rules.add(rule)

    @property
    def gfs(self):
        """This property will the output files (ie GTF/GFFs) of the wrapper."""
        return self.__gf_rules

    @property
    def runs(self):
        return self.configuration["programs"].get(self.toolname, dict()).get("runs", [])

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))

    @property
    def bams(self):
        return self.__bams

    @bams.setter
    def bams(self, bams: BamStats):
        assert isinstance(bams, list) and all(isinstance(bam, AtomicOperation) for bam in bams)
        assert all("bam" in bam.input for bam in bams), [bam.output for bam in bams]
        self.__bams = bams

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "{toolname}.done".format(toolname=self.toolname))

    @property
    def mono_threshold(self):

        return self.configuration["programs"].get(
            self.__toolname__, {}).get("monoexonic_abundance_threshold", 1)

    @property
    def multi_threshold(self):

        return self.configuration["programs"].get(
            self.__toolname__, {}).get("multiexonic_abundance_threshold", 0)

    @property
    def tag(self):
        if "abundance_tag" in self.configuration["programs"].get(
                self.__toolname__, {}):
            return self.configuration["programs"][self.__toolname__]["abundance_tag"]
        else:
            return self.__tag__


class FilterGF(AtomicOperation):

    """This atomic operation will remove low-abundance fragments using the filter_assemblies_by_quant script."""

    def __init_subclass__(cls, **kwargs):

        pass

    def __init__(self,
                 asm_run: ShortAssembler,
                 tag, monoexonic_threshold, multiexonic_threshold
                 ):
        super().__init__()
        self.configuration = asm_run.configuration
        self.input["gf"] = asm_run.output["gf"]
        # Get some parameters from the parent node ...
        self.toolname = asm_run.toolname
        self.sample = asm_run.sample
        self.run, self.alrun, self.suffix = asm_run.run, asm_run.alrun, asm_run.suffix
        self.label = asm_run.label
        # Now that we have the parameters, define the rest ...
        self.output["gf"] = os.path.splitext(self.input["gf"])[0] + ".filtered" + os.path.splitext(self.input["gf"])[1]
        self.output["link"] = self.link

        self.message = "Filtering low-abundance transcripts for: {input[gf]}".format(input=self.input)
        self.tag, self.mono, self.multi = tag, monoexonic_threshold, multiexonic_threshold

    @property
    def rulename(self):
        return "abundance_filter_{label}_{toolname}_{run}_{alrun}".format(
            toolname=self.toolname, run=self.run, alrun=self.alrun, label=self.sample.label)

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]  # The script is dependent on Mikado

    @property
    def is_small(self):
        return False

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):

        load, tag, mono, multi, input, output = self.load, self.tag, self.mono, self.multi, self.input, self.output
        cmd = "{load} filter_assemblies_by_quant.py -q {tag} -m {mono} -mu {multi} {input[gf]} {output[gf]} "
        link_src = self.link_src
        cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def link_src(self):
        return os.path.join("..", self.toolname,
                            "{sample}-{run}-{alrun}".format(sample=self.sample.label, run=self.run, alrun=self.alrun),
                            os.path.basename(self.output["gf"]))

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "2-assemblies"))

    @property
    def link(self):
        return os.path.join(self.outdir, "output", "{toolname}-{sample}-{run}-{alrun}.{suffix}".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run,
            suffix=self.suffix, alrun=self.alrun
        ))


class AsmStats(AtomicOperation):

    def __init__(self, asm_run: [ShortAssembler, FilterGF]):

        super().__init__()
        self.configuration = asm_run.configuration
        try:
            self.input["gf"] = asm_run.output["link"]
        except KeyError:
            raise KeyError((asm_run.rulename, asm_run.output))
        self.sample = asm_run.sample
        self.output["stats"] = self.input["gf"] + ".stats"
        self.message = "Computing assembly stats for: {input[gf]}".format(input=self.input)
        self.log = self.output["stats"] + ".log"
        self.__asm_run = asm_run

    @property
    def rulename(self):
        return "asm_gf_stats_{}".format(os.path.basename(os.path.splitext(self.input["gf"])[0]))

    @property
    def cmd(self):
        load, output, input = self.load, self.output, self.input
        log = self.log
        cmd = "{load} mikado util stats {input[gf]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["mikado"]

    @property
    def threads(self):
        return 1

    @property
    def label(self):
        return self.__asm_run.label

    @property
    def is_small(self):
        return True
