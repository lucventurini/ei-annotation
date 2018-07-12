import abc
from ...abstract import ShortSample, LongSample, AtomicOperation, EIWrapper
import os
from .bam import BamIndex, BamSort, BamStats
import glob


class IndexBuilder(AtomicOperation, metaclass=abc.ABCMeta):

    """Abstract class for building the index of a genome, given a tool"""

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError("Operation {} does not have a defined toolname!".format(cls.__name__))
        super().__init_subclass__()


    def __init__(self, configuration, outdir):
        super(IndexBuilder, self).__init__()
        self.configuration = configuration
        self.input = {"genome": self.genome}

        if self.configuration.get("reference", dict()).get("transcriptome", ""):
            self.input["ref_transcriptome"] = os.path.abspath(self.configuration["reference"]["transcriptome"])
        self.log = os.path.join(outdir, "index", "log", "{}.log".format(self.toolname))
        self.__threads = 1

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def rulename(self):
        return "{toolname}_index".format(toolname=self.toolname)

    @property
    @abc.abstractmethod
    def index(self):
        pass

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", self.toolname)

    @property
    def extra(self):
        return self.configuration["programs"].get(self.toolname, dict()).get("index", {}).get("extra", '')

    @property
    def resources(self):

        """This property will define the resources to be used for each rule."""

        # default = self.configuration["resources"]["default"]
        #
        # if self.step and self.step in self.configuration["resources"]:
        #     default.update(self.configuration["resources"][self.__step])
        default = dict()
        assert "programs" in self.configuration, self.rulename

        for key in "memory", "queue", "threads":
            default[key] = self.__retrieve_resource_from_programs(key)
            if key == "threads":
                self.__threads = default[key]
        assert self.threads is not None

        return default

    def __retrieve_resource_from_programs(self, resource):

        res = None
        if self.toolname in self.configuration["programs"]:
            if ("index" in self.configuration["programs"][self.toolname] and
                    resource in self.configuration["programs"][self.toolname]["index"]):
                res = self.configuration["programs"][self.toolname]["index"][resource]
            elif resource in self.configuration["programs"][self.toolname]:
                res = self.configuration["programs"][self.toolname][resource]
            else:
                res = None
        if res is None:
            res = self.configuration["programs"]["default"][resource]
        # assert res is not None, (resource,)
        if resource == "threads":
            assert res != 0 and res is not None
            self.__threads = res
        return res

    @property
    def threads(self):
        if self.__threads is None:
            self.__set_threads()
        assert self.__threads is not None
        return self.__threads


class IndexLinker(AtomicOperation, metaclass=abc.ABCMeta):

    """Abstract atomic operation for linking a pre-built genome index into place, if available."""

    """Abstract class for building the index of a genome, given a tool"""

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError("Operation {} does not have a defined toolname!".format(cls.__name__))
        super().__init_subclass__()

    def __init__(self, configuration, outdir):
        super(IndexLinker, self).__init__()
        self.configuration = configuration
        self.input = {"genome": self.genome}
        self.log = os.path.join(outdir, "index", "log", "{}.log".format(self.toolname))

    @property
    def index_folder(self):
        return self.configuration["programs"][self.toolname]["index"]["index_folder"]

    @property
    def index_name(self):
        return self.configuration["programs"][self.toolname]["index"]["index_name"]

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def rulename(self):
        return "{toolname}_index_link".format(toolname=self.toolname)

    @property
    @abc.abstractmethod
    def index(self):
        pass

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", self.toolname)

    @property
    @abc.abstractmethod
    def cmd(self):
        pass

    @property
    def threads(self):
        return 1


class ShortAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, indexer, sample, run):

        super(ShortAligner, self).__init__()
        self.__indexer = indexer
        self.__configuration, self.__sample, self.__run = None, None, None
        self.input = indexer.output
        self.configuration = indexer.configuration
        self.__sample = None
        self.sample = sample
        self.run = run
        self.input["read1"] = self.sample.read1
        if self.sample.read2 is not None:
            self.input["read2"] = self.sample.read2
        self.log = os.path.join(self.outdir, "{toolname}-{sample}-{run}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run
        ))

        self.message = "Aligning input with {toolname} (sample {sample} - run {run})".format(
            toolname=self.toolname,
            sample=self.sample.label,
            run=self.run
        )
        self.output["link"] = self.link

    @property
    def indexer(self):
        return self.__indexer

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def bamdir(self):
        return os.path.join(self.outdir, self.toolname,
                            "{sample}-{run}".format(sample=self.sample.label, run=self.run))

    @property
    def link(self):
        return os.path.join(self.outdir, "output", "{toolname}-{sample}-{run}.bam".format(
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
        runs = self.configuration["programs"].get(self.toolname, dict()).get("runs", [])
        if run not in range(len(runs)):
            raise ValueError
        self.__run = run

    @property
    def extra(self):
        try:
            return self.configuration["programs"][self.toolname]["runs"][self.run]
        except (IndexError, TypeError) as exc:
            raise TypeError(self.configuration["programs"][self.toolname])

    @property
    def ref_transcriptome(self):
        return self.configuration.get("reference", dict()).get("transcriptome", "")

    @property
    def index(self):
        return self.indexer.index

    @property
    @abc.abstractmethod
    def strand(self):
        pass

    @property
    @abc.abstractmethod
    def input_reads(self):
        pass

    @property
    def link_src(self):
        return os.path.relpath(self.output["bam"], start=os.path.dirname(self.output["link"]))

    @property
    def rulename(self):
        # We have to give a unique name to each alignment
        return "{toolname}-{sample}-{run}".format(sample=self.sample.label,
                                                  run=self.run, toolname=self.toolname)

    @property
    def sample(self):
        return self.__sample

    @sample.setter
    def sample(self, sample):
        if isinstance(sample, str):
            if sample not in self.configuration['short_reads']:
                raise KeyError(sample)
            sample = self.configuration['short_reads'][sample]
        elif isinstance(sample, ShortSample):
            if sample.label not in self.configuration['short_reads']:
                raise KeyError(sample.label)
        self.__sample = sample

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments"))

    @property
    def alrun(self):
        return "{toolname}-{sample.label}-{run}".format(toolname=self.toolname,
                                                        sample=self.sample, run=self.run)


class LongAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError()
        super().__init_subclass__()

    def __init__(self,
                 indexer: IndexBuilder,
                 sample,
                 run):

        super().__init__()
        self.__indexer = indexer
        self.__configuration, self.__sample, self.__run = None, None, None
        self.input = indexer.output
        self.configuration = indexer.configuration
        self.__sample = None
        self.sample = sample
        self.run = run
        # self.input = {"indexer": indexer.index}
        self.input["read1"] = self.sample.read1
        self.log = os.path.join(self.outdir, "{toolname}-{sample}-{run}.log".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run
        ))

        self.message = "Aligning input with {toolname} (sample {sample} - run {run})".format(
            toolname=self.toolname,
            sample=self.sample.label,
            run=self.run
        )

    @property
    def toolname(self):
        return self.__toolname__

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

    @property
    def indexer(self):
        return self.__indexer

    @property
    def link(self):
        return os.path.join(self.outdir, "output", "{toolname}-{sample}-{run}{suffix}".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run, suffix=self.suffix)
                            )

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "rnaseq", "4-long-alignments")

    @property
    def extra(self):
        return self.configuration["programs"][self.toolname]["runs"][self.run]

    @property
    def label(self):
        return "{toolname}-{sample}-{run}".format(
            toolname=self.toolname, sample=self.sample.label, run=self.run)


class ShortWrapper(EIWrapper, metaclass=abc.ABCMeta):

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError("Wrapper {} does not have a defined toolname!".format(cls.__name__))
        cls.__final_rulename__ = "{toolname}_flag".format(toolname=cls.__toolname__)
        super().__init_subclass__()

    def __init__(self, configuration, prepare_flag):
        self.__finalised = False
        self.__prepare_flag = prepare_flag
        super().__init__()
        self.add_node(prepare_flag.exit)
        self.__bam_rules = set()
        self.configuration = configuration
        self.__stats = []

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "{toolname}.done".format(toolname=self.toolname))

    def finalise(self):

        if self.__finalised:
            return
        new_bams = set()

        for bam in self.bams:
            sorter = BamSort(bam)
            self.add_edge(bam, sorter)
            indexer = BamIndex(sorter)
            self.add_edge(sorter, indexer)
            stater = BamStats(indexer)
            self.add_edge(indexer, stater)
            new_bams.add(stater)
            self.__stats.append(stater)
        self.__bam_rules = new_bams
        self.__add_flag_to_inputs()
        self.add_final_flag()
        self.__finalised = True

    def add_to_bams(self, rule):
        if not isinstance(rule, ShortAligner):
            raise TypeError
        if "link" not in rule.output:
            raise KeyError("Link not found for rule {}; output: {}; input: {}".format(
                rule.rulename,
                rule.output,
                rule.input))
        self.__bam_rules.add(rule)

    @property
    def prebuilt(self):
        if "index" not in self.configuration["programs"][self.toolname]:
            print("Index not found for {}".format(self.toolname))
            return False

        if (self.configuration["programs"][self.toolname]["index"].get("index_folder", None) and
                self.configuration["programs"][self.toolname]["index"].get("index_name", None)):
            return True

    @property
    def bams(self):
        """This property will the output files (ie BAMs) of the wrapper."""
        return self.__bam_rules

    @property
    def runs(self):
        return self.configuration["programs"].get(self.toolname, dict()).get("runs", [])

    def __add_flag_to_inputs(self):
        for rule in self.nodes:
            if rule == self.__prepare_flag.exit:
                continue
            rule.input["prep_flag"] = self.__prepare_flag.output["fai"]
            self.add_edge(self.__prepare_flag.exit, rule)

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
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments"))


class LongWrapper(EIWrapper, metaclass=abc.ABCMeta):

    def __init_subclass__(cls):

        if not hasattr(cls, "__toolname__"):
            raise NotImplementedError()
        cls.__final_rulename__ = "{toolname}_flag".format(toolname=cls.__toolname__)
        super().__init_subclass__()

    def __init__(self, prepare_flag):
        self.__prepare_flag = prepare_flag
        super().__init__()
        self.__gf_rules = set()
        self.configuration = prepare_flag.configuration
        self.__finalised = False

    def finalise(self):
        if self.__finalised is True:
            return

        new_gfs = set()

        for gf in self.gfs:
            stats = LongAlignerStats(gf)
            self.add_edge(gf, stats)
            new_gfs.add(stats)

        self.__gf_rules = new_gfs
        self.add_flag_to_inputs(self.__prepare_flag, "prep_flag", "fai")
        self.add_final_flag()
        self.__finalised = True

    @property
    def prebuilt(self):
        if "index" not in self.configuration["programs"][self.toolname]:
            print("Index not found for {}".format(self.toolname))
            return False

        if (self.configuration["programs"][self.toolname]["index"].get("index_folder", None) and
                self.configuration["programs"][self.toolname]["index"].get("index_name", None)):
            return True

    @property
    @abc.abstractmethod
    def indexer(self):
        pass

    @property
    def samples(self):
        return self.configuration["long_reads"]

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "4-long-alignments"))

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "{toolname}.done".format(toolname=self.toolname))

    def add_to_gfs(self, rule):
        if not isinstance(rule, LongAligner):
            raise TypeError
        if "link" not in rule.output:
            raise KeyError("Link not found for rule {}".format(rule.rulename))
        self.__gf_rules.add(rule)

    @property
    def runs(self):
        return self.configuration["programs"].get(self.toolname, dict()).get("runs", [])

    @property
    def gfs(self):
        return self.__gf_rules

    @property
    def toolname(self):
        return self.__toolname__


class LongAlignerStats(AtomicOperation):

    def __init__(self, aligner: LongAligner):
        super().__init__()
        self.configuration = aligner.configuration
        self.__aligner = aligner
        self.input = aligner.output
        self.outdir = aligner.outdir
        try:
            self.output = {"stats": os.path.splitext(self.input["link"])[0] + "stats"}
        except KeyError:
            raise KeyError(type(aligner), aligner.rulename)
        self.message = "Calculating statistics for: {input[link]}".format(input=self.input)
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return "long_alignment_stats_{toolname}-{sample}_{run}".format(
            toolname=self.__aligner.toolname,
            sample=self.__aligner.sample,
            run=self.__aligner.run
        )

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["gf"])

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[link]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def sample(self):
        return self.__aligner.sample

    @property
    def label(self):
        return self.__aligner.label
