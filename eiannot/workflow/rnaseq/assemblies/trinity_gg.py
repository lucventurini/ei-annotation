from .abstract import ShortAssembler
from ... import AtomicOperation, EIWrapper, ShortSample
from ..alignments.gmap import GmapIndex, gmap_intron_lengths
from ..alignments.abstract import ShortAligner
import functools
import subprocess as sp
import re
import os


@functools.lru_cache(maxsize=4, typed=True)
def trinity_version(loader):
    cmd = "{} Trinity --version && set -u".format(loader)
    output = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout.read().decode()
    try:
        version = [_ for _ in output.split("\n") if re.match("Trinity version:", _)][0]
    except IndexError as exc:
        print("Error in getTrinityVersion")
        raise IndexError(exc)
    version = re.sub("Trinity version: [^_]*_[rv]", "", version)
    return version


class TrinityGGWrapper(EIWrapper):

    def __init__(self, configuration, bams, outdir):
        super().__init__()
        self._outdir = outdir
        self.configuration = configuration

        runs = configuration["programs"]["hisat"]["runs"]

        if len(runs) > 0:
            indexer = GmapIndex(configuration, outdir)

            trinities = []
            mappers = []

            for bam in bams:
                for run in runs:
                    trinity = TrinityGG(bam, configuration, outdir, run)
                    mapper = TrinityGmap()


    # rule
    # shell: "touch {output}"
    # trinity_all:
    # input: expand(ASM_DIR + "/output/trinity-{run2}-{alrun}.gff", run2=TRINITY_RUNS, alrun=ALIGN_RUNS)
    # output: ASM_DIR + "/trinity.done"


class TrinityGG(ShortAssembler):

    def __init__(self, bam, configuration, outdir, run):
        super().__init__(bam, configuration, outdir, run)
        self.input["reference"] = self.genome

        self.output = {"transcripts": os.path.join(self._outdir, "Trinity.fasta")}

        self.message = "Using trinity in genome guided mode to assemble (run {run}): {input[bam]}".format(
            input=self.input, run=run
        )

    @property
    def loader(self):
        return ["trinity"]

    @property
    def rulename(self):
        # TODO: implement
        return ""

    @property
    def cmd(self):
        extra = self.extra
        strand = self.strand
        base = self.base
        load = self.load
        max_intron = self.max_intron
        threads = self.threads
        input = self.input
        cmd = "{load} "
        cmd += "Trinity --seqType=fq {strand} --output={outdir} {extra}"
        cmd += "--genome_guided_max_intron={max_intron}"
        cmd += "--CPU={threads} {base}={input[bam]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def memory(self):
        # TODO: the location of this parameter will most probably be changed
        memory = self.configuration.get("tgg", dict()).get("max_mem", 5000)
        return memory

    @property
    def memory_string(self):
        memory = self.memory
        if str(memory).isdigit() is True:
            memory = int(memory)
            if memory >=1000:
                memory = "{}G".format(int(memory / 1000))
        if self.__version.startswith("201"):
            memory = "--JM={}".format(memory)
        else:
            memory = "--max_memory={}".format(memory)
        return memory

    @property
    def base(self):
        version = self.__version
        memory = self.memory_string
        if version.startswith("201"):  # Old versions
            reads = self.input_reads
            cmd = "{reads} --genome={ref} {memory} --genome_guided_use_bam".format(
                reads=reads, ref=self.genome, memory=self.memory
            )
        else:
            cmd = "--full_cleanup {memory} --genome_guided_bam".format(memory=memory)
        return cmd

    @property
    def strand(self):

        # TODO: double check the various strandedness option
        if self.sample.strandedness == "fr-firststrand":
            return "--SS_lib_type=RF"
        elif self.sample.strandedness == "fr-secondstrand":
            return "--SS_lib_type=FR"
        elif self.sample.strandedness == "f":
            return "--SS_lib_type=F"
        elif self.sample.strandedness == "r":
            return "--SS_lib_type=R"
        else:
            return ""

    @property
    def __version(self):
        return trinity_version(self.load)

    @property
    def input_reads(self):
        assert isinstance(self.sample, ShortSample)
        if self.sample.paired:
            return "--left={} --right={}".format(self.sample.read1, self.sample.read2)
        else:
            return "--single={} ".format(self.sample.read1)


class TrinityGmap(ShortAssembler):

    def __init__(self, transcripts, run, configuration, outdir):

        super().__init__(bam=None, run=run, configuration=configuration, outdir=outdir)
        self.input['transcripts'] = transcripts
        self.message = "Mapping trinity transcripts to the genome: {input[transcripts]}"

    @property
    def rulename(self):
        # TODO: implement
        return ''

    @property
    def toolname(self):
        return "trinity"

    @property
    def loader(self):
        return ["gmap"]

    @property
    def coverage(self):
        return self.configuration["programs"][self.toolname]["coverage"]

    @property
    def identity(self):
        return self.configuration["programs"][self.toolname]["identity"]

    @property
    def paths(self):
        return self.configuration["programs"][self.toolname]["paths"]

    @property
    def cmd(self):
        load = self.load
        cmd = "{load}"
        index_dir = ''  # TODO: define
        db = ''  # TODO: define
        strand = self.strand
        min_intron = self.min_intron
        max_intron = gmap_intron_lengths(self.load, self.max_intron)
        cmd += "gmap --dir={index_dir} {strand} --db {db} --min-intronlength={min_intron}"
        cmd += "{max_intron}"
        identity = self.identity
        coverage = self.coverage
        paths = self.paths
        cmd += "--format=3 --min-trimmed-coverage={coverage} --min-identity={identity} -n {paths}"
        threads = self.threads
        input = self.input
        output = self.output
        log = self.log
        cmd += "-t {threads} {input[transcripts]} > {output[gf]} 2> {log}"
        link_src = ""
        gff = ""
        cmd += "&& ln -sf {link_src} {output[gf]} && touch -h {output[gf]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def strand(self):
        if self.sample.stranded:
            return "-z sense_filter"
        else:
            return ""
