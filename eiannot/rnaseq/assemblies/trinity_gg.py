from .abstract import ShortAssembler, ShortAssemblerWrapper
from ...abstract import AtomicOperation, EIWrapper, ShortSample
from ..alignments.gmap import GmapIndex, gmap_intron_lengths
from ..alignments.abstract import ShortAligner  # TODO: got to decide how to set the links
import functools
import subprocess as sp
import re
import os
import itertools


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


class TrinityGGWrapper(ShortAssemblerWrapper):

    def __init__(self, aln_wrapper):
        super().__init__(aln_wrapper)

        if len(self.runs) > 0 and len(bams) > 0:
            indexer = GmapIndex(configuration, self.outdir)

            trinities = []
            mappers = []

            for bam, run in itertools.product(bams, range(len(self.runs))):
                trinity = TrinityGG(bam, configuration, self.outdir, run)
                trinities.append(trinity)
                mapper = TrinityGmap(trinity, indexer)
                self.add_edge(trinity, mapper)
                self.add_edge(indexer, mapper)
                mappers.append(mapper)
                self.add_to_gf(mapper)

            flag = TrinityFlag(self.outdir, mappers)
            self.add_edges_from([mapper, flag] for mapper in mappers)

    @property
    def toolname(self):
        return "trinity"


class TrinityFlag(AtomicOperation):

    def __init__(self, outdir, gmaps):
        super().__init__()
        self.input = {"gffs": [gmap.output["gf"] for gmap in gmaps]}
        self.touch = True
        self.output = {"flag": os.path.join(outdir, "trinity_gg.done")}


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

    def __init__(self, trinitygg: TrinityGG, index: GmapIndex):

        run = trinitygg.run
        configuration = trinitygg.configuration
        outdir = trinitygg._outdir
        super().__init__(bam=None, run=run, configuration=configuration, outdir=outdir)
        self.input['transcripts'] = trinitygg.output["transcripts"]
        self.input.update(index.output)
        self.message = "Mapping trinity transcripts to the genome: {input[transcripts]}".format(
            input=self.input
        )
        self._trinitygg = trinitygg
        self._gmapdb = index
        self.log = os.path.join(self._outdir, "logs", "trinitygmap-{run}-{alrun}.log".format(
            run=run, alrun=trinitygg.alrun))

    @property
    def strand(self):
        if self.sample.stranded:
            return " -z sense_force "
        else:
            return ""

    @property
    def input_reads(self):
        return self.input["transcripts"]

    @property
    def suffix(self):
        return ".gff"

    @property
    def rulename(self):
        return 'trinity-gmap-{run}-{alrun}'.format(run=self.run, alrun=self._trinitygg.alrun)

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
        index_dir = os.path.dirname(self._gmapdb.out_prefix)
        db = os.path.basename(self._gmapdb.out_prefix)
        strand = self.strand
        min_intron = self.min_intron
        max_intron = gmap_intron_lengths(self.load, self.max_intron)
        cmd += "gmap --dir={index_dir} {strand} --db {db} --min-intronlength={min_intron} "
        cmd += "{max_intron} "
        identity = self.identity
        coverage = self.coverage
        paths = self.paths
        cmd += "--format=3 --min-trimmed-coverage={coverage} --min-identity={identity} -n {paths} "
        threads = self.threads
        input = self.input
        output = self.output
        log = self.log
        cmd += "-t {threads} {input[transcripts]} > {output[gf]} 2> {log} "
        link_src = self.link_src
        cmd += "&& ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def strand(self):
        if self.sample.stranded:
            return "-z sense_filter"
        else:
            return ""