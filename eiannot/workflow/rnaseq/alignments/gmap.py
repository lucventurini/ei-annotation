from ...abstract import AtomicOperation, EIWrapper, ShortSample
from .abstract import IndexBuilder, ShortAligner, ShortWrapper, LongWrapper, LongAligner
import os
import itertools
import functools
import subprocess


@functools.lru_cache(4, typed=True)
def gmap_intron_lengths(loader, max_intron):
    """Function to check the exact format of GMAP intron lengths"""
    cmd = "{} gmap --help".format(loader)
    output = subprocess.Popen(cmd, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.read().decode()
    if "--max-intronlength-middle" in output:
        return "--max-intronlength-middle={mi} --max-intronlength-ends={mi}".format(mi=max_intron)
    else:
        return "--intronlength={}".format(max_intron)


class GsnapWrapper(ShortWrapper):

    def __init__(self, configuration, outdir):

        # First, we have to build the index

        super().__init__(configuration)

        # Then we have to do all the alignments
        # Retrieve the running parameters
        if len(self.runs) > 0 and len(self.samples) > 0:
            # Start creating the parameters necessary for the run
            self.add_node(self.indexer)
            # Optionally build the reference splice catalogue
            gsnap_runs = []
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                hisat_run = GsnapAligner(configuration=configuration,
                                         index=self.indexer.out_prefix,
                                         sample=sample,
                                         outdir=outdir,
                                         run=run)
                gsnap_runs.append(hisat_run)
            self.add_edges_from([(self.indexer, run) for run in gsnap_runs])
            flag = GsnapFlag(outdir, [run.output["link"] for run in gsnap_runs])
            self.add_edges_from([(run, flag) for run in gsnap_runs])


class GsnapFlag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    def __init__(self, outdir, runs=[]):

        super().__init__()
        self.input["runs"] = []
        for number, run in enumerate(runs):
            self.input["run{}".format(number)] = run
        self.output["flag"] = os.path.join(outdir, "gsnap.done")
        self.touch = True

    @property
    def rulename(self):
        return "gsnap_flag"


class GmapIndex(IndexBuilder):

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        # TODO: probably the input should be the cleaned up genome
        self.output = {"index": os.path.join(self._outdir, "{}.sachildguide1024".format(self.species))}
        self.touch = True

    @property
    def toolname(self):
        return "gmap"

    @property
    def out_prefix(self):
        return os.path.abspath(os.path.dirname(self.output["index"]))

    @property
    def message(self):

        return "Building the GMAP index for {} in {}".format(self.input["genome"],
                                                             self._outdir)

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):

        load = None
        threads = self.threads
        input = self.input["genome"]
        outdir = os.path.dirname(os.path.dirname(self.output["index"]))
        species = self.species
        log = self.log
        cmd = "{load} gmap_build --dir={outdir} --db={species} {extra} {input[genome]} > {log} 2>&1".format(
            **locals()
        )
        return cmd

    @property
    def loader(self):
        return ["gmap"]

    @property
    def indexdir(self):
        return os.path.abspath(os.path.dirname(os.path.dirname(self.output["index"])))

    @property
    def dbname(self):
        return os.path.basename(os.path.dirname(self.output["index"]))


class GsnapAligner(ShortAligner):
    # infiles = lambda wildcards: tophatInput(wildcards.sample)  # Can use tophat function safely here

    def __init__(self,
                 indexer,
                 sample,
                 run):

        super().__init__(indexer=indexer, sample=sample, run=run)
        self.output["link"] = self.link
        self.output["bam"] = os.path.join(self.bamdir, "gsnap.bam")

    @property
    def compression_option(self):
        if self.sample.suffix == ".gz":
            return "--readFilesCommand zcat"
        elif self.sample.suffix in (".bz", ".bz2"):
            return "--readFilesCommand bzcat"
        else:
            return ""

    @property
    def input_reads(self):
        read1 = os.path.abspath(self.input["read1"])
        read2 = os.path.abspath(self.input["read2"])
        if read2:
            snippet = "{read1} {read2}"
        else:
            snippet = "{read1}"
        return snippet.format(**locals())

    @property
    def cmd(self):

        load = self.load
        index = os.path.dirname(self.input["index"])
        output = self.output
        threads = self.threads
        cmd = "{load}"
        indexdir = self.indexer.indexdir
        dbname = self.indexer.dbname
        cmd += "gsnap --dir={indexdir} --db={dbname} {extra} --novelsplicing=1 "
        min_intron, max_intron = self.min_intron, self.max_intron
        extra = self.extra
        bamdir = self.bamdir
        log = self.log
        cmd += "--localsplicedist={max_intron} --nthreads={threads} --format=sam --npaths=20"
        infiles = self.input_reads
        cmd += "{infiles} 2> {log}"
        cmd += "| samtools view -b -@ {threads} - > {output[bam]}"
        cwd = os.getcwd()
        cmd += " && cd {cwd} && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["gmap", "samtools"]

    @property
    def toolname(self):
        return "gsnap"

    @property
    def indexname(self):
        return "gmap"

    @property
    def strand(self):
        """GSNAP does not accept specifying the strand of reads, so this property returns an empty string."""
        return ""


class GmapLongReads(LongAligner):

    def __init__(self, indexer: GmapIndex, sample, run):

        super().__init__(indexer=indexer,
                         sample=sample, run=run)
        self.output = {"link": self.link,
                       "gf": os.path.join()}

        self.message = "Mapping long reads to the genome with gmap (sample: {sample.label} - run: {run})".format(
            sample=self.sample, run=self.run)

    @property
    def loader(self):
        return ["gmap"]

    @property
    def toolname(self):
        return "gmap"

    @property
    def suffix(self):
        return ".gff"

    @property
    def rulename(self):
        return "gmap_long_{sample.label}_{run}".format(sample=self.sample, run=self.run)

    @property
    def cmd(self):

        load = self.load
        max_intron = self.max_intron_cli
        threads = self.threads
        min_intron = self.min_intron
        index = self.indexer.indexdir
        dbname = self.indexer.dbname
        cmd = "{load} gmap --dir={index} --db {dbname} --min-intronlength={min_intron} {max_intron}"
        input, output, log = self.input, self.output, self.log
        cmd += "--format=3 {input[read1]} > {output[gf]} 2> {log} "
        link_src = os.path.relpath(self.link, start=os.path.dirname(self.output["gf"]))
        cmd += "&& ln -sf {link_src} {output[link]} && touch -h {output[link]}"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def max_intron_cli(self):
        return gmap_intron_lengths(self.load, self.max_intron)
