from .abstract import IndexBuilder, ShortAligner, ShortWrapper, LongWrapper, LongAligner, IndexLinker
import os
import itertools
import functools
import subprocess
import glob
import re


@functools.lru_cache(4, typed=True)
def gmap_intron_lengths(loader, max_intron, max_intron_middle=None):
    """Function to check the exact format of GMAP intron lengths"""
    cmd = "{} gmap --help".format(loader)
    output = subprocess.Popen(cmd, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.read().decode()
    if max_intron_middle is None:
        max_intron_middle = max_intron
    if "--max-intronlength-middle" in output:
        return "--max-intronlength-middle={mim} --max-intronlength-ends={mi}".format(mi=max_intron,
                                                                                     mim=max_intron_middle)
    else:
        return "--intronlength={}".format(max_intron)


class GsnapWrapper(ShortWrapper):

    __toolname__ = "gsnap"

    def __init__(self, configuration, prepare_flag):

        # First, we have to build the index

        super().__init__(configuration, prepare_flag=prepare_flag)

        # Then we have to do all the alignments
        # Retrieve the running parameters
        if len(self.runs) > 0 and len(self.samples) > 0:
            # Start creating the parameters necessary for the run
            # Optionally build the reference splice catalogue
            gsnap_runs = []
            indexer = self.indexer(self.configuration, self.outdir)
            self.add_node(indexer)
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                gsnap_run = GsnapAligner(indexer=indexer,
                                         sample=sample,
                                         run=run)
                self.add_to_bams(gsnap_run)
                gsnap_runs.append(gsnap_run)
            self.add_edges_from([(indexer, run) for run in gsnap_runs])

    @property
    def indexer(self):
        if self.prebuilt:
            return GmapLink
        else:
            return GmapIndex


class GmapIndex(IndexBuilder):

    __toolname__ = "gmap"

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        self.output = {"flag": os.path.join(self.outdir, "gmap_index.done")}
        self.touch = True

    @property
    def out_prefix(self):
        return os.path.abspath(os.path.dirname(self.output["index"]))

    @property
    def message(self):

        return "Building the GMAP index for {} in {}".format(self.input["genome"],
                                                             self.outdir)

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):

        load = self.load
        input = self.input
        outdir = os.path.dirname(self.output["flag"])
        species = self.species
        log = self.log
        extra = self.extra

        cmd = "{load} gmap_build --dir={outdir} --db={species} {extra} {input[genome]} > {log} 2>&1".format(**locals())

        return cmd

    @property
    def loader(self):
        return ["gmap"]

    @property
    def indexdir(self):
        return os.path.abspath(os.path.dirname(self.output["flag"]))

    @property
    def dbname(self):
        return self.species

    @property
    def index(self):
        return self.dbname


class GmapLink(IndexLinker):

    __toolname__ = "gmap"

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        self.input["index_folder"] = self.index_folder
        self.input["index_files"] = glob.glob(os.path.join(self.index_folder,
                                                           self.index_name,
                                                           "{index_name}*".format(index_name=self.index_name)))
        self.input["index_files"] += glob.glob(os.path.join(self.index_folder,
                                                            "{index_name}*".format(index_name=self.index_name)))
        self.output = {"flag": os.path.join(self.outdir, "gmap_index.done")}
        self.touch = True

    @property
    def cmd(self):
        outdir = self.outdir
        species = self.species

        cmd = "mkdir -p {outdir} && mkdir -p {outdir}/{species} && cd {outdir} ".format(**locals())
        for fname in glob.glob(os.path.join(
                os.path.abspath(self.index_folder),
                "{index_name}*".format(index_name=self.index_name))):
            if os.path.isdir(fname):
                continue
            link_src = os.path.relpath(os.path.abspath(fname), start=os.path.join(self.outdir))
            link_dest = self.species + re.sub(self.index_name, '', os.path.basename(fname))
            cmd += " && cp {link_src} {link_dest}".format(**locals())
        # Now inside the folder

        species = self.species
        cmd += " && cd {species} ".format(**locals())
        for fname in glob.glob(os.path.join(
                os.path.abspath(self.index_folder),
                self.index_name,
                "{index_name}*".format(index_name=self.index_name))):
            if os.path.isdir(fname):
                continue
            link_src = os.path.relpath(os.path.abspath(fname), start=os.path.join(self.outdir,
                                                                                  self.species))
            link_dest = self.species + re.sub(self.index_name, '', os.path.basename(fname))
            cmd += " && cp {link_src} {link_dest}".format(**locals())

        return cmd

    @property
    def indexdir(self):
        return os.path.abspath(os.path.dirname(self.output["flag"]))

    @property
    def loader(self):
        return []

    @property
    def dbname(self):
        return self.species

    @property
    def index(self):
        return self.dbname


class GsnapAligner(ShortAligner):
    # infiles = lambda wildcards: tophatInput(wildcards.sample)  # Can use tophat function safely here

    def __init__(self,
                 indexer,
                 sample,
                 run):

        super().__init__(indexer=indexer, sample=sample, run=run)
        self.output["link"] = self.link
        self.output["bam"] = os.path.join(self.bamdir, "gsnap.bam")
        # self.index = indexer.index

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
        # index = self.index
        output = self.output
        threads = self.threads
        cmd = "{load}"
        indexdir = self.indexer.indexdir
        dbname = self.indexer.dbname
        genome = self.genome
        cmd += "$(determine_gmap.py -s {genome}) --dir={indexdir} --db={dbname} {extra} --novelsplicing=1 "
        min_intron, max_intron = self.min_intron, self.max_intron
        extra = self.extra
        bamdir = self.bamdir
        log = self.log
        compression = self.compression
        cmd += " {compression} --localsplicedist={max_intron} --nthreads={threads} --format=sam --npaths=20 "
        infiles = self.input_reads
        cmd += " {infiles} 2> {log}"
        cmd += " | samtools view -b -@ {threads} - > {output[bam]}"
        cwd = os.getcwd()
        link_src = os.path.relpath(output["bam"], start=os.path.dirname(self.output["link"]))
        cmd += " && cd {cwd} && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["gmap", "samtools", "ei-annotation"]

    __toolname__ = "gsnap"

    @property
    def indexname(self):
        return "gmap"

    @property
    def strand(self):
        """GSNAP does not accept specifying the strand of reads, so this property returns an empty string."""
        return ""

    @property
    def compression(self):
        if self.sample.suffix == ".gz":
            return " --gunzip "
        elif self.sample.suffix == ".bz2":
            return "--bunzip2"
        else:
            return ""


class GmapLongReads(LongAligner):

    def __init__(self, indexer: GmapIndex, sample, run):

        super().__init__(indexer=indexer,
                         sample=sample, run=run)

        self.output = {
            "link": self.link,
            "gf": os.path.join(self.outdir,  "gmap", "{sample}-{run}", "gmap-{sample}-{run}{suffix}").format(
                sample=self.sample.label, run=self.run, suffix=self.suffix)
        }

        self.message = "Mapping long reads to the genome with gmap (sample: {sample.label} - run: {run})".format(
            sample=self.sample, run=self.run)

    @property
    def loader(self):
        return ["gmap", "ei-annotation"]

    __toolname__ = "gmap"

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
        genome = self.genome

        cmd = "{load} $(determine_gmap.py {genome}) --dir={index} --db {dbname} "
        cmd += " --min-intronlength={min_intron} {max_intron} "
        input, output, log = self.input, self.output, self.log
        coverage, identity, cross, strand = self.coverage, self.identity, self.cross_species, self.strand
        cmd += " {coverage} {identity} {cross} {strand}"
        threads = self.threads
        extra = self.extra
        cmd += " --format=2 {extra} --nthreads={threads} {input[read1]} > {output[gf]} 2> {log} "
        link_src = os.path.relpath(self.output["gf"], start=os.path.dirname(self.output["link"]))
        cmd += " && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def cross_species(self):
        if self.sample.type in ("pacbio", "ont", "ont-direct"):
            return " --cross-species "
        else:
            return ""

    @property
    def coverage(self):
        if "coverage" in self.configuration["programs"][self.toolname]:
            return "--min-trimmed-coverage={}".format(self.configuration["programs"][self.toolname]["coverage"])
        else:
            return ""

    @property
    def identity(self):
        if "identity" in self.configuration["programs"][self.toolname]:
            return "--min-identity={}".format(self.configuration["programs"][self.toolname]["identity"])
        else:
            return ""

    @property
    def strand(self):
        if self.sample.strandedness == "fr-firststrand":
            return " -z sense_force "
        elif self.sample.strandedness == "fr-secondstrand":
            return " -z antisense_force "
        else:
            return " "

    @property
    def max_intron_middle(self):
        return self.configuration["programs"][self.toolname].get("max_intron_middle", None)

    @property
    def max_intron_cli(self):
        return gmap_intron_lengths(self.load, self.max_intron, max_intron_middle=self.max_intron_middle)


class GmapLongWrapper(LongWrapper):

    def __init__(self, prepare_flag):

        # First, we have to build the index

        super().__init__(prepare_flag)

        # Then we have to do all the alignments
        # Retrieve the running parameters

        if len(self.runs) > 0 and len(self.samples) > 0:
            # Start creating the parameters necessary for the run
            indexer = self.indexer(self.configuration, self.outdir)
            self.add_node(indexer)
            # Optionally build the reference splice catalogue
            gmap_runs = []
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                gmap_run = GmapLongReads(indexer=indexer, sample=sample, run=run)
                gmap_runs.append(gmap_run)
                self.add_edge(indexer, gmap_run)
                self.add_to_gfs(gmap_run)

    __toolname__ = "gmap"

    @property
    def indexer(self):
        if self.prebuilt:
            return GmapLink
        else:
            return GmapIndex
