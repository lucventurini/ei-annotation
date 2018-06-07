from eiannot.workflow import AtomicOperation, EIWrapper, ShortSample, LongSample
from .abstract import IndexBuilder, ShortAligner
import os
import itertools


class StarIndex(IndexBuilder):

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)

        self.output = {"index": os.path.join(self._outdir, "SAindex")}
        self.touch = False
        self.message = "Indexing genome with star"

    @property
    def toolname(self):
        return "star"

    @property
    def loader(self):
        return ["star"]

    def cmd(self):
        load = self.load
        threads = self.threads
        input = self.input
        log = self.log
        align_dir = os.path.abspath(os.path.dirname(self.output["index"]))
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        if "ref_transcriptome" in self.input:
            trans = "--sjdbGTFfile {input[ref_transcriptome]}"
        else:
            trans = ""
        extra = self.__configuration["programs"]["star"]["index"]
        cmd = "{load}"
        cmd += "cd {align_dir} && "
        cmd += "STAR --runThreadN {threads} --runMode GenomeGenerate --genomeDir . {trans}"
        cmd += "--genomeFastaFiles {input[genome]} {extra} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class StarAligner(ShortAligner):

    def __init__(self, index, sample, run, configuration, outdir, ref_transcriptome=None):

        super(StarAligner, self).__init__(configuration=configuration,
                                          sample=sample, run=run,
                                          outdir=outdir, index=index, ref_transcriptome=ref_transcriptome)

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

        snippet = "--readFilesIn {read1}"
        if read2:
            snippet += " {read2}"
        snippet = snippet.format(**locals())
        return snippet

    @property
    def cmd(self):

        load = self.load
        index = os.path.dirname(self.input["index"])
        outdir = self.bamdir
        threads = self.threads
        cmd = "{load}"
        cmd += "cd {outdir} &&"
        cmd += "STAR --runThreadN {threads} --runMode alignReads --genomeDir {index}"
        rfc = self.compression_option
        infiles = self.input_reads
        cmd += "{rfc} --runMode alignReads {infiles}"
        cmd += "--outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM XS NM MD"
        cmd += "--outSAMstrandField intronMotif"
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += "--alignIntronMin {min_intron} --alignIntronMax {max_intron} --alignMatesGapMax {max_intron}"
        if self.ref_transcriptome:
            ref_transcriptome = self.ref_transcriptome
            cmd += "--sjdbGTFfile {ref_transcriptome}"
        extra = self.extra
        bamdir = self.bamdir
        log = self.log
        cwd = os.getcwd()  # TODO: this is probably *not* what we need
        cmd += "--outFileNamePrefix {bamdir} {extra} > {log} 2>&1"
        cmd += " && cd {cwd} && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["star"]

    @property
    def toolname(self):
        return "star"

    @property
    def strand(self):
        """STAR does not accept specifying the strand of reads, so this property returns an empty string."""
        return ""


class StarFlag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    def __init__(self, outdir, runs=[]):

        super().__init__()
        for number, run in enumerate(runs):
            self.input["run{}".format(number)] = run
        self.output["flag"] = os.path.join(outdir, "star.done")
        self.touch = True


class StarWorkflow(EIWrapper):

    def __init__(self, configuration, outdir):

        # First, we have to build the index

        super().__init__()

        # Then we have to do all the alignments
        # Retrieve the running parameters
        runs = configuration["programs"]["star"]["runs"]
        samples = configuration["short_reads"]["samples"]

        if len(runs) > 0:
            # Start creating the parameters necessary for the run
            indexer = StarIndex(configuration, outdir)
            self.add_node(indexer)
            # Optionally build the reference splice catalogue
            star_runs = []
            for sample, run in itertools.product(samples, range(len(runs))):
                hisat_run = StarAligner(configuration=configuration,
                                        index=indexer.output["index"],
                                        sample=sample,
                                        outdir=outdir,
                                        run=run)
                star_runs.append(hisat_run)
            self.add_edges_from([(indexer, run) for run in star_runs])
            flag = StarFlag(outdir, [run.output["link"] for run in star_runs])
            self.add_edges_from([(run, flag) for run in star_runs])
