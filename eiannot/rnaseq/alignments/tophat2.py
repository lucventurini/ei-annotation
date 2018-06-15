from ...abstract import AtomicOperation, EIWrapper, ShortSample, LongSample
from .abstract import IndexBuilder, ShortAligner, ShortWrapper
import os
import itertools


class TopHat2Index(IndexBuilder):

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)

        self.output = {"index": os.path.join(self._outdir, "{}_index.done".format(self.toolname))}
        self.log = os.path.join(self._outdir, "{}.index.log".format(self.toolname))
        self.touch = True
        self.message = "Indexing genome with {}".format(self.toolname)

    @property
    def out_prefix(self):
        return os.path.abspath(os.path.join(self._outdir[0], self.species))

    @property
    def toolname(self):
        return "tophat2"

    @property
    def loader(self):
        return ["tophat2"]

    @property
    def threads(self):
        """TopHat2 only builds in single core fashion."""
        return 1

    def cmd(self):
        load = self.load
        threads = 1
        input = self.input
        log = self.log
        align_dir = os.path.abspath(os.path.dirname(self.output["index"]))
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        extra = self.__configuration["programs"][self.toolname]["index"]
        cmd = "{load}"
        cmd += "bowtie2-build {input[genome]} {extra} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class TopHat2Aligner(ShortAligner):

    def __init__(self, index, sample, run, configuration, outdir, ref_transcriptome=None):

        super(TopHat2Aligner, self).__init__(configuration=configuration,
                                            sample=sample, run=run,
                                            outdir=outdir, index=index, ref_transcriptome=ref_transcriptome)
        self.output = {"bam": os.path.join(self.bamdir, "accepted_hits.bam"),
                       "link": self.link}

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
        cmd = "{load}"
        outdir = self.bamdir
        cmd += "tophat2 --num-threads={threads} --output-dir={outdir} --no-sort-bam"
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += "--min-intron-length={min_intron} --max-intron-length={max_intron} "
        if self.input.get("transcriptome", None):
            cmd += "--GTF={}".format(self.input.get("transcriptome"))
        extra = self.extra
        strand = self.strand
        infiles = self.input_reads
        threads = self.threads
        output = self.output
        index = self.index
        cmd += "{strand} {extra} {index} {infiles} > {log}"
        cmd += "| samtools view -b -@ {threads} - > {output[bam]}"
        link_src = self.link_src
        cmd += "&& ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def loader(self):
        return ["tophat2"]

    @property
    def toolname(self):
        return "tophat2"

    @property
    def strand(self):

        if self.sample.strandedness:
            return '--library-type={}'.format(self.sample)
        else:
            return ''


class TopHat2Flag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    def __init__(self, outdir, runs=[]):

        super().__init__()
        for number, run in enumerate(runs):
            self.input["run{}".format(number)] = run
        self.output["flag"] = os.path.join(outdir, "tophat2.done")
        self.touch = True


class TopHat2Wrapper(ShortWrapper):

    __indexer = TopHat2Index

    def __init__(self, configuration, prepare_flag):

        # First, we have to build the index

        super().__init__(configuration, prepare_flag)

        # Then we have to do all the alignments
        # Retrieve the running parameters

        if len(self.runs) > 0 and len(self.samples) > 0:
            # Start creating the parameters necessary for the run
            indexer = self.indexer(configuration, self.outdir)
            self.add_node(indexer)
            # Optionally build the reference splice catalogue
            star_runs = []
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                hisat_run = TopHat2Aligner(configuration=configuration,
                                           index=indexer.output["index"],
                                           sample=sample,
                                           outdir=self.outdir,
                                           run=run)
                star_runs.append(hisat_run)
            self.add_edges_from([(indexer, run) for run in star_runs])
            flag = TopHat2Flag(self.outdir, [run.output["link"] for run in star_runs])
            self.add_edges_from([(run, flag) for run in star_runs])
        self.finalise()

    @property
    def toolname(self):
        return "tophat2"

    @property
    def indexer(self):
        return TopHat2Index