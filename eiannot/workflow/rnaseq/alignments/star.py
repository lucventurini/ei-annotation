from eiannot.workflow import AtomicOperation, EIWrapper, ShortSample, LongSample
from .abstract import IndexBuilder, ShortAligner, LongAligner, ShortWrapper, LongWrapper
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

    @property
    def index(self):
        return self.output["index"]


class StarAligner(ShortAligner):

    def __init__(self, indexer, sample, run):

        super(StarAligner, self).__init__(indexer=indexer,
                                          sample=sample, run=run)

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


class StarWrapper(ShortWrapper):

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
                star_run = StarAligner(
                    indexer=indexer,
                    sample=sample, run=run)
                star_runs.append(star_run)
                self.add_to_bams(star_run)
            self.add_edges_from([(indexer, run) for run in star_runs])
            flag = StarFlag(self.outdir, [run.output["link"] for run in star_runs])
            self.add_edges_from([(run, flag) for run in star_runs])

    @property
    def toolname(self):
        return "star"

    @property
    def indexer(self):
        return StarIndex


class StarLong(LongAligner):

    def __init__(self, indexer: StarIndex, sample, run):
        if not isinstance(indexer, StarIndex):
            raise TypeError("Invalid indexer: {}".format(type(indexer)))

        super().__init__(indexer, sample, run)
        self.output = {"bam": self.bam}
        self.message = "Aligning long read sample {sample.label} with STAR, run {run}".format(
            sample=self.sample, run=self.run)

    @property
    def loader(self):
        return ["star"]

    @property
    def cmd(self):
        log = self.log  # TODO: implement
        index = os.path.dirname(self.input["index"])
        input_reads = self.input_reads
        threads = self.threads
        min_intron, max_intron = self.min_intron, self.max_intron
        ref_transcriptome = self.ref_transcriptome

        cmd = "{load} STARlong --runThreadN {threads} --runMode alignReads --outSAMattributes NH HI NM MD "
        cmd += " --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 "
        cmd += " --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 "
        cmd += "--scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 "
        cmd += " --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 "
        cmd += " --alignTranscriptsPerWindowNmax 10000 --genomeDir {index}"
        cmd += " {params.infiles} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif "
        cmd += " --alignIntronMin {min_intron} --alignIntronMax {max_intron} "
        cmd += " {ref_transcriptome} --outFileNamePrefix {params.outdir}/ > {log} 2>&1"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def input_reads(self):
        return "--readFilesIn {}".format(os.path.abspath(self.input["read1"]))

    @property
    def bam(self):
        return os.path.join(self.outdir, "lr_star", "{label}-{run}", "Aligned.out.bam").format(
            label=self.sample.label,
            run=self.run
        )

    @property
    def ref_transcriptome(self):
        return ""  # TODO: implement

    @property
    def rulename(self):
        return "star_long_{sample.label}_{run}".format(sample=self.sample, run=self.run)

    @property
    def toolname(self):
        return "star"

    @property
    def suffix(self):
        return ".bam"


class StarBam2Gtf(LongAligner):

    def __init__(self, aligner: StarLong):
        super().__init__(indexer=aligner.indexer, sample=aligner.sample, run=aligner.run)
        self.input = aligner.output
        self.output = {"gf": os.path.join(self.outdir, "lr_output", "lr_star-{sample}-{run}.gtf").format(
            sample=self.sample.label, run=self.run)}
        self.message = "Converting STAR long reads from BAM to GTF (sample: {sample.label} - run: {run})".format(
            sample=self.sample, run=self.run)

    @property
    def loader(self):
        return ["mikado"]

    @property
    def suffix(self):
        return ".gtf"

    @property
    def toolname(self):
        return "star"

    @property
    def rulename(self):
        return "starbam2gtf_{sample}_{run}".format(sample=self.sample.label, run=self.run)

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} bam2gtf.py {input[bam]} > {output[gf]} 2> {log}".format(**locals())
        return cmd


class StarLongWrapper(LongWrapper):

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
                star_run = StarLong(indexer=indexer, sample=sample, run=run)
                self.add_edge(indexer, star_run)
                bam2gtf_run = StarBam2Gtf(star_run)
                star_runs.append(bam2gtf_run)

            flag = StarLongFlag(self.outdir, [run.output["gf"] for run in star_runs])
            self.add_edges_from([(run, flag) for run in star_runs])

    @property
    def toolname(self):
        return "star_long"

    @property
    def indexer(self):
        return StarIndex


class StarLongFlag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    def __init__(self, outdir, runs=[]):

        super().__init__()
        for number, run in enumerate(runs):
            self.input["run{}".format(number)] = run
        self.output["flag"] = os.path.join(outdir, "star_long.done")
        self.touch = True
