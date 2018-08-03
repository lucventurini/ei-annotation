from .abstract import IndexBuilder, IndexLinker, ShortAligner, LongAligner, ShortWrapper, LongWrapper
import os
import itertools


class StarIndexLink(IndexLinker):

    __toolname__ = "star"

    def __init__(self, configuration, outdir):
        super().__init__(configuration, outdir)
        self.input["index_files"] = [os.path.join(self.index_folder, self.index_name, fname) for fname in
                                     ["Genome", "SA", "chrLength.txt", "chrNameLength.txt",
                                      "genomeParameters.txt", "Log.out",
                                      "SAindex", "chrName.txt", "chrStart.txt"]]

        self.output = {"index": os.path.join(self.outdir, "SAindex")}
        self.touch = False
        self.message = "Linking pre-indexed genome with star"

    @property
    def loader(self):
        return []

    @property
    def cmd(self):
        outdir = self.outdir
        cmd = "mkdir -p {outdir} && cd {outdir} ".format(**locals())
        for fname in self.input["index_files"]:
            rel_path = os.path.relpath(fname, self.outdir)
            new_path = os.path.basename(fname)
            cmd += " && ln -sf {rel_path} {new_path}".format(**locals())

        return cmd

    @property
    def index(self):
        return self.output["index"]


class StarIndex(IndexBuilder):

    __toolname__ = "star"

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)

        self.output = {"index": os.path.join(self.outdir, "SAindex")}
        self.touch = False
        self.message = "Indexing genome with star"

    @property
    def loader(self):
        return ["star"]

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        input = self.input
        align_dir = os.path.abspath(os.path.dirname(self.output["index"]))
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        if self.transcriptome is not None:
            ref_transcriptome = os.path.abspath(
                self.input["transcriptome"])
            trans = " --sjdbGTFfile {ref_transcriptome} ".format(ref_transcriptome=ref_transcriptome)
        else:
            trans = " "
        genome = os.path.relpath(self.input["genome"],
                                 start=align_dir)
        log = os.path.relpath(self.log, start=align_dir)
        extra = self.extra
        cmd = "{load} cd {align_dir} && "
        cmd += "STAR --runThreadN={threads} --runMode genomeGenerate --genomeDir . {trans}"
        cmd += "--genomeFastaFiles {genome} -- {extra} --outStd Log > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def index(self):
        return self.output["index"]


class StarAligner(ShortAligner):

    def __init__(self, indexer, sample, run):

        super(StarAligner, self).__init__(indexer=indexer,
                                          sample=sample, run=run)
        self.output["bam"] = os.path.join(self.bamdir, "Aligned.out.bam")

    @property
    def compression_option(self):
        if self.sample.suffix == ".gz":
            return "--readFilesCommand 'gzip -dc'"
        elif self.sample.suffix in (".bz", ".bz2"):
            return "--readFilesCommand 'bzip2 -dc'"
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
        outdir = self.bamdir
        index = os.path.relpath(os.path.dirname(self.input["index"]), start=outdir)
        threads = self.threads
        cmd = "{load}"
        cmd += " mkdir -p {outdir} && cd {outdir} &&"
        cmd += " STAR --runThreadN={threads} --genomeDir {index} "
        rfc = self.compression_option
        infiles = self.input_reads
        cmd += "{rfc} --runMode alignReads {infiles}"
        cmd += " --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM XS NM MD"
        cmd += " --outSAMstrandField intronMotif"
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += " --alignIntronMin {min_intron} --alignIntronMax {max_intron} --alignMatesGapMax {max_intron}"
        if self.transcriptome:
            ref_transcriptome = os.path.abspath(self.transcriptome)
            cmd += " --sjdbGTFfile {ref_transcriptome}"
        extra = self.extra
        # bamdir = self.bamdir
        log = os.path.relpath(self.log, start=outdir)
        link_src = self.link_src
        cmd += "  {extra} > {log} 2>&1"
        final_dir = os.path.dirname(os.path.abspath(self.configuration["outdir"]))
        cmd += " && cd {final_dir} && ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        output = self.output
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


class StarWrapper(ShortWrapper):

    __toolname__ = "star"

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

    @property
    def indexer(self):
        if self.prebuilt:
            return StarIndexLink
        else:
            return StarIndex


class StarLong(LongAligner):

    __toolname__ = "star_long"

    def __init__(self, indexer: StarIndex, sample, run):
        if not isinstance(indexer, (StarIndex, StarIndexLink)):
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

        outdir = os.path.dirname(self.output["bam"])
        try:
            index = os.path.relpath(os.path.dirname(self.input["index"]), start=outdir)
        except KeyError:
            raise KeyError(self.input)

        log = os.path.relpath(self.log, start=outdir)
        input_reads = self.input_reads
        threads = self.threads
        min_intron, max_intron = self.min_intron, self.max_intron
        ref_transcriptome = self.ref_transcriptome
        load = self.load

        extra = self.extra

        cmd = "{load} "
        cmd += " mkdir -p {outdir} && cd {outdir} && "
        cmd += "STARlong --runThreadN={threads} --runMode alignReads --outSAMattributes NH HI NM MD "
        cmd += " --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 "
        cmd += " --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 "
        cmd += "--scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 "
        cmd += " --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 "
        cmd += " --alignTranscriptsPerWindowNmax 10000 --genomeDir {index}"
        cmd += " {input_reads} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif "
        cmd += " --alignIntronMin {min_intron} --alignIntronMax {max_intron} "
        cmd += " {ref_transcriptome} {extra} > {log} 2>&1 "

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
    def suffix(self):
        return ".bam"


class StarBam2Gtf(LongAligner):

    __toolname__ = "star"

    def __init__(self, aligner: StarLong):
        super().__init__(indexer=aligner.indexer, sample=aligner.sample, run=aligner.run)
        self.input = aligner.output
        self.output = {"link": self.link,
                       "gf": os.path.splitext(aligner.output["bam"])[0] + self.suffix}
        self.message = "Converting STAR long reads from BAM to GTF (sample: {sample.label} - run: {run})".format(
            sample=self.sample, run=self.run)

    @property
    def loader(self):
        return ["mikado"]

    @property
    def suffix(self):
        return ".gtf"

    @property
    def rulename(self):
        return "starbam2gtf_{sample}_{run}".format(sample=self.sample.label, run=self.run)

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        link_source = os.path.relpath(self.output["gf"], start=os.path.dirname(self.link))
        link = self.link
        cmd = "{load} bam2gtf.py {input[bam]} > {output[gf]} 2> {log} && ln -sf {link_source} {output[link]}".format(
            **locals())
        return cmd


class StarLongWrapper(LongWrapper):

    __toolname__ = "star_long"

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
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                star_run = StarLong(indexer=indexer, sample=sample, run=run)
                self.add_edge(indexer, star_run)
                bam2gtf_run = StarBam2Gtf(star_run)
                self.add_to_gfs(bam2gtf_run)
                self.add_edge(star_run, bam2gtf_run)

    @property
    def indexer(self):
        if self.prebuilt:
            return StarIndexLink
        else:
            return StarIndex
