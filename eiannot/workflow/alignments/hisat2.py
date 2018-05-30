from .. import AtomicOperation, EIWrapper
import os
import itertools



class HisatWrapper(EIWrapper):

    def __init__(self, configuration, outdir):

        # First, we have to build the index

        indexer = HisatBuild(configuration, outdir)

        # Optionally build the reference splice catalogue


        # Then we have to do all the alignments

        alignments = []
        # Retrieve the running parameters
        runs = configuration["programs"]["hisat"]["runs"]
        samples = configuration["short_reads"]["samples"]

        if len(runs) > 0:
            # Start creating the parameters necessary for the run
            indexer = HisatBuild(configuration, outdir)
            self.add_node(indexer)
            for sample, run  in itertools.product(samples, range(len(runs))):






class HisatFlag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    pass



class HisatLink(AtomicOperation):

    @property
    def cmd(self):
        return "ln -s {input} {output}".format(input=self.input[0],
                                               output=self.output[0])


class HisatAligner(AtomicOperation):

    __loader = ["hisat", "samtools"]

    def __init__(self, configuration, outdir, sample, run):

        self.__sample = sample
        self._outdir = outdir
        self.__run = run
        self.__configuration = configuration
        self.output = {"bam": os.path.join(self._outdir, "hisat", "{sample}-{run,\d+}", "hisat.bam"),
                       "link": os.path.join(self._outdir, "output", "hisat-{sample}-{run,\d+}.bam")}

    @property
    def message(self):
        message = "Aligning input with hisat (sample {sample} - run {run})".format(
            sample=self.__sample,
            run=self.__run
        )
        return message

    @property
    def cmd(self):

        load = None  # TODO: maybe put this into the AtomicOperation class?

        cmd = "{load} ".format(load=load if load else "")
        cmd += "hisat2 -p {threads} "
        cmd += "--min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} "
        if self.input.get("transcriptome", None):
            cmd += "--known-splicesite-infile={}".format(self.input.get("transcriptome"))
        extra = self.__configuration["programs"][self.__run]["runs"]
        strand = ""  # TODO: implement
        infiles = ""  # TODO: implement
        cmd += "{strand} {extra} -x {index} {infiles} > {log} | samtools view -b -@ {threads} - > {output}".format(
            strand=strand, infiles=infiles,
            threads=threads, output=self.output[0])
        return cmd

        # rule align_hisat:
    # 	input:
    # 		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
    # 		index=rules.align_hisat_index.output
    # 	output:
    # 		bam=ALIGN_DIR+"/hisat/{sample}-{run,\d+}/hisat.bam",
    # 		link=ALIGN_DIR+"/output/hisat-{sample}-{run,\d+}.bam"
    # 	params:
    # 		indexdir=ALIGN_DIR+"/hisat/index/"+NAME,
    # 		load=loadPre(config, "hisat"),
    # 		load_samtools=loadPre(config, "samtools"),
    # 		link_src="../hisat/{sample}-{run}/hisat.bam",
    # 		extra=lambda wildcards: config["align_methods"]["hisat"][int(wildcards.run)],
    # 		ss_gen="hisat2_extract_splice_sites.py " + REF_TRANS + " > " + ALIGN_DIR + "/hisat/{sample}-{run}/splice_sites.txt &&" if REF_TRANS else "",
    # 		trans="--known-splicesite-infile=" + ALIGN_DIR + "/hisat/{sample}-{run}/splice_sites.txt" if REF_TRANS else "",
    # 		strand=lambda wildcards: hisatStrandOption(wildcards.sample),
    # 		infiles=lambda wildcards: hisatInput(wildcards.sample)
    # 	log: ALIGN_DIR+"/hisat-{sample}-{run}.log"
    # 	threads: THREADS
    # 	message: "Aligning input with hisat (sample {wildcards.sample} - run {wildcards.run})"
    #     	shell: "{params.load} {params.load_samtools} {params.ss_gen} hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {params.trans} {params.strand} {params.extra} -x {params.indexdir} --dta {params.infiles} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"
    #


class HisatExtractSplices(AtomicOperation):

    def __init__(self, configuration, outdir):

        self._outdir = outdir
        self.input = {"ref_trans": configuration.get("reference", dict()).get("transcriptome", "")}
        self.__configuration = configuration
        self.output = {"splices":  os.path.join(outdir, "rnaseq", "2-alignments", "index", "splice_sites.txt")}
        pass

    @property
    def cmd(self):
        # ss_gen="hisat2_extract_splice_sites.py " + REF_TRANS + " > " + ALIGN_DIR + "/hisat/{sample}-{run}/
        # splice_sites.txt &&" if REF_TRANS else "",
        load = self.load
        input = self.input
        output = self.output

        cmd = "{load} hisat2_extract_splice_sites.py {input[ref_trans]} > {output[splices]}".format(**locals())
        return cmd


class HisatBuild(AtomicOperation):

    def __init__(self, configuration, outdir):

        super().__init__()
        # TODO: probably the input should be the cleaned up genome
        self.input = {"genome": configuration["reference"]["genome"]}
        self._outdir = os.path.join(outdir, "rnaseq", "2-alignments", "index", "hisat")
        self.output = {"flag": os.path.join(self._outdir, "hisat_index.done")}
        self.log = os.path.join(outdir, "rnaseq", "2-alignments", "index", "log", "hisat.log")
        self.touch = True
        self.__configuration = configuration

    @property
    def out_prefix(self):
        return os.path.join(self._outdir[0], "genome")

    @property
    def message(self):

        return "Building the HISAT2 index for {} in {}".format(self.input["genome"],
                                                               self._outdir)

    @property
    def cmd(self):

        load = None
        threads = self.__configuration["threads"]  # TODO: better way to get threads
        input = self.input["genome"]
        out_prefix = self.out_prefix
        log = self.log
        cmd = "{load} hisat2-build {extra} -p {threads} {input[genome]} {out_prefix} > {log} 2>&1".format(
            **locals()
        )
        return cmd
