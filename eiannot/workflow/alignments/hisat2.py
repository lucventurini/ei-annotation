from .. import AtomicOperation, EIWrapper, ShortSample
from . import IndexBuilder, ShortAligner
import os
import itertools


class HisatWrapper(EIWrapper):

    def __init__(self, configuration, outdir):

        # First, we have to build the index

        super().__init__()

        # Then we have to do all the alignments
        # Retrieve the running parameters
        runs = configuration["programs"]["hisat"]["runs"]
        samples = configuration["short_reads"]["samples"]

        if len(runs) > 0:
            # Start creating the parameters necessary for the run
            indexer = HisatBuild(configuration, outdir)
            self.add_node(indexer)
            # Optionally build the reference splice catalogue
            splices = HisatExtractSplices(configuration, outdir)
            splice_out = splices.output["splices"]
            if splice_out:
                self.add_node(splices)
            hisat_runs = []
            for sample, run in itertools.product(samples, range(len(runs))):
                hisat_run = HisatAligner(configuration=configuration,
                                         index=indexer.out_prefix,
                                         sample=sample,
                                         outdir = outdir,
                                         run=run,
                                         ref_transcriptome=splice_out)
                hisat_runs.append(hisat_run)
            self.add_edges_from([(indexer, run) for run in hisat_runs])
            if splice_out:
                self.add_edges_from([(splices, run) for run in hisat_runs])
            flag = HisatFlag(outdir, [run.output["link"] for run in hisat_runs])
            self.add_edges_from([(run, flag) for run in hisat_runs])


class HisatFlag(AtomicOperation):

    # rule hisat_all:
    # 	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
    # 	output: ALIGN_DIR+"/hisat.done"
    # 	shell: "touch {output}"

    def __init__(self, outdir, runs=[]):

        super().__init__()
        for number, run in enumerate(runs):
            self.input["run{}".format(number)] = run
        self.output["flag"] = os.path.join(outdir, "hisat.done")
        self.touch = True


class HisatLink(AtomicOperation):

    @property
    def cmd(self):
        return "ln -s {input} {output}".format(input=self.input[0],
                                               output=self.output[0])

    @property
    def loader(self):
        return []


class HisatAligner(ShortAligner):

    def __init__(self, configuration, index, outdir, sample: ShortSample, run: int, ref_transcriptome=None):

        super().__init__(configuration=configuration, sample=sample,
                         run=run, index=index, ref_transcriptome=ref_transcriptome,
                         outdir=outdir)
        self.output = {"bam": os.path.join(self.bamdir,
                                           "hisat.bam"),
                       "link": self.link}

    @property
    def cmd(self):

        load = self.load
        cmd = "{load}"
        cmd += "hisat2 -p {threads} "
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += "--min-intronlen={min_intron} --max-intronlen={max_intron} "
        if self.input.get("transcriptome", None):
            cmd += "--known-splicesite-infile={}".format(self.input.get("transcriptome"))
        extra = self.extra
        strand = self.strand
        infiles = self.input_reads
        threads = self.threads
        output = self.output
        index = self.index
        cmd += "{strand} {extra} -x {index} {infiles} > {log}"
        cmd += "| samtools view -b -@ {threads} - > {output[bam]}"
        link_src = self.link_src
        cmd += "&& ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def strand(self):
        if self.sample.strandedness == "fr-firststrand":
            return "--rna-strandness=RF"
        elif self.sample.strandedness == "fr-secondstrand":
            return "--rna-strandness=FR"
        elif self.sample.strandedness == "f":
            return "--rna-strandness=F"
        elif self.sample.strandedness == "r":
            return "--rna-strandness=R"
        else:
            return ""

    @property
    def input_reads(self):

        read1, read2 = self.input["read1"], self.input.get("read2", None)

        if read2:
            return "-1 {read1} -2 {read2}".format(**locals())
        else:
            return "-U {read1}".format(**locals())

    @property
    def loader(self):
        return ["hisat", "samtools"]

    @property
    def toolname(self):
        return "hisat"


class HisatExtractSplices(AtomicOperation):

    __name__ = "hisat_extract_splices"

    def __init__(self, configuration, outdir):

        super().__init__()
        self._outdir = outdir
        self.input = {"ref_trans": configuration.get("reference", dict()).get("transcriptome", "")}
        self.__configuration = configuration
        if self.input["ref_trans"]:
            self.output = {"splices":  os.path.join(outdir, "rnaseq", "2-alignments", "index", "splice_sites.txt")}
        else:
            self.output = dict()

    @property
    def cmd(self):
        # ss_gen="hisat2_extract_splice_sites.py " + REF_TRANS + " > " + ALIGN_DIR + "/hisat/{sample}-{run}/
        # splice_sites.txt &&" if REF_TRANS else "",
        load = self.load
        input = self.input
        output = self.output

        cmd = "{load} hisat2_extract_splice_sites.py {input[ref_trans]} > {output[splices]}".format(**locals())
        return cmd

    @property
    def loader(self):
        return ["hisat"]

    @property
    def rulename(self):
        return self.__name__


class HisatBuild(IndexBuilder):

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        # TODO: probably the input should be the cleaned up genome
        self.output = {"flag": os.path.join(self._outdir, "hisat_index.done")}
        self.touch = True

    @property
    def toolname(self):
        return "hisat"

    @property
    def out_prefix(self):
        return os.path.abspath(os.path.join(self._outdir[0], "genome"))

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

    @property
    def loader(self):
        return ["hisat"]