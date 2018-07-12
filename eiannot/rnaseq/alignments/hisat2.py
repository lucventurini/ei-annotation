from ...abstract import AtomicOperation, ShortSample
from .abstract import IndexBuilder, ShortAligner, ShortWrapper, IndexLinker
import os
import itertools
import glob


class HisatWrapper(ShortWrapper):

    __toolname__ = "hisat2"

    def __init__(self, configuration, prepare_flag):

        # First, we have to build the index

        super().__init__(configuration, prepare_flag)

        # Then we have to do all the alignments
        # Retrieve the running parameters

        if len(self.runs) > 0 and len(self.samples) > 0:
            # Start creating the parameters necessary for the run
            if self.prebuilt is True:
                indexer = HisatLinker(configuration, self.outdir)
            else:
                indexer = self.indexer(configuration, self.outdir)
            self.add_node(indexer)
            assert len(indexer.output) == 1, indexer.output
            # Optionally build the reference splice catalogue
            splices = HisatExtractSplices(configuration, self.outdir)
            splice_out = splices.output.get("splices", None)
            if splice_out:
                self.add_node(splices)
            for sample, run in itertools.product(self.samples, range(len(self.runs))):
                hisat_run = HisatAligner(indexer, sample, run)
                # hisat_runs.append(hisat_run)
                self.add_to_bams(hisat_run)
            self.add_edges_from([(indexer, run) for run in self.bams])
            if splice_out:
                self.add_edges_from([(splices, run) for run in self.bams])
            assert len(indexer.output) == 1

    @property
    def indexer(self):
        return HisatBuild


class HisatLinker(IndexLinker):

    __toolname__ = "hisat2"

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        self.configuration = configuration
        self.input["index_folder"] = self.index_folder

        index_files = []

        for fname in glob.glob(os.path.join(self.index_folder,
                                                           "{index_name}*".format(index_name=self.index_name))):
            if fname.endswith("ht2") or fname.endswith("ht2l"):
                index_files.append(fname)
        assert len(index_files) > 0
        self.input["index_files"] = index_files
        self.output = {"flag": os.path.join(self.outdir, "hisat_index.done")}
        self.touch = True

    @property
    def cmd(self):
        outdir = self.outdir
        cmd = "mkdir -p {outdir} ".format(**locals())
        for fname in self.input["index_files"]:
            link_src = os.path.relpath(os.path.abspath(fname), start=self.outdir)

            link_dest = os.path.join(self.outdir,
                                     self.species + '.' + '.'.join(fname.split('.')[1:]))
            cmd += " && ln -sf {link_src} {link_dest}".format(**locals())

        return cmd

    @property
    def index(self):
        return os.path.abspath(os.path.join(self.outdir, self.species))

    @property
    def loader(self):
        return []


class HisatAligner(ShortAligner):

    __toolname__ = "hisat2"

    def __init__(self, indexer, sample: ShortSample, run: int):

        super().__init__(indexer=indexer, sample=sample, run=run)
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
        log = self.log
        cmd += "{strand} {extra} -x {index} {infiles} 2> {log}"
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
        return ["hisat2", "samtools"]


class HisatExtractSplices(AtomicOperation):

    __name__ = "hisat_extract_splices"

    def __init__(self, configuration, outdir):

        super().__init__()
        self._outdir = outdir
        self.input = {"ref_trans": configuration.get("reference", dict()).get("transcriptome", "")}
        self.configuration = configuration
        if self.input["ref_trans"]:
            self.output = {"splices":  os.path.join(outdir, "index", "splice_sites.txt")}
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
        return ["hisat2"]

    @property
    def rulename(self):
        return self.__name__


class HisatBuild(IndexBuilder):

    __toolname__ = "hisat2"

    def __init__(self, configuration, outdir):

        super().__init__(configuration, outdir)
        # TODO: probably the input should be the cleaned up genome
        self.output = {"flag": os.path.join(self.outdir, "hisat_index.done")}
        self.touch = True

    @property
    def out_prefix(self):
        return os.path.abspath(os.path.join(self.outdir, self.species))

    @property
    def message(self):

        return "Building the HISAT2 index for {} in {}".format(self.input["genome"],
                                                               self.outdir)

    @property
    def cmd(self):

        load = self.load
        threads = self.threads
        input = self.input
        out_prefix = self.out_prefix
        log = self.log
        extra = self.extra
        cmd = "{load} hisat2-build {extra} -p {threads} {input[genome]} {out_prefix} > {log} 2>&1".format(
            **locals()
        )
        return cmd

    @property
    def loader(self):
        return ["hisat2"]

    @property
    def index(self):
        return self.out_prefix
