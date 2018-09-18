from ...abstract import AtomicOperation
from ...preparation import FaidxGenome
import os


class BamSort(AtomicOperation):

    __name__ = "bam_sort"

    def __init__(self, bamrule):

        super().__init__()
        self.alrun = bamrule.alrun
        self.sample = bamrule.sample
        self.configuration = bamrule.configuration
        self.input = {"bam": bamrule.output["link"]}
        self.output = {"bam": os.path.splitext(self.input["bam"])[0] + ".sorted.bam"}
        output = self.output
        input = self.input
        load = self.load
        temp = os.path.join(os.path.dirname(self.input["bam"]), "..", "sort_{}".format(self.align_run))
        threads = self.threads
        self.message = "Using samtools to sort {input[bam]}".format(**locals())
        self.cmd = "{load} samtools sort -o {output[bam]} -O bam -m 1G -T {temp} -@ {threads} {input[bam]}".format(
            **locals()
        )

    @property
    def align_run(self):
        return os.path.splitext(os.path.basename(self.input["bam"]))[0]

    @property
    def rulename(self):
        return "{name}-{align_run}".format(name=self.__name__, align_run=self.align_run)

    @property
    def loader(self):
        return ["samtools"]


class BamIndex(AtomicOperation):

    __name__ = "bam_index"

    def __init__(self, bamrule):
        super().__init__()
        self.alrun = bamrule.alrun
        self.sample = bamrule.sample
        self.configuration = bamrule.configuration
        self.input = {"bam": bamrule.output["bam"]}
        self.output = {"index": self.input["bam"] + ".csi"}
        self.message = "Using samtools to index: {input[bam]}".format(input=self.input)
        self.log = os.path.join(self.input["bam"] + ".index.log")

    @property
    def loader(self):
        return ["samtools"]

    @property
    def align_run(self):
        return os.path.splitext(os.path.basename(self.input["bam"]))[0]

    @property
    def rulename(self):
        return "{name}-{align_run}".format(name=self.__name__, align_run=self.align_run)

    @property
    def cmd(self):
        load = self.load
        input = self.input
        log = self.log
        threads = self.threads
        cmd = "{load} samtools index -c -@ {threads} {input[bam]} >{log} 2>&1".format(**locals())
        return cmd


class BamStats(AtomicOperation):

    __name__ = "bam_stats"

    def __init__(self, bamrule: BamIndex):

        super().__init__()
        self.sample = bamrule.sample
        self.configuration = bamrule.configuration
        self.alrun = bamrule.alrun
        self.input = {"bam": bamrule.input["bam"], "index": bamrule.output["index"]}
        self.output = {"stats": bamrule.input["bam"] + ".stats"}
        input, output = self.input, self.output
        plot_dir = os.path.join(os.path.dirname(bamrule.input["bam"]), "plots", self.align_run, self.align_run)
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        self.message = "Using samtools to collect stats for: {input[bam]}".format(input=self.input)
        load = self.load
        self.cmd = "{load} samtools stats {input[bam]} > {output[stats]}"
        self.cmd += " && (plot-bamstats -p {plot_dir} {output[stats]} || touch {output[stats]})"
        self.cmd = self.cmd.format(**locals())

    @property
    def align_run(self):
        return os.path.splitext(os.path.basename(self.input["bam"]))[0]

    @property
    def loader(self):
        return ["samtools"]

    @property
    def rulename(self):
        return "{name}-{align_run}".format(name=self.__name__, align_run=self.align_run)

    @property
    def threads(self):
        return 1


class Bam2BigWig(AtomicOperation):

    __name__ = "bam_to_bigwig"

    """This class has the purpose of converting BAM RNASeq files into WIG files for Augustus.
    We will use deepTools to do the conversion."""

    def __init__(self, bamrule: BamIndex, strand=None):
        super().__init__()
        self.configuration = bamrule.configuration
        self.input = {"bam": bamrule.input["bam"], "index": bamrule.input["index"]}
        self.__strand = None
        self.strand = strand
        if self.strand:
            strand = "." + self.strand
        else:
            strand = ""

        self.output = {"bw": os.path.splitext(self.input["bam"])[0] + strand + ".bw"}
        self.sample = bamrule.sample

    @property
    def align_run(self):
        return os.path.splitext(os.path.basename(self.input["bam"]))[0]

    @property
    def rulename(self):
        return "{name}-{align_run}-{strand}".format(name=self.__name__,
                                                    align_run=self.align_run,
                                                    strand=self.strand)

    @property
    def loader(self):
        return ["deepTools"]

    @property
    def strand(self):
        return self.__strand

    @strand.setter
    def strand(self, strand):
        if strand in ("", None):
            self.__strand = None
        elif strand in ("forward", "reverse"):
            self.__strand = strand
        elif strand == "+":
            self.__strand = "forward"
        elif strand == "-":
            self.__strand = "reverse"
        else:
            raise ValueError("Invalid value for strand: {}".format(strand))

    @property
    def rna_sense(self):
        if not (self.strand and self.sample.stranded):
            return ""
        pref = " --filterRNAstrand "

        if self.sample.strandedness in ("fr-firststrand", "r"):
            # Normal one. Forward is forward, reverse is reverse
            return pref + self.strand
        elif self.sample.strandedness in ("f", "fr-secondstrand"):
            if self.strand == "forward":
                return pref + "reverse"
            else:
                return pref + "forward"

    @property
    def out_suffix(self):
        if not (self.strand and self.sample.stranded):
            return ""

        if self.sample.strandedness in ("fr-firststrand", "r"):
            return self.strand
        elif self.sample.strandedness in ("f", "fr-secondstrand"):
            if self.strand == "forward":
                return "reverse"
            else:
                return "forward"

    @property
    def cmd(self):

        load = self.load
        bam = self.input["bam"]
        rna_sense = self.rna_sense
        out = self.output["bw"]
        threads = self.threads
        cmd = "{load} "
        cmd += "bamCoverage -b {bam} --ignoreDuplicates {rna_sense} -o {out} -p {threads} --normalizeUsing BPM"
        cmd = cmd.format(**locals())
        return cmd


class MergeWigs(AtomicOperation):

    def __init__(self, bigwigs: [Bam2BigWig], faidx: FaidxGenome, strand=None):

        super().__init__()
        self.__strand = None
        self.strand = strand
        self.configuration = bigwigs[0].configuration
        self.input = {"bigwigs": [rule.output["bw"] for rule in bigwigs],
                      "fai": faidx.output["fai"] }
        name_prefix = "coverage"
        if self.strand is not None:
            name_prefix += "." + self.strand

        self.output = {"bigwig": os.path.join(self.outdir, name_prefix+ ".bw"),
                       "bedgraph": os.path.join(self.outdir, name_prefix + ".bg"),
                       "sorted_bedgraph": os.path.join(self.outdir, name_prefix + ".sorted.bg"),
                       "bw_list": os.path.join(self.outdir, name_prefix + ".txt"),
                       }
        self.add_to_params("chrom_list", os.path.join(self.outdir, "chrom_list.txt"))
        self.temps = ["bedgraph", "bw_list", "sorted_bedgraph"]

    @property
    def outdir(self):
        return os.path.join(os.path.join(self.configuration["outdir"], "rnaseq", "1-alignments", "output"))

    @property
    def is_small(self):
        return True

    @property
    def strand(self):
        return self.__strand

    @strand.setter
    def strand(self, strand):
        if strand not in (None, "forward", "reverse"):
            raise ValueError("Invalid strand: {}".format(strand))
        self.__strand = strand

    @property
    def rulename(self):
        name = "merge_wigs"
        if self.strand is not None:
            name += "_" + self.strand
        return name

    @property
    def max_value(self):
        return self.configuration.get("abinitio", {}).get("max_coverage", 50)

    @property
    def loader(self):
        return ["kent"]

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        bigwigs = "\\n".join(self.input["bigwigs"])
        bw_list = self.output["bw_list"]
        bedgraph = self.output["bedgraph"]
        out = self.output["bigwig"]
        fai = self.input["fai"]
        max_value = self.max_value
        chrom_list = self.get_from_params("chrom_list")
        cmd += """echo  "{bigwigs}" >  {bw_list} && bigWigMerge -inList -threshold=0 -clip={max_value} """
        cmd += """{bw_list} {bedgraph} && cut -f 1,2 {fai} > {chrom_list} && """
        sorted_bedgraph = self.output["sorted_bedgraph"]
        cmd += " sort -k1,1 -k2,2n {bedgraph} > {sorted_bedgraph} &&"
        cmd += " bedGraphToBigWig {sorted_bedgraph} {chrom_list} {out} && rm {chrom_list}"

        cmd = cmd.format(**locals())
        return cmd
