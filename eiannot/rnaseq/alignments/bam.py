from ...abstract import AtomicOperation
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

    # @property
    # def threads(self):
    #     return 1

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
    #
    # @property
    # def is_small(self):
    #     return True


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
