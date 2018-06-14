from .bam import BamStats, BamIndex, BamSort, AlnFlag
from .hisat2 import HisatWrapper
from .gmap import GsnapWrapper
from .tophat2 import TopHat2Wrapper
from .portcullis import PortcullisWrapper
from .star import StarWrapper, StarLongWrapper
from ...abstract import EIWrapper, EIWorfkflow, AtomicOperation
import os


class ShortAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarWrapper,
                "hisat": HisatWrapper,
                "tophat2": TopHat2Wrapper,
                "gsnap": GsnapWrapper}

    def __init__(self, prepare_wrapper: EIWrapper):

        super().__init__()
        stats = []
        bams = []
        self.configuration = prepare_wrapper.configuration
        prepare_flag = prepare_wrapper.output
        for wrapper in self.wrappers.values():
            instance = wrapper(self.configuration, bams, prepare_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            bams.extend(instance.bams)
            for bam in instance.bams:
                # TODO: here add the entrance from the "preparation" stage
                # self.add_edges_from()
                sorter = BamSort(bam)
                self.add_edge(bam, sorter)
                indexer = BamIndex(sorter)
                self.add_edge(sorter, indexer)
                stater = BamStats(indexer)
                self.add_edge(indexer, stater)
                stats.append(stater)

        final_flag = AlnFlag(stats)
        self.add_edges_from([(stat, final_flag) for stat in stats])
        if self.execute_portcullis is True:
            self.portcullis = PortcullisWrapper(self.configuration, bams)
            self.merge(self.portcullis)
            flag = self.portcullis.flag
            self.add_edges_from([stat, flag] for stat in stats)

    @property
    def execute_portcullis(self):

        # TODO: probably we should push this somewhere else
        return self.configuration["portcullis"]["execute"]

    @property
    def portcullis_junctions(self):
        return self.portcullis.output["bed"]


class LongAlignmentsWrapper(EIWrapper):

    wrappers = {"star": StarLongWrapper}

    def __init__(self, prepare_wrapper):

        super().__init__()
        stats = []
        gfs = []
        self.configuration = prepare_wrapper.configuration
        prepare_flag = prepare_wrapper.output

        for wrapper in self.wrappers.values():
            instance = wrapper(self.configuration, prepare_flag)
            instance.add_flag_to_inputs()
            self.merge([instance])
            gfs.extend(instance.gfs)
            for gf in gfs:
                stat = LongAlignerStats(gf)
                self.add_edge(gf, stat)
                stats.append(stat)

        flag = LongAlignersFlag(stats)
        self.add_edges_from([(stat, flag) for stat in stats])
        pass


class LongAlignerStats(AtomicOperation):

    def __init__(self, aligner):
        super().__init__()
        self.input = aligner.output
        self.outdir = aligner.outdir
        self.output = {"stats": os.path.splitext(self.input["gf"])[0] + ".stats"}
        self.message = "Calculating statistics for: {input[gf]}".format(input=self.input)
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["gf"])

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[gf]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class LongAlignersFlag(AtomicOperation):

    def __init__(self, stats_runs: [LongAlignerStats]):

        super().__init__()
        self.touch = True
        outdir = os.path.dirname(os.path.dirname(stats_runs[0].output["stats"]))
        self.output["flag"] = os.path.join(outdir, "long_reads.done")

    @property
    def rulename(self):
        return "aln_all"

    @property
    def loader(self):
        return []


    # rule lr_gmap:
    #     input:
    #         index = rules.gmap_index.output,
    #         reads = lambda wildcards: L_INPUT_MAP[wildcards.lsample]
    # output: link = ALIGN_DIR + "/lr_output/lr_gmap-{lsample}-{lrun}.gff",
    # gff = ALIGN_DIR + "/gmap/{lsample}-{lrun}/lr_gmap-{lsample}-{lrun}.gff"
    # params:
    #     load = loadPre(config, "gmap"),
    #     link_src = "../gmap/{lsample}-{lrun}/lr_gmap-{lsample}-{lrun}.gff",
    #     intron_length = gmap_intron_lengths(loadPre(config, "gmap"), MAX_INTRON)
    # log: ALIGN_DIR + "/gmap-{lsample}-{lrun}.log"
    # threads: THREADS
    # message: "Mapping long reads to the genome with gmap (sample: {wildcards.lsample} - run: {wildcards.lrun})"
    # shell: "{params.load} gmap --dir={ALIGN_DIR}/gmap/index --db={NAME} --min-intronlength={MIN_INTRON} {params.intron_length} --format=3 {input.reads} > {output.gff} 2> {log} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


