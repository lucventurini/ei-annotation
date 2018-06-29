from ...abstract import AtomicOperation, EIWrapper, ShortSample
from .workflow import ShortAlignmentsWrapper
import os
import re
import functools
import subprocess


@functools.lru_cache(maxsize=8, typed=True)
def portcullis_help(command, step):
    cmd = subprocess.Popen("{} portcullis {} --help".format(command, step),
                           shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    return cmd


@functools.lru_cache(maxsize=8, typed=True)
def portcullisStrandOption(strandedness, command, step):
    cmd = portcullis_help(command, step)
    if not any("strandedness" in _ for _ in cmd.split("\n")):
        return ""
    else:
        if strandedness == "fr-firststrand":
            return "--strandedness=firststrand"
        elif strandedness == "fr-secondstrand":
            return "--strandedness=secondstrand"
        else:
            return "--strandedness=unstranded"


class PortcullisWrapper(EIWrapper):

    # TODO: probably I need something here to start the runs
    def __init__(self, short_alignments: ShortAlignmentsWrapper):

        super().__init__(configuration=short_alignments.configuration)

        filters = []
        juncs = []
        if self.execute and short_alignments.bams:
            preps = []
            filters = []
            refprep = PortcullisPrepRef(self.configuration, outdir=self.outdir)
            refout = refprep.output["refbed"]
            if refout:
                self.add_node(refprep)
            else:
                refprep = None

            for bam in short_alignments.bams:
                prep = PortcullisPrep(bam, self.outdir)
                self.add_edge(short_alignments.exit, prep)
                preps.append(prep)

                junc = PortcullisJunc(prep)
                self.add_edge(prep, junc)
                juncs.append(junc)

                filt = PortcullisFilter(junc, portcullis_prep_ref=refprep)
                self.add_edge(junc, filt)
                if refprep is not None:
                    self.add_edge(refprep, filt)
                filters.append(filt)
            self.merger = PortcullisMerge(self.configuration, filters, self.outdir)
            self.add_edges_from([(filt, self.merger) for filt in filters])
            # self.failed_merger = PortcullisMergeFailed(juncs, self.merger)
            # self.add_edge(self.merger, self.failed_merger)
        self.add_final_flag(os.path.join(self.outdir, "portcullis.done"), rulename="portcullis_flag")

    @property
    def junctions(self):

        return self.merger.output["bed"]

    @property
    def junctions_task(self):
        return self.merger

    @property
    def toolname(self):
        return 'portcullis'

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "rnaseq", "3-portcullis")

    @property
    def execute(self):
        return self.configuration["programs"].get("portcullis", dict()).get("execute", True)

    # @property
    # def failed_junctions(self):
    #     return self.failed_merger.output["bed"]
    #
    # @property
    # def failed_junctions_task(self):
    #     return self.failed_merger


class PortcullisPrep(AtomicOperation):

    def __init__(self, bam, outdir):
        super().__init__()
        self.configuration = bam.configuration
        self.input = {"genome": self.genome, "bam": bam.input["bam"], "index": bam.input["index"]}
        self.sample = bam.sample
        self._outdir = outdir
        self.output = {"bai": os.path.join(
            self.sample_dir,
            "1-prep", "portcullis.sorted.alignments.bam.bai")}
        self.log = os.path.join(self._outdir, "log", "portcullis_{alrun}-prep.log".format(alrun=self.alrun))
        if not os.path.exists(os.path.dirname(self.log)):
            os.makedirs(os.path.dirname(self.log))
        self.message = "Using portcullis to prepare in folder {sample_dir}: {bam}".format(bam=self.input["bam"],
                                                                                          sample_dir=self.sample_dir)

    def get_sample(self):
        bam = os.path.basename(self.input["bam"])
        groups = re.search("(.*)-([^-]*)-([^-]*)", os.path.splitext(bam)[0]).groups()
        if groups is None:
            raise ValueError(bam)
        sample = "{groups[1]}".format(**locals())
        return self.configuration["short_reads"]["samples"][sample]

    @property
    def sample_dir(self):
        return os.path.join(self._outdir, self.alrun)

    @property
    def alrun(self):
        return os.path.splitext(os.path.basename(self.input["bam"]))[0]

    @property
    def rulename(self):
        return "portcullis-prep-{}".format(self.alrun)

    @property
    def loader(self):
        return ["portcullis"]

    @property
    def cmd(self):
        load = self.load
        cmd = "{load}"
        outdir = os.path.dirname(self.output["bai"])
        input = self.input
        log = self.log
        threads = self.threads
        cmd += "portcullis prep -o {outdir} -t {threads} {input[genome]} {input[bam]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class PortcullisPrepRef(AtomicOperation):

    def __init__(self, configuration, outdir):
        super().__init__()
        self.configuration = configuration
        self.input = {"reference": self.configuration["reference"].get("ref_transcriptome", "")}
        if not self.input["reference"]:
            self.output = {"refbed": ""}
            return
        self._outdir = outdir
        self.output = {"refbed": os.path.join(self._outdir, "ref_juncs.bed")}
        self.log = "Converting reference transcriptome to "

    @property
    def rulename(self):
        return "portcullis_prep_ref_transcriptome"

    @property
    def loader(self):
        return ["portcullis"]

    @property
    def cmd(self):

        load = self.load
        cmd = "{load}"
        input = self.input
        output = self.output
        cmd += "junctools convert -if gtf -of ebed {input[reference]} > {output[refbed]}"
        cmd = cmd.format(**locals())
        return cmd


class PortcullisJunc(AtomicOperation):

    def __init__(self, portcullis_prep: PortcullisPrep):

        super().__init__()
        self.configuration, self.sample = portcullis_prep.configuration, portcullis_prep.sample
        self._outdir = portcullis_prep._outdir
        self.input = {"bai": portcullis_prep.output["bai"]}
        self.alrun = portcullis_prep.alrun
        self.output = {"tab": os.path.join(self.sample_dir,
                                           "2-junc",
                                           "portcullis.junctions.tab"
                                           ),
                       "bed": os.path.join(self.sample_dir,
                                           "2-junc",
                                           "portcullis.junctions.bed"
                                           )}
        self.log = os.path.join(self._outdir, "log", "portcullis_{alrun}-junc.log".format(alrun=self.alrun))
        self.message = "Using portcullis to analyse potential junctions: {alrun}".format(
            alrun=self.alrun
        )

    @property
    def rulename(self):
        return "portcullis_junc_{}".format(self.alrun)

    @property
    def loader(self):
        return ["portcullis"]

    @property
    def strand(self):
        return portcullisStrandOption(self.sample.strandedness, self.load, "junc")

    @property
    def cmd(self):
        strand = self.strand
        prepdir = os.path.dirname(self.input["bai"])
        outdir = os.path.dirname(self.output["tab"])
        log = self.log
        load = self.load
        threads = self.threads
        alrun = self.alrun
        cmd = "{load}"
        cmd += " portcullis junc -o {outdir}/portcullis {strand} -t {threads} {prepdir} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def sample_dir(self):
        return os.path.join(self._outdir, self.alrun)


class PortcullisFilter(AtomicOperation):

    def __init__(self, portcullis_junc: PortcullisJunc, portcullis_prep_ref=None):

        super().__init__()

        self.sample, self.configuration = portcullis_junc.sample, portcullis_junc.configuration
        self.input = portcullis_junc.output
        self.alrun = portcullis_junc.alrun
        self.sample = portcullis_junc.sample
        self.sample_dir = portcullis_junc.sample_dir
        self._outdir = portcullis_junc._outdir
        self.output = {"bed_link": os.path.join(self._outdir,
                                                "output",
                                                "portcullis_{alrun}.pass.junctions.bed".format(alrun=self.alrun)),
                       "tab_link": os.path.join(self._outdir,
                                                "output",
                                                "portcullis_{alrun}.pass.junctions.tab".format(alrun=self.alrun))}
        self.message = "Using portcullis to filter invalid junctions: {alrun}".format(alrun=self.alrun)
        self.log = os.path.join(self._outdir, "log", "portcullis_{alrun}-filter.log".format(alrun=self.alrun))
        if portcullis_prep_ref is not None:
            assert isinstance(portcullis_prep_ref, PortcullisPrepRef)
            self.input["refbed"] = portcullis_prep_ref.output["refbed"]

    @property
    def rulename(self):
        return "portcullis_filter_{}".format(self.alrun)

    @property
    def loader(self):
        return ["portcullis"]

    @property
    def canonical_option(self):
        return self.configuration["programs"]["portcullis"].get("canonical_juncs", "OFF")

    @property
    def min_intron(self):
        return max(self.configuration["reference"]["min_intron"], 20)

    @property
    def max_intron(self):
        return self.configuration["reference"]["max_intron"]

    @property
    def cmd(self):

        load = self.load
        cmd = "{load}"
        outdir = os.path.join(self.sample_dir, "3-filt", "portcullis")
        threads = self.threads
        canonical = self.canonical_option
        cmd += " (portcullis filter -o {outdir} --canonical={canonical} "
        max_intron = self.max_intron
        if "refbed" in self.input:
            trans = "--reference {input[refbed]}"
        else:
            trans = ""
        prepdir = os.path.join(self.sample_dir, "1-prep").format(alrun=self.alrun)
        cmd += " --max_length={max_intron} {trans} --threads={threads} {prepdir} "
        input = self.input
        log = self.log
        cmd += " {input[tab]} > {log} 2>&1"

        # Now time for the links
        output = self.output
        bed_link_src = os.path.relpath(os.path.join(self.sample_dir, "3-filt", "portcullis.pass.junctions.bed"),
                                       start=os.path.dirname(self.output["bed_link"]))
        tab_link_src = os.path.relpath(os.path.join(self.sample_dir, "3-filt", "portcullis.pass.junctions.bed"),
                                       start=os.path.dirname(self.output["tab_link"]))

        # We need the unfiltered files because, if there are too few reads, portcullis filter will fail
        # as this is not recoverable (it's by design) we provide the unfiltered junctions instead
        bed_link_unfilt = os.path.relpath(self.input["bed"],
                                          start=os.path.dirname(self.output["bed_link"]))
        #
        # bed_link_unfilt = os.path.join("..", "portcullis_{alrun}", "2-junc", "{alrun}.junctions.bed").format(
        #     alrun=self.alrun)
        tab_link_unfilt = os.path.relpath(self.input["tab"],
                                          start=os.path.dirname(self.output["tab_link"]))

        cmd += " && ln -sf {bed_link_src} {output[bed_link]} && ln -sf {tab_link_src} {output[tab_link]} )"
        cmd += " || ln -sf {bed_link_unfilt} {output[bed_link]} && ln -sf {tab_link_unfilt} {output[tab_link]}"
        cmd += " && touch -h {output[tab_link]} && touch -h {output[bed_link]}"
        cmd = cmd.format(**locals())
        return cmd


class PortcullisMerge(AtomicOperation):

    def __init__(self, configuration, filters: [PortcullisFilter], outdir):
        super().__init__()
        assert hasattr(filters, "__iter__") and all(isinstance(_, PortcullisFilter) for _ in filters)
        # assert len(filters) > 0
        self.outdir = outdir
        if filters:
            self.configuration = filters[0].configuration
        else:
            self.configuration = configuration

        self.input = {"beds": [filter.output["bed_link"] for filter in filters],
                      "tabs": [filter.output["tab_link"] for filter in filters]}
        self.output = {"bed": os.path.join(self.outdir, "output", "portcullis.merged.bed"),
                       "tab": os.path.join(self.outdir, "output", "portcullis.merged.tab"),
                       "gff3": os.path.join(self.outdir, "output", "portcullis.merged.gff3"),}

        self.message = "Taking the union of portcullis results"
        self.log = os.path.join(self.outdir, "logs", "portcullis.merge.log")

    @property
    def loader(self):
        return ["portcullis"]

    @property
    def rulename(self):
        return "portcullis_merge"

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):

        beds = " ".join(self.input["beds"])
        tabs = " ".join(self.input["tabs"])
        output = self.output

        if len(self.input["beds"]) == 0:
            cmd = "touch {output[bed]} && touch {output[tab]}".format(**locals())
        elif len(self.input["beds"]) == 1:
            cmd = "cat {beds} > {output[bed]} && cat {tabs} > {output[tab]}".format(**locals())
        else:
            load = self.load
            prefix = "--prefix=portcullis_merged"
            log = self.log
            cmd = "{load}"
            cmd += "(junctools set {prefix} --output={output[tab]} --operator=mean union {tabs} > {log} 2>&1"
            cmd += " || touch {output[tab]}) && "
            cmd += "junctools convert -if portcullis -of ebed -o {output[bed]}  {output[tab]} && "
            cmd += "junctools convert -if portcullis -of igff -o {output[gff3]}  {output[tab]}"
        cmd = cmd.format(**locals())
        return cmd


# class PortcullisMergeFailed(AtomicOperation):
#
#     def __init__(self, juncs: [PortcullisJunc], merged: PortcullisMerge):
#
#         super().__init__()
#         self.configuration = merged.configuration
#         self.input["merged"] = merged.output["tab"]
#         self.input["juncs"] = [junc.output["tab"] for junc in juncs]
#         self.outdir = merged.outdir
#         self.output["tab"] = os.path.join(self.outdir, "portcullis.failed.tab"),
#         self.output["bed"] = os.path.join(self.outdir, "portcullis.failed.bed")
#         self.output["gff3"] = os.path.join(self.outdir, "portcullis.failed.gff3")
#         self.output["all"] = os.path.join(self.outdir, "portcullis.all.bed")
#         self.temps = ["all"]
#
#     @property
#     def loader(self):
#         return ["portcullis"]
#
#     @property
#     def rulename(self):
#         return "portcullis_extract_failed"
#
#     @property
#     def cmd(self):
#
#         load = self.load
#         input, output = self.input, self.output
#         tab_inputs = " ".join(self.input["juncs"])
#         prefix = "--prefix=portcullis_failed"
#
#         cmd = "{load} junctools set union --output {output[all]} {tab_inputs} && "
#         cmd += " junctools set subtract {prefix} --output={output[tab]} {output[all]} {input[merged]} && "
#         cmd += " junctools convert -if portcullis -of ebed -o {output[bed]} {output[tab]} && "
#         cmd += "junctools convert -if portcullis -of igff -o {output[gff3]}  {output[tab]}"
#         cmd = cmd.format(**locals())
#         return cmd
#
