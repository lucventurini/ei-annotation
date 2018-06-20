from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats
from ...abstract import AtomicOperation
import os
import itertools
import subprocess
import functools


@functools.lru_cache(maxsize=8, typed=True)
def get_class_location(load):
    cmd = "{load} which class_run.py".format(load=load)
    loc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    if not loc:
        raise OSError("CLASS2 not found!")
    return loc.rstrip()


class Class2Wrapper(ShortAssemblerWrapper):

    def __init__(self, aln_wrapper):

        super().__init__(aln_wrapper)
        outdir = ''  # TODO: implement!
        if len(self.runs) > 0 and len(self.bams) > 0:
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                class2 = Class2(bam, run)
                stat = AsmStats(class2)
                self.add_edge(class2, stat)
                self.add_to_gf(stat)
                continue

    @property
    def toolname(self):
        return "class2"


class Class2Flag(AtomicOperation):

    def __init__(self, class2s, outdir):
        super().__init__()
        self.input = {"gtfs": [class2.output["link"] for class2 in class2s]}
        self.touch = True
        self.output = {"flag": os.path.join(outdir, "class2.done")}


class Class2(ShortAssembler):

    def __init__(self, bam, run):

        super().__init__(bam, run)

    @property
    def toolname(self):
        return "class2"

    @property
    def loader(self):
        return ["class2", "mikado"]

    @property
    def strand(self):
        return ""

    @property
    def program(self):
        return get_class_location(self.load)

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        extra = self.extra
        if self.extra:
            extra = "-c \"{}\"".format(self.extra)
        else:
            extra = ""
        threads = self.threads
        outdir = self.gfdir
        program = self.program
        link_src = self.link_src
        input = self.input
        output = self.output
        log = self.log
        gtf = os.path.join(self.gfdir, "transcripts.gtf")
        cmd += "{program} --clean --force {extra} -p {threads} {input[bam]} > {output[gf]} 2> {log} &&"
        cmd += " ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def input_reads(self):
        return ""
