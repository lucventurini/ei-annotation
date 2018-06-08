from .abstract import ShortAssembler, ShortAssemblerWrapper
from ... import AtomicOperation, EIWrapper, ShortSample
import os
import itertools
import subprocess
import functools


@functools.lru_cache(maxsize=8, typed=True)
def get_class_location(load):
    cmd = "{load} && which class_run.py".format(load=load)
    loc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    if not loc:
        raise OSError("CLASS2 not found!")
    return loc


class Class2Wrapper(ShortAssemblerWrapper):

    def __init__(self, configuration, bams):

        super().__init__()
        self.configuration = configuration
        outdir = ''  # TODO: implement!
        if len(self.runs) > 0 and len(bams) > 0:
            scallops = []
            for bam, run in itertools.product(bams, range(len(self.runs))):
                scallop = Class2(bam, run, configuration, outdir)
                scallops.append(scallop)
                self.add_to_gf(scallop)
                continue
            flag = Class2Flag(scallops, outdir)
            self.add_edges_from([(scallop, flag) for scallop in scallops])

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

    def __init__(self, bam, run, configuration, outdir):

        super().__init__(bam, run, configuration, outdir)

    @property
    def toolname(self):
        return "class2"

    @property
    def loader(self):
        return ["class2"]

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
        threads = self.threads
        outdir = self.gfdir
        program = self.program
        link_src = self.link_src
        input = self.input
        output = self.output
        log = self.log
        gtf = os.path.join(self.gfdir, "transcripts.gtf")
        cmd += "{program} --clean --force -c \"{extra}\" -p {threads} {input[bam]} > {output[gf]} 2> {log} &&"
        cmd += " ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"
