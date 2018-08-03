from .abstract import ShortAssembler, ShortAssemblerWrapper, AsmStats, FilterGF
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

    __toolname__ = "class2"
    __tag__ = "Abundance"

    def __init__(self, aln_wrapper):

        super().__init__(aln_wrapper)
        if len(self.runs) > 0 and len(self.bams) > 0:
            for bam, run in itertools.product(self.bams, range(len(self.runs))):
                class2 = Class2(bam, run, create_link=False)
                class_filter = FilterGF(class2,
                                        tag=self.tag,
                                        monoexonic_threshold=self.mono_threshold,
                                        multiexonic_threshold=self.multi_threshold)
                self.add_edge(class2, class_filter)
                stat = AsmStats(class_filter)
                self.add_edge(class_filter, stat)
                self.add_to_gf(stat)
                continue


class Class2(ShortAssembler):

    def __init__(self, bam, run, create_link=True):

        super().__init__(bam, run, create_link=create_link)

    __toolname__ = "class2"

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
        cmd += "{program} --clean --force {extra} -p {threads} {input[bam]} > {output[gf]} 2> {log} "
        if self._create_link:
            cmd += " ln -sf {link_src} {output[link]} && touch -h {output[link]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def suffix(self):
        return "gtf"

    @property
    def input_reads(self):
        return ""
