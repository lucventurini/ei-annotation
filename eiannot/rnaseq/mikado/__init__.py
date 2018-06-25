from ...abstract import AtomicOperation, EIWrapper
import os
import networkx as nx
# import abc
from ..assemblies.workflow import AssemblyWrapper
from ..alignments.workflow import LongAlignmentsWrapper
from ..alignments.portcullis import PortcullisWrapper
from .prepare import MikadoPrepare, MikadoConfig
from .orfs import Prodigal, TransdecoderLongOrf, TransdecoderPred, OrfCaller
from .homology import MikadoHomologyWrapper
from ...preparation import FaidxGenome


__modes__ = ("permissive", "stringent", "nosplit", "split", "lenient")


class MikadoSerialise(AtomicOperation):

    def __init__(self,
                 prepare: MikadoPrepare,
                 homology: MikadoHomologyWrapper,
                 orfs: OrfCaller,
                 faidx: FaidxGenome,
                 portcullis: PortcullisWrapper):

        super().__init__()
        self.input = prepare.output
        self.input.update(orfs.output)
        self.input.update(homology.output)
        self.input.update(faidx.output)
        self.input["cfg"] = prepare.config
        self.configuration = prepare.configuration
        self.outdir = prepare.outdir
        self.output = {"db": os.path.join(self.outdir, "mikado.db")}
        self.log = os.path.join(self.outdir, "mikado_serialise.err")
        self.homology = homology
        self.message = "Running Mikado serialise to move numerous data sources into a single database"
        self.portcullis = portcullis

    @property
    def rulename(self):
        return "mikado_serialise"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def blast_xmls(self):
        pref = "--xml="
        xmls = self.homology.blast_xmls
        if len(xmls) == 0:
            return ""
        elif len(set([os.path.dirname(xml) for xml in xmls])) == 1:
            return "{pref}{fold}".format(fold=os.path.dirname(xmls[0]), **locals())
        else:
            return "{pref}{commas}".format(",".join(xmls), **locals())

    @property
    def blast_targets(self):
        if self.homology.blast_targets:
            return "--blast_targets={}".format(self.homology.blast_targets)
        else:
            return ""

    @property
    def orfs(self):
        return " --orfs={input[orfs]}".format(input=self.input)

    @property
    def junctions(self):
        if self.portcullis:
            junctions = self.portcullis.junctions
            return "--junctions={junctions}".format(**locals())
        else:
            return ""

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        blast_xmls = self.blast_xmls
        blast_targets = self.blast_targets
        junctions = self.junctions
        input = self.input
        orfs = self.orfs
        cmd += "mikado serialise {blast_xmls} {blast_targets} {junctions} {orfs} --transcripts={input[fa]} "
        cmd += "--genome_fai={input[fai]} --json-conf={input[cfg]} --force --start-method=spawn "
        outdir = self.outdir
        threads = self.threads
        log = self.log
        cmd += " -od {outdir} --procs={threads} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class MikadoPick(AtomicOperation):

    def __init__(self, serialise: MikadoSerialise):

        super().__init__()
        self.configuration = serialise.configuration
        self.input = serialise.output
        self.input["gtf"] = serialise.input["gtf"]
        self.input["cfg"] = serialise.input["cfg"]
        self.__serialise_dir = serialise.outdir
        self.output = {"loci": os.path.join(
            self.loci_dir, "mikado-{mode}.loci.gff3").format(**locals()),
                       "link": os.path.join(self.outdir, "mikado.loci.gff3")}
        self.log = os.path.join(self.outdir, "mikado-{mode}.pick.err").format(
            **locals()
        )
        self.message = "Running mikado picking stage in mode: {mode}".format(**locals())

    @property
    def mode(self):
        return self.configuration["mikado"]["pick"]["mode"]

    @mode.setter
    def mode(self, mode):
        assert mode in __modes__
        self.__mode = mode

    @property
    def loader(self):
        return ["mikado"]

    @property
    def rulename(self):
        return "mikado_pick_{mode}".format(mode=self.mode)

    @property
    def cmd(self):
        load = self.load

        cmd = "{load} "
        mode = self.mode
        threads = self.threads
        cmd += "mikado pick --source Mikado_{mode} --mode={mode} --procs={threads} "
        input = self.input
        cmd += "--start-method=spawn --json-conf={input[cfg]} "
        loci_out = os.path.basename(self.output["loci"])
        outdir = self.outdir
        log = self.log
        cmd += "-od {outdir} --loci_out {loci_out}  -lv INFO -db {input[db]} {input[gtf]} > {log} 2>&1 "
        cmd += " && cd {link_dir} && ln -s {link_src} {output[link]}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def loci_dir(self):
        return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado", "pick")

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado", "output")


class IndexMikado(AtomicOperation):

    def __init__(self, mikado: MikadoPick):
        super().__init__()
        self.input = mikado.output
        self.output = {"midx": self.input["link"] + ".midx"}
        self.log = os.path.join(os.path.dirname(self.input["link"]), "index_loci.log")

    @property
    def rulename(self):
        return "index_mikado_loci"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load
        input, log = self.input, self.log

        cmd = "{load} mikado compare -r {input[loci]} -l {log} --index"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.dirname(self.input["loci"])


class MikadoStats(AtomicOperation):

    def __init__(self, index: IndexMikado):

        super().__init__()
        self.input = index.output
        self.input.update(index.input)
        self.outdir = index.outdir
        self.configuration = index.configuration
        self.output = {"stats": os.path.splitext(self.input["link"])[0] + ".stats"}
        self.message = "Calculating statistics for Mikado run"
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return "mikado_stats"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["link"])

    @property
    def cmd(self):

        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[link]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


# class MikadoCollectStats(AtomicOperation):
#
#     def __init__(self, stats: [MikadoStats]):
#
#         assert len(stats) > 0 and all(isinstance(stat, MikadoStats) for stat in stats)
#         super().__init__()
#         self.input["stats"] = [stat.output["stats"] for stat in stats]
#         self.configuration = stats[0].configuration
#         assert "programs" in self.configuration, (stats[0].rulename, self.configuration)
#         self.outdir = os.path.dirname(os.path.dirname(stats[0].outdir))
#         self.output = {"stats": os.path.join(self.outdir, "pick", "comparison.stats")}
#         self.message = "Collecting statistics for Mikado in: {output[stats]}".format(output=self.output)
#
#     @property
#     def loader(self):
#         return ["mikado"]
#
#     @property
#     def rulename(self):
#         return "mikado_collect_stats"
#
#     @property
#     def threads(self):
#         return 1
#
#     @property
#     def cmd(self):
#         inputs = " ".join(self.input["stats"])
#         output = self.output
#         load = self.load
#         cmd = "{load} asm_collect.py {inputs} > {output[stats]}".format(**locals())
#         return cmd
