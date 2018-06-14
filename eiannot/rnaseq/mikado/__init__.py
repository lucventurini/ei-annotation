from ...abstract import AtomicOperation, EIWorfkflow, ShortSample, LongSample
import os
# import abc
from ..assemblies.workflow import AssemblyWrapper
from ..alignments.workflow import ShortAlignmentsWrapper, LongAlignmentsWrapper, PortcullisWrapper
from .prepare import MikadoPrepare, MikadoConfig
from .orfs import OrfCaller
from .homology import MikadoHomologyWrapper
from ...preparation import FaidxGenome


class Mikado(EIWorfkflow):

    def __init__(self,
                 assemblies: AssemblyWrapper,
                 long_alignments: LongAlignmentsWrapper,
                 portcullis: PortcullisWrapper):

        super().__init__()
        self.configuration = assemblies.configuration
        self.assemblies = assemblies
        self.long_alignments = long_alignments
        self.portcullis = portcullis

        self.configurer = MikadoConfig(portcullis_wrapper=self.portcullis,
                                       assemblies=self.assemblies,
                                       long_aln_wrapper=long_alignments)
        self.add_edges_from([step, self.configurer] for step in
                            [self.assemblies, self.long_alignments, self.portcullis])

        self.preparer = MikadoPrepare(self.configuration)
        self.add_edge(self.configuration, self.preparer)
        self.orfs = OrfCaller(self.preparer)
        self.add_edge(self.preparer, self.orfs)  # TODO: this will fail currently. I have to modify the abstract class
        self.homologies = MikadoHomologyWrapper(self.preparer)
        self.faidx_genome = FaidxGenome(self.configuration)
        self.serialiser = MikadoSerialise(prepare=self.preparer,
                                          homology=self.homologies,
                                          orfs=self.orfs,
                                          faidx=self.faidx_genome)
        mikados = []
        for mode in self.modes:
            picker = MikadoPick(self.serialiser, mode)
            self.add_edge(self.serialiser, picker)
            stats = MikadoStats(picker)
            self.add_edge(picker, stats)
            mikados.append(stats)
        collection = MikadoCollectStats(mikados)
        self.add_edges_from([mikado, collection] for mikado in mikados)

    @property
    def gfs(self):
        return self.assemblies.gfs

    @property
    def junctions(self):
        return self.portcullis.junctions

    @property
    def modes(self):
        # TODO: verify this is correct
        return self.configuration["mikado"]["modes"]


class MikadoSerialise(AtomicOperation):

    def __init__(self,
                 prepare: MikadoPrepare,
                 homology: MikadoHomologyWrapper,
                 orfs: OrfCaller,
                 faidx: FaidxGenome,
                 portcullis=PortcullisWrapper):

        super().__init__()
        self.input = prepare.output
        self.input.update(orfs.output)
        self.input.update(homology.output)
        self.input.update(faidx.output)
        self.input["cfg"] = prepare.config
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
            return "--blast-targets={}".format(self.homology.blast_targets)
        else:
            return ""

    @property
    def orfs(self):
        return " --orfs={input[orfs]}".format(input=self.input)

    @property
    def junctions(self):
        if self.portcullis:
            return "--junctions={portcullis.junctions}".format(portcullis=self.portcullis)
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

    def __init__(self, serialise: MikadoSerialise, mode):

        super().__init__()
        self.configuration = serialise.configuration
        self.input = serialise.output
        self.input["gtf"] = serialise.input["gtf"]
        self.outdir = serialise.outdir
        self.__mode = mode
        self.output = {"loci": os.path.join(
            self.outdir, "pick", "{mode}", "mikado-{mode}.loci.gff3").format(**locals())}
        self.log = os.path.join(os.path.basename(self.output["loci"]), "mikado-{mode}.pick.err").format(
            **locals()
        )
        self.message = "Running mikado picking stage in mode: {mode}".format(**locals())

    @property
    def mode(self):
        return self.__mode

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
        cmd += "--start-method=spawn --json-conf={input[cfg]}"
        loci_out = os.path.basename(self.output["loci"])
        outdir = self.outdir
        log = self.log
        cmd += "-od {outdir} --loci_out {loci_out}  -lv INFO -db {input[db]} {input[gtf]} > {log} 2>&1"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def loci_dir(self):
        return os.path.dirname(self.output["loci"])


class MikadoStats(AtomicOperation):

    def __init__(self, pick: MikadoPick):

        super().__init__()
        self.input = pick.output
        self.outdir = pick.outdir
        self.__mode = pick.mode
        self.output = {"stats": os.path.join(self.loci_dir, "mikado-{mode}.loci.stats")}
        self.message = "Calculating statistics for Mikado run in mode: {mode}".format(mode=self.mode)
        self.log = self.output["stats"] + ".log"

    @property
    def mode(self):
        return self.__mode

    @property
    def rulename(self):
        return "mikado_stats_{mode}".format(mode=self.mode)

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["loci"])

    @property
    def cmd(self):

        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[loci]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class MikadoCollectStats(AtomicOperation):

    def __init__(self, stats: [MikadoStats]):

        assert len(stats) > 0 and all(isinstance(stat, MikadoStats) for stat in stats)
        super().__init__()
        self.input["stats"] = [stat.output["stats"] for stat in stats]
        self.configuration = stats[0].configuration
        self.outdir = stats[0].outdir
        self.output = {"stats": os.path.join(self.outdir, "pick", "comparison.stats")}
        self.message = "Collecting statistics for Mikado in: {output[stats]}".format(output=self.output)

    @property
    def loader(self):
        return ["mikado"]

    @property
    def rulename(self):
        return "mikado_collect_stats"

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        inputs = " ".join(self.input["stats"])
        output = self.output
        load = self.load
        cmd = "{load} asm_collect.py {inputs} > {output[stats]}".format(**locals())
        return cmd
