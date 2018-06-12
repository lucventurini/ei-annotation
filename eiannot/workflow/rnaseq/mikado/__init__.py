from ... import AtomicOperation, EIWorfkflow, ShortSample, LongSample
import os
import abc
from .prepare import MikadoPrepare
from .orfs import OrfCaller
from .homology import HomologyWrapper  # TODO: implement
from ...preparation import FaidxGenome


class Mikado(EIWorfkflow):

    def __init__(self):

        super().__init__()



class MikadoSerialise(AtomicOperation):

    def __init__(self,
                 prepare: MikadoPrepare,
                 homology: HomologyWrapper,
                 orfs: OrfCaller,
                 faidx: FaidxGenome):

        super().__init__()
        self.input = prepare.output
        self.input.update(orfs.output)
        self.input.update(homology.output)
        self.input.update(faidx.output)
        self.input["cfg"] = prepare.config
        self.outdir = prepare.outdir
        self.output = {"db": os.path.join(self.outdir, "mikado.db")}
        self.log = os.path.join(self.outdir, "mikado_serialise.err")
        self.message = "Running Mikado serialise to move numerous data sources into a single database"

    @property
    def rulename(self):
        return "mikado_serialise"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def blast_xmls(self):
        # TODO: implement
        pass

    @property
    def blast_targets(self):
        # TODO: implement
        pass

    @property
    def orfs(self):
        return " --orfs={input[orfs]}".format(input=self.input)

    @property
    def cmd(self):
        load = self.load
        cmd = "{load} "
        blast_xmls = self.blast_xmls
        blast_targets = self.blast_targets
        input = self.input
        cmd += "mikado serialise {blast_xmls} {blast_targets} --start-method=spawn --transcripts={input[fa]} "
        orfs = self.orfs
        cmd += "--genome_fai={input[fai]} --json-conf={input[cfg]} --force {orfs} "
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
