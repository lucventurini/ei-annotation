from ..abstract import EIWrapper, AtomicOperation
from ..preparation import SanitizeGenome, SanitizeProteinBlastDB
import os
import glob


class ModelerWorkflow(EIWrapper):

    def __init__(self, sanitiser: SanitizeGenome,
                 sanitised_proteins: SanitizeProteinBlastDB):

        super().__init__()
        self.configuration = sanitiser.configuration
        if self.model_repeats is True:
            builder = BuildModelerDB(sanitiser)
            modeler = RepeatModeler(builder)
            # polisher = PolishRepeats(modeler)
            self.add_edges_from([(sanitiser, builder), (builder, modeler)])
            # TODO: how are we going to polish this?
            # self.add_edges_from([(modeler, polisher)])
            assert self.exit

    @property
    def model_repeats(self):
        return self.configuration.get("repeats", dict()).get("model", True)


class BuildModelerDB(AtomicOperation):

    def __init__(self, sanitiser: SanitizeGenome):

        super().__init__()
        self.configuration = sanitiser.configuration
        self.input = sanitiser.output
        self.output = {"db": "{db}.nog".format(db=self.dbname)}
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "build.log")

    @property
    def loader(self):
        return ["repeatmodeler"]

    @property
    def rulename(self):
        return "build_repeatmodeler_db"

    @property
    def cmd(self):

        name = self.dbname
        input, log = self.input, self.log
        load = self.load
        logdir = os.path.dirname(self.log)

        cmd = "{load} mkdir -p {logdir} && BuildDatabase -name {name} -engine ncbi {input[genome]} > {log} 2>&1"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "repeats", "database")

    @property
    def dbname(self):

        return os.path.join(self.outdir, "genome")

    @property
    def threads(self):
        return 1


class RepeatModeler(AtomicOperation):

    def __init__(self, builder: BuildModelerDB):
        super().__init__()
        self.input = builder.output
        self.output["families"] = os.path.join(self.outdir, "consensi.fa.classified")
        self.configuration = builder.configuration
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "modeler.log")


    @property
    def dbname(self):
        return os.path.splitext(self.input["db"])[0]

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "modeler")

    @property
    def rulename(self):
        return "genome_repeat_modeler"

    @property
    def loader(self):
        return ["repeatmodeler"]

    @property
    def cmd(self):

        dbname, log = self.dbname, self.log
        load = self.load
        threads = self.threads
        outdir = self.outdir
        log = os.path.abspath(self.log)
        logdir = os.path.dirname(self.log)
        outfile = os.path.basename(self.output["families"])
        # Remove failed runs
        for el in glob.glob(os.path.join(self.outdir, "RM*")):
            os.remove(el)

        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && cd {outdir} && "
        cmd += "(RepeatModeler -engine ncbi -pa {threads} -database {dbname} 2> {log} > {log} && "
        cmd += " ln -s RM*/{outfile} . ) || touch {outfile})"
        cmd = cmd.format(**locals())
        return cmd


class PolishRepeats(AtomicOperation):

    def __init__(self, modeler: RepeatModeler):

        super().__init__()
        self.configuration = modeler.configuration
        self.input = modeler.output

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "polish_modeler_repeats"


