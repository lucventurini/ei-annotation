from ..abstract import EIWrapper, AtomicOperation, Linker
from ..preparation import PrepareWrapper, SanitizeProteinBlastDB
import os


class BuildModellerDB(AtomicOperation):

    def __init__(self, sanitiser: PrepareWrapper):

        super().__init__()
        self.configuration = sanitiser.configuration
        self.input = sanitiser.output
        self.input["genome"] = self.genome
        self.output = {"db": "{db}.nog".format(db=self.dbname)}
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "build.log")

    @property
    def loader(self):
        return ["repeatmodeller"]

    @property
    def rulename(self):
        return "build_repeatmodeller_db"

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


class RepeatModeller(AtomicOperation):

    outfile = "consensi.fa.classified"

    def __init__(self, builder: BuildModellerDB):
        super().__init__()
        self.configuration = builder.configuration
        self.input = builder.output
        self.output["families"] = os.path.join(self.outdir, self.outfile)
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "modeller.log")

    @property
    def dbname(self):
        return os.path.splitext(self.input["db"])[0]

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "modeller")

    @property
    def rulename(self):
        return "genome_repeat_modeller"

    @property
    def loader(self):
        return ["repeatmodeller"]

    @property
    def cmd(self):

        dbname, log = self.dbname, self.log
        load = self.load
        threads = self.threads
        outdir = self.outdir
        # We have to take into account we will be in another folder
        dbname = os.path.relpath(dbname, start=outdir)
        log = os.path.abspath(self.log)
        logdir = os.path.dirname(self.log)
        outfile = os.path.basename(self.output["families"])
        # Remove failed runs
        # for el in glob.glob(os.path.join(self.outdir, "RM*")):
        #     os.remove(el)

        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && cd {outdir} && "
        cmd += "(RepeatModeler -engine ncbi -pa {threads} -database {dbname} 2> {log} > {log} && "
        cmd += " cp RM*/{outfile} . && rm -rf RM* || touch {outfile})"
        cmd = cmd.format(**locals())
        return cmd


class PolishRepeats(AtomicOperation):

    masked_file = RepeatModeller.outfile + ".masked"
    outfile = "modelled_repeats.fa"
    __rulename__ = "polish_modeller_repeats"

    def __init__(self, modeller: RepeatModeller, proteins: SanitizeProteinBlastDB):

        super().__init__()
        self.configuration = modeller.configuration
        self.input = modeller.output
        self.input.update(proteins.output)
        self.output["families"] = os.path.join(modeller.outdir, self.masked_file)
        self.masked_dir = os.path.dirname(modeller.output["families"])
        self.output["link"] = os.path.join(self.outdir, self.outfile)
        self.log = os.path.join(self.masked_dir, "polish.log")

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output")

    @property
    def loader(self):
        return ["repeatmasker"]

    @property
    def rulename(self):
        return self.__rulename__

    @property
    def cmd(self):

        load = self.load
        maskdir = self.masked_dir
        outdir = os.path.relpath(self.outdir, start=self.masked_dir)
        rm_library = os.path.abspath(self.input["db"])
        link_src = os.path.relpath(self.output["families"], start=self.masked_dir)
        threads = self.threads
        families = os.path.abspath(self.input["families"])
        link_dest = self.outfile
        log = self.log

        cmd = "{load} mkdir -p {maskdir} && cd {maskdir} && "
        cmd += "RepeatMasker  -s –no_is –nolow -x -dir . -lib {rm_library} "
        cmd += "-pa {threads} {families} 2> {log} > {log} && "
        cmd += "rm -rf RM_* && mkdir -p {outdir}  && cd {outdir} && "
        cmd += " ln -s {link_src} {link_dest} && touch -h {link_dest}"

        cmd = cmd.format(**locals())

        return cmd

    @property
    def polishing_models(self):
        return self.configuration.get("repeats", dict()).get("safe_pro", None)


class ModellerWorkflow(EIWrapper):

    __final_rulename__ = PolishRepeats.__rulename__

    def __init__(self, sanitiser: PrepareWrapper):

        super().__init__()
        self.configuration = sanitiser.configuration
        if self.model_repeats is True:
            builder = BuildModellerDB(sanitiser)
            modeller = RepeatModeller(builder)
            # polisher = PolishRepeats(modeller)
            self.add_edges_from([(sanitiser, builder), (builder, modeller)])
            if len(self.safe_proteins) > 0:
                proteins = SanitizeProteinBlastDB(self.configuration,
                                                  db="repeatsafe",
                                                  dbs=self.safe_proteins)
                self.polisher = PolishRepeats(modeller, proteins)
                self.add_edge(proteins, self.polisher)
            else:
                self.polisher = Linker(modeller.output["families"],
                                  os.path.join(self.outdir, PolishRepeats.outfile),
                                  "families", "families", "link_unpolished_repeats",
                                  self.configuration)
            self.add_edges_from([(modeller, self.polisher)])
            assert self.exit

    @property
    def model_repeats(self):
        return self.configuration.get("repeats", dict()).get("model", True)

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output")

    @property
    def flag_name(self):
        return self.polisher.output["link"]

    @property
    def safe_proteins(self):
        return self.configuration["repeats"].get("safe_proteins", [])