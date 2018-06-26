from ..abstract import EIWrapper, AtomicOperation, Linker
from ..preparation import PrepareWrapper, SanitizeProteinBlastDB, FaidxProtein
from .modeler import ModelerWorkflow
import os
# from .modeler import


class RepeatMasking(EIWrapper):

    def __init__(self, sanitised: PrepareWrapper):

        super().__init__()
        self.configuration = sanitised.configuration

        if self.execute is True:
            proteins = SanitizeProteinBlastDB(self.configuration)

            if self.model is True:
                modeler = ModelerWorkflow(sanitised, proteins)
                self.add_edges_from([(sanitised, modeler), (proteins, modeler)])
            else:
                modeler = None
            if self.retrieve_known is True:
                retriever = RetrieveLibraries(self.configuration)
                self.add_edge(sanitised, retriever)
            else:
                retriever = None
            library_creator = LibraryCreator(sanitiser=sanitised, retriever=retriever, modeler=modeler)
            [self.add_edge(_, library_creator) for _ in (retriever, modeler) if _ is not None]
            masker = Masker(sanitised, library_creator)
            self.add_edge(library_creator, masker)
        else:
            linker = Linker(sanitised.exit.genome, sanitised.exit.masked_genome,
                            "genome", "masked", "link_genome_to_masked")
            self.add_node(linker)
        assert self.exit

    @property
    def model(self):
        return self.configuration.get("repeats", dict()).get("model", True)

    @property
    def retrieve_known(self):
        return (self.configuration.get("repeats", dict()).get("species", None) is not None or
                self.configuration.get("repeats", dict()).get("clade", None) is not None)

    @property
    def execute(self):
        return self.model or self.retrieve_known


class RetrieveLibraries(AtomicOperation):

    def __init__(self, configuration):

        super().__init__()
        self.configuration = configuration
        self.input = {}
        self.output = {"libraries": os.path.join(self.outdir, "retrieved.fa")}
        self.log = os.path.join(self.outdir, "logs", "extract_libraries.log")

    @property
    def cmd(self):

        load = self.load
        order = self.order
        output = self.output
        outdir = self.outdir
        logdir = os.path.dirname(self.log)
        log = self.log

        cmd = "{load} mkdir -p {logdir} && mkdir -p {outdir} && "
        cmd += " queryRepeatDatabase.pl {order} > {output[libraries]} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def rulename(self):
        return "extract_rm_libraries"

    @property
    def loader(self):
        return ["repeatmasker"]

    @property
    def order(self):
        if self.configuration["repeats"]["species"]:
            return " -species {species}".format(species=self.configuration["repeats"]["species"])
        elif self.configuration["repeats"]["clade"]:
            return " -clade {clade}".format(clade=self.configuration["repeats"]["clade"])

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "repeats", "masker")


class LibraryCreator(AtomicOperation):

    def __init__(self,
                 sanitiser: PrepareWrapper,
                 retriever: RetrieveLibraries,
                 modeler: ModelerWorkflow
                 ):
        if retriever is None and modeler is None:
            raise ValueError("No input libraries!")
        super().__init__()
        self.configuration = sanitiser.configuration

        if retriever is not None:
            self.input["retrieved"] = retriever.output["libraries"]
        if modeler is not None:
            self.input["modeled"] = modeler.output["families"]
        self.output = {"libraries": os.path.join(self.outdir, "rm_library.fa")}

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "masker")

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "create_rm_library"

    @property
    def cmd(self):
        retrieved = self.input.get("retrieved", " ")
        modeled = self.input.get("modeled", " ")
        output = self.output
        cmd = "cat {retrieved} {modeled} > {output[libraries]}".format(**locals())
        return cmd


class Masker(AtomicOperation):

    def __init__(self,
                 sanitised: PrepareWrapper,
                 library_creator: [None]  # Modeller
                 ):
        super().__init__()
        self.configuration = sanitised.configuration
        self.input["rm_library"] = library_creator.output["libraries"]
        self.input["genome"] = self.genome
        self.output["masked"] = self.masked_genome

    @property
    def rulename(self):
        return "masker"

    @property
    def loader(self):
        return ["repeatmasker"]

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        input, threads = self.input, self.threads

        cmd = "{load} mkdir -p {outdir} && cd {outdir} && "
        cmd += "RepeatMasker -xsmall -dir . -lib {input[rm_library]} -pa {threads} {input[genome]} && "
        cmd += " ln -s genome.fa.masked genome.masked.fa"

        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "masker")
