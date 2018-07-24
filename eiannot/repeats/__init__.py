from ..abstract import EIWrapper, AtomicOperation, Linker, Toucher
from ..preparation import PrepareWrapper, SanitizeProteinBlastDB, FaidxGenome
from .modeler import ModelerWorkflow
import os
# from .modeler import


class RetrieveLibraries(AtomicOperation):

    def __init__(self, configuration):

        super().__init__()
        self.configuration = configuration
        self.input = {"mock": "."}
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
        self.output["genome"] = self.masked_genome
        self.output["masked"] = os.path.join(self.maskdir, "genome.fa.masked")
        self.output["mask_table"] = os.path.join(self.maskdir, "genome.fa.out")
        self.output["table"] = os.path.join(self.outdir, "repeats.table")
        self.log = os.path.join(self.maskdir, "repeat_masker.log")

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
        rm_library = os.path.basename(self.input["rm_library"])
        log = os.path.basename(self.log)
        outdir = os.path.abspath(self.outdir)
        maskdir = self.maskdir
        link_src = os.path.relpath(self.output["masked"], start=outdir)
        link_dest = os.path.basename(self.output["genome"])
        genome = os.path.relpath(os.path.abspath(self.genome), start=maskdir)
        gff_link_src = os.path.relpath(self.output["mask_table"], start=outdir)
        gff_link_dest = os.path.basename(self.output["table"])

        cmd = "{load} mkdir -p {outdir} && mkdir -p {maskdir} && cd {maskdir} && "
        cmd += "RepeatMasker -nolow -xsmall -dir -gff . -lib {rm_library} -pa {threads} {genome} 2> {log} > {log} && "
        cmd += "rm -rf RM_* && cd {outdir} && ln -s {link_src} {link_dest} && touch -h {link_dest} && "
        cmd += " ln -s {gff_link_src} {gff_link_dest} && touch -h {gff_link_dest}"

        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output")

    @property
    def maskdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "masker")


class FaidxMaskedGenome(FaidxGenome):

    def __init__(self, masker: Masker):

        super().__init__(sanitiser=masker)

    @property
    def rulename(self):
        return "faidx_masked_genome"

    @property
    def outdir(self):
        return os.path.dirname(self.masked_genome)


class RepeatMasking(EIWrapper):

    __final_rulename__ = "masker_done"

    def __init__(self, sanitised: PrepareWrapper):

        super().__init__()
        self.configuration = sanitised.configuration

        if self.execute is True:
            proteins = SanitizeProteinBlastDB(self.configuration)

            if self.model is True:
                modeler = ModelerWorkflow(sanitised)
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
            faidx = FaidxMaskedGenome(masker)
            self.add_edge(masker, faidx)

        else:
            linker = Linker(sanitised.exit.genome, sanitised.exit.masked_genome,
                            "genome", "genome", "link_genome_to_masked", self.configuration)
            linker.outdir = self.outdir
            self.add_node(linker)
            fai_linker = Linker(sanitised.fai.output["fai"], sanitised.exit.masked_genome + ".fai",
                                "fai", "fai", "link_genome_fai_to_masked", self.configuration)
            fai_linker.outdir = self.outdir
            self.add_node(fai_linker)
            self.add_edge(linker, fai_linker)

        self.add_final_flag()
        assert self.exit

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output")

    # @property
    # def output(self):
    #     return

    @property
    def model(self):
        return self.configuration.get("repeats", dict()).get("model", True)

    @property
    def retrieve_known(self):
        return (self.configuration.get("repeats", dict()).get("species", None) is not None or
                self.configuration.get("repeats", dict()).get("clade", None) is not None)

    @property
    def execute(self):

        return (self.model or self.retrieve_known) and self.configuration.get("repeats", dict()).get("execute", False)

    @property
    def masked_genome(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output", "genome.masked.fa")

    @property
    def fai(self):
        return os.path.join(self.configuration["outdir"], "repeats", "output", "genome.masked.fa.fai")

    @property
    def flag_name(self):
        return os.path.join(self.exit.outdir, "repeat_masking.done")