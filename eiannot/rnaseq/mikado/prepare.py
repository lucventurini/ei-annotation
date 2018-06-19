from ...abstract import AtomicOperation, ShortSample, LongSample
from ..alignments.portcullis import PortcullisWrapper
from ..assemblies.workflow import AssemblyWrapper
from ..assemblies.abstract import ShortAssembler
from ..alignments.workflow import LongAlignmentsWrapper  # , LongAligner
import os
import itertools


class MikadoConfig(AtomicOperation):

    def __init__(self,
                 portcullis_wrapper: PortcullisWrapper,
                 assemblies: AssemblyWrapper,
                 long_aln_wrapper: LongAlignmentsWrapper):
        super().__init__()
        self.configuration = portcullis_wrapper.configuration
        self.portcullis = portcullis_wrapper
        self.assemblies = assemblies
        self.long_aln_wrapper = long_aln_wrapper
        self.input["genome"] = self.genome
        self.input["asm_list"] = os.path.join(self.outdir, "models_list.txt")
        self.input["gfs"] = [gf.input["gf"] for gf in assemblies.gfs]
        self.input["gfs"].extend(gf.input["gf"] for gf in long_aln_wrapper.gfs)
        self.input["portcullis"] = self.portcullis.junctions
        self.input['long_flag'] = long_aln_wrapper.output["flag"]
        self.input['asm_flag'] = assemblies.output["flag"]
        self.__create_file_list()
        self.output = {"cfg": os.path.join(self.outdir, "mikado.yaml")}

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "rnaseq", "5-mikado")

    @property
    def threads(self):
        return 1

    @property
    def message(self):
        return "Creating Mikado configuration file."

    @property
    def loader(self):
        return ["mikado"]

    @property
    def log(self):
        return os.path.join(self.outdir, "mikado_config.log")

    @property
    def cmd(self):

        load = self.load
        cmd = "{load} "
        scoring_file = self.scoring_file
        if not os.path.exists(os.path.dirname(self.input["asm_list"])):
            os.makedirs(os.path.dirname(self.input["asm_list"]))
        file_list = self.input["asm_list"]

        cmd += "mikado configure --scoring={scoring_file} --list={input[asm_list]} "
        log = self.log
        input, output = self.input, self.output
        external = self.external
        cmd += "--reference={input[genome]} {external} {output[cfg]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    def __create_file_list(self):
        if not os.path.exists(os.path.dirname(self.input["asm_list"])):
            os.makedirs(os.path.dirname(self.input["asm_list"]))

        if not os.path.exists(self.input["asm_list"]):
            with open(self.input["asm_list"], mode="wt") as file_list:
                for gf in itertools.chain(self.assemblies.gfs, self.long_aln_wrapper.gfs):
                    # Write out the location of the file, and all other details
                    try:
                        line = [gf.input["link"], gf.label, gf.sample.stranded]
                    except KeyError:
                        raise KeyError((gf.rulename, gf.output))
                    print(*line, file=file_list, sep="\t")

    @property
    def scoring_file(self):
        return self.configuration["mikado"]["pick"]["scoring_file"]

    @property
    def junctions(self):
        if self.portcullis is not None:
            return "--junctions{portcullis.output[bed]}".format(portcullis=self.portcullis)
        else:
            return ""

    @property
    def rulename(self):
        return 'mikado_config'

    @property
    def external(self):
        # TODO: implement
        return ""


class MikadoPrepare(AtomicOperation):

    def __init__(self, config: MikadoConfig):
        super().__init__()
        self.outdir = config.outdir
        self.configuration = config.configuration
        assert 'programs' in config.configuration, config.configuration
        assert 'mikado' in config.configuration["programs"], config.configuration["programs"]
        self.input = config.output
        self.output = {"gtf": os.path.join(self.outdir, "mikado_prepared.gtf"),
                       "fa": os.path.join(self.outdir, "mikado_prepared.fasta")}
        self.message = "Preparing transcripts using Mikado"
        self.log = os.path.join(os.path.dirname(self.output["gtf"]), "mikado_prepare.log")

    @property
    def rulename(self):
        return "mikado_prepare"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        log = self.log
        threads = self.threads
        input, output = self.input, self.output
        load = self.load
        cmd = "{load} "
        outdir = self.outdir
        cmd += "mikado prepare --procs={threads} "
        cmd += "--json-conf={input[cfg]} -od {outdir} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def config(self):
        return self.input["cfg"]
