from ..alignments.portcullis import PortcullisWrapper
# from ..assemblies.workflow import AssemblyWrapper
from .abstract import MikadoOp
from ..alignments.__init__ import LongAlignmentsWrapper  # , LongAligner
import os


class MikadoConfig(MikadoOp):

    def __init__(self,
                 portcullis_wrapper: PortcullisWrapper,
                 assemblies,  #: AssemblyWrapper,
                 long_aln_wrapper: LongAlignmentsWrapper,
                 is_long: bool):
        super().__init__(is_long=is_long)
        if portcullis_wrapper:
            self.configuration = portcullis_wrapper.configuration
        elif long_aln_wrapper:
            self.configuration = long_aln_wrapper.configuration
        elif assemblies:
            self.configuration = assemblies.configuration
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
        return self.mikado_dir

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
                if not self.is_long:
                    for gf in self.assemblies.gfs:
                        try:
                            line = [gf.input["gf"], gf.label, gf.sample.stranded]
                        except KeyError:
                            raise KeyError((gf.rulename, gf.output))
                        print(*line, file=file_list, sep="\t")
                    score_add = self.long_bias_score
                else:
                    score_add = 0
                for gf in self.long_aln_wrapper.gfs:
                    try:
                        line = [gf.input["gf"], gf.label, gf.sample.stranded, score_add]
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
    def long_bias_score(self):
        return self.configuration["mikado"]["pick"].get("long_bias_score", 1)

    @property
    def _rulename(self):
        return 'mikado_config'

    @property
    def external(self):
        # TODO: implement
        return ""


class MikadoPrepare(MikadoOp):

    def __init__(self, config: MikadoConfig):
        super().__init__(is_long=config.is_long)
        self.configuration = config.configuration
        assert 'programs' in config.configuration, config.configuration
        assert 'mikado' in config.configuration["programs"], config.configuration["programs"]
        self.input = config.output
        self.output = {"gtf": os.path.join(self.outdir, "mikado_prepared.gtf"),
                       "fa": os.path.join(self.outdir, "mikado_prepared.fasta")}
        self.message = "Preparing transcripts using Mikado"
        self.log = os.path.join(os.path.dirname(self.output["gtf"]), "mikado_prepare.log")

    @property
    def _rulename(self):
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

    @property
    def outdir(self):
        return self.mikado_dir
