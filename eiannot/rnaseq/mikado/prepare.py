from ...abstract import AtomicOperation, ShortSample, LongSample
from ..alignments.portcullis import PortcullisWrapper
from ..assemblies.workflow import AssemblyWrapper
from ..assemblies.abstract import ShortAssembler
# from ..alignments.workflow import LongAlignmentsWrapper, LongAligner
import os


class MikadoConfig(AtomicOperation):

    def __init__(self,
                 portcullis_wrapper: PortcullisWrapper,
                 assemblies: AssemblyWrapper,
                 long_aln_wrapper):
        super().__init__()
        self.configuration = portcullis_wrapper.configuration
        self.portcullis = portcullis_wrapper
        self.assemblies = assemblies
        self.long_aln_wrapper = long_aln_wrapper
        self.input["genome"] = self.genome
        self.input["asm_list"] = os.path.join(self.outdir, "models_list.txt")
        self.input["gfs"] = [gf.output["gf"] for gf in assemblies.gfs]
        self.input["gfs"].extend([gf.output["gf"] for gf in long_aln_wrapper.gfs])  # TODO: implement
        self.input["portcullis"] = self.portcullis.output["bed"]
        self.__create_file_list()
        self.output = {"cfg": os.path.join(self.outdir, "mikado.yaml")}


    @property
    def outdir(self):

        return os.path.join(self.configuration["out_dir"],
                            "rnaseq", "4-mikado")

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
        os.makedirs(os.path.dirname(self.input["asm_list"]))
        file_list = self.input["asm_list"]

        cmd += "mikado configure --scoring={scoring_file} --list={input[asm_list]} "
        log = self.log
        input, output = self.input, self.output
        external = ''  # TODO: implement
        cmd += "--reference={input[genome]} --external={external} {output[cfg]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    def __create_file_list(self):
        with open(self.input["asm_list"], mode="wt") as file_list:
            for gf in self.assemblies.gfs:
                # this GFs are the *rules*
                assert isinstance(gf, ShortAssembler)
                sample = gf.sample
                assert isinstance(sample, ShortSample)
                print(gf.output["gf"], sample.label, sample.stranded, sep="\t", file=file_list)
            for gf in self.long_aln_wrapper:
                # assert isinstance(gf, LongAligner)
                sample = gf.sample
                assert isinstance(sample, LongSample)
                print(gf.output["gf"], sample.label, sample.stranded, sep="\t", file=file_list)

    @property
    def scoring_file(self):
        return self.configuration["mikado"]["pick"]["scoring_file"]

    @property
    def junctions(self):
        if self.portcullis is not None:
            return "--junctions{portcullis.output[bed]}".format(portcullis=self.portcullis)
        else:
            return ""


class MikadoPrepare(AtomicOperation):

    def __init__(self, config: MikadoConfig):
        super().__init__()
        self.outdir = config.outdir
        self.configuration = config.configuration
        self.input = config.output
        self.output = {"gtf": os.path.join(self.outdir, "mikado_prepared.gtf"),
                       "fa": os.path.join(self.outdir, "mikado_prepared.fasta")}
        self.message = "Preparing transcripts using Mikado"

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
        cmd = "{load} "
        outdir = self.outdir
        cmd += "mikado prepare --procs={threads} --fasta={input[ref]} "
        cmd += "--json-conf={input[cfg]} -od {outdir} > {log} 2> &1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def config(self):
        return self.input["cfg"]
