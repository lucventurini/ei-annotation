from .abstract import MikadoOp
from ...preparation import FaidxGenome
from ..alignments.portcullis import PortcullisWrapper
from .prepare import MikadoPrepare
from .orfs import OrfCaller
from .homology import MikadoHomologyWrapper
import os


class MikadoSerialise(MikadoOp):

    def __init__(self,
                 prepare: MikadoPrepare,
                 homology: MikadoHomologyWrapper,
                 orfs: OrfCaller,
                 # faidx: FaidxGenome,
                 portcullis: PortcullisWrapper):

        super().__init__(is_long=prepare.is_long)
        self.input = prepare.output
        self.homology = homology
        self.orfs = orfs
        self.portcullis = portcullis
        if self.orfs:
            self.input["orfs"] = self.orfs.output["orfs"]
        if homology:
            for key, val in self.homology.output.items():
                self.input["homology_{key}".format(key=key)] = val
        if self.portcullis:
            for key, val in self.portcullis.output.items():
                self.input["portcullis_{key}".format(key=key)] = val
        # self.input.update(faidx.output)
        self.input["cfg"] = prepare.config
        self.configuration = prepare.configuration
        self.outdir = prepare.outdir
        self.output = {"db": os.path.join(self.outdir, "mikado.db")}
        self.log = os.path.join(self.outdir, "mikado_serialise.err")
        self.message = "Running Mikado serialise to move numerous data sources into a single database"
        self.portcullis = portcullis

    @property
    def _rulename(self):
        return "mikado_serialise"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def blast_xmls(self):
        if not self.homology:
            return ""

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
    def orfs_cli(self):
        if not self.orfs:
            return ""
        else:
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
        orfs = self.orfs_cli
        fai = self.genome + ".fai"
        cmd += "mikado serialise {blast_xmls} {blast_targets} {junctions} {orfs} --transcripts={input[fa]} "
        cmd += "--genome_fai={fai} --json-conf={input[cfg]} --force --start-method=spawn "
        outdir = self.outdir
        threads = self.threads
        log = self.log
        cmd += " -od {outdir} --procs={threads} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd