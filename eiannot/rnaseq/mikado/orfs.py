import os
from .abstract import MikadoOp
import abc
from Bio.Data import CodonTable
from .prepare import MikadoPrepare


class OrfCaller(MikadoOp, metaclass=abc.ABCMeta):

    def __init__(self, prepare: MikadoPrepare):

        super().__init__(is_long=prepare.is_long)
        self.configuration = prepare.configuration
        self.outdir = os.path.join(self.mikado_dir, "orfs")
        self.input = prepare.output

    @property
    def minprot(self):

        return self.configuration["mikado"]["orfs"]["min_protein_length"]

    @abc.abstractmethod
    def genecode(self):
        pass

    @property
    def orf(self):
        return self.output["orfs"]


class Prodigal(OrfCaller):

    def __init__(self, prepare: MikadoPrepare):
        super().__init__(prepare)
        self.output = {"orfs": os.path.join(self.outdir, "transcripts.fasta.prodigal.gff3")}
        self.log = os.path.join(self.outdir, "prodigal.log")
        self.message="Running PRODIGAL on Mikado prepared transcripts: {input[fa]}".format(input=self.input)

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "mikado_prodigal"

    @property
    def loader(self):
        return ["prodigal"]

    @property
    def genecode(self):
        code = self.configuration.get("homology", dict()).get("genecode", "Standard")

        return CodonTable.generic_by_name[code].id

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        cmd = "{load} "
        cmd += " mkdir -p {outdir} && "
        log = self.log
        input, output = self.input, self.output
        genecode = self.genecode
        cmd += "prodigal -f gff -g {genecode} -i {input[fa]} -o {output[orfs]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd


class TransdecoderLongOrf(OrfCaller):

    def __init__(self, prepare: MikadoPrepare):
        super().__init__(prepare)
        self.output = {"orfs": os.path.join(self.outdir,
                                            os.path.basename(self.input["fa"]) + ".transdecoder_dir",
                                            "longest_orfs.gff3"),
                       "fa_link": os.path.join(self.outdir, os.path.basename(self.input["fa"]))}
        self.log = os.path.join(self.outdir, "transdecoder.longorf.log")
        self.message = "Running transdecoder longorf on Mikado prepared transcripts: {input[fa]}".format(
            input=self.input)

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "mikado_transdecoder_longorfs"

    @property
    def loader(self):
        return ["transdecoder"]

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        cmd = "{load} mkdir -p {outdir} && cd {outdir} && "
        genecode = self.genecode
        link_src = os.path.relpath(self.input["fa"], start=self.outdir)
        fa = os.path.basename(self.input["fa"])
        input_fa = os.path.relpath(self.input["fa"], outdir)
        cmd += "ln -sf {input_fa} {fa} && "
        checkpoint_dir = self.input["fa"] + ".transdecoder_dir.__checkpoints"
        cmd += " if [ -d {checkpoint_dir} ]; then rm -rf {checkpoint_dir}; fi && "
        minprot = self.minprot
        log = os.path.relpath(self.log, self.outdir)
        genecode = self.genecode
        cmd += "TransDecoder.LongOrfs -m {minprot} -t {fa} -G {genecode} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def genecode(self):
        """Brian Haas implemented his own version of the genetic codes within his Perl Library.
        As such, it will require time to link his names to the official NCBI codon tables.
        See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for a primer in the future"""
        # TODO: expand to other species
        return "Universal"


class TransdecoderPred(OrfCaller):

    def __init__(self,
                 prepare: MikadoPrepare,
                 long_orfs: TransdecoderLongOrf):
        super().__init__(prepare)
        self.input["orfs"] = long_orfs.output["orfs"]
        self.log = os.path.join(self.outdir, "transdecoder.predict.log")
        self.message="Running transdecoder predict on Mikado prepared transcripts: {input[fa]}".format(
            input=self.input
        )
        self.output = {"orfs": long_orfs.output["fa_link"] + ".transdecoder.bed",
                       "cds": long_orfs.output["fa_link"] + ".transdecoder.cds",
                       "gff": long_orfs.output["fa_link"] + ".transdecoder.gff3",
                       "pep": long_orfs.output["fa_link"] + ".transdecoder.pep"}

    @property
    def _rulename(self):
        return "mikado_transdecoder_pred"

    @property
    def loader(self):
        return ["transdecoder"]

    @property
    def cmd(self):
        # mikado = os.path.relpath(self.input["fa"], start=self.outdir)
        longorf = os.path.relpath(self.input["orfs"], start=self.outdir)
        fa = os.path.basename(self.input["fa"])
        log = self.log
        load = self.load
        outdir = self.outdir
        cmd = "{load} "
        log = os.path.relpath(self.log, self.outdir)
        genecode = self.genecode
        checkpoint_dir = self.input["fa"] + ".transdecoder_dir.__checkpoints"
        checkpoint_file = os.path.join(checkpoint_dir, "make_final_cds*ok")
        cmd += "cd {outdir} && if [ -e {checkpoint_file} ]; then rm -rf {checkpoint_dir}; fi &&"
        cmd += " TransDecoder.Predict -t {fa} -G {genecode} > {log} 2>&1"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def genecode(self):
        """Brian Haas implemented his own version of the genetic codes within his Perl Library.
                As such, it will require time to link his names to the official NCBI codon tables.
                See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for a primer in the future"""
        # TODO: expand to other species
        return "Universal"
