import os
from ...abstract import AtomicOperation
import abc
from Bio.Data import CodonTable
from .prepare import MikadoPrepare


class OrfCaller(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self, prepare: MikadoPrepare):

        super().__init__()
        self.configuration = prepare.configuration
        self.outdir = os.path.join(prepare.outdir, "orfs")
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
    def rulename(self):
        return "mikado_prodigal"

    @property
    def loader(self):
        return ["prodigal"]

    @property
    def genecode(self):
        return CodonTable.generic_by_name[self.configuration["homology"]["genecode"]].id

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
                                            "transcripts.fasta.transdecoder_dir",
                                            "longest_orfs.gff3")}
        self.log = os.path.join(self.outdir, "transdecoder.longorf.log")
        self.message = "Running transdecoder longorf on Mikado prepared transcripts: {input[fa]}".format(
            input=self.input)

    @property
    def threads(self):
        return 1

    @property
    def rulename(self):
        return "mikado_transdecoder_longorfs"

    @property
    def loader(self):
        return ["transdecoder"]

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        cmd = "{load} mkdir -p {outdir} && cd {outdir} &&"
        genecode = self.genecode
        link_src = os.path.relpath(self.input["fa"], start=self.outdir)
        fa = os.path.basename(self.input["fa"])
        cmd += "ln -sf {params.tr_in} {fa} && "
        minprot = self.minprot
        log = self.log
        genecode = self.genecode
        cmd += "TransDecoder.LongOrfs -m {minprot} -t {fa} --genetic_code {genecode} > {log} 2>&1"
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
        self.output = {"orfs": os.path.join(self.outdir, "transcripts.fasta.transdecoder.bed")}

    @property
    def rulename(self):
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

        cmd += "cd {outdir} && TransDecoder.Predict -t {fa} > {log} 2>&1"

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