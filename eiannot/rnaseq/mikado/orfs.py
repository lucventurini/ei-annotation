import os
from .abstract import MikadoOp
import abc
from Bio.Data import CodonTable
from .prepare import MikadoPrepare
from ...preparation import DiamondIndex, SanitizeProteinBlastDB, BlastxIndex
from ...abstract import AtomicOperation, EIWrapper


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


class GTCDS(OrfCaller):

    def __init__(self, prepare: MikadoPrepare):

        super().__init__(prepare)
        self.input["genome"] = self.genome
        self.output["gff3"] = os.path.join(self.outdir, "mikado_prepared.gt_cds.gff3")
        self.output["orfs"] = os.path.join(self.outdir, "mikado_prepared.gt_cds.trans.bed12")

    @property
    def loader(self):
        return ["ei-annot", "genometools", "mikado"]

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "mikado_gtcds"

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        log = self.log
        cmd = "{load}"
        cmd += "  mkdir -p {outdir} && "
        input, output = self.input, self.output
        cmd += """ cat {input[gtf]} | awk '$3!~"(CDS|UTR)"' | """
        cmd += " mikado util convert -if gtf -of gff3 - | "
        cmd += " gt gff3 -tidy -retainids -addids | gt cds -seqfile {input[genome]} -matchdesc | "
        cmd += "gff3_name_to_id.py - {output[gff3]} "
        cmd += " && mikado util convert -t -of bed12 {output[gff3]} {output[orfs]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def genecode(self):
        return None


class TransdecoderLongOrf(OrfCaller):

    def __init__(self, prepare: MikadoPrepare):
        super().__init__(prepare)
        self.output = {"orfs": os.path.join(self.outdir,
                                            os.path.basename(self.input["fa"]) + ".transdecoder_dir",
                                            "longest_orfs.gff3"),
                       "pep": os.path.join(self.outdir,
                                            os.path.basename(self.input["fa"]) + ".transdecoder_dir",
                                            "longest_orfs.pep"),
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
        checkpoint_dir = os.path.basename(self.input["fa"]) + ".transdecoder_dir.__checkpoints"
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


class TransdecoderOrfDiamond(OrfCaller):

    __toolname__ = "diamond"

    def __init__(self, prepare: MikadoPrepare, long_orfs: TransdecoderLongOrf, diamond: DiamondIndex):

        super().__init__(prepare)
        self.input["pep"] = long_orfs.output["pep"]
        self.folder = os.path.dirname(self.input["pep"])
        self.output["txt"] = os.path.join(os.path.dirname(self.input["pep"]), "longest_orfs.blastp.txt")
        self.input["db"] = diamond.output["db"]

    @property
    def _rulename(self):
        return "mikado_transdecoder_blastp"

    @property
    def loader(self):
        return [self.__toolname__]

    @property
    def log(self):
        return os.path.join(self.folder, "blastp.log")

    @property
    def max_target_seqs(self):
        return self.configuration["mikado"]["homology"]["max_target_seqs"]

    @property
    def evalue(self):
        return self.configuration["mikado"]["homology"]["evalue"]

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        log = self.log
        input, output = self.input, self.output

        cmd = "{load} "
        cmd += "diamond blastp --threads {threads} --outfmt 6"
        cmd += " --out {output[txt]} --db {input[db]} --salltitles --query {input[pep]} --sensitive "
        evalue, max_targets = self.evalue, self.max_target_seqs
        cmd += " --max-target-seqs {max_targets} --evalue {evalue} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def genecode(self):
        return None


class TransdecoderOrfBlast(OrfCaller):

    __toolname__ = "blast"

    def __init__(self, prepare: MikadoPrepare, long_orfs: TransdecoderLongOrf, diamond: BlastxIndex):

        super().__init__(prepare)
        self.input["pep"] = long_orfs.output["pep"]
        self.folder = os.path.dirname(self.input["pep"])
        self.output["txt"] = os.path.join(os.path.dirname(self.input["pep"]), "longest_orfs.blastp.txt")
        self.input["db"] = diamond.output["db"]

    @property
    def _rulename(self):
        return "mikado_transdecoder_blastp"

    @property
    def loader(self):
        return [self.__toolname__]

    @property
    def log(self):
        return os.path.join(self.folder, "blastp.log")

    @property
    def max_target_seqs(self):
        return self.configuration["mikado"]["homology"]["max_target_seqs"]

    @property
    def evalue(self):
        return self.configuration["mikado"]["homology"]["evalue"]

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        log = self.log
        input, output = self.input, self.output

        cmd = "{load} "
        cmd += "blastp -num_threads {threads} -outfmt 6"
        cmd += " -out {output[txt]} -db {input[db]} -query {input[pep]} "
        evalue, max_targets = self.evalue, self.max_target_seqs
        cmd += " -max_target_seqs {max_targets} -evalue {evalue} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def genecode(self):
        return None


class TransdecoderPred(OrfCaller):

    def __init__(self,
                 prepare: MikadoPrepare,
                 long_orfs: TransdecoderLongOrf,
                 homology: TransdecoderOrfBlast):
        super().__init__(prepare)
        self.input["orfs"] = long_orfs.output["orfs"]
        if homology is not None:
            self.input["txt"] = homology.output["txt"]
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
        input, output = self.input, self.output
        cmd = "{load} "
        log = os.path.relpath(self.log, self.outdir)
        genecode = self.genecode
        checkpoint_dir = os.path.basename(self.input["fa"]) + ".transdecoder_dir.__checkpoints"
        checkpoint_file = os.path.join(checkpoint_dir, "make_final_cds*ok")
        cmd += "cd {outdir} && if [ -e {checkpoint_file} ]; then rm -rf {checkpoint_dir}; fi &&"
        cmd += " TransDecoder.Predict -t {fa} -G {genecode} "
        if self.input.get("txt", None) is not None:
            cmd += " --retain_blastp_hits {input[txt]} "
        cmd += " > {log} 2>&1"
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


class TransDecoderWrapper(EIWrapper):

    def __init__(self, preparer: MikadoPrepare):
        self.preparer = preparer
        super().__init__()
        self.configuration = preparer.configuration
        self.outdir = preparer.outdir
        longorf = TransdecoderLongOrf(self.preparer)
        self.add_edge(self.preparer, longorf)
        self.sanitizer = SanitizeProteinBlastDB(self.configuration, dbs=self.dbs)
        if self.execute:
            if self.program == "blastx":  # This will actually be blastp
                indexer = BlastxIndex(self.sanitizer)
                homology = TransdecoderOrfBlast(self.preparer, longorf, indexer)
            else:
                indexer = DiamondIndex(self.sanitizer)
                homology = TransdecoderOrfDiamond(self.preparer, longorf, indexer)
            self.add_edge(indexer, homology)
            self.add_edge(longorf, homology)
        else:
            homology = None
        self.orfs = TransdecoderPred(self.preparer, longorf, homology)
        self.add_edge(longorf, self.orfs)

        if homology is not None:
            self.add_edge(homology, self.orfs)

    @property
    def execute(self):
        return len(self.sanitizer.protein_dbs) > 0

    @property
    def dbs(self):
        return self.configuration["mikado"]["homology"]["prot_dbs"]

    @property
    def program(self):
        return self.configuration["mikado"]["homology"]["program"]

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "homology.done")
