import os
from ...preparation import DiamondIndex
from .abstract import MikadoOp
from .pick import MikadoStats


class MikadoSequenceExtractor(MikadoOp):

    def __init__(self, stats: MikadoStats):

        super().__init__(stats)
        self.input["loci"] = stats.input["link"]
        self.input["genome"] = self.genome
        self.__outdir = stats.outdir
        self.output = {"transcripts": os.path.join(self.outdir, "mikado.loci.transcripts.fasta"),
                       "proteins": os.path.join(self.outdir, "mikado.loci.proteins.fasta"),
                       "cds": os.path.join(self.outdir, "mikado.loci.cds.fasta")}

    @property
    def outdir(self):
        return self.__outdir

    @property
    def loader(self):
        return ["gffread"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        genome = self.genome
        outdir = self.outdir
        cmd = "{load} gffread -g {genome} -w {output[transcripts]}  "
        cmd += " -x {output[cds]} -y {output[proteins]} "
        cmd += " {input[loci]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "extract_mikado_fasta"

    @property
    def is_small(self):
        return True


class MikadoSequenceIndexer(MikadoOp):

    def __init__(self, extractor: MikadoSequenceExtractor):
        super().__init__(extractor)
        self.__outdir = extractor.outdir
        self.output = {"fai": self.input["transcripts"] + ".fai"}
        self.log = os.path.join(self.outdir, "mikado.faidx.log")

    @property
    def outdir(self):
        return self.__outdir

    @property
    def _rulename(self):
        return "faidx_mikado_transcripts"

    @property
    def loader(self):
        return ["samtools"]

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} samtools faidx {input[transcripts]} > {log} 2>&1".format(**locals())

        return cmd

    @property
    def is_small(self):
        return True


class MikadoDiamondIndex(DiamondIndex):

    def __init__(self, extractor: MikadoSequenceExtractor):
        self.is_long = extractor.is_long
        super().__init__(extractor, key="proteins", rule_suffix="fln{}".format("_long" if self.is_long else ""))


class SelfDiamondP(MikadoOp):

    def __init__(self, index: MikadoDiamondIndex):

        super().__init__(is_long=index.is_long)
        self.configuration = index.configuration
        self.input = index.output
        self.input["query"] = index.input["db"]
        self.output["blast_txt"] = os.path.join(
            self.outdir,
            os.path.splitext(os.path.basename(self.input["query"]))[0] + ".dmnd.txt")
        self.__loader = index.loader
        self.log = os.path.join(self.outdir, "diamond_proteins.log")

    @property
    def outdir(self):
        return os.path.dirname(self.input["query"])

    @property
    def fmt(self):
        line = "6 qseqid sseqid pident qstart qend sstart send "
        line += "qlen slen length nident mismatch positive gapopen gaps evalue bitscore"
        return line

    @property
    def loader(self):
        return self.__loader

    @property
    def cmd(self):
        load, input, output = self.load, self.input, self.output
        threads = self.threads
        log = self.log
        fmt = self.fmt
        cmd = "{load} "
        cmd += "diamond blastp --threads {threads} --outfmt {fmt} --compress 0 "
        cmd += " --out {output[blast_txt]} --db {input[db]} --query {input[query]} --sensitive "
        evalue = 1
        cmd += " --evalue {evalue} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def _rulename(self):
        return "self_diamond_blastp"


class ConvertMikadoToBed12(MikadoOp):

    def __init__(self, mikado: MikadoStats):

        super().__init__(is_long=mikado.is_long)
        self.configuration = mikado.configuration
        self.input = mikado.output
        self.input.update(mikado.input)
        self.output = {"bed12": os.path.splitext(mikado.input["link"])[0] + ".bed12"}
        self.log = os.path.join(self.mikado_dir, "logs", "convert_to_bed12.log")

    @property
    def _rulename(self):
        return "convert_mikado_to_bed12"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        log = self.log
        cmd = "{load} mikado util convert -of bed12 {input[loci]} {output[bed12]} 2> {log} > {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def is_small(self):
        return True