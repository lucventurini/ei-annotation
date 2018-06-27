from ..abstract import EIWrapper, AtomicOperation
from ..preparation import SanitizeProteinBlastDB, FaidxProtein
from ..repeats.workflow import RepeatMasking
import os


class ExonerateProteinWrapper(EIWrapper):

    def __init__(self,
                 masker: RepeatMasking):

        super().__init__()
        self.configuration = masker.configuration
        sanitised = SanitizeProteinBlastDB(self.configuration)

        if sanitised.protein_dbs:
            faidx = FaidxProtein(sanitised)
            self.add_edge(sanitised, faidx)
            chunk_proteins = ChunkProteins(sanitised)
            chunks = [Exonerate(chunk_proteins, chunk, masker) for chunk in range(1, chunk_proteins.chunks + 1)]
            self.add_edges_from([(chunk_proteins, chunk) for chunk in chunks])
            self.add_edges_from([(masker, chunk) for chunk in chunks])
            collapsed = CollapseExonerate(chunks)
            self.add_edges_from([(chunk, collapsed) for chunk in chunks])
            convert = ConvertExonerate(collapsed, faidx)
            self.add_edge(faidx, convert)
            self.add_edge(collapsed, convert)


class ChunkProteins(AtomicOperation):

    def __init__(self, sanitised: SanitizeProteinBlastDB):

        super().__init__()
        self.configuration = sanitised.configuration
        self.input = sanitised.output  # input[db] is our file
        self.output["flag"] = os.path.join(os.path.dirname(self.outdir), "chunking.done")
        self.output["chunks"] = [os.path.join(self.outdir, "chunk_{}.fasta".format(
            str(chunk).zfill(3))) for chunk in range(1, self.chunks + 1)]
        self.log = os.path.join(os.path.dirname(self.outdir), "logs", "chunking.log")

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"],
                            "proteins",
                            "chunks")

    @property
    def rulename(self):
        return "chunk_proteins_for_exonerate"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load

        outdir = self.outdir
        chunks = self.chunks
        input, output, log = self.input, self.output, self.log
        logdir = os.path.dirname(self.log)
        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && "
        cmd += " split_fasta.py -m {chunks} {input[db]} {outdir}/chunk 2> {log} > {log}"
        cmd += " && touch {output[flag]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def chunks(self):
        return self.configuration["homology"]["protein_chunks"]


class Exonerate(AtomicOperation):

    def __init__(self,
                 chunks: ChunkProteins,
                 chunk,
                 masked: RepeatMasking
                 ):
        super().__init__()
        self._chunks = chunks
        self.configuration = chunks.configuration
        self.__chunk = None
        self.chunk = chunk
        self.input["flag"] = chunks.output["flag"]
        self.input["genome"] = masked.output["masked"]
        assert self.input["fasta"] in chunks.output["chunks"]
        self.output["txt"] = os.path.join(self.outdir, "{chunk}.exonerate.txt").format(chunk=self.chunk)

    @property
    def chunk(self):
        return self.__chunk

    @chunk.setter
    def chunk(self, chunk):

        assert isinstance(chunk, int)
        fasta = os.path.join(self._chunks.outdir, "chunk_{}.fasta".format(str(chunk).zfill(3)))
        assert fasta in self._chunks.output["chunks"], (fasta, self._chunks.output["chunks"])
        self.input["fasta"] = fasta
        self.__chunk = chunk

    @property
    def rulename(self):
        return "align_protein_{chunk}".format(chunk=self.chunk)

    @property
    def loader(self):
        return ['exonerate']

    @property
    def cmd(self):

        load = self.load
        cmd = "{load}"
        fasta = self.input["fasta"]
        cmd += " exonerate --model protein2genome --showtargetgff yes --showvulgar yes "
        min_intron, max_intron = self.min_intron, self.max_intron
        cmd += " --softmaskquery yes --softmasktarget yes --bestn 10  --minintron {min_intron} "
        cmd += " --maxintron {max_intron} --percent 30 --score 50 --geneseed 50 --showalignment no "
        cmd += " --query {input[fasta]} --target {input[genome]} "
        input, output, log = self.input, self.output, self.log
        cmd += " -ryo '>%qi\\tlength=%ql\\talnlen=%qal\\tscore=%s\\tpercentage=%pi\\nTarget>%ti\\tlength=%tl\\talnlen=%tal\\n' "
        cmd += " > {output[txt]} 2> {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")


class CollapseExonerate(AtomicOperation):

    def __init__(self, runs: [Exonerate]):
        super().__init__()
        self.configuration = runs[0].configuration
        self.input["chunks"] = [run.output["txt"] for run in runs]
        self.output["collapsed"] = os.path.join(self.outdir, "exonerate.txt")

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "collapse_exonerate"

    @property
    def cmd(self):
        inputs = " ".join(self.input["chunks"])
        output = self.output
        return "cat {inputs} > {output[collapsed]}".format(**locals())

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")


class ConvertExonerate(AtomicOperation):

    def __init__(self, collapsed: CollapseExonerate, fai: FaidxProtein):

        super().__init__()
        self.configuration = collapsed.configuration
        self.input = collapsed.output
        self.input["proteins"] = fai.input["db"]
        self.input["fai"] = fai.output['fai']
        self.output = {'gff3': os.path.join(self.outdir, "exonerate.gff3")}
        self.log = os.path.join(self.outdir, "convert.log")

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "exonerate_to_gff3"

    @property
    def identity(self):
        return self.configuration["homology"]["identity"]

    @property
    def coverage(self):
        return self.configuration["homology"]["coverage"]

    @property
    def cmd(self):

        identity, coverage = self.identity, self.coverage
        input, output, log = self.input, self.output, self.log

        cmd = "exonerate2gff.pl --in {input[collapsed]} --minIdentity {identity} --minCoverage {coverage} "
        cmd += "--fasta {input[proteins]} > {output[gff3]} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")