from ..abstract import AtomicOperation, EIWrapper
import abc
import os
from ..repeats import RepeatMasking
from ..preparation import SanitizeProteinBlastDB
from ..rnaseq.alignments.portcullis import PortcullisWrapper


def _get_value(conf, dbname, value):
    if value in conf["homology"]["prot_dbs"].get(dbname, {}):
        return conf["homology"]["prot_dbs"][dbname][value]
    elif value in conf["homology"]:
        return conf["homology"][value]
    else:
        return None


class ChunkProteins(AtomicOperation):

    def __init__(self, sanitised: SanitizeProteinBlastDB):

        super().__init__()
        self.configuration = sanitised.configuration
        self.input = sanitised.output  # input[db] is our file
        self.dbname = sanitised.dbname
        self.output["flag"] = os.path.join(os.path.dirname(self.outdir),
                                           "{dbname}_chunking.done".format(dbname=self.dbname))
        self.output["chunks"] = [os.path.join(self.outdir, "{dbname}_{cid}.fasta".format(
            cid=str(chunk).zfill(3), dbname=sanitised.dbname)) for chunk in range(1, self.chunks + 1)]
        self.log = os.path.join(os.path.dirname(self.outdir), "logs",
                                "{dbname}_chunking.log".format(dbname=self.dbname))

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"],
                            "proteins",
                            "chunks")

    @property
    def rulename(self):
        return "chunk_proteins_for_alignment_{dbname}".format(dbname=self.dbname)

    @property
    def loader(self):
        return ["mikado", "genometools"]

    @property
    def cmd(self):
        load = self.load

        outdir = self.outdir
        chunks = self.chunks
        input, output, log = self.input, self.output, self.log
        logdir = os.path.dirname(self.log)
        dbname = self.dbname
        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && "
        cmd += " split_fasta.py -m {chunks} {input[db]} {outdir}/{dbname} 2> {log} > {log}"
        cmd += " && touch {output[flag]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def chunks(self):
        if not _get_value(self.configuration, self.dbname, "chunks"):
            return self.configuration["homology"]["protein_chunks"]
        else:
            return _get_value(self.configuration, self.dbname, "chunks")

    @property
    def is_small(self):
        return True



class ProteinChunkAligner(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self,
                 chunks: ChunkProteins,
                 chunk,
                 masked: RepeatMasking):
        super().__init__()
        self.dbname = chunks.dbname
        self._chunks = chunks
        self.configuration = chunks.configuration
        self.__chunk = None
        self.chunk = chunk
        self.input["flag"] = chunks.output["flag"]
        self.input["genome"] = self.masked_genome
        self.input["repeat_flag"] = masked.exit.output["flag"]
        assert self.input["fasta"] in chunks.output["chunks"]

    @property
    def rulename(self):
        return "{toolname}_{dbname}_{chunk}".format(toolname=self.toolname,
                                                    dbname=self.dbname,
                                                    chunk=str(self.chunk).zfill(3))

    @property
    def _coverage_value(self):
        return _get_value(self.configuration, self.dbname, "coverage")

    @property
    def _identity_value(self):
        return _get_value(self.configuration, self.dbname, "identity")

    @property
    def toolname(self):
        return self.__toolname__

    @property
    def chunk(self):
        return self.__chunk

    @chunk.setter
    def chunk(self, chunk):
        assert isinstance(chunk, int)
        fasta = os.path.join(self._chunks.outdir, "{dbname}_{chunk}.fasta".format(
            dbname=self.dbname, chunk=str(chunk).zfill(3)))
        assert fasta in self._chunks.output["chunks"], (fasta, self._chunks.output["chunks"])
        self.input["fasta"] = fasta
        self.__chunk = chunk

    @property
    @abc.abstractmethod
    def cmd(self):
        pass

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")

    @property
    def logdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "logs")


class FilterAlignments(AtomicOperation):

    __rulename__ = "filter_alignments"

    def __init__(self,
                 collapse,  #: CollapseGTH,
                 portcullis: PortcullisWrapper,
                 masker: RepeatMasking):

        super().__init__()
        self.dbname = collapse.dbname
        self.configuration = collapse.configuration
        self.input = collapse.output
        if portcullis.merger.input["beds"]:
            self.input["junctions"] = portcullis.junctions
        self.input["fai"] = masker.fai
        self.input["genome"] = self.masked_genome
        self.log = os.path.join(os.path.dirname(self.outdir), "logs",
                                "filter_alignments_{dbname}.log".format(dbname=self.dbname))
        self.output["gff3"] = os.path.join(self.outdir, "{dbname}.filtered.gff3".format(dbname=self.dbname))

    @property
    def rulename(self):
        return self.__rulename__ + "_{dbname}".format(dbname=self.dbname)

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]  # TODO: this will need to be updated with the proper details

    @property
    def min_intron(self):
        return _get_value(self.configuration, self.dbname, "min_intron")

    @property
    def max_intron_middle(self):
        return _get_value(self.configuration, self.dbname, "max_intron_middle")

    @property
    def max_intron_ends(self):
        return _get_value(self.configuration, self.dbname, "max_intron_ends")

    @property
    def coverage(self):
        cov = _get_value(self.configuration, self.dbname, "coverage")
        if cov is not None:
            return cov
        else:
            return 0

    @property
    def identity(self):
        iden = _get_value(self.configuration, self.dbname, "identity")
        if iden is not None:
            return iden
        else:
            return 0

    @property
    def cmd_coverage(self):
        if self.coverage:
            return " -mincov {} ".format(self.coverage)
        else:
            return ""

    @property
    def cmd_identity(self):
        if self.identity:
            return " -minid {} ".format(self.identity)
        else:
            return ""

    @property
    def threads(self):
        return 1

    @property
    def source(self):
        identity = self.identity or 0
        coverage = self.coverage or 0
        return "{dbname}.{cov}cov{iden}id".format(dbname=self.dbname,
                                                  cov=int(round(coverage)),
                                                  iden=int(round(identity)))

    @property
    def cmd(self):

        load = self.load
        mini = self.min_intron
        maxe, maxm = self.max_intron_ends, self.max_intron_middle
        genome = self.masked_genome
        outdir = self.outdir
        logdir = os.path.dirname(self.log)
        min_coverage, min_identity = self.cmd_coverage, self.cmd_identity
        input, output, log = self.input, self.output, self.log
        source = self.source
        if "junctions" in self.input:
            junctions = "-j {input[junctions]}".format(**locals())
        else:
            junctions = ""

        cmd = "{load} mkdir -p {outdir} && mkdir -p {logdir} && "
        cmd += " filter_exonerate.py -s {source} -minI {mini} -maxE {maxe} -maxM {maxm} {junctions} -g {genome} "
        cmd += " {min_identity} {min_coverage} {input[gff3]} {output[gff3]} 2> {log} > {log}"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "output")

    @property
    def is_small(self):
        return True


class ProteinWrapper(EIWrapper, metaclass=abc.ABCMeta):

    __final_rulename__ = "proteins_done"

    def __init__(self, masker: RepeatMasking, portcullis: PortcullisWrapper):
        super().__init__(configuration=masker.configuration)
        self.masker = masker
        self.portcullis = portcullis
        self.__proteins = []
        if self.execute:
            for db in self.dbs:
                self.execute_protein(db)
        self.add_final_flag()
        assert self.exit

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "proteins.done")

    @property
    def execute(self):
        exe = self.configuration["homology"].get("execute", True)
        return exe and self.dbs

    @property
    def dbs(self):
        return self.configuration["homology"]["prot_dbs"]

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "output")

    @property
    def use_exonerate(self):
        return self.configuration["homology"].get("use_exonerate", False)

    @abc.abstractmethod
    def execute_protein(self, db):
        pass

    @property
    def proteins(self):
        return self.__proteins

    def _add_to_proteins(self, protein: FilterAlignments):
        assert isinstance(protein, FilterAlignments)
        self.__proteins.append(protein)
