from ..abstract import AtomicOperation
from .abstract import ProteinChunkAligner, FilterAlignments, ProteinWrapper, _get_value
import functools
from .chunking import ChunkProteins
from ..repeats import RepeatMasking
import re
import subprocess as sp
import os
from ..preparation import FaidxProtein, SanitizeProteinBlastDB
from ..rnaseq.alignments.portcullis import PortcullisWrapper


@functools.lru_cache(maxsize=4, typed=True)
def exonerate_multithread(loader):
    cmd = "{} exonerate --help && set -u".format(loader)
    output = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout.read().decode()
    try:
        _ = [_ for _ in output.split("\n") if re.match("--cores", _)][0]
    except IndexError:
        return False
    return True


class Fasta2Esi(AtomicOperation):

    def __init__(self, masked: RepeatMasking):

        super().__init__()
        self.masked = masked
        self.configuration = masked.configuration
        self.input.update(masked.output)
        self.input["genome"] = self.masked_genome
        self.output["genome"] = os.path.join(self.outdir, "genome.masked.esi")
        self.log = os.path.join(self.outdir, "fasta2esi.log")

    @property
    def loader(self):
        return ["exonerate"]

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "indices", "exonerate")

    @property
    def cmd(self):

        esd = os.path.join(os.path.dirname(self.output["genome"]),
                           "genome.masked.esd")
        load = self.load
        outdir = self.outdir
        input, output, log = self.input, self.output, self.log
        memory = self.resources["memory"]
        cmd = "{load} mkdir -p {outdir} && fasta2esd --softmask yes {input[genome]} {esd} 2> {log} > {log} && "
        cmd += " esd2esi --translate yes --memorylimit {memory}  {esd} {output[genome]} >> {log} 2>> {log}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def rulename(self):
        return "fasta_to_esi"


class Exonerate(ProteinChunkAligner):

    __toolname__ = "exonerate"

    def __init__(self,
                 chunks: ChunkProteins,
                 chunk,
                 esi: Fasta2Esi
                 ):
        super().__init__(chunks, chunk, esi.masked)
        self.input["genome"] = esi.output["genome"]
        self.output["txt"] = os.path.join(self.outdir, "{dbname}_{chunk}.exonerate.txt").format(
            chunk=self.chunk, dbname=self.dbname)
        self.log = os.path.join(self.logdir, "exonerate.{dbname}_{chunk}.log".format(
            dbname=self.dbname, chunk=self.chunk))

    @property
    def loader(self):
        return ['exonerate', "ei-annotation"]

    @property
    def min_intron_cli(self):
        if _get_value(self.configuration, self.dbname, "min_intron"):
            return _get_value(self.configuration, self.dbname, "min_intron")
        else:
            return self.min_intron

    @property
    def max_intron_cli(self):
        if _get_value(self.configuration, self.dbname, "max_intron_middle"):
            return _get_value(self.configuration, self.dbname, "max_intron_middle")
        else:
            return self.max_intron

    @property
    def geneseed(self):
        return _get_value(self.configuration, self.dbname, "geneseed") or 250

    @property
    def cmd(self):

        load = self.load
        outdir, logdir = self.outdir, self.logdir
        cmd = "{load} mkdir -p {logdir} && mkdir -p {outdir} && set +u && "
        # Now get the port and get the server
        min_intron, max_intron = self.min_intron_cli, self.max_intron_cli
        identity = self._identity_value
        threads = self.threads
        cmd += " exonerate_wrapper.py -ir {min_intron} {max_intron} --identity {identity} -t {threads}"
        memory = max(self.resources["memory"] * 0.95, self.resources["memory"] - 500)   # Do not give all memory
        geneseed = self.geneseed
        cmd += " -M {memory} --geneseed {geneseed} "
        input, output, log = self.input, self.output, self.log
        server_log = os.path.join(self.logdir,
                                  "exonerate.{dbname}_{chunk}.server.log".format(
                                      dbname=self.dbname, chunk=self.chunk)
                                  )
        cmd += "  --serverlog {server_log} --log {log} {input[genome]} {input[fasta]} {output[txt]} >> {log} 2>> {log}"
        cmd = cmd.format(**locals())
        return cmd


class CollapseExonerate(AtomicOperation):

    def __init__(self, runs: [Exonerate]):
        super().__init__()
        assert len(set(_.dbname for _ in runs)) == 1
        self.dbname = runs[0].dbname
        self.configuration = runs[0].configuration
        self.input["chunks"] = [run.output["txt"] for run in runs]
        self.output["collapsed"] = os.path.join(self.outdir, "{dbname}.exonerate.txt".format(dbname=self.dbname))

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "collapse_exonerate_{dbname}".format(dbname=self.dbname)

    @property
    def cmd(self):
        inputs = " ".join(self.input["chunks"])
        output = self.output
        outdir = self.outdir
        return "mkdir -p {outdir} && cat {inputs} > {output[collapsed]}".format(**locals())

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")


class ConvertExonerate(AtomicOperation):

    def __init__(self, collapsed: CollapseExonerate, fai: FaidxProtein):

        super().__init__()
        self.dbname = collapsed.dbname
        assert collapsed.dbname == fai.dbname  # Sanity check
        self.configuration = collapsed.configuration
        self.input = collapsed.output
        self.input["proteins"] = fai.input["db"]
        self.input["fai"] = fai.output['fai']
        self.output = {'gff3': os.path.join(self.outdir, "{dbname}.exonerate.gff3".format(dbname=self.dbname))}
        self.log = os.path.join(self.outdir, "{dbname}.convert.log".format(dbname=self.dbname))

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "exonerate_to_gff3_{}".format(self.dbname)

    @property
    def _identity_value(self):
        return _get_value(self.configuration, self.dbname, "identity")

    @property
    def _coverage_value(self):
        return _get_value(self.configuration, self.dbname, "coverage")

    @property
    def coverage(self):
        if self._coverage_value:
            return "--minCoverage {coverage}".format(coverage=self._coverage_value)
        else:
            return ""

    @property
    def identity(self):
        if self._identity_value:
            return "--minIdentity {identity}".format(identity=self._identity_value)
        else:
            return ""

    @property
    def cmd(self):

        identity, coverage = self.identity, self.coverage
        input, output, log = self.input, self.output, self.log

        cmd = "exonerate2gff.pl --in {input[collapsed]} {identity} {coverage} "
        cmd += " --fasta {input[proteins]} > {output[gff3]} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], "proteins", "alignments")


class ExonerateProteinWrapper(ProteinWrapper):

    def __init__(self,
                 masker: RepeatMasking, portcullis: PortcullisWrapper):

        super().__init__(masker, portcullis)

    def execute_protein(self, db):
        sanitised = SanitizeProteinBlastDB(self.configuration, db)
        faidx = FaidxProtein(sanitised)
        self.add_edge(sanitised, faidx)
        chunk_proteins = ChunkProteins(sanitised)
        fasta2esi = Fasta2Esi(masked=self.masker)
        self.add_edge(self.masker, fasta2esi)
        chunks = [Exonerate(chunk_proteins, chunk, fasta2esi) for chunk in range(1, chunk_proteins.chunks + 1)]
        self.add_edges_from([(chunk_proteins, chunk) for chunk in chunks])
        self.add_edges_from([(fasta2esi, chunk) for chunk in chunks])
        collapsed = CollapseExonerate(chunks)
        self.add_edges_from([(chunk, collapsed) for chunk in chunks])
        convert = ConvertExonerate(collapsed, faidx)
        self.add_edge(faidx, convert)
        self.add_edge(collapsed, convert)
        filterer = FilterAlignments(collapse=convert,
                                    portcullis=self.portcullis,
                                    masker=self.masker)
        self.add_edge(self.masker, filterer)
        self.add_edge(self.portcullis, filterer)
        self.add_edge(convert, filterer)
