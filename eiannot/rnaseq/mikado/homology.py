import os
import abc
from ...abstract import AtomicOperation, EIWrapper
from .prepare import MikadoPrepare
from ...preparation import DiamondIndex, BlastxIndex, SanitizeProteinBlastDB
from .abstract import MikadoOp


class SplitMikadoPrepareFasta(MikadoOp):

    def __init__(self, prepare: MikadoPrepare):
        super().__init__(is_long=prepare.is_long)
        self.configuration = prepare.configuration
        self.input = prepare.output
        self.outdir = self.mikado_dir
        self.output["split_flag"] = os.path.join(self.outdir, 'homology', 'split.done')
        self.output["split"] = ["{outprefix}_{chunk}.fasta".format(
            outprefix=self.outprefix, chunk=str(_).zfill(3)) for _ in range(1, self.chunks + 1)]

    @property
    def message(self):
        return "Splitting fasta: {input[fa]}".format(input=self.input)

    @property
    def chunks(self):
        return self.configuration["mikado_homology"]["chunks"]

    @property
    def loader(self):
        return ["mikado"]

    @property
    def _rulename(self):
        return "mikado_prepare_split_fa"

    @property
    def log(self):
        return os.path.join(os.path.dirname(self.output["split_flag"]), "split.log")

    @property
    def outprefix(self):
        return os.path.join(self.fasta_dir, "chunk")

    @property
    def cmd(self):
        load = self.load
        log = self.log
        input, output = self.input, self.output
        chunks = self.chunks
        outprefix = self.outprefix

        cmd = "{load} split_fasta.py -m {chunks} {input[fa]} {outprefix} && touch {output[split_flag]}".format(
            **locals()
        )

        cmd = cmd.format(**locals())
        return cmd

    @property
    def fasta_dir(self):
        return os.path.join(self.mikado_dir, "mikado_homology", "fastas")


class MikadoHomology(MikadoOp, metaclass=abc.ABCMeta):

    def __init__(self,
                 chunk_id,
                 split: SplitMikadoPrepareFasta,
                 index: AtomicOperation):  # We are going to create a Diamond/Blastx index somewhere else

        if "db" not in index.output:
            raise KeyError("Missing input from index rule {rule}".format(rule=index.rulename))

        super().__init__(is_long=split.is_long)
        self.configuration = split.configuration
        self.input = split.output
        self.input["db"] = index.output["db"]
        self.__chunk_id = chunk_id
        chunk_id = self.chunk_id
        del self.input["split"]
        self.input["chunk"] = os.path.join(self.fasta_dir,
                                           "chunk_{chunk_id}.fasta".format(**locals()))
        self.output = {"xml": os.path.join(self.outdir,
                                           "xmls",
                                           "chunk_{chunk_id}.xml.gz".format(**locals()))}

    @property
    def chunk_id(self):
        return str(self.__chunk_id).zfill(3)

    @property
    def max_target_seqs(self):
        return self.configuration["mikado_homology"]["max_target_seqs"]

    @property
    def evalue(self):
        return self.configuration["mikado_homology"]["evalue"]

    @property
    def log(self):
        logdir = os.path.join(self.outdir, "mikado_homology", "logs")
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        name = self.loader[0]
        return os.path.join(logdir, "{name}_{chunk_id}.log".format(chunk_id=self.chunk_id, name=name))

    @property
    def _rulename(self):
        name = self.loader[0]
        return "mikado_{name}_chunk_{chunk_id}".format(chunk_id=self.chunk_id, name=name)

    @property
    def outdir(self):
        return os.path.join(self.mikado_dir, "mikado_homology")

    @property
    def fasta_dir(self):
        return os.path.join(self.mikado_dir, "mikado_homology", "fastas")


class MikadoDiamond(MikadoHomology):

    def __init__(self,
                 chunk_id: int,
                 split: SplitMikadoPrepareFasta,
                 diamond_index: DiamondIndex):

        super().__init__(chunk_id, split, diamond_index)

    @property
    def loader(self):
        return ["diamond"]

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        log = self.log
        input, output = self.input, self.output

        cmd = "{load} "
        cmd += "diamond blastx --threads {threads} --outfmt xml --compress 1"
        cmd += " --out {output[xml]} --db {input[db]} --salltitles --query {input[chunk]} --sensitive "
        evalue, max_targets = self.evalue, self.max_target_seqs
        cmd += " --max-target-seqs {max_targets} --evalue {evalue} > {log} 2>&1"

        cmd = cmd.format(**locals())
        return cmd


class MikadoBlastx(MikadoHomology):

    def __init__(self,
                 chunk_id: int,
                 split: SplitMikadoPrepareFasta,
                 blast_index: BlastxIndex):
        super().__init__(chunk_id, split, blast_index)

    @property
    def loader(self):
        return ["blast"]

    @property
    def cmd(self):

        load = self.load
        log = self.log
        cmd = "{load} "
        # uncompressed = os.path.splitext(self.output["xml"])[0]
        threads = self.threads
        db = os.path.splitext(self.input["db"])
        input, output = self.input, self.output
        cmd += " blastx -num_threads {threads} -outfmt 5 -query {input[chunk]} "
        cmd += " -db {db} "
        evalue, max_targets = self.evalue, self.max_target_seqs
        cmd += " -evalue {evalue} -max_target_seqs {max_targets} 2> {log} | gzip -c - > {output[xml]}"
        cmd = cmd.format(**locals())
        return cmd


class MikadoHomologyFlag(MikadoOp):

    def __init__(self, homologies: [MikadoHomology]):

        assert len(homologies) and all(isinstance(homology, MikadoHomology) for homology in homologies)
        super().__init__(is_long=homologies[0].is_long)
        self.configuration = homologies[0].configuration
        self.outdir = homologies[0].outdir
        self.touch = True
        self.input = {"xmls": [homology.output["xml"] for homology in homologies]}
        self.output = {"flag": os.path.join(self.outdir, "mikado_homology", "homology.done")}

    @property
    def loader(self):
        return []

    @property
    def _rulename(self):
        return "homology_all"


class MikadoHomologyWrapper(EIWrapper):

    __final_rulename__ = "mikado_homology_flag"

    def __init__(self, preparer: MikadoPrepare):

        super().__init__()
        self.configuration = preparer.configuration
        self.outdir = preparer.outdir
        self.sanitizer = SanitizeProteinBlastDB(self.configuration, dbs=self.dbs)
        if self.execute is True:
            split = SplitMikadoPrepareFasta(preparer)
            if self.program == "blastx":
                indexer, executer = BlastxIndex, MikadoBlastx
            elif self.program == "diamond":
                indexer, executer = DiamondIndex, MikadoDiamond
            else:
                raise KeyError("Unrecognised homology assessment program")

            indexer = indexer(self.sanitizer)
            assert indexer.configuration["programs"]
            self.add_edge(self.sanitizer, indexer)
            executers = []
            for chunk in range(1, self.chunks + 1):
                exe = executer(chunk, split, indexer)
                self.add_edge(indexer, exe)
                self.add_edge(split, exe)
                executers.append(exe)

            self.flag = MikadoHomologyFlag(executers)
            self.add_edges_from([(exe, self.flag) for exe in executers])
        if preparer.is_long:
            self.__final_rulename__ += "_long"

    @property
    def dbs(self):
        return self.configuration["mikado_homology"]["prot_dbs"]

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "homology.done")

    @property
    def chunks(self):
        return self.configuration["mikado_homology"]["chunks"]

    @property
    def program(self):
        return self.configuration["mikado_homology"]["program"]

    @property
    def execute(self):
        return self.configuration["mikado_homology"]["execute"] and len(self.sanitizer.protein_dbs) > 0

    @property
    def blast_xmls(self):
        if self.execute:
            return self.flag.input["xmls"]
        else:
            return []

    @property
    def blast_targets(self):
        if self.execute:
            return self.sanitizer.output["db"]
        else:
            return None
