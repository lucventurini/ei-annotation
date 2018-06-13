import os
import abc
from ... import AtomicOperation, EIWrapper
from .prepare import MikadoPrepare
from ...preparation import DiamondIndex, BlastxIndex, SanitizeProteinBlastDB


class SplitMikadoFasta(AtomicOperation):

    def __init__(self, prepare: MikadoPrepare):
        super().__init__()
        self.configuration = prepare.configuration
        self.input = prepare.output
        self.outdir = prepare.outdir
        self.output = {"split": os.path.join(self.outdir, 'homology', 'split.done')}

    @property
    def message(self):
        return "Splitting fasta: {input[fa]}".format(input=self.input)

    @property
    def chunks(self):
        return self.configuration["mikado"]["blastx"]["chunks"]

    @property
    def loader(self):
        return ["mikado"]

    @property
    def rulename(self):
        return "mikado_split_fa"

    @property
    def log(self):
        return os.path.join(os.path.dirname(self.output["split"]), "split.log")

    @property
    def outprefix(self):
        return os.path.join(os.path.dirname(self.output["split"]), "fastas", "chunk")

    @property
    def cmd(self):
        load = self.load
        log = self.log
        input, output = self.input, self.output
        chunks = self.chunks
        outprefix = self.outprefix

        cmd = "{load} split_fasta.py -m {chunks} {input[fa]} {outprefix} && touch -h {output[split]}".format(
            **locals()
        )

        cmd = cmd.format(**locals())
        return cmd

class MikadoHomology(AtomicOperation, metaclass=abc.ABCMeta):

    def __init__(self,
                 chunk_id,
                 split: SplitMikadoFasta,
                 index: AtomicOperation): # We are going to create a Diamond/Blastx index somewhere else

        if "db" not in index.output:
            raise KeyError("Missing input from index rule {rule}".format(rule=index.rulename))

        super().__init__()
        self.configuration = split.configuration
        self.outdir = split.outdir
        self.input = split.output
        self.input["db"] = index.output["db"]
        self.__chunk_id = chunk_id
        chunk_id = self.chunk_id
        self.input["chunk"] = os.path.join(self.outdir,
                                           "homology",
                                           "fastas",
                                           "chunk_{chunk_id}.fasta".format(**locals()))
        self.output = {"xml": os.path.join(self.outdir,
                                          "homology",
                                          "xmls",
                                          "chunk_{chunk_id}.xml.gz".format(**locals()))}

    @property
    def chunk_id(self):
        return str(self.__chunk_id).zfill(3)

    @property
    def max_target_seqs(self):
        return self.configuration["homology"]["max_target_seqs"]

    @property
    def evalue(self):
        return self.configuration["homology"]["evalue"]

    @property
    def log(self):
        logdir = os.path.join(self.outdir, "homology", "logs")
        os.makedirs(logdir)
        name = self.loader[0]
        return os.path.join(logdir, "{name}_{chunk_id}.log".format(chunk_id=self.chunk_id, name=name))

    @property
    def rulename(self):
        name = self.loader[0]
        return "mikado_{name}_chunk_{chunk_id}".format(chunk_id=self.chunk_id, name=name)


class MikadoDiamond(MikadoHomology):

    def __init__(self,
                 chunk_id: int,
                 split: SplitMikadoFasta,
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
                 split: SplitMikadoFasta,
                 blast_index: BlastxIndex):
        super().__init__(chunk_id, split, blast_index)

    @property
    def rulename(self):
        return "mikado_blast_chunk_{}".format(self.chunk_id)

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


class MikadoHomologyFlag(AtomicOperation):

    def __init__(self, homologies: [MikadoHomology]):

        assert len(homologies) and all(isinstance(homology, MikadoHomology) for homology in homologies)
        super().__init__()
        self.configuration = homologies[0].configuration
        self.outdir = homologies[0].outdir
        self.touch = True
        self.input = {"xmls": [homology.output["xml"] for homology in homologies]}
        self.output = {"flag": os.path.join(self.outdir, "homology", "homology.done")}

    @property
    def loader(self):
        return []

    @property
    def rulename(self):
        return "homology_all"


class MikadoHomologyWrapper(EIWrapper):

    def __init__(self, preparer: MikadoPrepare):

        super().__init__()
        self.configuration = preparer.configuration
        self.sanitizer = SanitizeProteinBlastDB(self.configuration)
        if self.execute is True:
            split = SplitMikadoFasta(preparer)
            if self.program == "blastx":
                indexer, executer = BlastxIndex, MikadoBlastx
            elif self.program == "diamond":
                indexer, executer = DiamondIndex, MikadoDiamond
            else:
                raise KeyError("Unrecognised homology assessment program")

            indexer = indexer(self.sanitizer)
            self.add_edge(self.sanitizer, indexer)
            executers = []
            for chunk in range(1, self.chunks + 1):
                exe = executer(chunk, split, indexer)
                self.add_edge(indexer, exe)
                executers.append(exe)

            self.flag = MikadoHomologyFlag(executers)
            self.add_edges_from([(exe, self.flag) for exe in executers])

    @property
    def chunks(self):
        return self.configuration["homology"]["chunks"]

    @property
    def program(self):
        return self.configuration["homology"]["program"]

    @property
    def execute(self):
        return self.configuration["homology"]["execute"] and len(self.sanitizer.protein_dbs) > 0

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
