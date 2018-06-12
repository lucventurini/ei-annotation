import os
import abc
from ... import AtomicOperation
from .prepare import SplitMikadoFasta


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
                 diamond_index: AtomicOperation):

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
                 blast_index: AtomicOperation  # This will have to be the BLASTX index
                 ):
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