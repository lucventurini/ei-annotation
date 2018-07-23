from ..abstract import AtomicOperation
from ..preparation import SanitizeProteinBlastDB
import os


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
