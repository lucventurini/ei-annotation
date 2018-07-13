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


