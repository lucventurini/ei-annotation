from ..abstract import AtomicOperation, EIWrapper
import os


class SanitizeGenome(AtomicOperation):

    def __init__(self, configuration, genome):

        super().__init__()
        self.configuration = configuration
        self.input["genome"] = genome
        self.output["genome"] = self.genome
        self.log = os.path.join(os.path.dirname(self.genome), 'sanitize.log')
        self.message = 'Sanitizing the genome FASTA file before any other operation.'

    @property
    def loader(self):
        return ["ei-annotation"]

    @property
    def rulename(self):
        return "sanitize_genome"

    @property
    def cmd(self):

        load = self.load
        input, output = self.input, self.output
        log = self.log

        cmd = "{load} sanitize_blast_db.py -o {output[genome]} {input[genome]} 2> {log} > {log}"
        cmd = cmd.format(**locals())
        return cmd


class FaidxGenome(AtomicOperation):

    def __init__(self, sanitised: SanitizeGenome):

        super().__init__()
        self.input = sanitised.output
        self.output["fai"] = self.input["genome"] + ".fai"
        self.configuration = sanitised.configuration
        self.log = os.path.join(os.path.dirname(self.input["genome"]), "faidx.log")
        self.message = "Indexing the sanitised genome with SAMtools"

    @property
    def loader(self):
        return ["samtools"]

    @property
    def rulename(self):
        return "faidx_genome"

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        cmd = "{load} samtools index {input[genome]} > {log} 2>&1".format(load=self.load, input=self.input,
                                                                          log=self.log)
        return cmd


class SanitizeProteinBlastDB(AtomicOperation):

    def __init__(self, configuration):
        super().__init__()
        self.configuration = configuration
        self.input["db"] = self.protein_dbs
        self.output["db"] = os.path.join(self.outdir, "homologyDB.fa")
        self.log = ""

    @property
    def loader(self):
        return ["ei-annotation"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        log = self.log

        cmd = "{load} sanitize_blast_db.py -o {output[db]} {input[db]} 2> {log} > {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "inputs", "proteins")

    @property
    def protein_dbs(self):
        return self.configuration["homology"]["prot_db"]

    @property
    def rulename(self):
        return "sanitize_protein_db"


class DiamondIndex(AtomicOperation):

    def __init__(self, sanitizer: SanitizeProteinBlastDB):

        super().__init__()
        self.input = sanitizer.output
        self.output["db"] = os.path.splitext(self.input["db"])[0] + ".dmnd"
        self.outdir = os.path.dirname(self.output["db"])
        self.log = os.path.join(self.outdir, "diamond.index.log")
        self.message = "Making DIAMOND protein database for: {input[db]}".format(input=self.input)

    @property
    def loader(self):
        return ["diamond"]

    @property
    def rulename(self):
        return "diamond_protein_index"

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        db = os.path.splitext(self.output["db"])[0]
        input, log = self.input, self.log
        cmd = "{load} diamond makedb --threads {threads} --in {input[db]} --db {db} > {log} 2>{log}".format(**locals())
        return cmd


class BlastxIndex(AtomicOperation):

    def __init__(self, sanitizer: SanitizeProteinBlastDB):
        super().__init__()
        self.input = sanitizer.output
        self.output["db"] = os.path.splitext(self.input["db"])[0] + ".pog"
        self.outdir = os.path.dirname(self.output["db"])
        self.log = os.path.join(self.outdir, "blast.index.log")
        self.message = "Making BLASTX protein database for: {input[db]}".format(input=self.input)

    @property
    def loader(self):
        return ["blast"]

    @property
    def rulename(self):
        return "protein_blast_index"

    @property
    def cmd(self):
        load = self.load
        db = os.path.splitext(self.output["db"])[0]
        log, input = self.log, self.input
        cmd = "{load} makeblastb -in {input[db]} -out {db} -dbtype prot -parse_seqids > {log} 2>&1".format(**locals())
        return cmd
