from .. import AtomicOperation, EIWrapper
import os


class SanitizeGenome(AtomicOperation):

    def __init__(self, configuration, genome):

        super().__init__()
        self.configuration = configuration
        self.input["genome"] = genome
        self.output["genome"] = self.genome
        self.log = ""  # TODO: implement this

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


class SanitizeBlastDB(AtomicOperation):

    def __init__(self, configuration, db):
        super().__init__()
        self.configuration = configuration
        self.input["db"] = db
        self.output["db"] = ''  # TODO: implement this
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
