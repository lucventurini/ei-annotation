from .. import AtomicOperation, EIWrapper


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
