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
    def threads(self):
        return 1

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

    def __init__(self, sanitiser, configuration=None):

        super().__init__()
        if sanitiser is None:
            self.configuration = configuration
            self.input["genome"] = self.genome

        else:
            self.input = sanitiser.output
            self.configuration = sanitiser.configuration

        self.output["fai"] = self.input["genome"] + ".fai"
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
        cmd = "{load} samtools faidx {input[genome]} > {log} 2>&1".format(load=self.load, input=self.input,
                                                                          log=self.log)
        return cmd


class SanitizeProteinBlastDB(AtomicOperation):

    def __init__(self, configuration, db=None, dbs=None):
        super().__init__()
        self.__dbname = db
        self.__dbs = dbs
        self.configuration = configuration
        self.input["db"] = self.protein_dbs
        self.output["db"] = os.path.join(self.outdir, "{dbname}.fa".format(dbname=self.dbname))
        self.log = os.path.join(os.path.dirname(self.output["db"]), "sanitize.log")

    @property
    def dbname(self):
        if self.__dbname is not None:
            return self.__dbname
        else:
            return "homologyDB"

    @dbname.setter
    def dbname(self, dbname):
        if not isinstance(dbname, str):
            raise TypeError(type(dbname))
        self.__dbname = dbname

    @property
    def loader(self):
        return ["ei-annotation", "genometools"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        log = self.log
        dbs = " ".join(self.input["db"])

        cmd = "{load} sanitize_blast_db.py {dbs} 2> {log} |"
        cmd += " gt seqtransform -addstopaminos -width 60 > {output[db]} 2> {log}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def outdir(self):

        return os.path.join(self.configuration["outdir"], "inputs", "proteins")

    @property
    def protein_dbs(self):
        if self.__dbs is not None:
            return self.__dbs
        else:
            assert self.dbname in self.configuration["homology"]["prot_db"]
            assert 'fastas' in self.configuration["homology"]["prot_db"][self.dbname]
            return self.configuration["homology"]["prot_db"][self.dbname]['fastas']

    @property
    def rulename(self):
        return "sanitize_{dbname}_protein_db".format(dbname=self.dbname)

    @property
    def threads(self):
        return 1


class FaidxProtein(AtomicOperation):

    def __init__(self, sanitised: SanitizeProteinBlastDB):

        super().__init__()
        self.configuration = sanitised.configuration
        self.dbname = sanitised.dbname
        self.input = sanitised.output
        self.output["fai"] = sanitised.output["db"] + ".fai"

    @property
    def loader(self):
        return ["samtools"]

    @property
    def threads(self):
        return 1

    @property
    def rulename(self):
        return "faidx_proteins_{dbname}".format(dbname=self.dbname)

    @property
    def cmd(self):
        load = self.load
        input = self.input

        return "{load} samtools faidx {input[db]}".format(**locals())


class PrepareWrapper(EIWrapper):

    __final_rulename__ = "prepare_flag"

    def __init__(self, configuration, genome):
        super().__init__(configuration)
        self.sanitizer = SanitizeGenome(configuration, genome)
        self.faidx = FaidxGenome(self.sanitizer)
        self.add_edge(self.sanitizer, self.faidx)

    @property
    def fai(self):
        return self.faidx

    @property
    def flag_name(self):
        return os.path.join(os.path.dirname(self.faidx.output["fai"]), "prepare.done")


class DiamondIndex(AtomicOperation):

    def __init__(self, sanitizer: SanitizeProteinBlastDB):

        super().__init__()
        self.input = sanitizer.output
        self.configuration = sanitizer.configuration
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
        self.configuration = sanitizer.configuration
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
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        db = os.path.splitext(self.output["db"])[0]
        log, input = self.log, self.input
        cmd = "{load} makeblastb -in {input[db]} -out {db} -dbtype prot -parse_seqids > {log} 2>&1".format(**locals())
        return cmd
