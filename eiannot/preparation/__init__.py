from ..abstract import AtomicOperation, EIWrapper, Linker
import os


class PrepareWrapper(EIWrapper):

    __final_rulename__ = "prepare_flag"

    def __init__(self, configuration, genome):
        super().__init__(configuration)
        self.sanitizer = SanitizeGenome(configuration, genome)
        self.faidx = FaidxGenome(self.sanitizer)
        self.add_edge(self.sanitizer, self.faidx)
        if self.sanitizer._raw_transcriptome is not None:
            if self.sanitizer._raw_transcriptome.endswith("gtf"):
                prepare_transcriptome = Linker(self.sanitizer._raw_transcriptome,
                                               self.sanitizer.transcriptome,
                                               "transcriptome",
                                               "transcriptome",
                                               "link_transcriptome_annotation",
                                               configuration
                                               )
            else:
                prepare_transcriptome = ConvertTranscriptome(configuration)
            self.add_node(prepare_transcriptome)
            self.add_edge(self.sanitizer, prepare_transcriptome)

        self.add_final_flag()

    @property
    def fai(self):
        return self.faidx

    @property
    def flag_name(self):
        return os.path.join(os.path.dirname(self.faidx.output["fai"]), "prepare.done")


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

        cmd = "{load} sanitize_sequence_db.py -o {output[genome]} {input[genome]} 2> {log} > {log}"
        cmd = cmd.format(**locals())
        return cmd


class ConvertTranscriptome(AtomicOperation):

    def __init__(self, configuration):

        super().__init__()
        self.configuration = configuration
        self.input["transcriptome"] = self._raw_transcriptome
        self.output["transcriptome"] = self.transcriptome

    @property
    def loader(self):
        return ["mikado"]

    @property
    def rulename(self):
        return "prepare_transcriptome_annotation"

    @property
    def cmd(self):

        load, input, output = self.load, self.input, self.output

        cmd = "{load} mikado util convert -of gtf {input[transcriptome]} {output[transcriptome]}"
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

    @property
    def is_small(self):
        return True


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
        out_tmp = os.path.join(os.path.dirname(output["db"]),
                               "{}.temp".format(os.path.basename(output["db"])))

        # We have to use a temporary file because otherwise gt will pollute the folder
        cmd = "{load} sanitize_sequence_db.py -cstop {dbs} 2> {log} > {out_tmp} && "
        cmd += " gt seqtransform -addstopaminos -width 60 {out_tmp} > {output[db]} 2> {log} && "
        cmd += " rm {out_tmp}*"
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
            assert self.dbname in self.configuration["homology"]["prot_dbs"]
            assert 'fastas' in self.configuration["homology"]["prot_dbs"][self.dbname]
            return self.configuration["homology"]["prot_dbs"][self.dbname]['fastas']

    @property
    def rulename(self):
        return "sanitize_{dbname}_protein_db".format(dbname=self.dbname)

    @property
    def threads(self):
        return 1

    @property
    def is_small(self):
        return True


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

    @property
    def is_small(self):
        return True


class DiamondIndex(AtomicOperation):

    def __init__(self, sanitizer, key="db", rule_suffix=""):

        super().__init__()
        self.__rule_suffix = ""
        self._rule_suffix = rule_suffix
        self.input = sanitizer.output
        assert key in sanitizer.output
        if key != "db":
            del self.input[key]  # Otherwise there are repeated inputs
        self.input["db"] = sanitizer.output[key]
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
        rule = "diamond_protein_index"
        if self._rule_suffix == '':
            return rule
        else:
            return rule + "_" + self._rule_suffix

    @property
    def _rule_suffix(self):
        return self.__rule_suffix

    @_rule_suffix.setter
    def _rule_suffix(self, rule):
        if rule is None:
            rule = ""
        if not isinstance(rule, str):
            raise TypeError()
        self.__rule_suffix = rule

    @property
    def cmd(self):
        load = self.load
        threads = self.threads
        db = os.path.splitext(self.output["db"])[0]
        input, log = self.input, self.log
        cmd = "{load} diamond makedb --threads {threads} --in {input[db]} --db {db} > {log} 2>{log}".format(**locals())
        return cmd


class BlastxIndex(AtomicOperation):

    def __init__(self, sanitizer: SanitizeProteinBlastDB, key="db", rule_suffix=""):
        super().__init__()
        self.__rule_suffix = ""
        self._rule_suffix = rule_suffix
        assert key in sanitizer.output  # Sanity check!
        self.input = sanitizer.output
        self.input["db"] = sanitizer.output[key]
        if key != "db":
            del self.input[key]  # Otherwise there are repeated inputs
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
        rule = "protein_blast_index_"
        if self._rule_suffix != '':
            return rule
        else:
            return rule + "_" + self._rule_suffix

    @property
    def _rule_suffix(self):
        return self.__rule_suffix

    @_rule_suffix.setter
    def _rule_suffix(self, rule):
        if rule is None:
            rule = ""
        if not isinstance(rule, str):
            raise TypeError()
        self.__rule_suffix = rule

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        db = os.path.splitext(self.output["db"])[0]
        log, input = self.log, self.input
        cmd = "{load} makeblastdb -in {input[db]} -out {db} -dbtype prot -parse_seqids > {log} 2>&1".format(**locals())
        return cmd
