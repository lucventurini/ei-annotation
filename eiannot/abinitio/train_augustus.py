from .fln import FlnWrapper
from ..abstract import AtomicOperation, EIWrapper
from ..repeats import RepeatMasking
import os


class TrainAugustusWrapper(EIWrapper):

    __final_rulename__ = "augustus_training_done"

    def __init__(self, fln_wrapper: FlnWrapper, rm_wrapper: RepeatMasking):
        super().__init__(configuration=fln_wrapper.configuration)

        self.__is_long = fln_wrapper.is_long
        self.__training_candidates = fln_wrapper.training_candidates
        prepare_aug_folder = PrepareAugConfig(configuration=self.configuration, is_long=self.is_long)
        self.add_edge(fln_wrapper, prepare_aug_folder)
        self.add_edge(rm_wrapper, prepare_aug_folder)
        self.add_node(prepare_aug_folder)
        if self.training_candidates is None or self.augustus_species is not None:
            pass
        else:
            # First: convert candidates into GeneBank
            converter = ConvertToGb(fln_wrapper, prepare_aug_folder)
            self.add_edge(prepare_aug_folder, converter)
            trainer = TrainAugustus(converter, prepare_aug_folder)
            self.add_edge(converter, trainer)
            self.add_edge(prepare_aug_folder, trainer)
        self.add_final_flag()

    @property
    def training_candidates(self):
        return self.__training_candidates

    @property
    def outdir(self):
        return os.path.join(self._root_dir,
                            "abinitio", "2-training{long}".format(long="-long-reads" if self.is_long else ""))

    @property
    def is_long(self):
        return self.__is_long

    @property
    def augustus_species(self):
        return self.configuration.get("training", {}).get(
            "model_species", None
        )

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "augustus_training.done")


class PrepareAugConfig(AtomicOperation):

    def __init__(self, configuration, is_long):

        super().__init__()
        self.is_long = is_long
        self.configuration = configuration
        self.input["genome"] = self.genome
        self.output["generic_config"] = os.path.join(self.outdir,
                                                     "config",
                                                     "species",
                                                     "generic",
                                                     "generic_metapars.cfg")
        if self.augustus_species is not None:
            self.output["species_config"] = os.path.join(
                self.outdir,
                "config",
                "species",
                self.augustus_species,
                "{species}_parameters.cfg".format(species=self.augustus_species)
                )

    @property
    def loader(self):
        return ["augustus"]

    @property
    def rulename(self):
        return "prepare_augustus_config_folder"

    @property
    def outdir(self):
        return os.path.join(self._root_dir,
                            "abinitio", "2-training{long}".format(long="-long-reads" if self.is_long else ""))

    @property
    def augustus_species(self):
        return self.configuration.get("training", {}).get(
            "model_species", None
        )

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        cmd = "{load}"
        cmd += "mkdir -p {outdir} && cd {outdir} && "
        cmd += "mkdir -p config/species && cp -r ${{AUGUSTUS_CONFIG_PATH}}/species/generic config/species/ && "
        cmd += "cp -r ${{AUGUSTUS_CONFIG_PATH}}/model config/ "
        if self.augustus_species is not None:
            species = self.augustus_species
            cmd += "&& cp -r ${{AUGUSTUS_CONFIG_PATH}}/species/{species} config/species/"

        cmd = cmd.format(**locals())
        return cmd


class ConvertToGb(AtomicOperation):

    def __init__(self, fln_wrapper: FlnWrapper, preparer: PrepareAugConfig):

        super().__init__()
        self.configuration = fln_wrapper.configuration
        self.is_long = fln_wrapper.is_long
        self.flank = fln_wrapper.fln_filter.flank
        self.outdir = preparer.outdir
        self.input["Training"] = fln_wrapper.training_candidates
        self.input["Testing"] = fln_wrapper.testing_candidates
        self.input["genome"] = self.masked_genome  # Double check with Gemy and David!
        self.output["Training"] = os.path.join(self.outdir, "train.gb")
        self.output["Testing"] = os.path.join(self.outdir, "test.gb")
        self.log = os.path.join(self.outdir, "convert.log")

    @property
    def loader(self):
        return ["augustus"]

    @property
    def threads(self):
        return 1

    @property
    def is_small(self):
        return True

    @property
    def rulename(self):
        return "convert_{long}candidates_to_gb".format(long="long_" if self.is_long else "")

    @property
    def cmd(self):
        load = self.load
        flank = self.flank
        outdir = self.outdir
        out_training = self.output["Training"]
        out_testing = self.output["Testing"]
        genome = self.input["genome"]
        training = self.input["Training"]
        testing = self.input["Testing"]
        log = self.log
        cmd = "{load} mkdir -p {outdir} && "
        cmd += "gff2gbSmallDNA.pl {training} {genome} {flank} {out_training} --softmasked > {log} 2> {log} &&"
        cmd += "gff2gbSmallDNA.pl {testing} {genome} {flank} {out_testing} --softmasked > {log} 2> {log}"
        cmd = cmd.format(**locals())
        return cmd


class TrainAugustus(AtomicOperation):

    def __init__(self, converter: ConvertToGb, config_preparer: PrepareAugConfig):

        super().__init__()
        self.configuration = converter.configuration
        self.outdir = converter.outdir
        self.is_long = converter.is_long
        self.input.update(converter.output)
        self.input.update(config_preparer.output)
        self.log = os.path.join(self.outdir, "training.log")
        self.output["flag"] = os.path.join(self.outdir, "training.done")

    @property
    def rulename(self):
        return "train_augustus{long}".format(long="_long" if self.is_long else "")

    @property
    def loader(self):
        return ["augustus"]

    @property
    def quick(self):
        """Boolean in training/quick. Default: true.
        If False, EIannot will perform a full optimize_augustus run.
        Otherwise, it will just perform a single etraining."""
        return self.configuration.get("training", {}).get(
            "quick", True)

    @property
    def threads(self):
        if not self.quick:
            return super(self).threads
        else:
            return 1

    @property
    def is_small(self):
        return not self.quick

    @property
    def cmd(self):
        load = self.load
        species = self.species
        outdir = self.outdir
        training_src = os.path.relpath(os.path.abspath(self.input["Training"]), start=self.outdir)
        log = os.path.abspath(self.log)
        flag = os.path.basename(self.output["flag"])

        cmd = "{load} "
        cmd += "mkdir -p {outdir} && cd {outdir} && "
        cmd += "ln -sf {training_src} training.gb && "
        cmd += "new_species.pl --AUGUSTUS_CONFIG_PATH=$(pwd)/config --species={species} 2> new_species.log >&2 && "
        cmd += " sed -i '/^UTR/s/.*/UTR\ton/' config/species/{species}/{species}_parameters.cfg && "
        cmd += "AUGUSTUS_CONFIG_PATH=$(pwd)/config && "
        if self.quick is True:
            cmd += "(etraining --species={species} training.gb 2>&1 > {log} "
        else:
            cmd += "(optimize_augustus.pl --species={species} --kfold={threads} --cpus={threads} --UTR=on "
            cmd += "training.gb 2>&1 > {log}"
        cmd += " && touch {flag})"

        cmd = cmd.format(**locals())

        return cmd
