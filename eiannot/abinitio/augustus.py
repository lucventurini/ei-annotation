from ..abstract import AtomicOperation, EIWrapper
from .train_augustus import TrainAugustusWrapper
from ..rnaseq.mikado import Mikado
from ..repeats import RepeatMasking
from ..preparation import FaidxGenome
from .converters import ConvertToHints
from .fln import FlnWrapper
from ..proteins import ProteinWrapper
import pkg_resources
from .fasta_split import FastaMSplitter
from .abstract import augustus_root_dir, AugustusMethod
import os


class AugustusWrapper(EIWrapper):

    def __init__(self,
                 faidx: FaidxGenome,
                 mikado: Mikado,
                 mikado_long: Mikado,
                 rmasker: RepeatMasking,
                 proteins: ProteinWrapper
                 ):

        # First we have to do FLN

        super().__init__(faidx.configuration)
        self.fai = faidx
        self.rmasker = rmasker
        self.proteins = proteins
        if mikado_long:
            self.mikado_long = mikado_long
            self.fln_long = FlnWrapper(mikado_long, self.fai)
            self.add_edge(mikado_long, self.fln_long)

        else:
            self.mikado_long = None
            self.fln_long = None
        # assert self.fln_long.entries

        if mikado:
            self.mikado = mikado
            self.fln = FlnWrapper(mikado, self.fai)
            self.add_edge(mikado, self.fln)
        else:
            self.mikado = None
            self.fln = None

        self.trainer = TrainAugustusWrapper(fln_wrapper=self.fln, rm_wrapper=self.rmasker)
        self.fasta_splitter = FastaMSplitter(faidx)

        dir = None
        fnames = []

        for run in self.runs:

            runner = AugustusRunWrapper(fasta_splitter=self.fasta_splitter,
                                        mikado=self.fln,
                                        mikado_long=self.fln_long,
                                        rmasker=self.rmasker,
                                        proteins=self.proteins,
                                        trainer=self.trainer,
                                        run=run)
            if dir is None:
                dir = self._augustus_dir

    @property
    def _augustus_dir(self):
        return os.path.join(self._root_dir, augustus_root_dir)

    @property
    def flag_name(self):
        # TODO implement
        pass

    @property
    def runs(self):
        # TODO implement
        """This property returns the runs to be executed according to the configuration"""
        return []


class AugustusRunWrapper(EIWrapper):

    __final_rulename__ = "augustus_run"

    def __init__(self,
                 fasta_splitter: FastaMSplitter,
                 mikado: FlnWrapper,
                 mikado_long: FlnWrapper,
                 rmasker: RepeatMasking,
                 proteins: ProteinWrapper,
                 trainer:TrainAugustusWrapper,
                 run: int):

        super().__init__(configuration=mikado.configuration)
        self.__run = run
        self.trainer = trainer
        self.converter = ConvertToHints(mikado=mikado,
                                        repeats=rmasker,
                                        mikado_long=mikado_long,
                                        run=run,
                                        proteins=proteins)
        self.add_edges_from([(_, self.converter) for _ in (mikado, mikado_long, rmasker, proteins)])
        self.fasta_splitter = fasta_splitter
        self.runs = [AugustusChunk(self.fasta_splitter, self.trainer, self.converter, chunk, run)
                     for chunk in self.fasta_splitter.chunks]

        self.add_edges_from([(self.trainer, run) for run in self.runs])
        self.add_edges_from([(self.converter, run) for run in self.runs])

    @property
    def augustus_run(self):
        return self.__run

    @property
    def final_rule(self):
        return "{}_{}".format(self.__final_rulename__, self.augustus_run)

    @property
    def flag_name(self):
        # TODO implement
        pass


class AugustusChunk(AtomicOperation):

    def __init__(self,
                 fastam: FastaMSplitter,
                 trained: TrainAugustusWrapper,
                 hints: ConvertToHints,
                 chunk: int,
                 aug_run: int):

        super().__init__()
        self.input.update(fastam.output)
        self.input["hints"]
        self.__chunk = chunk
        self.__run = aug_run

    @property
    def run(self):
        return self.__run

    @property
    def chunk(self):
        return self.__chunk

    @property
    def rulename(self):
        return "run_augustus_run{}_{}".format(self.run, self.chunk)

    @property
    def species(self):

        if self.configuration.get("training", {}).get("train", True):
            return self.input["trained"]
        else:
            model = self.configuration.get("training", {}).get("model_species", None)
            assert model is not None
            return model

    @property
    def loader(self):
        return ["ei-annotation", "augustus", "mikado"]

    @property
    def extrinsic(self):
        if self.configuration["abinitio"]["extrinsic"] is None:
            return pkg_resources.resource_filename("eiannot.configuration", "extrinsic.ei_augustus_generic.cfg")
        else:
            return self.configuration["abinitio"]["extrinsic"]

    @property
    def cmd(self):

        cmd = "{load}"
        load = self.load
        species = self.species
        extrinsic = self.extrinsic
        input, output = self.input
        log = self.log
        chunk = self.chunk
        if input["hints"] is not None:
            hints = "--hintsfile={input[hints]}"
        else:
            hints = ""
        cmd += "execute_augustus.py {input[genome]} {input[chunk_db]} {chunk} {output[gtf]} "
        # Execute augustus will take as last argument the augustus command string
        cmd += "\"augustus --species {species}  --UTR=on --extrinsicCfgFile={extrinsic} "
        cmd += "--stopCodonExcludedFromCDS=true --genemodel=partial "
        cmd += " --alternatives-from-evidence=true {hints} --noInFrameStop=true "
        cmd += " --allow_hinted_splicesites=atac --errfile={log} \""
        cmd = cmd.format(**locals())
        return cmd


class JoinChunks(AugustusMethod):

    def __init__(self, chunks: [AugustusChunk]):

        if len(chunks) == 0:
            raise ValueError("No chunk provided")

        super().__init__(chunks[0].configuration)
        self.__run = chunks[0].run

        self.input["chunks"] = [chunk.output["gtf"] for chunk in chunks]
        self.input["list"] = os.path.join(self._augustus_dir, "output", "file_list.{}.txt".format(self.run))
        self.__create_file_list()

        self.output["gff"] = os.path.join(self._augustus_dir, "output",
                                          "joined_{}.gtf".format(self.run))

    @property
    def loader(self):
        return ['augustus']

    @property
    def cmd(self):
        # TODO: at the moment I am just using joingenes directly, giving in a list of files;
        # TODO: however, it might be more efficient to write a wrapper, like I did for Augustus itself, especially
        # TODO: because different chromosomes can certainly be analysed in parallel.
        load = self.load
        input, output = self.input, self.output
        minoverlap = self.minoverlap  # This is to ensure that we do not miss anything that might be on the margins
        cmd = "{load} joingenes --inputfile={input.list} -o {output.gff} -e {minoverlap} -m eukaryote -a"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def run(self):
        return self.__run

    @property
    def rulename(self):
        return "join_genes_augustus_run{}".format(self.run)

    def __create_file_list(self):
        if not os.path.exists(os.path.dirname(self.input["list"])):
            os.makedirs(os.path.dirname(self.input["list"]))

        with open(self.input["list"], "wt") as out:
            for gt in self.input["chunks"]:
                print(gt, 1, sep="\t", file=out)

    @property
    def __subfolder(self):
        return "output"

    @property
    def minoverlap(self):
        # TODO put this into the configuration
        return self.configuration.get("abinitio", dict()).get("minsize", 5 * 10 ** 5)
