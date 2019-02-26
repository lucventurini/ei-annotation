from .abstract import EIWorfkflow
from .preparation import PrepareWrapper
from .repeats import RepeatMasking
from .preparation.prepare import parse_samplesheet
from .rnaseq.alignments import ShortAlignmentsWrapper, LongAlignmentsWrapper
from .rnaseq.alignments.portcullis import PortcullisWrapper
from .rnaseq.mikado import Mikado
from .rnaseq.assemblies import AssemblyWrapper
from .proteins import ExonerateProteinWrapper, GTHProteinWrapper
from .abinitio.augustus import AugustusWrapper
from .abinitio.converters import ConvertToHints
import os
import yaml


class AnnotationWorklow(EIWorfkflow):

    def __init__(self, configuration, genome, samplesheet):

        # First thing: prepare the reads. This is *outside* the snakemake.
        self.configuration = parse_samplesheet(samplesheet, configuration)
        super().__init__(self.configuration)
        with open(".configuration.yaml", "wt") as out:
            yaml.dump(self.configuration, out)
        self.prepare = PrepareWrapper(self.configuration, genome)
        self.merge([self.prepare])
        faid = [_ for _ in self if _.rulename == "faidx_genome"].pop()
        assert list(faid.input.keys()) == ["genome"], faid.input
        # Second thing: prepare the genome
        self.short_wrapper = ShortAlignmentsWrapper(self.prepare)

        faid = [_ for _ in self if _.rulename == "faidx_genome"].pop()
        assert list(faid.input.keys()) == ["genome"], faid.input
        self.portcullis = PortcullisWrapper(self.short_wrapper)
        self.assemblies = AssemblyWrapper(self.short_wrapper)
        self.merge([self.short_wrapper, self.portcullis, self.assemblies])

        self.long_wrapper = LongAlignmentsWrapper(self.prepare)
        faid = [_ for _ in self if _.rulename == "faidx_genome"].pop()
        assert list(faid.input.keys()) == ["genome"], faid.input
        if self.long_wrapper.gfs:
            self.mikado_long = Mikado(assemblies=self.assemblies,
                                      long_alignments=self.long_wrapper,
                                      portcullis=self.portcullis,
                                      only_long=True)
            assert self.mikado_long.exit
            if not self.mikado_long.stats:
                raise ValueError("No stats for mikado, number of gfs: {}".format(len(self.long_wrapper.gfs)))
            self.add_edge(self.long_wrapper, self.mikado_long)
            self.add_edge(self.portcullis, self.mikado_long)

        self.mikado = Mikado(assemblies=self.assemblies,
                             long_alignments=self.long_wrapper,
                             portcullis=self.portcullis,
                             only_long=False)
        self.add_edge(self.assemblies, self.mikado)
        self.add_edge(self.long_wrapper, self.mikado)
        self.add_edge(self.portcullis, self.mikado)

        faid = [_ for _ in self if _.rulename == "faidx_genome"].pop()
        assert list(faid.input.keys()) == ["genome"], faid.input
        self.repeats = RepeatMasking(self.prepare)

        if self.use_exonerate:
            self.protein_alignments = ExonerateProteinWrapper(self.repeats, self.portcullis)
        else:
            self.protein_alignments = GTHProteinWrapper(self.repeats, self.portcullis)

        self.merge([self.repeats, self.protein_alignments])

        # self.augustus = AugustusWrapper(mikado=self.mikado, mikado_long=self.mikado_long,
        #                                 faidx=faid, proteins=self.protein_alignments,
        #                                 rmasker=self.repeats)

        self.add_final_flag()
        # faid = [_ for _ in self if _.rulename == "faidx_genome"].pop()
        # assert list(faid.input.keys()) == ["genome"], faid.input

    @property
    def flag_name(self):
        return os.path.join(self.configuration["outdir"], "all.done")

    @property
    def use_exonerate(self):
        return self.configuration["homology"].get("use_exonerate", True)

    __final_rulename__ = "all"
