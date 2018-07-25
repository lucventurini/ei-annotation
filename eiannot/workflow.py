from .abstract import EIWorfkflow
from .preparation import PrepareWrapper
from .repeats import RepeatMasking
from .preparation.prepare import parse_samplesheet
from .rnaseq.alignments import ShortAlignmentsWrapper, LongAlignmentsWrapper
from .rnaseq.alignments.portcullis import PortcullisWrapper
from .rnaseq.mikado import Mikado
from .rnaseq.assemblies import AssemblyWrapper
from .proteins import ExonerateProteinWrapper, GTHProteinWrapper
from .abinitio.fln import FlnWrapper
import os


class AnnotationWorklow(EIWorfkflow):

    def __init__(self, configuration, genome, samplesheet):

        # First thing: prepare the reads. This is *outside* the snakemake.
        self.configuration = parse_samplesheet(samplesheet, configuration)
        super().__init__(self.configuration)
        self.prepare = PrepareWrapper(self.configuration, genome)
        self.merge([self.prepare])
        # Second thing: prepare the genome
        self.short_wrapper = ShortAlignmentsWrapper(self.prepare)

        self.portcullis = PortcullisWrapper(self.short_wrapper)
        self.assemblies = AssemblyWrapper(self.short_wrapper)
        self.merge([self.short_wrapper, self.portcullis, self.assemblies])

        self.long_wrapper = LongAlignmentsWrapper(self.prepare)
        if self.long_wrapper.gfs:
            self.mikado_long = Mikado(assemblies=self.assemblies,
                                      long_alignments=self.long_wrapper,
                                      portcullis=self.portcullis,
                                      only_long=True)
            assert self.mikado_long.exit
            if not self.mikado_long.stats:
                raise ValueError("No stats for mikado, number of gfs: {}".format(len(self.long_wrapper.gfs)))
            self.merge([self.mikado_long])
            self.fln_long = FlnWrapper(self.mikado_long)
            # assert self.fln_long.entries
            self.add_edge(self.mikado_long, self.fln_long)
        else:
            self.fln_long = None

        if self.long_wrapper.gfs or self.assemblies.gfs:
            self.mikado = Mikado(assemblies=self.assemblies,
                                 long_alignments=self.long_wrapper,
                                 portcullis=self.portcullis,
                                 only_long=False)
            self.merge([self.mikado])
            if self.mikado.stats:
                self.fln = FlnWrapper(self.mikado)
                self.merge([self.fln])
                self.add_edge(self.mikado, self.fln)

        self.repeats = RepeatMasking(self.prepare)

        if self.use_exonerate:
            self.protein_alignments = ExonerateProteinWrapper(self.repeats, self.portcullis)
        else:
            self.protein_alignments = GTHProteinWrapper(self.repeats, self.portcullis)

        self.merge([self.repeats, self.protein_alignments])

        self.add_final_flag()

    @property
    def flag_name(self):
        return os.path.join(self.configuration["outdir"], "all.done")

    @property
    def use_exonerate(self):
        return self.configuration["homology"].get("use_exonerate", True)

    __final_rulename__ = "all"
