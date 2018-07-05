from .abstract import EIWorfkflow, FinalFlag
from .preparation import PrepareWrapper, FaidxProtein, SanitizeProteinBlastDB
from .repeats.workflow import RepeatMasking
from .preparation.prepare import parse_samplesheet
from .rnaseq.alignments.workflow import ShortAlignmentsWrapper, LongAlignmentsWrapper
from .rnaseq.alignments.portcullis import PortcullisWrapper
from .rnaseq.mikado.workflow import Mikado
from .rnaseq.assemblies.workflow import AssemblyWrapper
from .proteins.workflow import ExonerateProteinWrapper
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
            assert self.fln_long.entries
            self.add_edge(self.mikado_long, self.fln_long)
        else:
            self.fln_long = None

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
        self.protein_alignments = ExonerateProteinWrapper(self.repeats, self.portcullis)
        self.merge([self.repeats, self.protein_alignments])




        flag = os.path.join(self.configuration["outdir"], "all.done")
        self.add_final_flag(flag)
