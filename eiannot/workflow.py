from .abstract import EIWorfkflow, FinalFlag
from .preparation import PrepareWrapper
from .preparation.prepare import parse_samplesheet
from .rnaseq.alignments.workflow import ShortAlignmentsWrapper, LongAlignmentsWrapper
from .rnaseq.alignments.portcullis import PortcullisWrapper
from .rnaseq.mikado import Mikado
from .rnaseq.assemblies.workflow import AssemblyWrapper
import os


class AnnotationWorklow(EIWorfkflow):

    def __init__(self, configuration, genome, samplesheet):

        # First thing: prepare the reads. This is *outside* the snakemake.
        configuration = parse_samplesheet(samplesheet, configuration)
        super().__init__(configuration)
        self.prepare = PrepareWrapper(self.configuration, genome)
        self.merge([self.prepare])
        # Second thing: prepare the genome
        self.short_wrapper = ShortAlignmentsWrapper(self.prepare)
        self.long_wrapper = LongAlignmentsWrapper(self.prepare)
        self.portcullis = PortcullisWrapper(self.short_wrapper)
        self.merge([self.short_wrapper, self.long_wrapper, self.portcullis])
        self.assemblies = AssemblyWrapper(self.short_wrapper)
        self.merge([self.assemblies])
        self.mikado = Mikado(assemblies=self.assemblies, long_alignments=self.long_wrapper, portcullis=self.portcullis)
        self.merge([self.mikado])
        flag = os.path.join(self.configuration["outdir"], "all.done")
        self.add_final_flag(flag)
