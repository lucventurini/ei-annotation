from .abstract import EIWorfkflow
from .preparation import PrepareWrapper
from .rnaseq.alignments.workflow import ShortAlignmentsWrapper, LongAlignmentsWrapper
from .rnaseq.assemblies.workflow import AssemblyWrapper


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
        self.assemblies = AssemblyWrapper(self.short_wrapper)
