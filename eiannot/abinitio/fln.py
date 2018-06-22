from ..abstract import EIWrapper, AtomicOperation
from ..rnaseq.mikado import IndexMikado, Mikado


class FlnWrapper(EIWrapper):

    def __init__(self, mikado: Mikado):

        super().__init__()

        if mikado.picker:
            # Start extracting the data ...
            sequence_extractor = None



class MikadoSequenceExtractor(AtomicOperation):

    def __init__(self, picker):
        # rule extract_mikado_fasta:
        #   input:
        #     mikado=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
        #     genome=REF
        #   output: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.cdna.fasta")
        #   params:
        #     load=loadPreCmd(config.get("load", dict()).get("gffread", None))
        #   shell: """{params.load} gffread -g {input.genome} -w {output} -C {input.mikado}"""

        pass

    @property
    def loader(self):
        return ["gffread"]

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        genome = self.genome
        outdir = self.outdir

        cmd = "{load} gffread -g {genome} -w {output[transcripts]} -C {input[loci]}"

    @property
    def outdir(self):
        pass



class FLN(AtomicOperation):

    def __init__(self):
        pass


class ConvertModelsToGb(AtomicOperation):

    def __init__(self):

        """Task to convert models from FLN into GeneBank format"""
    #       rule convert_training_to_gb:
    #     input:
    #       gff=rules.extract_training_gff.output,
    #       fa=rules.sanitize_reference.output  # We might use the SM reference as well
    #     output: os.path.join(OU, "training_models.gb")
    #     params:
    #       load=loadPreCmd(config.get("load", dict()).get("augustus", None)),
    #       flank=2000,  # Move into configuration
    #       keep_mask=  # Potentially we could put "--softmasked" here
    #     shell: """{params.load} gff2gbSmallDNA.pl {input.gff} {input.softmasked} {params.flank} > {output}"""

    @property
    def rulename(self):
        return "convert_fln_to_gb"

    @property
    def loader(self):
        return ["augustus"]
