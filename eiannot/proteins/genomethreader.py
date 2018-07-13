from .chunking import ChunkProteins
from ..abstract import AtomicOperation


class GTH(AtomicOperation):

    ___toolname__ = "gth"

    def __init__(self):
        super().__init__()
        self.input["genome"] = self.genome

    @property
    def loader(self):
        return ["gth"]

    @property
    def toolname(self):
        return self.___toolname__

    @property
    def species(self):
        if "species" in self.configuration["programs"][self.toolname]:
            return " -species {} ".format(self.configuration["programs"][self.toolname]["species"])
        else:
            return " "

    @property
    def coverage(self):
        if "coverage" in self.configuration["programs"][self.toolname]:
            return " --gcmincoverage {} ".format(self.configuration["programs"][self.toolname]["coverage"])
        else:
            return " "

    @property
    def extra(self):
        return self.configuration["programs"].get(self.toolname, {}).get("extra", " ")

    @property
    def gcintron(self):
        return "-gcmaxgapwidth={}".format(self.max_intron)

    @property
    def cmd(self):

        # gth -intermediate  -introncutout -first 10 -species rice -gff3out -gcmincoverage 50 -paralogs
        # -o eitest/proteins/alignments/chunk_003.gff3 -force
        # -genomic eitest/inputs/reference/genome.fa -protein eitest/proteins/chunks/chunk_003.stops.fasta

        load = self.load
        cmd = "{load} gth -intermediate  -introncutout -first 10 "
        species = self.species
        extra = self.extra
        cmd += " {species} -gff3out {coverage} {extra} -paralogs "
        input, output = self.input, self.output
        cmd += " -o {output[gff3]} -force -genomic {input[genome]} -protein {input[chunk]}"

        cmd = cmd.format(**locals())
        return cmd
