import os
from eicore.external_process.snakemake_helper import loadPreCmd


def strandedness_calculator():

  pass


rule chrom_sizes:
  input: REF + ".fai"
  output: REF + "sizes.txt"
  shell: """cut -f 1,2 {input} > {output}"""


rule bam2wig_unstranded:
  input:
    bam=,
    chrom_file=rules.chrom_sizes.output
  output: touch(os.path.join(outdir, "Alignments", "{sample}.wig.done"))
  params:
    strandedness = strandedness_calculator,
    prefix=os.path.join(outdir, "Alignments", "{sample}")
    load=loadPreCmd(config["load"]["rseqc"])
  shell: """{params.load} bam2wig.py {params.strandedness} -s {input.chrom_file} -i {input.bam} -o {params.prefix}"""