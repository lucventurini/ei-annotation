import os
from eicore.external_process.snakemake_helper import loadPreCmd


def strandedness_calculator(sample):

  stranded = SAMPLE_MAP[sample]

  if stranded == "r":  # "fr-firststrand"
    return "-d '1-+,1+-,2++,2--'"
  elif stranded == "f":  # "fr-secondstrand"
    return "-d '1++,1--,2++,2--'"
  else:
    return ""


rule chrom_sizes:
  input: REF + ".fai"
  output: REF + "sizes.txt"
  shell: """cut -f 1,2 {input} > {output}"""


rule bam2wig:
  input:
    bam=,  # I have to decide on a folder structure first
    chrom_file=rules.chrom_sizes.output
  output:
    plus=os.path.join(outdir, "Alignments", "{sample}_Forward.wig"),
    minus=os.path.join(outdir, "Alignments", "{sample}_Reverse.wig"),
  params:
    strandedness = strandedness_calculator,
    prefix=os.path.join(outdir, "Alignments", "{sample}")
    load=loadPreCmd(config["load"]["rseqc"])
  run:
    shell("""{params.load} bam2wig.py {params.strandedness} -s {input.chrom_file} -i {input.bam} -o {params.prefix}""")
    if not params.strandedness:
      shell("mv {params.prefix}.wig {params.prefix}_Forward.wig")
      touch("{params.prefix}_Reverse.wig")
