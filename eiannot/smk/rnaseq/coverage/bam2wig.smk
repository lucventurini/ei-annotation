import os
from eicore.external_process.snakemake_helper import loadPreCmd


def get_sample_bam(wildcards):
  return samples[wildcards.sample]["bam"]


rule bam2wig:
  input:
    bam=get_sample_bam
  output:
    # Four files, one each for the four possible combinations of strand/mapping
    os.path.join(outdir, "{sample}.done")
  params:
    minquality=config.get("foo", 30)
    load=loadPreCmd(config["load"]["augustus"])
  ### -F 3844: exclude reads unmapped, not primary alignment, failing QC, optical dups, suppl. alignment
  ### flag 16: read on the reverse strand
  ### flag 64: read 1
  ### flag 128: read 2
  shell: """{params.load} \\
  bam2wig <(samtools view -b -q {params.minquality} -F 3844 -F 16 -f 64 {input.bam}) > {wildcards.sample}.read1_+.wig &&
  bam2wig <(samtools view -b -q {params.minquality} -F 3844 -f 16 -f 64 {input.bam}) > {wildcards.sample}.read1_-.wig &&
  bam2wig <(samtools view -b -q {params.minquality} -F 3844 -F 16 -f 128 {input.bam}) > {wildcards.sample}.read2_+.wig &&
  bam2wig <(samtools view -b -q {params.minquality} -F 3844 -f 16 -f 128 {input.bam}) > {wildcards.sample}.read2_-.wig &&
  touch {output}"""
