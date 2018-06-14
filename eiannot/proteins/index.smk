import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule index_protein:
  input: rules.sanitize_proteins.output
  output: os.path.join(outdir, "FASTAs", "{protein}.fasta.fai")
  params:
    load=loadPreCmd(config["load"]["samtools"])
  shell: """{params.load} samtools faidx {input}"""