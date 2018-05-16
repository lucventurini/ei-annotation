from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule sanitize_reference:
  input: REF
  output: os.path.join(outdir, "0-Reference", "genome.fa")
  params:
    load=loadPreCmd(config.get("load", dict()).get("mikado", None))
  shell: """{params.load} sanitize_blast_db.py -o {output} {input}"""


rule samtools_index_reference:
  input: rules.sanitize_reference.output
  output: rules.sanitize_reference.output + ".fai"
  params:
    load=loadPreCmd(config.get("load", dict()).get("samtools", None))
  shell: """{params.load} samtools faidx {input}"""