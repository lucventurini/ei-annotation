import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule sanitize_proteins:
  input: protein_to_file
  output: os.path.join(outdir, "FASTAs", "{protein}.fasta")
  params:
    load=loadPreCmd(config["load"]["mikado"])
  shell: """{params.load} sanitize_blast_db.py --out {output} {input}"""