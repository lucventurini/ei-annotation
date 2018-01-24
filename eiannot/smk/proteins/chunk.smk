import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule chunk_proteins:
  input: os.path.join(outdir, "FASTAs", "{protein}.fasta")  # rules.sanitize_proteins.output
  output: os.path.join(outdir, "FASTAs", "Chunks", "{protein}", "split.done")
  params:
    load=loadPreCmd(config["load"]["mikado"]),
    outdir=os.path.join(outdir, "FASTAs", "Chunks", "{protein}/{protein}"),
    chunks=protein_chunks,
  shell: """{params.load} mkdir -p {params.outdir} && split_fasta.py -m {params.chunks} {input} {params.outdir} && touch {output}"""