import os
from eicore.external_process.snakemake_helper import loadPreCmd

print("PROTEINS keys:", PROTEINS.keys())
print("ARRAY:", protein_chunk_array)

rule collapse_protein_alignment:
  input: expand(os.path.join(outdir, "Alignments", "Chunks", "{protein}", "{protein_chunk}.exonerate.txt"), protein_chunk=protein_chunk_array, protein=PROTEINS.keys())
  output: os.path.join(outdir, "Alignments", "{protein}.exonerate")
  shell: """cat {input} > {output}"""