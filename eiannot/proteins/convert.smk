import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule convert_protein_alignment:
  input:
    protein=os.path.join(outdir, "Alignments", "{protein}.exonerate"),
    protein_fai=os.path.join(outdir, "FASTAs", "{protein}.fasta.fai")   # rules.index_protein.output
  output: os.path.join(outdir, "Alignments", "{protein}.gff3")
  log: os.path.join(outdir, "Alignments", "{protein}.convert.log")
  params:
    coverage=protein_coverage,
    identity=protein_identity,
  threads: 1
  shell: """exonerate2gff.pl --in {input} --minIdentity {params.identity} --minCoverage {params.coverage} --fasta_count $(cat {input.protein_fai} | wc -l) > {output} 2> {log}"""