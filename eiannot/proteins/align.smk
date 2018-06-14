import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule align_protein:
  input:
    protein=rules.chunk_proteins.output,
    genome=REF
  output: os.path.join(outdir, "Alignments", "Chunks", "{protein}", "{protein_chunk}.exonerate.txt")
  params:
    load=loadPreCmd(config["load"]["exonerate"]),
    min_intron=MIN_INTRON,
    max_intron=MAX_INTRON,
    outdir=os.path.join(outdir, "Alignment", "{protein}"),
    fasta=os.path.join(outdir, "FASTAs", "Chunks", "{protein}", "{protein}_{protein_chunk}.fasta")
  log: os.path.join(outdir, "Alignments", "Chunks", "{protein}", "{protein_chunk}.exonerate.log")
  shell: """{params.load} exonerate --model protein2genome --showtargetgff yes --showvulgar yes \
  --softmaskquery yes --softmasktarget yes --bestn 10  --minintron {params.min_intron} --maxintron {params.max_intron} \
  --percent 30 --score 50 --geneseed 50 --showalignment no --query {params.fasta} --target {input.genome} --ryo '>%qi\tlength=%ql\talnlen=%qal\tscore=%s\tpercentage=%pi\nTarget>%ti\tlength=%tl\talnlen=%tal\n' \
  > {output} 2> {log}"""