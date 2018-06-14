import os
import sys
from eicore.external_process.snakemake_helper import loadPreCmd


#Use join_aug_pred.pl and the predictionStart / predictionEnd arguments to make the problem tractable

# Copying from panYeast

# Run Augustus and produce a gff file
rule augustus_gff:
  input:
    asm = lambda wildcards: "%s/asm.%s.fa" % (__GIVEN_ASM_OUTDIR__, wildcards.asm)
  output:
    gff = "%s/augustus_gff.{asm}.gff" % __AUGUSTUS_OUTDIR__
  threads: 4
  params:
    augustus_species = lambda wildcards: config["data"][wildcards.asm]["augustus_species"] if "augustus_species" in config["data"][wildcards.asm] else tconfig["augustus_species"],
    augustus_params  = tconfig["augustus_params"],
    rule_outdir = __AUGUSTUS_OUTDIR__
  conda: "%s/conda_envs/augustus.yaml" % __PIPELINE_COMPONENTS__
  benchmark: "%s/augustus_gff.{asm}" % __LOGS_OUTDIR__
  shell: """
    augustusspecies="{params.augustus_species}"
    if [ -e '{params.augustus_species}' ]; then
      cp -r '{params.augustus_species}' "$AUGUSTUS_CONFIG_PATH/species/"
      augustusspecies=`basename $augustusspecies`
    fi
    which augustus
    mkdir -p {params.rule_outdir}
    augustus {params.augustus_params} \
             --gff3=on \
             --genemodel=complete \
             --strand=both  \
             --species=$augustusspecies \
             {input.asm} \
      > {output.gff}
  """

###############################################################################
# Rename the proteins in the GFF file

rule augustus_gff_sample:
  input:
    gff = lambda wildcards: "%s/augustus_gff.%s.gff" % (__AUGUSTUS_OUTDIR__, wildcards.asm)
  output:
    gff = "%s/augustus.{asm}.gff" % __AUGUSTUS_OUTDIR__
  params:
    geneid_prefix = lambda wildcards: wildcards.asm
  shell: """
    sed -e "s/\([= ]\)\(g[0-9]\+\)/\\1{params.geneid_prefix}|\\2/g" {input.gff} > {output.gff}
  """


rule augustus_prot_fasta:
  input:
    nt_fasta = lambda wildcards: "%s/transcripts.%s.fa" % (__TRANS_OUTDIR__, wildcards.asm)
  output:
    prot_fasta = "%s/augustus.{asm}.prots.fa" % __AUGUSTUS_OUTDIR__
  run:
    from utils import seq as seq
    with open(output.prot_fasta, 'w') as fo:
      for (fa_header, fa_seq) in zip(*seq.read_fasta(input.nt_fasta)):
        fo.write(">%s\n" % fa_header)
        fo.write("%s\n" % seq.nt_translate(fa_seq))
      #efor
    #ewith
#  shell: """
#    gffread -y {output.prot_fasta} -g {input.asm} {input.gff}
#  """