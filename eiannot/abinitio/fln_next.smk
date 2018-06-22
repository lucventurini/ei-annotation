import os
import subprocess
from eicore.external_process.snakemake_helper import loadPreCmd


rule fln_all:
	input: expand(os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln_results", "pt_seqs"), chunk_id=CHUNK_ARRAY)
	output:
	  flag=os.path.join(outdir, "Hints", "FLN", "fln.done"),
	  table=os.path.join(outdir, "Hints", "FLN", "fln.table.txt")
	shell: """cat {input} | awk 'NR==1 || $1!="Query_id"'   > {output.table} && touch {output.flag}"""


rule filter_fln:
  input:
    table=rules.fln_all.output.table,
    bed12=rules.convert_mikado_to_bed12.output
  output:
    table=os.path.join(outdir, "Hints", "FLN", "fln_filtered.table.txt"),
    list=os.path.join(outdir, "Hints", "FLN", "fln_filtered.list.txt")
  params:
    load=loadPreCmd(config.get("load", {}).get("mikado", None), config.get("load", {}).get("ei-annot", None)),
    out=os.path.join(outdir, "Hints", "FLN", "fln_filtered"),
    mikado=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci")
  shell: """{params.load} filter_fln.py {input.table} {params.mikado} {params.out}"""

rule extract_gold_gff:
  input:
    list=rules.filter_fln.output.list,
    mikado=rules.index_mikado.input.mikado
  output:
    gff3=os.path.join(outdir, "Hints", "FLN", "gold_models.gff3"),
    list=os.path.join(outdir, "Hints", "FLN", "gold_models.txt")
  params:
    load=loadPreCmd(config.get("load", {}).get("mikado", None))
  shell: """{params.load} awk '{{if ($4=="Gold") {{OFS="\\t"; print $1,$2}}}}' {input.list} > {output.list} &&
  mikado util grep {output.list} {input.mikado} {output.gff3}"""

rule extract_silver_gff:
  input:
    list=rules.filter_fln.output.list,
    mikado=rules.index_mikado.input.mikado
  output:
    gff3=os.path.join(outdir, "Hints", "FLN", "silver_models.gff3"),
    list=os.path.join(outdir, "Hints", "FLN", "silver_models.txt")
  params:
    load=loadPreCmd(config.get("load", {}).get("mikado", None))
  shell: """{params.load} awk '{{if ($4=="Silver") {{OFS="\\t"; print $1,$2}}}}' {input.list} > {output.list} &&
  mikado util grep {output.list} {input.mikado} {output.gff3}"""


rule extract_bronze_gff:
  input:
    list=rules.filter_fln.output.list,
    mikado=rules.index_mikado.input.mikado
  output:
    gff3=os.path.join(outdir, "Hints", "FLN", "bronze_models.gff3"),
    list=os.path.join(outdir, "Hints", "FLN", "bronze_models.txt")
  params:
    load=loadPreCmd(config.get("load", {}).get("mikado", None))
  shell: """{params.load} awk '{{if ($4=="Bronze") {{OFS="\\t"; print $1,$2}}}}' {input.list} > {output.list} &&
  mikado util grep {output.list} {input.mikado} {output.gff3}"""


rule extract_training_gff:
  input:
    list=rules.filter_fln.output.list,
    mikado=rules.index_mikado.input.mikado
  output:
    gff3=os.path.join(outdir, "Hints", "FLN", "training_models.gff3"),
    list=os.path.join(outdir, "Hints", "FLN", "training_models.txt")
  params:
    load=loadPreCmd(config.get("load", {}).get("mikado", None))
  shell: """{params.load} awk '{{if ($3=="True") {{OFS="\\t"; print $1,$2}}}}' {input.list} > {output.list} &&
  mikado util grep {output.list} {input.mikado} {output.gff3}"""


rule fln:
  input:
    split=rules.split_mikado.output
  output:
    ptseq=touch(os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln_results", "pt_seqs")),
    flag=touch(os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln.done"))
  params:
    taxon=config["fln"]["taxon"],
    dbs=config["fln"]["dbs"],
    fasta=lambda wildcards: os.path.join(outdir, "Hints", "FLN", "fastas", "chunk_%s.fasta" % wildcards.chunk_id),
    outdir=os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}"),
    load=loadPreCmd(config.get("load", dict()).get("full_lengther_next", None)),
  threads: THREADS
  log: os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln.log")
  shell: """{params.load} mkdir -p {params.outdir} && ln -rst {params.outdir} {params.fasta} &&
  cd {params.outdir} && mv $(basename {params.fasta}) chunk.fasta &&
  full_lengther_next -w {threads} -g {params.taxon} -a {params.dbs} -f chunk.fasta 2>&1 > fln.log"""