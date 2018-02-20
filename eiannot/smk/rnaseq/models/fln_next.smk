import os
import subprocess
from eicore.external_process.snakemake_helper import loadPreCmd


rule extract_mikado_fasta:
  # A problem here is that we should select ONE mikado output. Instead Daijin returns multiple
  input:
    mikado=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
    genome=REF
  output: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.cds.fasta")
  params:
    load=loadPreCmd(config["load"]["gffread"])
  shell: """{params.load} gffread -g {input.ref} -w {output} {input.mikado}"""

rule convert_mikado_to_bed12:
  input:
    os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
  output: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.bed12")
  params:
    load=loadPreCmd(config["load"]["mikado"])
  shell: """{params.load} mikado util convert -of bed12 {input.mikado} {output}"""


rule fln_all:
	input: expand(os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln_results", "pt_seqs"), chunk_id=CHUNK_ARRAY)
	output:
	  flag=os.path.join(outdir, "Hints", "FLN", "fln.done"),
	  table=os.path.join(outdir, "Hints", "FLN", "fln.table.txt")
	shell: "cat {input} | sed '2,/^Query_id/d' > {output.table} && touch {output.flag}"


rule filter_fln:
  input:
    table=rules.fln_all.output.table,
    bed12=rules.convert_mikado_to_bed12.output
  output: os.path.join(outdir, "Hints", "FLN", "fln.filtered_table.txt")
  params:
    load=loadPreCmd(config["load"]["mikado"])
  shell: """{params.load}  """




rule split_mikado:
	input: tr=rules.extract_mikado_fasta.output
	output: os.path.join(outdir, "Hints", "FLN", "fastas", "split.done")
	params:
		outdir=os.path.join(outdir, "Hints", "FLN", "fastas"),
		chunks=config["fln"]["chunks"],
		load=loadPreCmd(config["load"]["mikado"])
	threads: 1
	message: "Splitting fasta: {input.tr}"
	shell: """{params.load} mkdir -p {params.outdir} &&
	split_fasta.py -m {params.chunks} {input.tr} {params.outdir} && touch {output}"""

rule fln:
  input:
    split=rules.split_fa.output
  output:
    ptseq=os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln_results", "pt_seqs"),
    flag=os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln.done"),
  params:
    taxon=config["fln"]["taxon"],
    outdir=os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}"),
    load=loadPreCmd(config["load"]["full_lengther_next"]),
  threads: THREADS
  log: os.path.join(outdir, "Hints", "FLN", "Chunks", "chunk-{chunk_id}", "fln.log")
  shell: """{params.load} mkdir -p {params.outdir} && ln -rst {params.outdir} {input} &&
  cd {params.outdir} && mv $(basename {input}) chunk.fasta &&
  full_lengther_next -w {threads} -g {params.taxon} -f chunks.fasta 2>&1 > fln.log &&
  touch {params.outdir}/fln.done"""