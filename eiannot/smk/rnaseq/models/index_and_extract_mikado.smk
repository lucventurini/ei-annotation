import os
from eicore.external_process.snakemake_helper import loadPreCmd

# A problem here is that we should select ONE mikado output. Instead Daijin returns multiple options in this regard

CHUNK_ARRAY = [str(_).zfill(3) for _ in range(1, config["fln"]["chunks"] + 1)]

rule index_mikado:
  input:
    mikado=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3")
  output:
    midx=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3.midx")
  params:
    load=loadPreCmd(config.get("load", dict()).get("mikado", None))
  log: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "index.log")
  shell: """{params.load} mikado compare -r {input} -l {log} --index"""

rule extract_mikado_fasta:
  input:
    mikado=os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
    genome=REF
  output: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.cdna.fasta")
  params:
    load=loadPreCmd(config.get("load", dict()).get("gffread", None))
  shell: """{params.load} gffread -g {input.genome} -w {output} -C {input.mikado}"""

rule convert_mikado_to_bed12:
  input:
    os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
  output: os.path.join(outdir, "Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.bed12")
  params:
    load=loadPreCmd(config.get("load", dict()).get("mikado", None))
  shell: """{params.load} mikado util convert -of bed12 {input.mikado} {output}"""

rule split_mikado:
	input: tr=rules.extract_mikado_fasta.output
	output: os.path.join(outdir, "Hints", "FLN", "fastas", "split.done")
	params:
		out=os.path.join(outdir, "Hints", "FLN", "fastas", "chunk"),
		outdir=os.path.join(outdir, "Hints", "FLN", "fastas"),
		chunks=config["fln"]["chunks"],
		load=loadPreCmd(config.get("load", dict()).get("mikado", None))
	threads: 1
	message: "Splitting fasta: {input.tr}"
	shell: """{params.load} mkdir -p {params.outdir} &&
	split_fasta.py -m {params.chunks} {input.tr} {params.out} && touch {output}"""