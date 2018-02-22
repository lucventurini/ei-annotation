from eicore.external_process.snakemake_helper import loadPreCmd
import os

rule asm_class:
	input:
		bam=os.path.join(ALIGN_DIR, "output", "{alrun}.sorted.bam"),
		align=rules.align_all.output,
		ref=REF
	output:
		link=os.path.join(ASM_DIR, "output", "class-{run2,\d+}-{alrun}.gtf"),
		gtf=os.path.join(ASM_DIR, "class-{run2,\d+}-{alrun}", "class-{run2}-{alrun}.gtf")
	params:
		outdir=os.path.join(ASM_DIR, "class-{run2}-{alrun}"),
		load=loadPreCmd(config.get("load", dict()).get("class", None)),
		extra=lambda wildcards: config.get("asm_methods", dict()).get("class", [""]*int(wildcards.run2))[int(wildcards.run2)],
		link_src=os.path.join("..", "class-{run2}-{alrun}", "class-{run2}-{alrun}.gtf")
	log: os.path.join(ASM_DIR, "class-{run2}-{alrun}.log")
	threads: THREADS
	message: "Using class to assemble (run {wildcards.run2}): {input.bam}"
	shell: "{params.load} {CLASS} --clean --force -c \"{params.extra}\" -p {threads} {input.bam} > {output.gtf} 2> {log} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"

rule class_all:
	input: expand(ASM_DIR+"/output/class-{run2}-{alrun}.gtf", run2=CLASS_RUNS, alrun=ALIGN_RUNS)
	output: touch(os.path.join(ASM_DIR, "class.done")
