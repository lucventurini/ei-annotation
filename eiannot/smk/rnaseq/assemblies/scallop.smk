from eicore.external_process.snakemake_helper import loadPreCmd
import os

rule scallop_all:
	input: expand(os.path.join(ASM_DIR, "output", "scallop-{run2}-{alrun}.gtf"), run2=SCALLOP_RUNS, alrun=ALIGN_RUNS)
	output: touch(os.path.join(ASM_DIR, "scallop.done"))


rule asm_scallop:
	input:
		bam=os.path.join(ALIGN_DIR, "output", "{alrun}.sorted.bam"),
		align=rules.align_all.output,
	output:
		gtf=os.path.join(ASM_DIR, "output", "cufflinks-{run2,\d+}-{alrun}.gtf")
	params:
		outdir=os.path.join(ASM_DIR, "scallop-{run2}-{alrun}"),
		gtf=os.path.join(ASM_DIR, "scallop-{run2}-{alrun}", "transcripts.gtf"),
		link_src=os.path.join("..", "scallop-{run2}-{alrun}", "transcripts.gtf"),
		load=loadPreCmd(config.get("load", dict()).get("scallop", None)),
		extra=lambda wildcards: config.get("asm_methods", dict()).get("scallop", [""]*len(int(wildcards.run2)+1)[int(wildcards.run2)],
		strand=lambda wildcards: tophatStrandOption(extractSample(wildcards.alrun), SAMPLE_MAP)
	log: os.path.join(ASM_DIR, "scallop-{run2}-{alrun}.log")
	threads: 1
	message: "Using Scallop to assemble (run {wildcards.run2}): {input.bam}"
	shell: """{params.load} mkdir -p {params.outdir} &&
	 scallop -i {input.bam} -o {params.gtf} {params.strand} {params.extra} > {log} 2>&1 && ln -sf {params.link_src} {output.gtf} && touch -h {output.gtf}"""
