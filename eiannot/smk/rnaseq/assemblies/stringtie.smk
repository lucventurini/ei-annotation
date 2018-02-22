from eicore.external_process.snakemake_helper import loadPreCmd
import os

rule asm_stringtie:
	input:
		bam=os.path.join(ALIGN_DIR, "output", "{alrun}.sorted.bam"),
		align=rules.align_all.output
	output:
		link=os.path.join(ASM_DIR, "output", "stringtie-{run2,\d+}-{alrun}.gtf"),
		gtf=os.path.join(ASM_DIR, "stringtie-{run2,\d+}-{alrun}", "stringtie-{run2}-{alrun}.gtf")
	params:
		load=loadPreCmd(config.get("load", dict()).get("stringtie", None)),
		extra=lambda wildcards: config.get("asm_methods", dict()).get("stringtie", [""]*int(wildcards.run2))[int(wildcards.run2)],
		gtf=os.path.join(ASM_DIR, "stringtie-{run2}-{alrun}", "stringtie-{run2}-{alrun}.gtf"),
		trans="-G " + REF_TRANS if REF_TRANS else "",
		link_src=os.path.join("..", "stringtie-{run2}-{alrun}", "stringtie-{run2}-{alrun}.gtf")
	log: os.path.join(ASM_DIR, "stringtie-{run2}-{alrun}.log")
	threads: THREADS
	message: "Using stringtie to assemble (run {wildcards.run2}): {input.bam}"
	shell: "{params.load} stringtie {input.bam} -l Stringtie_{wildcards.run2}_{wildcards.alrun} -f 0.05 -m 200 {params.extra} {params.trans} -o {params.gtf} -p {threads} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"

rule stringtie_all:
	input: expand(ASM_DIR+"/output/stringtie-{run2}-{alrun}.gtf", run2=STRINGTIE_RUNS, alrun=ALIGN_RUNS)
	output: touch(os.path.join(ASM_DIR, "stringtie.done"))
