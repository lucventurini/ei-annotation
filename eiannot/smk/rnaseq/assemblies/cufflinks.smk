from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import tophatStrandOption


rule cufflinks_all:
	input: expand(os.path.join(ASM_DIR, "output", "cufflinks-{run2}-{alrun}.gtf"), run2=CUFFLINKS_RUNS, alrun=ALIGN_RUNS)
	output: touch(os.path.join(ASM_DIR, "cufflinks.done"))


rule asm_cufflinks:
	input:
		bam=os.path.join(ALIGN_DIR, "output", "{alrun}.sorted.bam"),
		align=rules.align_all.output,
		ref=rules.sanitize_reference.output
	output:
		gtf=os.path.join(ASM_DIR, "output", "cufflinks-{run2,\d+}-{alrun}.gtf")
	params:
		outdir=os.path.join(ASM_DIR, "cufflinks-{run2}-{alrun}"),
		gtf=os.path.join(ASM_DIR, "cufflinks-{run2}-{alrun}", "transcripts.gtf"),
		link_src=os.path.join("..", "cufflinks-{run2}-{alrun}", "transcripts.gtf"),
		load=loadPreCmd(config.get("load", dict()).get("cufflinks", None)),
		extra=lambda wildcards: config.get("asm_methods", dict()).get("cufflinks", [""]*len(int(wildcards.run2)+1)[int(wildcards.run2)],
		trans="--GTF-guide=" + REF_TRANS if REF_TRANS else "",
		strand=lambda wildcards: tophatStrandOption(extractSample(wildcards.alrun), SAMPLE_MAP)
	log: os.path.join(ASM_DIR, "cufflinks-{run2}-{alrun}.log")
	threads: THREADS
	message: "Using cufflinks to assemble (run {wildcards.run2}): {input.bam}"
	shell: "{params.load} cufflinks --output-dir={params.outdir} --num-threads={threads} {params.trans} {params.strand} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --no-update-check {params.extra} {input.bam} > {log} 2>&1 && ln -sf {params.link_src} {output.gtf} && touch -h {output.gtf}"
