from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import tophatStrandOption, tophatInput
import os


rule tophat_all:
	input: expand(ALIGN_DIR+"/output/tophat-{sample}-{run}.bam", sample=SAMPLES, run=TOPHAT_RUNS)
	output: touch(os.path.join(ALIGN_DIR, "tophat.done"))


rule align_tophat_index:
	input: ref=rules.sanitize_reference.output
	output: os.path.join(ALIGN_DIR, "tophat", "index", NAME + ".4.bt2")
	params:
		idxdir=os.path.join(ALIGN_DIR, "tophat", "index" + NAME),
		load=loadPreCmd(config.get("load", dict()).get("tophat", None))
	log: os.path.join(ALIGN_DIR, "tophat.index.log")
	threads: 1
	message: "Indexing genome with tophat"
	shell: "{params.load} bowtie2-build {input.ref} {params.idxdir} > {log} 2>&1"


rule align_tophat:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.align_tophat_index.output
	output:
		link=os.path.join(ALIGN_DIR, "output", "tophat-{sample}-{run,\d+}.bam")
	params:
		outdir=os.path.join(ALIGN_DIR, "tophat", "{sample}-{run}"),
		bam=os.path.join(ALIGN_DIR, "tophat", "{sample}-{run}", "accepted_hits.bam"),
		indexdir=os.path.join(ALIGN_DIR, "tophat", "index", NAME),
		load=loadPreCmd(config.get("load", dict()).get("tophat", None)),
		extra=lambda wildcards: config.get("align_methods", dict()).get("tophat", [""]*int(wildcards.run))[int(wildcards.run)],
		link_src=os.path.join("..", "tophat", "{sample}-{run}", "accepted_hits.bam"),
		strand=lambda wildcards: tophatStrandOption(wildcards.sample, SAMPLE_MAP),
		infiles=lambda wildcards: tophatInput(wildcards.sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP),
		#trans="--GTF=" + REF_TRANS if REF_TRANS else ""
		trans=""
	log: os.path.join(ALIGN_DIR, "tophat-{sample}-{run}.log")
	threads: THREADS
	message: "Aligning RNAseq data with tophat (sample - {wildcards.sample} - run {wildcards.run})"
	shell: "{params.load} tophat2 --output-dir={params.outdir} --no-sort-bam --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} {params.strand} {params.trans} {params.extra} {params.indexdir} {params.infiles} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"

