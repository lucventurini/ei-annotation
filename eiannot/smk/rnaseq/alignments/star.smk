from eicore.external_process.snakemake_helper import loadPreCmd
from eicore.smk.rnaseq import starCompressionOption, starInput
import os


rule star_all:
	input: expand(os.path.join(ALIGN_DIR, "output", "star-{sample}-{run}.bam"), sample=SAMPLES, run=STAR_RUNS)
	output: touch(os.path.join(ALIGN_DIR, "star.done")


rule align_star_index:
	input: rules.sanitize_reference.output
	output: os.path.join(ALIGN_DIR, "star", "index", "SAindex")
	params:
		indexdir=os.path.join(ALIGN_DIR_FULL, "star", "index"),
		load=loadPreCmd(config.get("load", dict()).get("star", None)),
		trans="--sjdbGTFfile " + os.path.abspath(REF_TRANS) if REF_TRANS else "",
		extra=config.get("extra", dict()).get("star_index", "")
	log: os.path.join(ALIGN_DIR_FULL, "star.index.log")
	threads: THREADS
	message: "Indexing genome with star"
	shell: "{params.load} cd {ALIGN_DIR_FULL}/star && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} {params.trans} --genomeFastaFiles {input} {params.extra} > {log} 2>&1 && cd {CWD}"


rule align_star:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.align_star_index.output
	output:
		bam=os.path.join(ALIGN_DIR, "star", "{sample}-{run,\d+}", "Aligned.out.bam"),
		link=os.path.join(ALIGN_DIR, "output", "star-{sample}-{run,\d+}.bam")
	params:
		outdir=os.path.join(ALIGN_DIR_FULL, "star", "{sample}-{run}"),
		indexdir=os.path.join(ALIGN_DIR_FULL, "star", "index"),
		load=loadPreCmd(config.get("load", dict()).get("star", None)),
		extra=lambda wildcards: config.get("align_methods", dict()).get("star", [""]*int(wildcards.run))[int(wildcards.run)],
		link_src=os.path.join("..", "star", "{sample}-{run}", "Aligned.out.bam"),
		trans="--sjdbGTFfile " + os.path.abspath(REF_TRANS) if REF_TRANS else "",
		rfc=lambda wildcards: starCompressionOption(wildcards.sample, EXT_MAP),
		infiles=lambda wildcards: starInput(wildcards.sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP)
	log: os.path.join(ALIGN_DIR_FULL, "star-{sample}-{run}.log")
	threads: int(THREADS)
	message: "Aligning input with star (sample {wildcards.sample} - run {wildcards.run})"
	shell: "{params.load} cd {params.outdir}; STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} {params.rfc} {params.infiles} --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM XS NM MD --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} {params.trans} --alignMatesGapMax {MAX_INTRON} --outFileNamePrefix {params.outdir}/ {params.extra} > {log} 2>&1 && cd {CWD} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"

