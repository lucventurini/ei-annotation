from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import hisatStrandOption, hisatInput
import os


rule align_hisat_index:
	input: REF
	output: os.path.join(ALIGN_DIR, "hisat", "index", NAME + ".done")
	params: load=loadPreCmd(config.get("load", dict()).get("hisat", None))
	log: os.path.join(ALIGN_DIR, "hisat.index.log")
	threads: THREADS
	message: "Indexing genome with hisat"
	shell: "{params.load} hisat2-build -p {threads} {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1 && touch {output}"


rule align_hisat:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.align_hisat_index.output
	output:
		bam=os.path.join(ALIGN_DIR, "hisat", "{sample}-{run,\d+}", "hisat.bam"),
		link=os.path.join(ALIGN_DIR, "output", "hisat-{sample}-{run,\d+}.bam")
	params:
		indexdir=os.path.join(ALIGN_DIR, "hisat", "index", NAME),
		load=loadPreCmd(config.get("load", dict()).get("hisat", None), config.get("load", dict()).get("samtools", None)),
		link_src=os.path.join("..", "hisat", "{sample}-{run}", "hisat.bam"),
		extra=lambda wildcards: config["align_methods"]["hisat"][int(wildcards.run)],
		ss_gen="hisat2_extract_splice_sites.py " + REF_TRANS + " > " + os.path.join(ALIGN_DIR, "hisat", "{sample}-{run}", "splice_sites.txt") + "&&" if REF_TRANS else "",
		trans="--known-splicesite-infile=" + os.path.join(ALIGN_DIR, "hisat", "{sample}-{run}", "splice_sites.txt") if REF_TRANS else "",
		strand=lambda wildcards: hisatStrandOption(wildcards.sample),
		infiles=lambda wildcards: hisatInput(wildcards.sample)
	log: ALIGN_DIR+"/hisat-{sample}-{run}.log"
	threads: THREADS
	message: "Aligning input with hisat (sample {wildcards.sample} - run {wildcards.run})"
    	shell: "{params.load} {params.load_samtools} {params.ss_gen} hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {params.trans} {params.strand} {params.extra} -x {params.indexdir} --dta-cufflinks {params.infiles} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule hisat_all:
	input: expand(ALIGN_DIR+"/output/hisat-{sample}-{run}.bam", sample=SAMPLES, run=HISAT_RUNS)
	output: ALIGN_DIR+"/hisat.done"
	shell: "touch {output}"