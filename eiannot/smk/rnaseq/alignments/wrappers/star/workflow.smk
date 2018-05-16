from eiannot.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import starCompressionOption, starInput
import os

wrapper_dir=os.path.dirname(os.path.abspath(__name__))


rule align_star_index:
	input: rules.sanitize_reference.output
	output: os.path.join(config["aln_globals"]["align_dir"], "star", "index", "SAindex")
	log: os.path.join(config["aln_globals"]["align_dir"], "star.index.log")
	threads: THREADS
	message: "Indexing genome with star"
	wrapper: "file://" + os.path.abspath(wrapper_dir, "index")


rule star_all:
	input: expand(os.path.join(config["aln_globals"]["align_dir"], "output", "star-{sample}-{run}.bam"), sample=SAMPLES, run=STAR_RUNS)
	output: touch(os.path.join(config["aln_globals"]["align_dir"], "star.done")


rule align_star:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.align_star_index.output
	output:
		bam=os.path.join(config["aln_globals"]["align_dir"], "star", "{sample}-{run,\d+}", "Aligned.out.bam"),
		link=os.path.join(config["aln_globals"]["align_dir"], "output", "star-{sample}-{run,\d+}.bam")
	log: os.path.join(ALIGN_DIR_FULL, "star-{sample}-{run}.log")
	threads: int(THREADS)
	message: "Aligning input with star (sample {wildcards.sample} - run {wildcards.run})"
	wrapper: "file://" + os.path.join(wrapper_dir, "align")


rule lr_star_all:
	input: expand(os.path.join(config["aln_globals"]["align_dir"], "lr_output", "lr_star-{lsample}-{lrun}.gtf", lsample=L_SAMPLES, lrun=LR_STAR_RUNS)
	output: touch(os.path.join(config["aln_globals"]["align_dir"], "lr_star.done"))

rule align_lr_star:
	input:
		index=rules.align_star_index.output,
		reads=lambda wildcards: L_INPUT_MAP[wildcards.lsample]  # TODO: change this
	output:
	  bam=os.path.join(config["aln_globals"]["align_dir"], "lr_star", "{lsample}-{lrun}", "Aligned.out.bam"),
	log: os.path.join(config["aln_globals"]["align_dir"], "lr_star-{lsample}-{lrun}.log")
	threads: THREADS
	message: "Mapping long reads to the genome with star (sample: {wildcards.lsample} - run: {wildcards.lrun})"
	wrapper: "file://" + os.path.join(wrapper_dir, "lr_align")

rule starbam2gtf:
	input: rules.align_lr_star.output.bam
	output: gtf=os.path.join(config["aln_globals"]["align_dir"], "lr_output", "lr_star-{lsample}-{lrun}.gtf")
	params: load=loadPreCmd(config.get("load", {}).get("mikado", ""))
	message: "Converting STAR long reads from BAM to GTF (sample: {wildcards.lsample} - run: {wildcards.lrun})"
	shell: "{params.load} bam2gtf.py {input} > {output.gtf}"