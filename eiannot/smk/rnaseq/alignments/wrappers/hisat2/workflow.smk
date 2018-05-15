from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule hisat_all:
  # Change samples and runs
	input: expand(os.path.join(config["aln_globals"]["align_dir"], "output", "hisat-{sample}-{run}.bam"), sample=SAMPLES, run=HISAT_RUNS)
	output: touch(os.path.join(os.path.join(config["aln_globals"]["align_dir"], "hisat2.done"))

rule align_hisat_index:
	input: rules.sanitize_reference.output
	# Find the correct name for NAME
	output: os.path.join(config["aln_globals"]["align_dir"], "hisat", "index", NAME + ".done")
	log: os.path.join(config["aln_globals"]["align_dir"], "hisat.index.log")
	threads: THREADS
	message: "Indexing genome with hisat"
	wrapper: "file://" + os.path.abspath(os.path.join(os.path.dirname(__name__), "index"))


rule align_hisat:
	input:
	  # TODO: change the way we handle this, across the RNA-Seq
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.align_hisat_index.output
	output:
		bam=os.path.join(config["aln_globals"]["align_dir"], "hisat", "{sample}-{run,\d+}", "hisat.bam"),
		link=os.path.join(config["aln_globals"]["align_dir"], "output", "hisat-{sample}-{run,\d+}.bam")
	log: ALIGN_DIR+"/hisat-{sample}-{run}.log"
	threads: THREADS
	message: "Aligning input with hisat (sample {wildcards.sample} - run {wildcards.run})"
	wrapper: "file://" + os.path.abspath(os.path.join(os.path.dirname(__name__), "align"))

