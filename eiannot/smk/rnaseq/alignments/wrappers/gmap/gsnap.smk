from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import tophatInput
import os


rule align_gsnap:
	input:
		r1=lambda wildcards: INPUT_1_MAP[wildcards.sample],
		index=rules.gmap_index.output
	output:
		bam=os.path.join(ALIGN_DIR, "gsnap", "{sample}-{run,\d+}", "gsnap.bam"),
		link=os.path.join(ALIGN_DIR, "output", "gsnap-{sample}-{run,\d+}.bam")
	params:
		load=loadPreCmd(config.get("load", dict()).get("gmap", None), config.get("load", dict()).get("samtools", None)),
		extra=lambda wildcards: config.get("align_methods", dict()).get("gsnap", [""]* int(wildcards.run))[int(wildcards.run)],
		link_src=os.path.join("..", "gsnap", "{sample}-{run}", "gsnap.bam"),
		# Can use tophat function safely here
		infiles=lambda wildcards: tophatInput(wildcards.sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP))
	log: os.path.join(ALIGN_DIR, "gsnap-{sample}-{run}.log")
	threads: THREADS
	message: "Aligning RNAseq with gsnap (sample {wildcards.sample} - run {wildcards.run})"
	shell: "{params.load} gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {params.extra} --novelsplicing=1 --localsplicedist={MAX_INTRON} --nthreads={threads} --format=sam --npaths=20 {params.infiles} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule gsnap_all:
	input: expand(os.path.join(ALIGN_DIR, "output", "gsnap-{sample}-{run}.bam"), sample=SAMPLES, run=GSNAP_RUNS)
	output: touch(os.path.join(ALIGN_DIR, "gsnap.done"))
