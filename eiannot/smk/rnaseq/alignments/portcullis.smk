from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import portcullisStrandOption
import os


rule portcullis_merge:
	input:
	  beds=expand(os.path.join(PORTCULLIS_DIR, "output", "portcullis_{aln_method}.pass.junctions.bed"), aln_method=ALIGN_RUNS),
	  tabs=expand(os.path.join(PORTCULLIS_DIR, "output", "portcullis_{aln_method}.pass.junctions.tab"), aln_method=ALIGN_RUNS)
	output:
	  bed=touch(PORTCULLIS_DIR, "output", "portcullis.merged.bed"),
	  tab=touch(PORTCULLIS_DIR, "output", "portcullis.merged.tab"),
	params: load=loadPreCmd(config.get("load", dict()).get("portcullis", None))
	log: os.path.join(PORTCULLIS_DIR, "output", "portcullis.merged.log")
	threads: 1
	message: "Taking intersection of portcullis results"
	run:
		if RUN_PORTCULLIS:
		  if len(input) == 1:
		    shell("cat {input.beds} > {output.bed}")
		    shell("cat {input.tabs} > {output.tab}")
		  else:
		    shell("{params.load} junctools set --prefix=portcullis_merged --output={output.tab} --operator=mean union {input.tabs} > {log} || touch {output.tab}")
		    shell("{params.load} junctools convert -if portcullis -of ebed -o {output.bed}  {output.tab} ")
		else:
			shell("touch {output}")

rule portcullis_prep:
	input:
	  ref=rules.sanitize_reference.output,
		aln_done=rules.align_all.output
	output: os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "1-prep", "portcullis.sorted.alignments.bam.bai")
	params:
		outdir=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "1-prep"),
		load=loadPreCmd(config.get("load", dict()).get("portcullis", None)),
		files=lambda wildcards: PORTCULLIS_IN[wildcards.aln_method],
	log: os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "prep.log")
	threads: THREADS
	message: "Using portcullis to prepare: {wildcards.aln_method}"
	shell: "{params.load} portcullis prep -o {params.outdir} -t {threads} {input.ref} {params.files} > {log} 2>&1"


rule portcullis_junc:
	input:
		bai=rules.portcullis_prep.output
	output: os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "2-junc", "{aln_method}.junctions.tab")
	params:
		prepdir=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "1-prep"),
		outdir=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "2-junc"),
		load=loadPreCmd(config.get("load", dict()).get("portcullis", None)),
		strand=lambda wildcards: portcullisStrandOption(wildcards.aln_method, config, "junc", SAMPLE_MAP)
	log: os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}-junc.log")
	threads: THREADS
	message: "Using portcullis to analyse potential junctions: {wildcards.aln_method}"
	shell: "{params.load} portcullis junc -o {params.outdir}/{wildcards.aln_method} {params.strand} -t {threads} {params.prepdir} > {log} 2>&1"

rule portcullis_filter:
	input: rules.portcullis_junc.output
	output:
		link=os.path.join(PORTCULLIS_DIR, "output", "portcullis_{aln_method}.pass.junctions.bed")
	params:
		outdir=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "3-filt"),
		prepdir=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "1-prep"),
		load=loadPreCmd(config.get("load", dict()).get("portcullis", None)),
		bed=os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "3-filt", "{aln_method}.pass.junctions.bed"),
		ss_gen="mkdir -p " + os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "3-filt") + " && gtf2bed.py " + REF_TRANS + " > " + os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "3-filt", "ref_juncs.bed") + " && " if REF_TRANS else "",
		trans="--reference=" + os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}", "3-filt", "ref_juncs.bed") if REF_TRANS else "",
		link_src=os.path.join("..", "portcullis_{aln_method}", "3-filt", "{aln_method}.pass.junctions.bed"),
		link_unfilt=os.path.join("..", "portcullis_{aln_method}", "2-junc", "{aln_method}.junctions.bed")
	log: os.path.join(PORTCULLIS_DIR, "portcullis_{aln_method}-filter.log")
	threads: THREADS
	message: "Using portcullis to filter invalid junctions: {wildcards.aln_method}"
	shell: "{params.load} {params.ss_gen} portcullis filter -o {params.outdir}/{wildcards.aln_method} --canonical={CANONICAL_JUNCS} --max_length={MAX_INTRON} {params.trans} --threads={threads} {params.prepdir} {input} > {log} 2>&1 && ln -sf {params.link_src} {output.link} || ln -sf {params.link_unfilt} {output.link} && touch -h {output.link}"
