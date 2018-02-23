from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule mikado_prepare:
	input:
		ref=rules.sanitize_reference.output,
		fai=rules.samtools_index_reference.output
	output:
		gtf=os.path.join(MIKADO_DIR, "mikado_prepared.gtf"),
		fa=os.path.join(MIKADO_DIR, "mikado_prepared.fasta")
	params:
		load=loadPreCmd(config.get("load", dict()).get("mikado", None)),
		cfg=CFG
	threads: THREADS
	message: "Preparing transcripts using mikado"
	shell: "{params.load} mikado prepare --start-method=spawn --procs={threads} --fasta={input.ref} --json-conf={params.cfg} -od {MIKADO_DIR} 2>&1"