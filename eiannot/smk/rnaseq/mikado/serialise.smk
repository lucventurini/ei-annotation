from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule mikado_serialise:
	input:
		blast=rules.blast_all.output,
		orfs=orf_out,
		fai=rules.samtools_index_reference.output,
		transcripts=rules.mikado_prepare.output.fa
	output:
	  db=os.path.join(MIKADO_DIR, "mikado.db")
	log: os.path.join(MIKADO_DIR, "mikado_serialise.err")
	params:
	  cfg=CFG,
		blast="--xml=" + os.path.join(BLAST_DIR, "xmls") if len(BLASTX_TARGET) > 0 else "",
		load=loadPreCmd(config.get("load", dict()).get("mikado", None)),
		blast_target="--blast_targets=" + os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa") if len(BLASTX_TARGET) > 0 else ""
	threads: THREADS
	message: "Running Mikado serialise to move numerous data sources into a single database"
	run:
	    if config.get("mikado", dict()).get("use_prodigal") is True:
	        params.orfs = "--orfs="+str(rules.prodigal.output)
	    elif config.get("transdecoder", dict()).get("execute", True) is True:
	        params.orfs = "--orfs="+str(rules.transdecoder_pred.output)
	    else:
	        params.orfs = ""
	    shell("{params.load} mikado serialise {params.blast} {params.blast_target} --start-method=spawn --transcripts={input.transcripts} --genome_fai={input.fai} --json-conf={params.cfg} --force {params.orfs} -od {MIKADO_DIR} --procs={threads} > {log} 2>&1")