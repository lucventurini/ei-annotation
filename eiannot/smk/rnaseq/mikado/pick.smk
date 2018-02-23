from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule mikado_collect_stats:
	input:
	    mikado=expand(os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.stats"), mode=MIKADO_MODES)
	output: os.path.join(MIKADO_DIR, "pick", "comparison.stats")
	params:
	    load=loadPreCmd(config.get("load", dict()).get("mikado", None))
	threads: 1
	message: "Collecting mikado statistics"
	shell: "{params.load} {ASM_COLLECT} {input.mikado} > {output}"


rule mikado_pick:
	input:
		gtf=rules.mikado_prepare.output.gtf,
		db=rules.mikado_serialise.output
	output:
		loci=os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.gff3")
	log: os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.pick.err")
	params:
		cfg=CFG,
		load=loadPreCmd(config.get("load", dict()).get("mikado", None)),
		outdir=os.path.join(MIKADO_DIR, "pick", "{mode}")
	threads: THREADS
	message: "Running mikado picking stage"
	shell: "{params.load} mikado pick --source Mikado_{wildcards.mode} --mode={wildcards.mode} --procs={threads} --start-method=spawn --json-conf={params.cfg} -od {params.outdir} --loci_out mikado-{wildcards.mode}.loci.gff3 -lv INFO {input.gtf} -db {input.db} > {log} 2>&1"


rule mikado_stats:
    input:
        rules.mikado_pick.output.loci
    output:
        stats=os.path.join(MIKADO_DIR, "pick", "{mode}", "mikado-{mode}.loci.stats")
    params: load=loadPreCmd(config.get("load", dict()).get("mikado", None))
    shell: "{params.load} mikado util stats {input} > {output}"
