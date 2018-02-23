from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import portcullisStrandOption
import os

rule mikado_cfg:
	input:
		asm=rules.asm_all.output,
		lr=rules.lreads_all.output,
		portcullis=rules.portcullis_merge.output,
		cfg=CFG,
		ref=rules.sanitize_reference.output
	output:
		mikado=os.path.join(OUT_DIR, "mikado.yaml")
	params:
		load=loadPreCmd(config.get("load", dict()).get("mikado", None)),
		scoring=config["mikado"]["pick"]["scoring_file"],
		junctions="--junctions={}".format(rules.portcullis_merge.output.bed)
	log: os.path.join(OUT_DIR, "mikado.yaml.log")
	threads: 1
	message: "Creating Mikado configuration file"
	shell: "{params.load} mikado configure --gff={MIKADO_IN_STR} --labels={MIKADO_LABEL_STR} --strand-specific-assemblies={MIKADO_SS_STR} {params.junctions} --scoring {params.scoring} --reference={input.ref} --external={input.cfg} {output} 2> {log}"