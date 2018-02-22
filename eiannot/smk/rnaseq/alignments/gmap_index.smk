from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule gmap_index:
  input: rules.sanitize_reference.output
  output: os.path.join(ALIGN_DIR, "gmap", "index", NAME, NAME+".sachildguide1024")
	params:
	  load=loadPreCmd(config.get("load", dict()).get("gmap", None)),
	  dir=os.path.join(ALIGN_DIR, "gmap", "index")
	threads: 1
  log: os.path.join(ALIGN_DIR, "gmap.index.log")
  message: "Indexing genome with gmap"
  shell: "{params.load} gmap_build --dir={params.dir} --db={NAME} {input} > {log} 2>&1"