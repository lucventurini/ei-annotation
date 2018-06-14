import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule extract:
  input:  # TODO
  output: os.path.join(outdir, "reference.rm.fa")
  params:
    load=loadPreCmd(config["load"]["repeatmasker"]),
    species=config["repeats"]["species"]
  shell: """{params.load} queryRepeatDatabase.pl --species {params.species} > {output}"""