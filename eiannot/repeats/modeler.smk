## Part of ei-annot. This SnakeMake will execute RepeatModeller on the selected genome.

import os
from eicore.external_process.snakemake_helper import loadPreCmd


rule modeler:
  input: os.path.join(outdir, "RModel", "genome.nhr")
  output: os.path.join(outdir, "RModel", "genome-families.fa")
  params:
    load=loadPreCmd(config["load"]["repeatmodeler"]),
    name=os.path.join(outdir, "RModel", "genome"),
    dir=os.path.join(outdir, "RModel")
  threads: THREADS
  shell: """{params.load} cd {params.dir} && RepeatModeler -pa {threads} -engine ncbi -database {params.name} && ln -s RM*/genome-families.fa ."""