## Part of ei-annot. This SnakeMake will create the database for RepeatModeller.
## TODO: at the moment the engine is fixed to NCBI. I think it is appropriate given that it is endorsed
## Directly by the creators of Repeat(Masker|Modeler), but we might want to make it more flexible

import os
from eicore.external_process.snakemake_helper import loadPreCmd

rule build_rm_db:
  input: REF
  output: os.path.join(outdir, "RModel", "genome.nhr")
  params:
    load=loadPreCmd(config["load"]["exonerate"]),
    name=os.path.join(outdir, "RModel", "genome"),
    refdir=os.path.join(outdir, "RModel")
  shell: """{params.load} mkdir -p {params.refdir} && BuildDatabase -name {params.name} -engine ncbi {input}"""


