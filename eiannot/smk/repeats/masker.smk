import os
from eicore.external_process.snakemake_helper import loadPreCmd

def get_input():
  pass


rule masker:
  input:
    genome=REF,
    library=os.path.join(outdir, "reference_and_rmodeler.lib")
  output: os.path.join(outdir, "genome.fa.masked")
  params:
    load=loadPreCmd(config["load"]["repeatmasker"])
  threads: THREADS
