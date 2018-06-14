import os
from eicore.external_process.snakemake_helper import loadPreCmd

rule polish_rm:
  input: os.path.join(outdir, "RModel", "genome-families.fa")
  output: os.path.join(outdir, "RModel", "genome-families.polished.fa")
