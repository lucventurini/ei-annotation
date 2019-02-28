import os
from eiannot import load_pre_cmd

rule polish_rm:
  input: os.path.join(outdir, "RModel", "genome-families.fa")
  output: os.path.join(outdir, "RModel", "genome-families.polished.fa")
