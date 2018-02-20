import os
import sys
from eicore.external_process.snakemake_helper import loadPreCmd

### We should have two paths: one for training, one for not training

if config.get("augustus", dict()).get("train", False) is True:
  


rule train