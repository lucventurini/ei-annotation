import os
import sys
from eicore.external_process.snakemake_helper import loadPreCmd

### We should have two paths: one for training, one for not training

if config.get("augustus", dict()).get("train", False) is True:

  rule convert_training_to_gb:
    input:
      gff=rules.extract_training_gff.output,
      fa=rules.sanitize_reference.output  # We might use the SM reference as well
    output: os.path.join(OU, "training_models.gb")
    params:
      load=loadPreCmd(config.get("load", dict()).get("augustus", None)),
      flank=2000,  # Move into configuration
      keep_mask=  # Potentially we could put "--softmasked" here
    shell: """{params.load} gff2gbSmallDNA.pl {input.gff} {input.softmasked} {params.flank} > {output}"""

  rule split_training:
    input: rules.convert_models_gb.output
    output:
      train=,
      test=
    params:
      load=loadPreCmd(config.get("load", dict()).get("augustus", None)),
      size=
    log:
    shell: """{params.load} randomSplit {input} {params.size} 2> {log} > {log}"""

  rule train_augustus:
    input: rules.split_training.output.train
    output:
    params:
      load=loadPreCmd(config.get("load", dict()).get("augustus", None)),
    log:
    shell: """{params.load} """
