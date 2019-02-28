import os
from eiannot import load_pre_cmd


rule merge_wigs_plus:
  input: rules.bam2wig.output.plus
  output: os.path.join(outdir, "Alignments", "final.forward.wig")
  params:
    load=""
  shell: """{params.load} merge_wigs.py -o {params.output} {params.input}"""


rule merge_wigs_plus:
  input: rules.bam2wig.output.minus
  output: os.path.join(outdir, "Alignments", "final.reverse.wig")
  params:
    load=""
  shell: """{params.load} merge_wigs.py -o {params.output} {params.input}"""