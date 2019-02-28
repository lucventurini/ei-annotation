import os
from eiannot import load_pre_cmd


rule wig2hints:
  input: rules.combine_wigs.output
  output: os.path.join(outdir, "Alignments", "hints.ep.gff")
  params:
    load=loadPreCmd(config["load"]["augustus"]),
    priority=
  shell: """{params.load} bigWigToWig {input} /dev/stdout | wig2hints.pl --src=W --pri=3 --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --radius=4.5 > {output}"""