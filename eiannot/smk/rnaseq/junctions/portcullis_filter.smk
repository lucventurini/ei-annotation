## The idea behind this filter is pretty simple:
## Take the "pass" file from Portcullis, remove any intron greater than 10kbps,
## put anything with score == 1 as "gold", and the rest as "silver"


import os
import subprocess
from eicore.external_process.snakemake_helper import loadPreCmd
MAX_LENGTH = 10000  # Ideally this should be in the config!

rule filter:
  input: os.path.join(out_dir, "Daijin", "4-portcullis", "output", "portcullis.merged.tab")
  output:
    gold=os.path.join(outdir, "portcullis.gold.bed"),
    silver=os.path.join(outdir, "portcullis.silved.bed"),
    tmp=temp(os.path.join(outdir, "temp.gff3"))
  params:
    load=loadPreCmd(config["load"]["portcullis"])
  run:
    called=subprocess.call("{params.load} junctools convert -of igff -o {output.tmp} {input}", shell=True)
    if called > 0:
      raise OSError
    with open(output.tmp, "rt") as pgff, open(output.gold, "wt") as gold, open(output.silver) as silver:
      for line in pgff:
        line = line.strip()
        fields = line.split("\t")
        if fields[0][0] == "#":  # Comment
          print(line, file=gold)
          print(line, file=silver)
          continue
        try:
          start, end = int(fields[3]), int(fields[4])
        except TypeError, ValueError as exc:
          continue
        length = end - start + 1
        if length > MAX_LENGTH:
          continue
        score = float(fields[5])
        if score == 1:
          print(line, file=gold)
        else:
          print(line, file=silver)
