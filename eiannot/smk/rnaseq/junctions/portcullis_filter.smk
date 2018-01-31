## The idea behind this filter is pretty simple:
## Take the "pass" file from Portcullis, remove any intron greater than 10kbps,
## put anything with score == 1 as "gold", and the rest as "silver"


import os
MAX_LENGTH = 10000  # Ideally this should be in the config!


rule filter:
  input: "portcullis.pass.bed"
  output:
    gold="portcullis.gold.bed",
    silver="portcullis.silved.bed"
  run:
    with open(input, "rt") as pbed, open(output[0], "wt") as gold, open(output[1], "wt") as silver:
      for line in pbed:
        fields = line.strip().split("\t")
        if len(fields) != 12:
          print(line, end='', file=gold)
          print(line, end='', file=silver)
          continue
        score = float(fields[4])
        length = int(fields[7]) - int(fields[6])
        if length >= MAX_LENGTH:
          continue
        if score == 1:
          print(line, end='', file=gold)
        else:
          print(line, end='', file=silver)
