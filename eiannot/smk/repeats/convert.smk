import os
from eicore.external_process.snakemake_helper import loadPreCmd

rule convert:
  input: rules.masker.output
  output: os.path.join(outdir, "repeats.gff3")
  run:
    ## Following advice from http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.IncorporateRepeats
    src = "RepeatMasker"
    feature = "nonexonpart"
    with open("{input}") as inp, open("{output}", "wt") as out:
      for iii, line in enumerate(inp):
        if iii < 3:
          continue
        fields = line.rstrip().split()
        chrom = fields[4]
        start, end = fields[5:6]
        print(chrom, src, feature, start, end, 0, ".", ".", "src={}".format(src), file=out, sep="\t")
