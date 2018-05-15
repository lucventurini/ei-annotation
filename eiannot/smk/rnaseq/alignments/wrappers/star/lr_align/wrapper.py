from snakemake.shell import shell
from eicore.external_process.snakemake_helper import loadPreCmd
import os

config = snakemake.config


load = loadPreCmd(config.get("load", {}).get("star", ""))
outdir = os.path.abspath(os.path.dirname(snakemake.output.bam))

if not os.path.exists(outdir):
    os.makedirs(outdir)
elif not os.path.isdir(outdir):
    raise TypeError("{} is not a directory".format(outdir))

indexdir = os.path.abspath(os.path.dirname(snakemake.input.index))

ref_trans = config["aln_globals"]["ref_trans"]
if ref_trans:
    trans = "--sjdbGTFfile " + os.path.abspath(ref_trans)
else:
    trans = ""

# TODO: change this
infiles= starLongInput(snakemake.wildcards.lsample)

shell("{load}"
      "STARlong STARlong --runThreadN {snakemake.threads}"
      "--runMode alignReads --outSAMattributes NH HI NM MD"
      "--readNameSeparator space --outFilterMultimapScoreRange 1 "
      "--outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4"
      "--scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 "
      "--alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 "
      "--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 "
      "--genomeDir {indexdir}"
      "{infiles} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif"
      "--alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON}"
      "{trans} --outFileNamePrefix {params.outdir}/ > {snakemake.log} 2>&1"
      )

