# import snakemake as snm
from snakemake.shell import shell
from eicore.external_process.snakemake_helper import loadPreCmd
import os
from eicore.smk.rnaseq import seSample


def starCompressionOption(read_file):
    if read_file.endswith(".gz"):
        return "--readFilesCommand zcat"
    elif read_file.endswith(".bz") or read_file.endswith(".bz2"):
        return "--readFilesCommand bzcat"
    else:
        return ""

def starInput(sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP):
    if not seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
        return "--readFilesIn " + os.path.abspath(INPUT_1_MAP[sample]) + " " + os.path.abspath(INPUT_2_MAP[sample])
    else:
        return "--readFilesIn " + os.path.abspath(INPUT_1_MAP[sample])

config = snakemake.config
params = snakemake.params

load = loadPreCmd(config.get("load", dict()).get("star", None))
indexdir = params.get("indexdir")
outdir = os.path.join(config["aln_globals"]["align_dir"],
                      "star", "{}-{}".format(snakemake.wildcards.get("sample"), snakemake.wildcards.get("run")))

rfc = starCompressionOption(snakemake.input.reads[1])
ref_trans = config["aln_globals"]["ref_trans"]
# TODO: change this
infiles=starInput(snakemake.wildcards.sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP)

run_num = int(snakemake.wildcards.get("run"))


extra = config.get("align_methods", {}).get("star", [""] * run_num)[run_num]

if "align_methods" in config:
    if "star" in config["align_methods"]:
        extra = config["align_methods"]["star"][int(snakemake.wildards.get("run")]
    else:
        extra = ""
else:
    extra = ""

link_src = os.path.join("..", "star", "{}-{}".format(snakemake.wildcards.get("sample"),
                                                     snakemake.wildcards.get("run")), "Aligned.out.bam")

if not os.path.exists(outdir):
    os.makedirs(outdir)
elif not os.path.isdir(outdir):
    raise TypeError("{} should be a directory, not a file".format(outdir))

shell("{load} cd {outdir};"
      "STAR --runThreadN {snakemake.threads}"
      "{rfc}"
      "--runMode alignReads"
      "--genomeDir {indexdir}"
      "--outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM XS NM MD "
      "--outSAMstrandField intronMotif"
      "--alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON}"
      "--alignMatesGapMax {MAX_INTRON}"
      "--outFileNamePrefix {outdir}/"
      "{extra}"
      "{params.infiles} {params.trans} {params.extra}"
      "> {log} 2> &1"
      )

# Link the output
shell("ln -sf {link_src} {snakemake.output.link}")
# Touch output
shell("touch -h {snakemake.output.link}")
