from snakemake.shell import shell
import os
from eicore.external_process.snakemake_helper import loadPreCmd


config = snakemake.config

align_dir = config["aln_globals"]["align_dir"]
star_dir = os.path.join(align_dir, "star")
indexdir=os.path.join(star_dir, "index")
ref_trans = config["aln_globals"]["ref_trans"]
trans = "--sjdbGTFfile " + os.path.abspath(ref_trans) if ref_trans else ""
extra = config.get("extra", dict()).get("star_index", "")

load = loadPreCmd(config.get("load", dict()).get("star", None))

if not os.path.exists(indexdir):
    os.makedirs(indexdir)
elif not os.path.isdir(indexdir):
    raise ValueError("{} should be a directory.".format(indexdir))

shell("{params.load} "
      "cd {star_dir} && STAR "
      "--runThreadN {snakemake.threads} "
      "--runMode genomeGenerate "
      "--genomeDir {indexdir}"
      "{trans} {extra}"
      "--genomeFastaFiles {snakemake.input} > {snakemake.log} 2>&1")