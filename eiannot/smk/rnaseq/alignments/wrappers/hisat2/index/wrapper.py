import snakemake as snm
from snakemake.shell import shell
from eicore.external_process.snakemake_helper import loadPreCmd
import os


config = snakemake.config
load = loadPreCmd(config.get("load", dict()).get("hisat", None))
outdir = os.path.abspath(os.path.dirname(snakemake.output))
# TODO: get the right name
out_prefix = os.path.join(outdir, NAME)

shell("{load}"
      "hisat2-build -p {snakemake.threads}"
      "{snakemake.input} {out_prefix} > {snakemake.log} 2>&1")

shell("touch -f {snakemake.output}")