from snakemake.shell import shell
import os


config = snakemake.config

outdir = os.path.dirname(snakemake.output.gtf)
load = loadPreCmd(config.get("load", dict()).get("class", None))
run = int(snakemake.wildcards.run2)
alrun = snakemake.wildcards.alrun
sample = snakemake.wildcards.sample
extra = config.get("asm_methods", dict()).get("class", [""] * run)[run]
link_src = os.path.join("..", "class-{run2}-{alrun}", "class-{run2}-{alrun}.gtf")



