import snakemake as snm
from snakemake.shell import shell
from eicore.external_process.snakemake_helper import loadPreCmd
import os
from eiannot.smk.rnaseq import hisatStrandOption, hisatInput


def hisatStrandOption(sample, SAMPLE_MAP):
    if SAMPLE_MAP[sample] == "fr-firststrand":
        return "--rna-strandness=RF"
    elif SAMPLE_MAP[sample] == "fr-secondstrand":
        return "--rna-strandness=FR"
    elif SAMPLE_MAP[sample] == "f":
        return "--rna-strandness=F"
    elif SAMPLE_MAP[sample] == "r":
        return "--rna-strandness=R"
    else:
        return ""


def hisatInput(sample, SAMPLE_MAP, INPUT_1_MAP, INPUT_2_MAP):
    if not seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
        return "-1 " + INPUT_1_MAP[sample] + " -2 " + INPUT_2_MAP[sample]
    else:
        return "-U " + INPUT_1_MAP[sample]

config = snakemake.config
load = loadPreCmd(config.get("load", dict()).get("hisat", None))
load_samtools = loadPreCmd(config.get("load", dict()).get("samtools", None))

indexdir = os.path.dirname(snakemake.input.index)
sample = snakemake.wildcards.sample
run = snakemake.wildcards.run

link_src = os.path.join("..", "hisat", "{}-{}".format(sample, run), "hisat.bam")
extra = config["align_methods"]["hisat"][int(run)]
ref_trans = config["aln_globals"]["ref_trans"]
if ref_trans:
    # Extract the reference splicing sites
    ss_gen = "hisat2_extract_splice_sites.py {ref_trans}  > "
    ss_gen += os.path.join(config["aln_globals"]["align_dir"], "hisat", "{sample}-{run}", "splice_sites.txt")
    shell("{load} {load_samtools} {ss_gen}")
    # Now add it to the align call
    trans = "--known-splicesite-infile={}".format(os.path.join(
        config["aln_globals"]["align_dir"], "hisat", "{}-{}".format(sample, run), "splice_sites.txt"))
else:
    trans = ""




strand=lambda wildcards: hisatStrandOption(wildcards.sample),
infiles=lambda wildcards: hisatInput(wildcards.sample)

if ss_gen:



shell("{load} {load_samtools}")

shell: "{params.load} {params.load_samtools} {params.ss_gen} hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {params.trans} {params.strand} {params.extra} -x {params.indexdir} --dta-cufflinks {params.infiles} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"