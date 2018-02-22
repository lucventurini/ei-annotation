import functools
import subprocess
import re
import os
from eicore.external_process.snakemake_helper import loadPreCmd


# Generic functions


def seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
    s = SAMPLE_MAP[sample]
    if s == "f" or s == "r":
        if not INPUT_2_MAP[sample]:
            return True
        else:
            i=0
            # Should throw an error here
    return False


def extractSample(align_run):
    parts = align_run.split("-")
    return parts[1]


# TopHat specific functions


def tophatStrandOption(sample, SAMPLE_MAP):
    if SAMPLE_MAP[sample] == "f":
        return "--library-type=fr-secondstrand"
    elif SAMPLE_MAP[sample] == "r":
        return "--library-type=fr-firststrand"
    else:
        return "--library-type=" + SAMPLE_MAP[sample]


def tophatInput(sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP):
    if not seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
        return INPUT_1_MAP[sample] + " " + INPUT_2_MAP[sample]
    else:
        return INPUT_1_MAP[sample]


# Trinity specific functions

def trinityInput(sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP):
    if not seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
        return "--left={} --right={}".format(INPUT_1_MAP[sample], INPUT_2_MAP[sample])
    else:
        return "--single={} ".format(INPUT_1_MAP[sample])


@functools.lru_cache(maxsize=4, typed=True)
def getTrinityVersion(config):
    cmd = loadPreCmd(config.get("load", dict()).get("trinity", None))
    cmd += "Trinity --version && set -u"
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    version = [_ for _ in output.split("\n") if re.match("Trinity version:", _)][0]
    version = re.sub("Trinity version: [^_]*_(r|v)", "", version)
    return version


def trinityParameters(config, sample, REF, TGG_MAX_MEM):
    version = getTrinityVersion(config)
    if str(TGG_MAX_MEM).isdigit() is True:
        TGG_MAX_MEM = int(TGG_MAX_MEM)
        if TGG_MAX_MEM >= 1000:
            TGG_MAX_MEM = "{}G".format(int(TGG_MAX_MEM/1000))
    if version.startswith("201"):   # Old versions
        reads = trinityInput(sample)
        memory = "--JM={}".format(TGG_MAX_MEM)
        cmd = "{reads} --genome={ref} {memory} --genome_guided_use_bam".format(reads=reads, ref=REF, memory=memory)
    else:
        memory="--max_memory={}".format(TGG_MAX_MEM)
        cmd = "--full_cleanup {memory} --genome_guided_bam".format(memory=memory)
    return cmd


def trinityStrandOption(sample, SAMPLE_MAP):
    if SAMPLE_MAP[sample] == "fr-firststrand":
        return "--SS_lib_type=RF"
    elif SAMPLE_MAP[sample] == "fr-secondstrand":
        return "--SS_lib_type=FR"
    elif SAMPLE_MAP[sample] == "f":
        return "--SS_lib_type=F"
    elif SAMPLE_MAP[sample] == "r":
        return "--SS_lib_type=R"
    else:
        return ""


def gmap_intron_lengths(config, MAX_INTRON):
    cmd = "{} gmap --help".format(loadPreCmd(config.get("load", dict()).get("gmap", "")))
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.read().decode()
    if "--max-intronlength-middle" in output:
        return "--max-intronlength-middle={mi} --max-intronlength-ends={mi}".format(mi=MAX_INTRON)
    else:
        return "--intronlength={}".format(MAX_INTRON)

# STAR specific functions

def starCompressionOption(sample, EXT_MAP):
    if EXT_MAP[sample] == ".gz":
        return "--readFilesCommand zcat"
    elif EXT_MAP[sample] == ".bz":
        return "--readFilesCommand bzcat"
    else:
        return ""

def starInput(sample, INPUT_1_MAP, INPUT_2_MAP, SAMPLE_MAP):
    if not seSample(sample, SAMPLE_MAP, INPUT_2_MAP):
        return "--readFilesIn " + os.path.abspath(INPUT_1_MAP[sample]) + " " + os.path.abspath(INPUT_2_MAP[sample])
    else:
        return "--readFilesIn " + os.path.abspath(INPUT_1_MAP[sample])


def starLongInput(sample, L_INPUT_MAP):
    return "--readFilesIn " + os.path.abspath(L_INPUT_MAP[sample])

# HISAT2 specific options


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

# Portcullis options

@functools.lru_cache(maxsize=8, typed=True)
def portcullisHelp(config, step):
    command = loadPreCmd(config.get("load", dict()).get("portcullis", None))
    cmd = subprocess.Popen("{} portcullis {} --help".format(command, step),
                           shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    return cmd


def portcullisStrandOption(run, config, step, SAMPLE_MAP):
    parts=run.split("-")
    sample=parts[1]
    cmd = portcullisHelp(config, step)
    if not any("strandedness" in _ for _ in cmd.split("\n")):
        return ""
    else:
        if SAMPLE_MAP[sample] == "fr-firststrand":
            return "--strandedness=firststrand"
        elif SAMPLE_MAP[sample] == "fr-secondstrand":
            return "--strandedness=secondstrand"
        else:
            return "--strandedness=unstranded"