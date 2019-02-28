#!/usr/bin/env python3

"""This script has the function of linking all the reads inside the sample sheet,
and to create the sample objects to be put inside the configuration."""


from ..abstract import ShortSample, LongSample
import os
import csv
import re


def parse_samplesheet(samplesheet, configuration):

    """This function will parse the samplesheet, create the sample objects, and put them in the configuration itself."""

    configuration["long_reads"] = dict()
    configuration["short_reads"] = dict()
    outdir = os.path.join(configuration["outdir"], "inputs", "reads")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif os.path.exists(outdir) and not os.path.isdir(outdir):
        raise OSError("Read directory is not a directory at all: {}".format(outdir))

    if samplesheet is None:
        return configuration

    with open(samplesheet) as sheet:
        for line in csv.reader(sheet, delimiter="\t"):
            if line[0].lstrip().startswith("#"):
                continue  # Ignore comments!
            label, read1, read2, type, strandedness = line
            label = re.sub("\s", "_", label)  # Remove spaces!
            if type in ("illumina", "short"):
                sample = ShortSample(read1, read2, label, outdir, strandedness=strandedness)
                tag = "short_reads"
            else:
                sample = LongSample(read1, label, outdir, strandedness, type)
                tag = "long_reads"
            if label in configuration[tag]:  # Double label!
                raise KeyError(
                    "{short_tag} read label {label} was found at least twice in the sample sheet. Please recheck it.".format(
                        short_tag=tag.split("_")[0].capitalize(),
                        **locals()))
            configuration[tag][label] = sample

    return configuration
