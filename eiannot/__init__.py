import snakemake
import argparse
from .preparation.prepare import parse_samplesheet
# from .workflow import workflow  # TODO: this will be where we create the total workflow

def parser():
    parser = argparse.ArgumentParser()
    # TODO: implement the CLI parser

    return parser


def create_snakemake(args: argparse.ArgumentParser):

    configuration = args.configuration  # TODO: this must be following the JSON schema

    configuration = parse_samplesheet(configuration)

    workflow = ''