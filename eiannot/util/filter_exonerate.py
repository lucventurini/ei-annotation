#!/usr/bin/env python3

import sys
import argparse
from Mikado.parsers import to_gff
from Mikado.parsers.bed12 import Bed12Parser
from Mikado.transcripts import Transcript
import pyfaidx


__doc__ = """Script to filter exonerate GFFs, by doing the following:
- truncate alignments with terminal introns over a certain size
- remove alignments if:
    - any of their internal introns is outside of the specified range
    - non-canonical junctions are present, which are not confirmed by selected verified junctions
    - 

"""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-j", "--junctions", default=None, type=Bed12Parser)
    parser.add_argument("-ir", "--intron-range", dest="intron_range",
                        required=True, type=tuple)
    parser.add_argument("-g", "--genome", required=True, type=pyfaidx.Fasta)
    parser.add_argument("gff", type=to_gff)
    parser.add_argument("out", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()

    current = None

    introns = set()
    if args.junctions is not None:
        for bed in args.junctions:
            introns.add((bed.chrom, bed.thick_start, bed.thick_end))


    for line in args.gff:
        if line.header is True:
            continue
        elif line.feature == "match":
            if current:
                evaluate(current, args, introns)
            current = Transcript(line)
        elif line.feature == "match_part":
            current.add_exon(line)

    evaluate(current, args, introns)