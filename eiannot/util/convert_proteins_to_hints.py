#!/usr/bin/env python3

import sys
import argparse
from Mikado.parsers import to_gff, GFF
from Mikado.transcripts import Transcript


__doc__ = """Converter of exonerate protein alignments into augustus GFF hints."""


def print_protein(transcript: Transcript, args: argparse.Namespace):

    if not isinstance(transcript, Transcript):
        return  # Uninitialised

    transcript.finalize()
    source = args.source or transcript.source
    # src=P;grp=ENSP00000013807.4.m1;pri=4;
    attributes = "src={args.src};grp={transcript.id};pri={args.priority}".format(**locals())
    exons = sorted(transcript.exons)
    introns = sorted(transcript.introns)
    for index in range(len(transcript.exons)):
        exon = exons[index]
        if index in (0, len(transcript.exons) -1):
            feature = "CDSpart"
        else:
            feature = "CDS"
        line = [transcript.chrom, source, feature,
                exon[0], exon[1], transcript.score,
                transcript.strand, ".", attributes]
        print(*line, sep="\t", file=args.out)
        if index < len(introns):
            intron = introns[index]
            line = [transcript.chrom, source, "intron",
                    intron[0], intron[1], transcript.score,
                    transcript.strand, ".", attributes]
            print(*line, sep="\t", file=args.out)
    return


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-p", "--priority", type=int, default=4)
    parser.add_argument("-s", "--source", default=None, help="Alternative source for the input data.")
    parser.add_argument("-src", default="P",
                        help="Value for the \"src\" field in the hints file. Default: %(default)s")
    parser.add_argument("proteins", type=to_gff)
    parser.add_argument("out", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    current = None

    for line in args.proteins:
        assert isinstance(line, GFF.GffLine)
        if line.header is True:
            continue
        elif line.is_transcript:
            print_protein(current, args)
            current = Transcript(line)
            continue
        elif line.is_exon:
            current.add_exon(line)
            continue
        else:
            continue

    print_protein(current, args)
    return


if __name__ == "__main__":
    main()
