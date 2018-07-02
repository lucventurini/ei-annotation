#!/usr/bin/env python3

import sys
import argparse
from Mikado.parsers import to_gff
from Mikado.parsers.bed12 import Bed12Parser
from Mikado.transcripts import TranscriptChecker, Transcript
import pyfaidx
import operator
from collections import OrderedDict as odict


__doc__ = """Script to filter exonerate GFFs, by doing the following:
- truncate alignments with terminal introns over a certain size
- remove alignments if:
    - any of their internal introns is outside of the specified range
    - non-canonical junctions are present, which are not confirmed by selected verified junctions
    - 

"""


def _is_intron_valid(intron: (int, int), minI: int, maxI: int, verified_introns, canonical_junctions):

    is_verified = (intron in verified_introns or intron in canonical_junctions)

    return is_verified and minI <= intron[1] - intron[0] + 1 <= maxI


def evaluate(transcript: Transcript, args, verified_introns: set):

    if transcript.monoexonic is True:
        # No check is needed
        return transcript

    assert isinstance(args.genome, pyfaidx.Fasta)

    transcript = TranscriptChecker(transcript,
                                   args.genome[transcript.chrom][transcript.start - 1:transcript.end])
    transcript.finalize()

    t_introns = sorted(transcript.introns, key=operator.itemgetter(0))
    minI, maxE, maxM = args.minI, args.maxE, args.maxM

    if len(t_introns) == 2:
        if not _is_intron_valid(t_introns[0], minI, maxM, verified_introns, transcript.canonical_junctions):
            exon_lengths = dict((_[1] - _[0] + 1, _) for _ in transcript.exons)
            if max(exon_lengths) / transcript.cdna_length >= 0.8:
                # TODO: flexible fraction maybe?
                exon = exon_lengths[min(exon_lengths)]
                transcript.unfinalize()
                transcript.remove_exon(exon)
                transcript.finalize()
                return transcript
            else:
                return None
        else:
            return transcript
    else:
        # Check internal introns. If any of them is wrong, return None
        if (len(t_introns) > 2 and
                any(not _is_intron_valid(intron, minI, maxM, verified_introns, transcript.canonical_junctions)
                for intron in t_introns[1:-1])):
            return None
        if not _is_intron_valid(t_introns[0], minI, maxE, verified_introns, transcript.canonical_junctions):
            transcript.unfinalize()
            transcript.remove_exon(transcript.exons[0])
        if not _is_intron_valid(t_introns[-1], minI, maxE, verified_introns, transcript.canonical_junctions):
            transcript.unfinalize()
            transcript.remove_exon(transcript.exons[-1])
        transcript.finalize()
        return transcript


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-j", "--junctions", default=None, type=Bed12Parser)
    parser.add_argument("-min", "--min-intron-length", dest="minI",
                        type=int, required=True)
    parser.add_argument("-maxM", "--max-intron-length-middle", dest="maxM",
                        type=int, required=True)
    parser.add_argument("-maxE", "--max-intronlength-ends", dest="maxE", type=int,
                        default=None)
    parser.add_argument("-g", "--genome", required=True, type=pyfaidx.Fasta)
    parser.add_argument("gff", type=to_gff)
    parser.add_argument("out", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()

    current = None

    if not args.maxE:
        args.maxE = args.maxM

    if args.minI >= min(args.maxE, args.maxM):
        raise ValueError("The minimum intron length must be lower than maximum intron lengths!{} vs {}, {}".format(
            args.minI, args.maxE, args.maxM
        ))

    introns = set()
    if args.junctions is not None:
        for bed in args.junctions:
            introns.add((bed.chrom, bed.thick_start, bed.thick_end))

    for line in args.gff:
        if line.header is True:
            continue
        elif line.feature == "match":
            if current:
                current = evaluate(current, args, introns)
                if current:
                    print(current, file=args.out)
            current = Transcript(line, intron_range=[args.minI, max(args.maxE, args.maxM)])
        elif line.feature == "match_part":
            current.add_exon(line)

    current = evaluate(current, args, introns)
    if current:
        print(current, file=args.out)


if __name__ == "__main__":
    main()
