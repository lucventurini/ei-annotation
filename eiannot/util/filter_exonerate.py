#!/usr/bin/env python3

import sys
import argparse
from Mikado.parsers import to_gff
from Mikado.parsers.bed12 import Bed12Parser
from Mikado.transcripts import TranscriptChecker, Transcript
import pyfaidx
import operator
from Mikado.exceptions import IncorrectStrandError
import re


__doc__ = """Script to filter exonerate GFFs, by doing the following:
- truncate alignments with terminal introns over a certain size
- remove alignments if:
    - any of their internal introns is outside of the specified range
    - non-canonical junctions are present, which are not confirmed by selected verified junctions
    - 

"""


def _is_intron_valid(intron: (int, int), intron_index, minI: int, maxI: int, verified_introns, canonical_junctions):

    is_verified = (intron in verified_introns or intron_index in canonical_junctions)
    is_correct_length = (minI <= intron[1] - intron[0] + 1 <= maxI)

    return is_verified and is_correct_length


def evaluate(transcript: Transcript, args, verified_introns: set):

    if not transcript:
        print("Null transcript")
        return None

    transcript.finalize()
    if "Note" in transcript.attributes and "note" not in transcript.attributes:
        transcript.attributes["note"] = transcript.attributes["Note"]
        del transcript.attributes["Note"]

    if (args.min_coverage or args.min_identity) and ("note" in transcript.attributes or "Note" in transcript.attributes):
        note = transcript.attributes["note"].split("|")
        transcript.note = dict((re.search("^([^:]*):(.*)", no).groups()) for no in note)
        # try:
        #     transcript.note = dict((_[0], _[1]) for _ in transcript.attributes["note"])
        # except IndexError:
        #     raise IndexError(transcript.attributes["note"])
        if not ("cov" in transcript.note and "id" in transcript.note):
            # No coverage and identity!
            # print("No cov/iden for {}".format(transcript.id))
            return None
        transcript.note["cov"] = float(transcript.note["cov"])
        transcript.note["id"] = float(transcript.note["id"])
        if transcript.note["cov"] < args.min_coverage or transcript.note["id"] < args.min_identity:
            return None
    else:
        pass

    if transcript.monoexonic is True:
        # No further check is needed
        return transcript

    assert isinstance(args.genome, pyfaidx.Fasta)

    check_transcript = TranscriptChecker(transcript,
                                         args.genome[transcript.chrom][transcript.start -1 :transcript.end])
    check_transcript.finalize()
    try:
        check_transcript.check_strand()
    except IncorrectStrandError:
        # No correct strand found.
        pass

    t_introns = sorted(check_transcript.introns, key=operator.itemgetter(0))
    minI, maxE, maxM = args.minI, args.maxE, args.maxM

    # num = len(check_transcript.canonical_junctions)
    if len(t_introns) == 1:
        if not _is_intron_valid(t_introns[0], 0, minI, maxM, verified_introns, check_transcript.canonical_junctions):
            exon_lengths = dict((_[1] - _[0] + 1, _) for _ in transcript.exons)
            if max(exon_lengths) / transcript.cdna_length >= 0.8:  # TODO: maybe this has to be configurable?
                exon = exon_lengths[min(exon_lengths)]
                # Create new transcript object
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
                any(not _is_intron_valid(intron, pos, minI, maxM, verified_introns, check_transcript.canonical_junctions)
                for pos, intron in enumerate(t_introns[1:-1], 1))):
            return None
        is_start_correct = _is_intron_valid(t_introns[0], 0, minI, maxE, verified_introns,
                                            check_transcript.canonical_junctions)
        is_end_correct = _is_intron_valid(t_introns[-1], len(t_introns) -1, minI, maxE, verified_introns,
                                          check_transcript.canonical_junctions)

        if not (is_start_correct and is_end_correct):
            transcript.unfinalize()
            if not is_start_correct:
                transcript.remove_exon(transcript.exons[0])
            if not is_end_correct:
                transcript.remove_exon(transcript.exons[-1])
            transcript.unfinalize()

        transcript.finalize()
        return transcript


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-j", "--junctions", default=None, type=Bed12Parser)
    parser.add_argument("-minI", "--min-intron-length", dest="minI",
                        type=int, required=True)
    parser.add_argument("-maxM", "--max-intron-length-middle", dest="maxM",
                        type=int, required=True)
    parser.add_argument("-maxE", "--max-intronlength-ends", dest="maxE", type=int,
                        default=None)
    parser.add_argument("-minid", "--min-identity", dest="min_identity",
                        type=float, default=50)
    parser.add_argument("-mincov", "--min-coverage", dest="min_coverage",
                        type=float, default=80)
    parser.add_argument("-s", "--source", required=True, type=str)
    parser.add_argument("-g", "--genome", required=True, type=pyfaidx.Fasta)
    parser.add_argument("gff", type=to_gff)
    parser.add_argument("out", type=argparse.FileType("wt"), default=sys.stdout, nargs="?")
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

    print("##gff-version\t3", file=args.out)

    for line in args.gff:
        if line.header is True:
            continue
        elif line.feature == "match":
            if current:
                current = evaluate(current, args, introns)
                if current:
                    current.source = args.source
                    print(current.format("gff3"), file=args.out)
                    print("###", file=args.out)
            current = Transcript(line, intron_range=[args.minI, max(args.maxE, args.maxM)])
        elif line.feature == "match_part":
            current.add_exon(line)

    current = evaluate(current, args, introns)
    if current:
        current.source = args.source
        print(current.format("gff3"), file=args.out)


if __name__ == "__main__":
    main()
