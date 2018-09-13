#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import Mikado
from . import build_pos_index


def prepare_rows(transcript: Mikado.transcripts.Transcript,
                 score, source, category, print_cds=False):

    lines = []

    for row in Mikado.transcripts.transcript.create_lines_cds(transcript,
                                                              with_introns=True,
                                                              to_gtf=False, all_orfs=False, transcriptomic=False):
        assert isinstance(row, Mikado.parsers.GFF.GffLine)

        if row.is_transcript:
            continue
        elif row.feature == "intron":
            feature = "intron"

        elif row.feature == "CDS":
            if print_cds is False:
                continue
            elif row.start in (transcript.combined_cds_start, transcript.combined_cds_end):
                feature = "CDSpart"
            elif row.end in (transcript.combined_cds_start, transcript.combined_cds_end):
                feature = "CDSpart"
            else:
                feature = "CDS"
        elif row.feature == "exon":
            if row.start == transcript.start or row.end == transcript.end:
                feature = "exonpart"
            else:
                feature = "exon"
        else:
            continue

        line = [transcript.chrom, "mikado_{}".format(category), feature,
                row.start, row.end, ".", row.strand, ".",
                "group={id};source={source};pri={score}".format(source=source, id=transcript.id, score=score)]

        lines.append("\t".join(str(_) for _ in line))
    return lines


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-gs", "--gold-score", dest="gold_score", type=int)
    parser.add_argument("--gold-source", dest="gold_source", default="M")
    parser.add_argument("-ss", "--silver-score", dest="silver_score", type=int)
    parser.add_argument("--silver-source", dest="silver_source", default="M")
    parser.add_argument("-bs", "--bronze-score", dest="bronze_score", type=int)
    parser.add_argument("--bronze-source", dest="bronze_source", default="E")
    parser.add_argument("-as", "--all-score", dest="all_score", type=int)
    parser.add_argument("--all-source", dest="all_source", default="E")
    parser.add_argument("--cds", action="store_true", default=False)
    parser.add_argument("loci")
    parser.add_argument("table")
    parser.add_argument("outfile", type=argparse.FileType("wt"), default=sys.stdout, nargs="?")
    args = parser.parse_args()

    table = pd.read_csv(args.table, delimiter="\t", index_col=[0])
    indexer, positions, gene_positions, genes = build_pos_index(args.loci)
    grouped = table.groupby("parent")
    for parent in grouped.groups:
        gene = genes[parent]
        group = grouped.get_group(parent)
        for tid in group.index:
            row = group.loc[tid]
            transcript = gene[tid]
            lines = prepare_rows(transcript, category="all",
                                 score=args.all_score, source=args.all_source,
                                 print_cds=args.cds)
            print(*lines, sep="\n", file=args.outfile)
            if row.category.isna():
                continue
            elif row.category.astype(str) == "Gold":
                source, score, category = args.gold_source, args.gold_score, "gold"
            elif row.category.astype(str) == "Silver":
                source, score, category = args.silver_source, args.silver_score, "silver"
            elif row.category.astype(str) == "Bronze":
                source, score, category = args.bronze_source, args.bronze_score, "bronze"
            else:
                raise ValueError("Unrecognised category: {}".format(row.category.astype(str)))
            lines = prepare_rows(transcript, category=category,
                                 score=score, source=source,
                                 print_cds=args.cds)
            print(*lines, sep="\n")


main()
