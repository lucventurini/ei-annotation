#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import Mikado
import Mikado.transcripts
from Mikado.parsers.GFF import GffLine
from Mikado.parsers import to_gff
from Mikado.scales.gene_dict import GeneDict
from Mikado.scales.compare import check_index, create_index
from Mikado.utilities.log_utils import create_null_logger, create_default_logger

# from eiannot.util import build_pos_index
import os
import re


def prepare_rows(transcript: Mikado.transcripts.Transcript,
                 score, source, category, print_cds=False):

    lines = []

    for row in Mikado.transcripts.transcript.create_lines_cds(transcript,
                                                              with_introns=True,
                                                              to_gtf=False, all_orfs=False, transcriptomic=False):
        row = GffLine(row)

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

    if args.loci.endswith("midx"):
        midx = args.loci
    else:
        midx = os.path.join(os.path.dirname(args.loci),
                            re.sub("\.(gff3*|bed12)", ".gff3.midx", os.path.basename(args.loci)))
        if not os.path.exists(midx):
            create_index(to_gff(args.loci), create_null_logger(), midx, ref_gff=True,
                         exclude_utr=False, protein_coding=False)

    try:
        check_index(args.loci, create_null_logger())
    except Mikado.exceptions.CorruptIndex:
        os.remove(midx)
        create_index(to_gff(args.loci), create_null_logger(), midx, ref_gff=True,
                     exclude_utr=False, protein_coding=False)

    genes = GeneDict(midx)
    table = pd.read_csv(args.table, delimiter="\t", index_col=[0])

    grouped = table.groupby("parent")
    for parent in grouped.groups:
        gene = genes[parent]
        group = grouped.get_group(parent)
        for tidx in group.index:
            row = group.loc[tidx]
            try:
                tid = row.tid
            except AttributeError:
                raise AttributeError(row)
            try:
                transcript = gene[tid]
            except KeyError:
                raise KeyError(group.index, row.tid)
            lines = prepare_rows(transcript, category="all",
                                 score=args.all_score, source=args.all_source,
                                 print_cds=args.cds)
            print(*lines, sep="\n", file=args.outfile)
            if pd.isna(row.Category):
                continue
            elif row.Category == "Gold":
                source, score, category = args.gold_source, args.gold_score, "gold"
            elif row.Category == "Silver":
                source, score, category = args.silver_source, args.silver_score, "silver"
            elif row.Category == "Bronze":
                source, score, category = args.bronze_source, args.bronze_score, "bronze"
            else:
                raise ValueError("Unrecognised category: {}".format(row.Category))
            lines = prepare_rows(transcript, category=category,
                                 score=score, source=source,
                                 print_cds=args.cds)
            print(*lines, sep="\n", file=args.outfile)


main()
