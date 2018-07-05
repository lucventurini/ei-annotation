#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import Mikado
import itertools
import sqlite3
import collections
import magic
from Mikado.exceptions import CorruptIndex
from functools import partial


def build_pos_index(index_name):

    positions = dict()

    wizard = magic.Magic(mime=True)

    if wizard.from_file("{0}".format(index_name)) == b"application/gzip":
        raise CorruptIndex("Invalid index file")
    try:
        conn = sqlite3.connect("{0}".format(index_name))
        cursor = conn.cursor()
        tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
        if sorted(tables) != sorted([("positions",), ("genes",)]):
            raise CorruptIndex("Invalid database file")
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid database file")

    genes = dict()
    try:
        for counter, obj in enumerate(cursor.execute("SELECT * from positions")):
            chrom, start, end, gid = obj
            if chrom not in positions:
                positions[chrom] = collections.defaultdict(list)
            positions[chrom][(start, end)].append(gid)
            genes[gid] = (chrom, start, end)
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid index file. Rebuilding.")

    indexer = collections.defaultdict(list).fromkeys(positions)
    for chrom in positions:
        indexer[chrom] = Mikado.utilities.intervaltree.IntervalTree.from_tuples(
            positions[chrom].keys()
        )

    return indexer, positions, genes


def remove_genes_with_overlaps(training_candidates: pd.DataFrame,
                               indexer: dict, positions: dict, genes: dict, flank: int):

    to_remove = set()

    for num in range(training_candidates.shape[0]):
        gid = training_candidates["parent"].iloc[num]  # Get the gene name
        found = set()
        chrom, start, end = genes[gid]
        for interval in indexer[chrom].find(start-flank, end+flank):
            found = set.union(found, positions[chrom][(interval[0], interval[1])])
        found.remove(gid)
        if len(found) > 0:
            to_remove.add(gid)
    training_candidates = training_candidates[~
                              training_candidates[training_candidates.columns[0]].isin(to_remove)]
    return training_candidates


def determine_category(row, gold, silver, bronze):

    if row[0] in gold:
        return "Gold"
    elif row[0] in silver:
        return "Silver"
    elif row[0] in bronze:
        return "Bronze"
    else:
        return "NA"


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("fln")
    parser.add_argument("mikado")
    parser.add_argument("out")
    args = parser.parse_args()

    args.flank = abs(args.flank)
    bed12 = args.mikado + ".bed12"
    metrics = args.mikado + ".metrics.tsv"
    midx = args.mikado + ".gff3.midx"

    fln = pd.read_csv(args.fln, sep="\t")
    transcripts = [Mikado.transcripts.Transcript(_) for _ in Mikado.parsers.bed12.Bed12Parser(bed12)
                   if _.header is False]
    orfs = list(itertools.chain(
        *[[(t.id, _.thick_start, _.thick_end) for _ in t.get_internal_orf_beds()] for t in transcripts]))
    orf_df = pd.DataFrame(orfs, columns=[fln.columns[0], "mikado_orf_start", "mikado_orf_end"])
    metrics = pd.read_csv(metrics, sep="\t")
    metrics["transcripts_per_gene"] = metrics.groupby("parent")["parent"].transform("count")
    merged = pd.merge(pd.merge(fln, orf_df, on=fln.columns[0]),
                              
                              metrics, left_on=fln.columns[0], right_on=metrics.columns[0])

    # Now that we have a complete table, it is easy to filter out

    # Now we have to exclude genes that are within 1000bps of another gene
    # Load from the MIDX the positions

    indexer, positions, genes = build_pos_index(midx)

    # This instead are the gold/silver/bronze categories for the training

    gold = merged[(
        (merged.Status.str.match("complete", case=False)) &
        ((merged["ORF_start"] == merged["mikado_orf_start"]) & (merged["ORF_end"] == merged["mikado_orf_end"])) &
        ((merged.combined_cds_length == merged.selected_cds_length) & (merged.selected_cds_length>=300) &
         (merged.five_utr_num == 2) & (merged.three_utr_num == 1))
    )]  # [[merged.columns[0], "parent"]]

    training_candidates = gold[(
         (gold.combined_cds_fraction >= 0.5) &
         (gold.transcripts_per_gene == 1) & (gold.source_score == gold.source_score.max()))
    ][[gold.columns[0], "parent"]]
    training_candidates = remove_genes_with_overlaps(training_candidates, indexer, positions, genes, args.flank)

    silver = merged[(
        (~merged["parent"].isin(gold["parent"])) &
        ((merged.combined_cds_length == merged.selected_cds_length) & (merged.selected_cds_length>=300) &
         (merged.five_utr_num == 2) & (merged.three_utr_num == 1))
    )][[merged.columns[0], "parent"]]

    bronze = merged[((~merged["parent"].isin(silver["parent"])) &
                     (~merged["parent"].isin(gold["parent"])) &
                     (merged["max_intron_length"] <= 50000))][[merged.columns[0], "parent"]]

    merged["Training"] = merged[merged.columns[0]].isin(training_candidates[training_candidates.columns[0]])
    cat = partial(determine_category,
                  gold=set(gold.iloc[:, 0].unique()),
                  silver=set(silver.iloc[:, 0].unique()),
                  bronze=set(bronze.iloc[:, 0].unique()))

    if merged.shape[0] > 0:
        merged["Category"] = merged.apply(cat, axis=1)
    else:
        merged["Category"] = pd.Series()

    # Now write out the CSVs ...
    merged.to_csv(args.out + ".table.txt", sep="\t", index=False, header=True)
    merged[[merged.columns[0], "parent", "Training", "Category"]].to_csv(
        args.out + ".list.txt", sep="\t", index=False, header=True
    )
    return


main()
