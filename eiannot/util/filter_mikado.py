#!/usr/bin/env python3

import argparse
import pandas as pd
import Mikado
import Mikado.transcripts
import networkx as nx
import itertools
from eiannot.util import build_pos_index
from functools import partial
try:
    import ujson as json
except ImportError:
    import json


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


def perc(string):
    string = float(string)
    assert 0 <= string <= 100
    return string


def load_blast_db(args):
    blast_db = pd.read_csv(args.blast, delimiter="\t",
                           names=["qseqid",
                                  'sseqid',
                                  'pident',
                                  'qstart',
                                  'qend',
                                  'sstart',
                                  'send',
                                  'qlen',
                                  'slen',
                                  'length',
                                  'nident',
                                  'mismatch',
                                  'positive',
                                  'gapopen',
                                  'gaps',
                                  'evalue',
                                  'bitscore'], index_col=["qseqid"], comment="#")

    blast_db = blast_db[blast_db.evalue <= args.evalue]
    blast_db["subject_coverage"] = 100 * (blast_db.send - blast_db.sstart + 1) / blast_db.slen
    blast_db["query_coverage"] = 100 * (blast_db.qend - blast_db.qstart + 1) / blast_db.qlen
    blast_db["coverage"] = blast_db[["subject_coverage", "query_coverage"]].apply(min, axis=1)
    blast_db["distance"] = pd.DataFrame({"send_dist": blast_db["slen"] - blast_db["send"],
                                         "qend_dist": blast_db["qlen"] - blast_db["qend"],
                                         "qstart": blast_db.qstart, "sstart": blast_db.sstart}).apply(max, axis=1)
    if 0 <= blast_db.pident.min() <= blast_db.pident.max() <= 1:
        if args.identity > 1:
            args.identity /= 100

    blast_db["Complete"] = (blast_db.coverage >= args.coverage) & (blast_db.distance <= args.distance) & \
                           (blast_db.pident >= args.identity)

    # Select ONLY THE BEST VALUE FOR EACH MATCH
    blast_db = blast_db.sort_values(
        ["Complete", "evalue", "coverage", "pident", "distance"], ascending=False).groupby(by=blast_db.index).head(1)

    blast_db = blast_db.assign(Status="Fragment")
    blast_db.loc[blast_db.Complete == True, "Status"] = "Complete"
    return blast_db


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--max-intron", dest="max_intron", type=int, default=10000)
    parser.add_argument("-cov", "--coverage", type=perc, default=80)
    parser.add_argument("-id", "--identity", type=perc, default=80)
    parser.add_argument("-e", "--evalue", type=float, default=10**-6)
    parser.add_argument("-d", "--distance", type=int, default=20)  # TODO change name
    parser.add_argument("--max-training", dest="max_training", default=2000, type=int)
    parser.add_argument("-cs", "--cross-testing", dest="cross_testing", default=0.2, type=float,
                        help="Fraction of training candidates to be used for cross testing.")
    # parser.add_argument("fln")
    parser.add_argument("mikado")
    parser.add_argument("blast")
    parser.add_argument("out_prefix")
    args = parser.parse_args()

    args.flank = abs(args.flank)
    blast_db = load_blast_db(args)
    metrics = pd.read_csv(args.mikado + ".metrics.tsv", sep="\t", index_col=0)
    bed12 = args.mikado + ".bed12"
    midx = args.mikado + ".bed12.midx"

    metrics["transcripts_per_gene"] = metrics.groupby("parent")["parent"].transform("count")
    merged = pd.merge(metrics, blast_db, left_index=True, right_index=True, validate="one_to_one", how="outer")
    merged.loc[merged.Status.isnan(), ]

    # Purge everything which has a long intron
    merged = merged[merged.max_intron_length <= args.max_intron]
    # Now that we have a complete table, it is easy to filter out

    # Now we have to exclude genes that are within 1000bps of another gene
    # Load from the MIDX the positions

    indexer, positions, gene_positions, genes = build_pos_index(midx)

    # This instead are the gold/silver/bronze categories for the training

    gold = merged[(
        (merged.Status.str.contains("complete", case=False)) &
         # Mikado and FLN account for the end of the ORF differently. This ensures that they mean the same codon.
        ((merged.combined_cds_length == merged.selected_cds_length) &
         (merged.three_utr_num_complete <= 1)) & (merged.five_utr_num_complete <= 2) &
        (merged.three_utr_num <= 2) & (merged.five_utr_num <= 3)
    )]

    training_candidates = gold[(
         (gold.combined_cds_fraction >= 0.5) &
         (gold.transcripts_per_gene == 1)  # & (gold.source_score == gold.source_score.max())  # Not sure about this
    )][[gold.columns[0], "parent"]]
    training_candidates = remove_genes_with_overlaps(training_candidates,
                                                     indexer,
                                                     positions,
                                                     gene_positions,
                                                     args.flank)



    

    nodes = blast_db[(blast_db.qseqid.isin(training_candidates[gold.columns[0]])) &
                         (blast_db.qseqid != blast_db.sseqid) &
                         (blast_db.sseqid.isin(training_candidates[gold.columns[0]])) &
                         (blast_db.coverage >= args.coverage) & (blast_db.pident >= args.identity)]

    # Now we are going to use NetworKit, as networkx would just make us cry
    node_index = pd.Series(nodes.qseqid.append(nodes.sseqid).unique())

    graph = nx.Graph()
    graph.add_edges_from(zip(nodes.sseqid, nodes.qseqid))

    to_keep = [list(_)[0] for _ in nx.connected_components(graph)]
    training_candidates = training_candidates[(training_candidates[gold.columns[0]].isin(to_keep)) |
                                              (~training_candidates[gold.columns[0]].isin(node_index))]

    # Now select only X candidates
    # Sample would fail if asked to select an X > df.shape[0], so we enforce the minimum of the two
    training_candidates = training_candidates.sample(min(args.max_training,
                                                         training_candidates.shape[0]))
    training_final = training_candidates.sample(frac=args.cross_testing)
    testing_final = training_candidates[~training_candidates.index.isin(training_final.index)]

    silver = merged[
        (merged.index.str.endswith(".1")) &
        (
        (~merged["parent"].isin(gold["parent"])) &
        ((merged.combined_cds_length == merged.selected_cds_length) &
         (merged.selected_cds_length >= 900) &
          (merged.three_utr_num_complete <= 1)) & (merged.five_utr_num_complete <= 2) &
         (merged.three_utr_num <= 2) & (merged.five_utr_num <= 3)
    )][[merged.columns[0], "parent"]]

    bronze = merged[(
        (merged.index.str.endswith(".1")) &
        (~merged["parent"].isin(silver["parent"])) &
        (~merged["parent"].isin(gold["parent"]))
    )][[merged.columns[0], "parent"]]

    merged["Training"] = merged[merged.columns[0]].isin(training_final[training_final.columns[0]])
    merged["Testing"] = merged[merged.columns[0]].isin(testing_final[testing_final.columns[0]])
    cat = partial(determine_category,
                  gold=set(gold.iloc[:, 0].unique()),
                  silver=set(silver.iloc[:, 0].unique()),
                  bronze=set(bronze.iloc[:, 0].unique()))

    if merged.shape[0] > 0:
        merged["Category"] = merged.apply(cat, axis=1)
    else:
        merged["Category"] = pd.Series()

    # Now write out the CSVs ...
    merged.to_csv(args.out_prefix + ".table.txt", sep="\t", index=False, header=True)

    # Now write out the GFFs

    with open("{args.out_prefix}.Training.gff3".format(**locals()), "wt") as training:
        print("##gff-version\t3", file=training)
        for gene in merged[merged.Training == True].parent.astype(str):
            print(genes[gene].format("gff3"), file=training)

    with open("{args.out_prefix}.Testing.gff3".format(**locals()), "wt") as testing:
        print("##gff-version\t3", file=testing)
        for gene in merged[merged.Training == True].parent.astype(str):
            print(genes[gene].format("gff3"), file=testing)

    for category in ("Gold", "Silver", "Bronze"):
        with open("{args.out_prefix}.{category}.gff3".format(**locals()), "wt") as out_gff:
            print("##gff-version\t3", file=out_gff)
            for gene in merged[merged.Category == category].parent.astype(str):
                print(genes[gene].format("gff3"), file=out_gff)

    return


main()
