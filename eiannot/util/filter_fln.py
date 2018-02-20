#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import Mikado
import itertools


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("fln")
    parser.add_argument("bed12")
    parser.add_argument("out")
    args = parser.parse_args()

    fln = pd.read_csv(args.fln, sep="\t")

    complete = fln[(fln.Status.str.match("complete", case=False)) & (~fln.Status.str.match("putative", case=False))]
    transcripts = [Mikado.transcripts.Transcript(_) for _ in Mikado.parsers.bed12.Bed12Parser(args.bed12)
                   if _.header is False]
    orfs = list(itertools.chain(
        *[[(t.id, _.thick_start, _.thick_end) for _ in t.get_internal_orf_beds()] for t in transcripts]))
    orf_df = pd.DataFrame(orfs, columns=[complete.columns[0], "mikado_orf_start", "mikado_orf_end"])

    merged = pd.merge(complete, orf_df, on=complete.columns[0])
    correct = merged[
        (merged["ORF_start"] == merged["mikado_orf_start"]) & (merged["ORF_end"] == merged["mikado_orf_end"])]

    assert isinstance(correct, pd.DataFrame)
    correct.to_csv(path_or_buf=args.out, sep="\t", header=True, index=False, index_label=False)

    return


main()


