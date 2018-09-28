#!/usr/bin/env python3

__doc__ = """Quick script to filter portcullis tab junctions into proper Augustus junctions."""


import pandas as pd
import argparse
import subprocess as sp
import tempfile
import os


def is_threshold(arg):

    err = "Invalid threshold, it must be a valid float number between 0 and 1"
    try:
        arg = float(arg)
    except ValueError:
        raise ValueError(err)
    if not 0 <= arg <= 1:
        raise ValueError(err)
    return arg


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("tab")
    parser.add_argument("-s", "--suffix", nargs=2, default=["gold", "silver"],
                        help="Suffices for the pass/fail. Default: gold, silver")
    parser.add_argument("--sources", nargs=2, default=["E", "E"],
                        help="Tag for the src field in the GFF. Default: %(default)s")
    parser.add_argument("-p", "--priorities", nargs=2, type=int, default=[6, 4],
                        help="Priorities for pass/fail. Default: 6, 4")
    parser.add_argument("-mi", "--max-intron-length", default=10000, type=int, dest="mi")
    parser.add_argument("-t", "--threshold", nargs=1, type=is_threshold, default=1)
    parser.add_argument("prefix")
    args = parser.parse_args()

    pass_gff3 = "{args.prefix}.gold.{source}{priority}.gff3".format(source=args.sources[0],
                                                                    priority=args.priorities[0],
                                                                    **locals())
    fail_gff3 = "{args.prefix}.silver.{source}{priority}.gff3".format(source=args.sources[1],
                                                                      priority=args.priorities[1],
                                                                      **locals())

    if not os.path.exists(os.path.dirname(args.prefix)):
        os.makedirs(os.path.dirname(args.prefix))

    with tempfile.NamedTemporaryFile(mode="wt") as temp_pass:
        with tempfile.NamedTemporaryFile(mode="wt") as temp_fail:
            tab = pd.read_csv(args.tab, delimiter="\t")
            pass_df = tab[(tab.score >= args.threshold) & (tab.size <= args.mi)]
            fail_df = tab[(tab.score < args.threshold) | (tab.size > args.mi)]
            pass_df.to_csv(temp_pass, sep="\t", index=False)
            temp_pass.flush()
            fail_df.to_csv(temp_fail, sep="\t", index=False)
            temp_fail.flush()
            for gff, temp, pri in zip([pass_gff3, fail_gff3], [temp_pass, temp_fail], args.priorities):
                with open(gff, "wt") as out:
                    cmd = "junctools convert -if portcullis -of igff {temp.name}".format(**locals())
                    cmd = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
                    for line in cmd.stdout:
                        line = line.decode().rstrip().rstrip(";")
                        line += ";priority={pri}".format(**locals())
                        print(line, file=out)
    return


if __name__ == "__main__":
    main()
