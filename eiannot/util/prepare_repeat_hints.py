#!/usr/bin/env python3

import pandas as pd
import argparse
import subprocess as sp
import sys
import tempfile
import os


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("rmout")
    parser.add_argument("out", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()

    temp = tempfile.NamedTemporaryFile(mode="wt",
                                       dir=os.path.dirname(args.out.name) if args.out != sys.stdout else
                                       os.path.dirname(args.rmout))

    header = []
    with open(args.rmout) as rmout:
        for pos, line in enumerate(rmout):
            print(line, end='', file=temp)
            if pos == 2:
                break

    rmout = pd.read_csv(args.rmout,
                        delimiter="\s+",
                        skipinitialspace=True,
                        skiprows=3,
                        names=["score", "div", "del", "ins", "query", "qstart", "qend", "qleft", "strand",
                               "repeat", "family", "tstart", "tend", "tleft", "ID"])

    rmout[~rmout.family.isin(("Simple_repeat", "Low_complexity"))].to_csv(temp,
                                                                          header=False,
                                                                          sep=" ",
                                                                          index=False)

    temp.flush()
    if args.out != sys.stdout:
        end = " > {args.out.name}".format(**locals())
    else:
        end = ""

    outcode = sp.call("rmOutToGff3.pl {temp.name} {end}".format(**locals()), shell=True)
    if outcode != 0:
        raise OSError("rmOutToGff3.pl failed")

    temp.close()

    return


if __name__ == "__main__":
    main()
