#!/usr/bin/env python3

import pyBigWig
import argparse
import sys
import gzip


__doc__ = """This little script will use pyBigWig to convert a BigWig file to a variableStep compressed wig."""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-z", "--gzip", default=False, action="store_true")
    parser.add_argument("--clip", default=50, type=int, help="Maximum value for wig hints.")
    parser.add_argument("bw", type=pyBigWig.open)
    parser.add_argument("out", nargs="?", default=sys.stdout)
    args = parser.parse_args()

    if args.out != sys.stdout:
        if args.gzip is True:
            args.out = gzip.open(args.out, "wt")
        else:
            args.out = open(args.out, "wt")
    else:
        pass

    for chrom in sorted(args.bw.chroms()):
        print("variableStep chrom={}".format(chrom), file=args.out)
        for interval in args.bw.intervals(chrom):
            start, end, value = interval
            for base in range(start, end):
                print(base, min(value, args.clip), file=args.out)

    args.bw.close()
    args.out.close()


if __name__ == "__main__":
    main()