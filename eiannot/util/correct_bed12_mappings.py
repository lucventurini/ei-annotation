#!/usr/bin/env python3

import Mikado
import sys
import argparse
import gzip
import bz2
from collections import defaultdict
import operator


__doc__ = """Little script to correct the issue of multiple mappings with the same name in BED12 MiniMap2 files."""


def get_handles(args):

    bed_to_close = True
    if args.bed in (None, "-"):
        args.bed = sys.stdin
        bed_to_close = False
    elif args.bed.endswith("gz"):
        args.bed = gzip.open(args.bed, "rt")
    elif args.bed.endswith("bz", "bz2"):
        args.bed = bz2.open(args.bed, "rt")
    else:
        args.bed = open(args.bed, "rt")

    out_to_close = True
    if args.out in (None, "-"):
        args.out = sys.stdout
        out_to_close = False
    elif args.out.endswith("gz"):
        args.out = gzip.open(args.out, "wt")
    elif args.out.endswith("bz", "bz2"):
        args.out = bz2.open(args.out, "wt")
    else:
        args.out = open(args.bed, "wt")

    return args, bed_to_close, out_to_close


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-s", "--sort", default=False, action="store_true")
    parser.add_argument("bed", default=None, nargs="?")
    parser.add_argument("out", default=None, nargs="?")
    args = parser.parse_args()

    # Open handles as needed
    args, bed_to_close, out_to_close = get_handles(args)

    names = defaultdict(list)
    for line in args.bed:
        record = Mikado.parsers.bed12.BED12(line)
        if record.header is True or record.invalid is True:
            continue
        names[record.id].append(record)

    for name in names:
        records = sorted(names[name], key=lambda bed: (-bed.score, bed.chrom, bed.start, bed.end))
        for path, record in enumerate(records, 1):
            record.name = "{record.name}.path{path}".format(**locals())
            print(record, file=args.out)

    if bed_to_close:
        args.bed.close()
    if out_to_close:
        args.out.close()

    return 0


if __name__ == "__main__":
    main()
