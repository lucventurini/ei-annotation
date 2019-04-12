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

    names = dict()

    for line in args.bed:
        _fields = line.rstrip().split("\t")
        if ";" in _fields[3]:
            _splitted = _fields[3].split(";")
            name, _other = _splitted[0], _splitted[1:]

        else:
            name, _other = _fields[3], ""
        names.__setitem__(name, names.get(name, 0) + 1)
        _fields[3] = "{name}.path{counter}{addi}".format(
            counter=names.get(name),
            addi=";" + _other if _other else "",
            **locals()
        )
        print(*_fields, sep="\t", file=args.out)

    if bed_to_close:
        args.bed.close()
    if out_to_close:
        args.out.close()

    return 0


if __name__ == "__main__":
    main()
