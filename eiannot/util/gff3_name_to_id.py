#!/usr/bin/env python3


import argparse
import Mikado
import sys


__doc__ = """This script has the purpose of changing the Name of a transcript back to its ID. This is necessary for
correcting the mangling that happens with """


def main():

    tid2name = dict()

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("gff", default=sys.stdin)
    parser.add_argument("out", default=sys.stdout, nargs="?", type=argparse.FileType("wt"))
    args = parser.parse_args()

    if args.gff in ("-", "stdin"):
        args.gff = Mikado.parsers.to_gff(sys.stdin, input_format="gff3")
    else:
        args.gff = Mikado.parsers.to_gff(args.gff, input_format="gff3")

    for line in args.gff:
        if line.is_transcript:
            tid2name[line.id] = line.name
            line.id = line.name
        elif line.is_exon:
            nparents = []
            for parent in line.parent:
                assert parent in tid2name
                nparents.append(tid2name[parent])
            line.parent = nparents
        print(line, file=args.out)
        continue

    return 0


main()
