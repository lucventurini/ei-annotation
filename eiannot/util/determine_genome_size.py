#!/usr/bin/env python3
import pyfaidx
import argparse
import math


__doc__ = """Quick script to determine the size of the genome."""


def main():

    limit = 2 ** 32 - 1
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-G", action="store_true", default=False,
                        help="Print out the size as gigabases?")
    parser.add_argument("genome")
    args = parser.parse_args()

    fa = pyfaidx.Fasta(args.genome)
    total = sum(len(fa.records[_]) for _ in fa.records)
    if args.G is False:
        print(total)
    else:
        print("{}G".format(
            math.ceil(total * 100 / 10 ** 9) / 100
        ))

    return


if __name__ == "__main__":
    main()
