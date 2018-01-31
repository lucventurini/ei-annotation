import os
from scipy import sparse
import argparse
import numpy
import re
from collections import OrderedDict


def parse_wig(wig_name, chroms):

    matr = sparse.lil_matrix((len(chroms), max(chroms.values())))
    chrom_keys = list(chroms.keys())

    with open(wig_name) as wig:
        curr_chrom = None
        for line in wig:
            if "Step" in line:
                chrom = re.search(".*Step chrom=(\S*)(?:\s|$).*", line).groups()[0]
                assert chrom in chroms
                curr_chrom = chrom_keys.index(chrom)
            else:
                pos, val = line.rstrip().split()
                pos = int(pos)
                val = float(val)
                try:
                    matr[curr_chrom, pos] = val
                except IndexError as exc:
                    raise IndexError("{}. curr_chrom, pos: {}, {}. Shape: {}".format(exc, curr_chrom, pos, matr.shape))

    return matr


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("chrom_sizes")
    parser.add_argument("wig", nargs="+")
    args = parser.parse_args()

    chroms = OrderedDict()

    with open(args.chrom_sizes) as _:
        for line in _:
            chrom, size = line.rstrip().split()[:2]
            size = int(size)
            assert chrom not in chroms
            chroms[chrom] = size

    matr = sparse.lil_matrix((len(chroms), max(chroms.values())))
    for wig in args.wig:
        matr += parse_wig(wig, chroms)

    with open(args.out, "wt") as out:
        for cidx, chrom in enumerate(chroms.keys()):
            print("variableStep chrom={}".format(chrom), file=out)
            for pos in matr[cidx].nonzero()[1]:
                print(pos, matr[cidx, pos], sep="\t", file=out)

main()