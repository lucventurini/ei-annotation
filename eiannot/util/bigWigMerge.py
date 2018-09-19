#!/usr/bin/env python3


import pyBigWig
import argparse
import sys
import numpy as np
import multiprocessing as mp
import pyfaidx


def populate_array(array, chrom, bws, normalizer, clip):

    for index, bw in enumerate(bws):
        bw = pyBigWig.open(bw)
        for interval in bw.intervals(chrom):
            start, end, value = interval
            if interval == 0:
                continue
            array[index, start:end] = value

    normalized = normalizer(array, axis=0)
    normalized = np.clip(normalized, 0, clip)
    return normalized


def main():

    normalizers = {"mean": np.mean,
                   "median": np.median,
                   "max": np.max,
                   "sum": np.sum}

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--clip", type=float, default=50)
    parser.add_argument("--threshold", type=float, default=0,
                        help="Minimum value to report (non inclusive). Default: greater than 0")
    parser.add_argument("-n", "--normalization", choices=["max", "median", "mean", "sum"],
                        default="median",
                        help="Normalization to use. Default: %(default)s")
    parser.add_argument("--inList",
                        default=False,
                        action="store_true")
    parser.add_argument("-g", "--genome", required=True, type=pyfaidx.Fasta)
    parser.add_argument("-o", "--out", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("-p", "--processes", type=int, default=1)
    parser.add_argument("bw", nargs="+", help="Input BigWig files")
    args = parser.parse_args()

    if args.inList is True:
        assert len(args.bw) == 1
        args.bw = [line.rstrip() for line in open(args.bw[0])]

    zeros = dict(
        (chrom,
         np.zeros([len(args.bw), len(args.genome[chrom])])
         ) for chrom in args.genome.keys())

    pool = mp.Pool(processes=args.processes)

    for chrom in args.genome.keys():
        zeros[chrom] = pool.apply_async(populate_array,
            args=(np.zeros([len(args.bw), len(args.genome[chrom])]),
                  chrom,
                  args.bw,
                  normalizers[args.normalization],
                  args.clip))

    pool.close()
    pool.join()
    for chrom in args.genome.keys():
        ar = zeros[chrom].get()
        if not np.where(ar > 0):
            # Skip empty chromosomes
            continue
        print("variableStep chrom={chrom}".format(**locals()), file=args.out)
        for index in np.where(ar > args.threshold)[0]:
            print(index, ar[index], file=args.out)
            continue
        continue

    return


if __name__ == "__main__":
    main()



