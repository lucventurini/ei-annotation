import os
from scipy import sparse
import argparse
import numpy as np
import re
from collections import OrderedDict
import pandas
import pyBigWig
import time
# from sklearn.preprocessing import scale


def parse_wig(sample, chroms, tpm=True):

    # "Sample" will contain three things:
    # "mapped": number of mapped reads, for normalisation
    # "prefix" of the wig files

    def __parse(wig, chrom_keys, matr):
        for chrom in chrom_keys:
            if chrom not in wig.chroms():
                # print("Chrom {} not found in {}".format(chrom, wig_name))
                continue
            matr[chrom] = np.nan_to_num(np.array(wig.values(chrom, 0, chroms[chrom])))
        return matr

    chrom_keys = list(chroms.keys())

    print(time.ctime(), "Starting with", sample["prefix"], "total of", sample["mapped"], "reads mapped")

    pos_matr = dict()
    [pos_matr.setdefault(chrom, None) for chrom in chrom_keys]
    neg_matr = None

    if sample["stranded"] is True:
        neg_matr = pos_matr.copy()
        pos_name = "{}_Forward.bw".format(sample["prefix"])
        with pyBigWig.open(pos_name, "rb") as _pos:
            pos_matr = __parse(_pos, chrom_keys, pos_matr)
        neg_name = "{}_Reverse.bw".format(sample["prefix"])
        with pyBigWig.open(neg_name, "rb") as _neg:
            neg_matr = __parse(_neg, chrom_keys, neg_matr)
    else:
        pos_name = "{}.bw".format(sample["prefix"])
        with pyBigWig.open(pos_name, "rb") as _pos:
            pos_matr = __parse(_pos, chrom_keys, pos_matr)

    if tpm is True:
        # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
        # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
        # Divide the RPK values by the “per million” scaling factor. This gives you TPM.

        nonzero = sum(len(pos_matr[_].nonzero()[0]) for _ in chroms if pos_matr[_] is not None) / 10 **3
        if neg_matr:
            nonzero += sum(len(neg_matr[_].nonzero()[0]) for _ in chroms if pos_matr[_] is not None) / 10 ** 3

        print("Nonzero for", sample["prefix"], ":", nonzero, "; mapped:", sample["mapped"])

        denominator = (nonzero / 10**3) * (sample["mapped"] / 10 ** 6)

        [pos_matr.__setitem__(_,  pos_matr[_] / denominator) for _ in chroms if pos_matr[_] is not None]
        if neg_matr:
            [neg_matr.__setitem__(_, neg_matr[_] / denominator) for _ in chroms if neg_matr[_] is not None]

    print(time.ctime(), "Finished parsing", sample["prefix"])

    return pos_matr, neg_matr


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-c", "--chrom_sizes", required=True)
    parser.add_argument("--raw", default=False, action="store_true",
                        help="""By default, this script will calculate a pseudo-TPM count for each BAM before
                        calculating the total coverage. With this flag, raw values will be used instead.""")
    parser.add_argument("-m", "--multiplier", type=int, default=1)
    parser.add_argument("-s", "--samplesheet", required=True)
    args = parser.parse_args()

    chroms = OrderedDict()

    print(time.ctime(), "Started")

    ## Parse the samplesheet to get the information
    samples = dict()

    pattern = re.compile(".*reads mapped:\s*([0-9]*).*")
    with open(args.samplesheet) as _:
        for line in _:
            name, prefix, bam, stranded = line.strip().split()
            assert os.path.exists(bam)
            assert os.path.exists("{}.stats".format(bam)), "{}.stats".format(bam)
            with open("{}.stats".format(bam)) as bstats:
                for line in bstats:
                    if "reads mapped:" in line:
                        mapped = int(re.search(pattern, line).groups()[0])
                        break
            stranded = (stranded == "True")
            samples[name] = {"prefix": prefix, "mapped": mapped, "stranded": stranded}

    with open(args.chrom_sizes) as _:
        for line in _:
            chrom, size = line.rstrip().split()[:2]
            size = int(size)
            assert chrom not in chroms
            chroms[chrom] = size

    pos_matr = dict()
    [pos_matr.setdefault(chrom, np.zeros(chroms[chrom])) for chrom in chroms]

    neg_matr = dict()
    [neg_matr.setdefault(chrom, np.zeros(chroms[chrom])) for chrom in chroms]
    for sample in samples:
        pos_nmatr, neg_nmatr = parse_wig(samples[sample], chroms, tpm=not args.raw)
        for chrom in chroms:
            if pos_nmatr[chrom] is not None:
                pos_matr[chrom] += pos_nmatr[chrom]
            if neg_nmatr is not None and neg_nmatr[chrom] is not None:
                neg_matr[chrom] += neg_nmatr[chrom]
        print(time.ctime(), "Finished adding", sample)

    print(time.ctime(), "Starting writing out the final file")

    with pyBigWig.open(re.sub("\.bw$", "", args.out) + ".+.bw", "wb") as out:
        out.addHeader(list(chroms.items()))
        max_cov = -1
        for chrom in chroms.keys():
            if len(pos_matr[chrom].nonzero()) > 0:
                ar = pos_matr[chrom] * args.multiplier
                max_cov = max(max_cov, ar.max())
                out.addEntries(chrom, 0,
                                   values=pos_matr[chrom] * args.multiplier,
                                   span=chroms[chrom], step=1)
            continue
        pass
    print("Maximum positive coverage:", max_cov)

    if any(len(neg_matr[chrom].nonzero()[0]) > 0 for chrom in chroms):
        max_cov = -1
        with pyBigWig.open(re.sub("\.bw$", "", args.out) + ".-.bw", "wb") as out:
            out.addHeader(list(chroms.items()))
            for chrom in chroms.keys():
                if len(pos_matr[chrom].nonzero()) > 0:
                    ar = neg_matr[chrom] * args.multiplier
                    max_cov = max(max_cov, ar.max())
                    out.addEntries(chrom, 0,
                                   values=ar,
                                   span=chroms[chrom], step=1)
        print("Maximum negative coverage:", max_cov)


    print(time.ctime(), "Finished")


main()