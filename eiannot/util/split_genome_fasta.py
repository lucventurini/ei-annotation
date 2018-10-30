#!/usr/bin/env python3

import argparse
import os
import pyfaidx
from eiannot.library.chromclass import Chrom, Chunk, AugBase
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy.orm.session import Session
from collections import defaultdict
import numpy as np
from math import ceil


__doc__ = """Script to split FASTA sequences in a fixed number of multiple files.
This imitates the splitMfasta.pl script from Augustus, but with a key difference - the number of files is
pre-determined. Moreover, the overlap percentage can be adapted until a minimum overlap size is determined."""


def positive(string):

    string = int(string)
    if string < 1:
        raise ValueError("Only positive values are acceptable!")
    return string


def percentage(string):

    err = "Percentages must be values between 0 (inclusive) and 1 (exclusive)"

    try:
        string = float(string)
        if not (0 <= string < 1):
            raise ValueError(err)
        return string
    except ValueError:
        raise ValueError(err)


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-n", "--num-chunks", required=True, type=positive, dest="num_chunks")
    parser.add_argument("-ms", "--minsize", required=True, type=positive)
    parser.add_argument("-minO", "--min-overlap",
                        required=True, help="Minimum overlap length.",
                        type=positive, dest="min_overlap")
    parser.add_argument("fasta", type=pyfaidx.Fasta, help="Input FASTA file.")
    parser.add_argument("db")
    args = parser.parse_args()

    # We cannot have a minimum size that is lower than the overlap
    args.minsize = max(args.minsize, args.min_overlap / 0.99)

    # Create the tables
    if os.path.exists(args.db):
        os.remove(args.db)
    engine = create_engine("sqlite:///{args.db}".format(**locals()))
    assert isinstance(engine, Engine)
    AugBase.metadata.create_all(engine)
    session = Session(bind=engine, autocommit=True, autoflush=True)
    session.begin()
    # We are doing it as an explicit for cycle to save computation time (cycles are expensive)
    chroms = []

    for chrom in args.fasta.keys():
        l = len(args.fasta[chrom])
        chroms.append(Chrom(chrom, length=l))
    session.bulk_save_objects(chroms)

    chroms = [_ for _ in session.query(Chrom)]

    lengths = np.array([(_.length, _.name) for _ in chroms], dtype=[("length", np.int), ("name", np.object)])
    lengths.sort(order=["length"])
    lengths = lengths[::-1]

    # We are going to presume that the maximum overlap is 99%. So if there is no way of getting reasonably within
    # distance of the number of chunks like that, then the parameters are wrong.
    # The point here is that we are establishing a very high boundary for the overlap (99%)

    # Now we have to determine the chunks. Importantly, there can be *less* chunks than the number specified.
    # N = sum(Ni), where:
    # Ni: number of chunks for the scaffold i
    # N: total number of chunks
    # L = sum(Wi) + sum(Wi(Ni + NiO - O) - Wi), where
    # Wi: chunk size for a chromosome/scaffold of length Li. Wi can be:
    #    - W (standard chunk size), if Li > W * 2 - O
    #    - Li, otherwise
    # N: number of chunks
    # O: overlap (as percentage of W)

    # Let us verify the boundaries ...

    # The most basic control is to verify that, under the maximum overlap percentage, the
    # minimum overlap and minimum chunk size are respected.

    correct = set()

    # We are guaranteed to find some solutions to the problem. The issue is, which one is the best?

    # Increment the chunk and the overlap by fraction of the genome
    lsum = int(lengths["length"].sum())
    gp = genome_percentile = int(min(max(1000, int(lsum * 0.01)), int(args.minsize)))

    possibilities = []

    for chunk_size in range(int(args.minsize), int(lengths["length"].max()), int(gp)):
        ovl_perc = max(100, int(chunk_size * 0.01))
        for overlap_size in range(args.min_overlap, chunk_size - ovl_perc, ovl_perc):
            below = lengths[lengths['length'] < chunk_size * 2 - overlap_size]['length']
            _exhausted = int(ceil(below.sum() / chunk_size))
            num_chunks = _exhausted + int(ceil(
                (lsum - below.sum() - overlap_size) / (chunk_size - overlap_size)
            ))
            if num_chunks <= 0 or num_chunks > args.num_chunks or chunk_size * num_chunks < lsum:
                continue
            # print(num_chunks, chunk_size, overlap_size)
            assert chunk_size * num_chunks >= lsum, (chunk_size * num_chunks, lsum)
            total = below.sum() + (num_chunks - _exhausted) * chunk_size
            possibilities.append((num_chunks, chunk_size, overlap_size, total))

    # possibilities = dict((key, val) for key, val in possibilities.items() if key <= args.num_chunks)

    possibilities = np.array(possibilities, dtype=[("num_chunks", np.int),
                                                   ('chunk_size', np.int), ('overlap', np.int), ("total", np.int)
                                                   ])

    poss = possibilities.copy()
    poss['num_chunks'] = abs(args.num_chunks - poss['num_chunks'])

    # There might be more than one possibility for the best case. We have to find the best of the best, which
    # we define as the case with the least amount of overlapping regions

    best = None
    multiplier = 1.5

    while best is None:
        __new = possibilities[possibilities["total"] <= lsum * multiplier][
            poss[poss["total"] <= lsum * multiplier].argsort(
                order=['num_chunks', 'chunk_size', 'total', 'overlap'])]
        if __new.shape[0] == 0:
            print(multiplier)
            multiplier += .5
            continue
        else:
            best = __new[0]

    num_chunks, chunk_size, overlap, total = best

    # print(rm, remainder, overlap/chunk_size, total, chunk_size, overlap)
    assert chunk_size >= args.minsize
    assert overlap >= args.min_overlap, (chunk_size, overlap)

    print(num_chunks, chunk_size, overlap, total)

    chunks = defaultdict(list)
    chunk_id = 0
    chunks = []
    current_chunk = []

    for (length, name) in lengths:
        # We might have more than one chunk here, depending on the size of the chromosome.
        # In short, we want the following:
        # - cluster up the end of the chromosomes if we have less than the overlap remaining (thus reducing the
        # number of chunks)
        # - in the case where the end would be *beyond* the chromosome end but the remaining distance is *more* than
        # the overlap, store this as the beginning chunk, then move to the next chromosome
        # - otherwise, just store the overlap and continue
        start = 0
        end = 0
        chrom = session.query(Chrom).filter(Chrom.name == name).one()
        assert chrom.chrom_id is not None, chrom.name

        while end < chrom.length:
            end = start + chunk_size
            if end > chrom.length:
                current_chunk.append((chunk_id, chrom.chrom_id, start, chrom.length))
                break
            if end + chunk_size - overlap > chrom.length or (start == 0 and
                                                             chunk_size * 2 - overlap >= chrom.length):
                end = chrom.length
            current_chunk.append((chunk_id, chrom.chrom_id, start, end))
            chunk_id += 1
            chunks.extend(current_chunk)
            current_chunk = []
            start = end - overlap

    chunks.extend(current_chunk)
    session.autocommit = True
    session.bulk_save_objects([Chunk(*chunk) for chunk in chunks])
    session.flush()
    session.commit()
    session.flush()
    session.close()
    engine.dispose()

    # print("Final check")
    # assert chunk_id <= args.num_chunks, chunk_id
    return


if __name__ == "__main__":
    main()
