#!/usr/bin/env python3

import argparse
import os
import pyfaidx
from eiannot.library.chromclass import Chrom, Chunk, AugBase
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy.orm.session import Session
from collections import defaultdict
import numpy as np


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
    parser.add_argument("-minOP", required=False, default=0.1, type=percentage,
                        help="Maximum overlap in percentage of the chunk. Default: 10%.")
    parser.add_argument("-maxOP", required=False, default=0.7, type=percentage,
                        help="Maximum overlap in percentage of the chunk. Default: 70%.")
    parser.add_argument("-minO", "--min-overlap",
                        required=True, help="Minimum overlap length.",
                        type=positive, dest="min_overlap")
    # parser.add_argument("-o", "--overlap", required=True, type=percentage)
    parser.add_argument("fasta", type=pyfaidx.Fasta, help="Input FASTA file.")
    parser.add_argument("db")
    args = parser.parse_args()

    args.overlap = round(args.min_overlap * 1.0 / args.minsize, 2)
    if args.overlap >= args.maxOP:
        # raise ValueError("The overlap cannot be larger than the minimum size. Changing.")
        # This is a mistake. Correcting.
        args.minsize = args.min_overlap / args.maxOP

    # Create the tables
    if os.path.exists(args.db):
        os.remove(args.db)
    engine = create_engine("sqlite:///{args.db}".format(**locals()))
    assert isinstance(engine, Engine)
    AugBase.metadata.create_all(engine)
    session = Session(bind=engine, autocommit=True, autoflush=True)
    session.begin()
    # We are doing it as an explicit for cycle to save computation time (cycles are expensive)
    total_length = 0
    chroms = []
    for chrom in args.fasta.keys():
        l = len(args.fasta[chrom])
        chroms.append(Chrom(chrom, l))
        total_length += l
    session.bulk_save_objects(chroms)

    # Now we have to determine the chunks. Importantly, there can be *less* chunks than the number specified.
    # L = w * ( n - n*o + o )
    # L: total length
    # n: number of chunks
    # o: overlap percentage, 0 <= o < 1
    # w: chunk width
    # => w = L / (n - n*o + o)

    total = 0
    min_size = args.minsize * 2 - args.min_overlap

    chroms = [chrom for chrom in session.query(Chrom).order_by(Chrom.length.desc())]

    chunk_size = args.minsize
    rm = args.num_chunks

    correct = set()

    overdoing = 0.05
    __orig_chunks = args.num_chunks

    args.num_chunks = min(args.num_chunks,
                          (total_length * (1 - overdoing) / args.minsize - args.minOP) * 1 / (1 - args.minOP))

    while args.num_chunks < __orig_chunks * 0.9:
        args.num_chunks = min(args.num_chunks,
                             (total_length * (1 - overdoing) / args.minsize - args.minOP ) * 1 / (1 - args.minOP))

        overdoing += 0.01
        if overdoing > 1:
            raise ValueError("I cannot optimize!")

    maximum_chunk_size = total_length / (args.num_chunks + args.minOP - args.num_chunks * args.minOP)

    while chunk_size < maximum_chunk_size:
        chunk_size = max(chunk_size + 1, int(round(chunk_size * 1.01)))
        remainder = 0
        current_chunk = 0
        _exhausted_chunks = 0
        for chrom in chroms:
            if chrom.length > min_size:
                if current_chunk > 0:
                    remainder += current_chunk
                    _exhausted_chunks += 1
                    current_chunk = 0
                continue
            else:
                remainder += chrom.length
                current_chunk += chrom.length
                if current_chunk >= min_size:
                    current_chunk = 0
                    _exhausted_chunks += 1
        rm = remaining_chunks = args.num_chunks - _exhausted_chunks
        # If the overlap is the minimum, then the covered distance is maximum
        max_total = remainder + chunk_size * (rm - rm * args.minOP - args.minOP)
        if max_total < total_length:
            # print(remainder, chunk_size, args.maxOP, max_total, total_length)
            continue  # We will *never* reach the proper size
        else:
            # Now the name of the game is to find the minimum overlap which,
            # given this chunk size, will be bound and be above the total
            diff = float("inf")
            for ol in np.arange(args.minOP, args.maxOP, 0.01):
                if ol * chunk_size < args.min_overlap:
                    continue
                total = remainder + chunk_size * (rm - rm * ol - ol)
                if abs(total_length - total) / total_length <= overdoing:
                    correct.add((chunk_size, ol))

    if len(correct) == 0:
        # Nothing has been found that respects the parameters! Aborting
        raise OSError(
            """The parameters inputed for the program are unsatisfiable;
            I cannot find a suitable way to partition the genome. Please increase the maximum overlap or
            reduce the number of chunks.""")

    correct = {chunk_size * ol: (chunk_size, ol) for chunk_size, ol in correct}
    chunk_size, args.overlap = correct[min(correct.keys())]

    # print(rm, remainder, overlap/chunk_size, total, chunk_size, overlap)
    assert chunk_size >= args.minsize
    assert args.minOP <= args.overlap <= args.maxOP, args.overlap
    overlap = int(round(chunk_size * args.overlap))
    assert overlap >= args.min_overlap, (chunk_size, args.overlap, overlap, args.min_overlap)

    chunks = defaultdict(list)
    chunk_id = 0
    chunks = []
    current_chunk = []

    for chrom in chroms:
        # We might have more than one chunk here, depending on the size of the chromosome.
        # In short, we want the following:
        # - cluster up the end of the chromosomes if we have less than the overlap remaining (thus reducing the
        # number of chunks)
        # - in the case where the end would be *beyond* the chromosome end but the remaining distance is *more* than
        # the overlap, store this as the beginning chunk, then move to the next chromosome
        # - otherwise, just store the overlap and continue
        start = 0
        end = 0

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
