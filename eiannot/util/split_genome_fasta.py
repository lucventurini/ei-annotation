#!/usr/bin/env python3

from itertools import chain, repeat
import argparse
import os
import textwrap
import pyfaidx
from math import log, ceil, floor
from itertools import zip_longest
from eiannot.library.chromclass import Chrom, Chunk, AugBase
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy.orm.session import Session
from collections import defaultdict
import numbers


__doc__ = """Script to split FASTA sequences in a fixed number of multiple files.
This imitates the splitMfasta.pl script from Augustus, but with a key difference - the number of files is
pre-determined. Moreover, """






def positive(string):

    string = int(string)
    if string < 1:
        raise ValueError("Only positive values are acceptable!")
    return string


def fixed_grouper(number, iterable, padvalue=None):
    # return [iterable[x:x + n] for x in range(0, len(MyList), n)]
    return zip_longest(*[chain(iterable, repeat(padvalue, number - 1))] * number)


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
                        required=False, default=None,
                        type=positive, dest="min_overlap")
    parser.add_argument("-o", "--overlap", required=True, type=percentage)
    parser.add_argument("fasta", type=pyfaidx.Fasta, help="Input FASTA file.")
    parser.add_argument("db")
    args = parser.parse_args()

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

    chrom_cache = dict()

    for obj in session.query(Chrom):
        chrom_cache[obj.name] = obj.chrom_id

    # Now we have to determine the chunks. Importantly, there can be *less* chunks than the number specified.
    print(args.overlap)
    print(total_length)
    print(args.num_chunks)

    # L = w * ( n - n*o + o )
    # L: total length
    # n: number of chunks
    # o: overlap percentage, 0 <= o < 1
    # w: chunk width
    # => w = L / (n - n*o + o)

    chunk_size = 0
    overlap = 0
    if args.min_overlap is None:
        args.min_overlap = int(args.minsize * args.overlap)

    while chunk_size < args.minsize:  # and overlap < args.min_overlap:
        if args.overlap >= 1:
            raise ValueError("It is impossible to create a chunk distribution with the required parameters. Aborting")

        chunk_size = int(total_length / (args.num_chunks - args.num_chunks * args.overlap + args.overlap))
        overlap = int(round(chunk_size * args.overlap, 0))
        args.overlap += 0.01

    overlap = max(overlap, args.min_overlap)

    assert chunk_size >= args.minsize
    assert overlap >= args.min_overlap

    print(chunk_size, overlap, ((chunk_size + overlap) * args.num_chunks <= total_length))

    # Rules:
    # - each genome will have a chunk of size at least equal to minsize at the beginning and end
    # - each section of the genome should be covered at least three times, unless it's at the beginning or the end

    chunks = defaultdict(list)
    chunk_id = 0
    chunks = []
    current_chunk = []
    for chrom in args.fasta.keys():
        # We might have more than one chunk here, depending on the size of the chromosome.
        # In short, we want the following:
        # - cluster up the end of the chromosomes if we have less than the overlap remaining (thus reducing the
        # number of chunks)
        # - in the case where the end would be *beyond* the chromosome end but the remaining distance is *more* than
        # the overlap, store this as the beginning chunk, then move to the next chromosome
        # - otherwise, just store the overlap and continue
        start = 0
        end = 0
        print("Starting with", chrom)
        while end < len(args.fasta[chrom]):
            end = start + chunk_size
            if end > len(args.fasta[chrom]):
                current_chunk.append((chunk_id, chrom_cache[chrom], start, len(args.fasta[chrom])))
                break
            if end + overlap > len(args.fasta[chrom]):
                end = len(args.fasta[chrom])
            current_chunk.append((chunk_id, chrom_cache[chrom], start, end))
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
