#!/usr/bin/env python


from eiannot.library.chromclass import Chrom, Chunk, AugBase
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy.orm.session import Session
import argparse
import os
import subprocess
import Bio.SeqIO
import Bio.bgzf
from functools import partial
import tempfile
import gzip


__doc__ = """"""


def positive(string):

    string = int(string)
    if string < 1:
        raise ValueError("Only positive values are acceptable!")
    return string


def main():

    fasta_index = partial(Bio.SeqIO.index, format="fasta")

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("genome", type=fasta_index)
    parser.add_argument("db")
    parser.add_argument("chunk", type=positive)
    parser.add_argument("out")
    parser.add_argument("augustus")
    args = parser.parse_args()

    assert os.path.exists(args.db)
    engine = create_engine("sqlite:///{args.db}".format(**locals()))
    session = Session(bind=engine, autocommit=True, autoflush=True)
    # Check that the database is initialised properly
    assert set.intersection(set(engine.table_names()), set(AugBase.metadata.tables.keys())) == set(engine.table_names())
    # Now query the chunk ID

    with gzip.open(os.path.splitext(args.out)[0] + "_" + str(args.chunk) + ".gtf.gz", mode="wb") as out:
        for chunk in session.query(Chunk).filter(Chunk.chunk_id == args.chunk):
            chrom_file = Bio.bgzf.open(tempfile.mktemp(suffix=".fa.gz", dir=os.getcwd()), mode="wt")
            chrom, start, end = chunk.chrom, chunk.start, chunk.end
            # Using SeqIO because it is much faster than
            Bio.SeqIO.write(args.genome[chrom], chrom_file, "fasta")
            chrom_file.flush()
            assert os.path.exists(chrom_file._handle.name)
            assert os.stat(chrom_file._handle.name).st_size > 0
            # We will have to use bgzip to decompress on the fly
            command = "{args.augustus} --predictionStart={start} --predictionEnd={end} "
            command += "<(bgzip -cd {chrom_file._handle.name})"
            command = command.format(**locals())
            print(command)
            subprocess.call(command, shell=True, stdout=out, executable="/bin/bash")
            os.remove(chrom_file._handle.name)

    return


main()
