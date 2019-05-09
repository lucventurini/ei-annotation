#!/usr/bin/env python


from eiannot.library.chromclass import Chrom, Chunk, AugBase
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy.orm.session import Session
import argparse
import os
import subprocess
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.bgzf
import pysam
from functools import partial
import tempfile
import re
import gzip
import pybedtools
from Mikado.parsers.GTF import GtfLine


__doc__ = """"""


def positive(string):

    string = int(string)
    if string < 0:
        raise ValueError("Only positive values are acceptable!")
    return string


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("genome", type=pysam.FastaFile)
    parser.add_argument("db")
    parser.add_argument("hints")
    parser.add_argument("chunk", type=positive, help="Number of chunks to divide the genome into.")
    parser.add_argument("out")
    parser.add_argument("augustus", help="Command line to be used for this parallelisation.")
    args = parser.parse_args()

    assert os.path.exists(args.db)
    engine = create_engine("sqlite:///{args.db}".format(**locals()))
    session = Session(bind=engine, autocommit=True, autoflush=True)
    # Check that the database is initialised properly
    assert set.intersection(set(engine.table_names()), set(AugBase.metadata.tables.keys())) == set(engine.table_names())
    # Now query the chunk ID

    if args.out.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(args.out, mode="wt") as out:
        hints = pybedtools.BedTool(args.hints).as_intervalfile()
        for chunk in session.query(Chunk).filter(Chunk.chunk_id == args.chunk):
            chrom_file = Bio.bgzf.open(tempfile.mktemp(suffix=".fa.gz", dir=os.path.dirname(args.out)), mode="wt")
            hints_file = re.sub(r"\.fa.gz", ".gff", chrom_file._handle.name)
            # hints_file = tempfile.mktemp(suffix=".gff", dir=os.path.dirname(args.out))
            chrom, start, end = chunk.chrom, chunk.start, chunk.end
            interval = pybedtools.Interval(chrom, start, end)
            # Using SeqIO because it is much faster than PyFaidx for this purpose
            seq = Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(args.genome.fetch(chrom, start, end)), description="", id=chrom
            )

            Bio.SeqIO.write(seq, chrom_file, "fasta")
            with open(hints_file, "wt") as hint_handle:
                for hint in hints.all_hits(interval):
                    if hint.start < start or hint.end > end:  # Exclude partial overlaps
                        continue
                    hint.start -= start
                    hint.end -= start
                    print(hint, file=hint_handle, end='')
            chrom_file.close()
            assert os.path.exists(chrom_file._handle.name)
            assert os.stat(chrom_file._handle.name).st_size > 0
            # We will have to use bgzip to decompress on the fly
            command = "{args.augustus} --hintsfile={hints_file} {chrom_file._handle.name}"
            command = command.format(**locals())
            print(command)
            aug = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, executable="bash",
                                   universal_newlines=True)
            for line in aug.stdout:
                if line[0] == "#":
                    print(line, file=out, end="")
                    continue
                else:
                    fields = line.rstrip().split("\t")
                    fields[3] = int(fields[3]) + start
                    fields[4] = int(fields[4]) + start
                    print(*fields, file=out, sep="\t")

            aug.terminate()
            aug.communicate()
            if aug.returncode:
                raise OSError("Something went wrong during the augustus run. Please inspect the logs.")
            os.remove(chrom_file._handle.name)
            os.remove(hints_file)

    return


main()
