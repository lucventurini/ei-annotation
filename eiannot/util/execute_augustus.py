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
import io
import gzip
from Mikado.parsers.GTF import GtfLine


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
        for chunk in session.query(Chunk).filter(Chunk.chunk_id == args.chunk):
            chrom_file = Bio.bgzf.open(tempfile.mktemp(suffix=".fa.gz", dir=os.getcwd()), mode="wt")
            chrom, start, end = chunk.chrom, chunk.start, chunk.end
            # Using SeqIO because it is much faster than PyFaidx for this purpose
            Bio.SeqIO.write(args.genome[chrom], chrom_file, "fasta")
            chrom_file.flush()
            assert os.path.exists(chrom_file._handle.name)
            assert os.stat(chrom_file._handle.name).st_size > 0
            # We will have to use bgzip to decompress on the fly
            command = "{args.augustus} --predictionStart={start} --predictionEnd={end} "
            command += "<(bgzip -cd {chrom_file._handle.name})"
            command = command.format(**locals())
            print(command)
            aug = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
            for line in io.TextIOWrapper(aug.stdout):
                if line[0] == "#":
                    print(line, file=out)
                else:
                    # we have to change the line so that the gene/transcript names are correct
                    current_gene = None
                    fields = line.rstrip().split("\t")
                    if len(fields) != 9:  # Invalid line. Continue
                        print("#", line, file=out)
                    elif fields[2] in ("gene", "transcript"):  # This is an invalid GTF line ...
                        assert ";" not in fields[-1]
                        fields[-1] = "{chrom}_{start}_{end}-{f}".format(f=fields[-1], **locals())
                        print(*fields, file=out, sep="\t")
                    else:
                        line = GtfLine(line)
                        if "transcript_id" in line.attributes:
                            line.transcript = "{chrom}_{start}_{end}-{line.transcript}".format(**locals())
                        if "gene_id" in line.attributes:
                            line.gene = "{chrom}_{start}_{end}-{line.gene}".format(**locals())
                        print(line, file=out)

            os.remove(chrom_file._handle.name)

    return


main()
