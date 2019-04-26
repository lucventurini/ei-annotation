#!/usr/bin/env python3

import pandas as pd
import argparse
import subprocess as sp
import sys
import tempfile
import os
import Mikado
import io


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-p", "--priority", default=1, type=int)
    parser.add_argument("-s", "--source", default="RM")
    parser.add_argument("rmout")
    parser.add_argument("out", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()

    temp = tempfile.NamedTemporaryFile(mode="wt",
                                       dir=os.path.dirname(args.out.name) if args.out != sys.stdout else
                                       os.path.dirname(args.rmout))

    header = []
    with open(args.rmout) as rmout:
        for pos, line in enumerate(rmout):
            print(line, end='', file=temp)
            if pos == 2:
                break

    rmout = pd.read_csv(args.rmout,
                        delimiter="\s+",
                        skipinitialspace=True,
                        skiprows=3,
                        names=["score", "div", "del", "ins", "query", "qstart", "qend", "qleft", "strand",
                               "repeat", "family", "tstart", "tend", "tleft", "ID"])

    rmout[~rmout.family.isin(("Simple_repeat", "Low_complexity"))].to_csv(temp,
                                                                          header=False,
                                                                          sep=" ",
                                                                          index=False)

    temp.flush()
    pipe = sp.Popen("rmOutToGff3.pl {temp.name}".format(**locals()), shell=True, stdout=sp.PIPE)

    for line in Mikado.parsers.to_gff(pipe.stdout, input_format="gff3"):
        if line.header is True:
            print(line, file=args.out)
            continue
        else:
            assert isinstance(line, Mikado.parsers.GFF.GffLine)
            line.feature = "nonexonpart"
            line.add_attribute("source", args.source)
            line.add_attribute("priority", args.priority)
            group = line.attributes.get("target", line.attributes.get("Target"))
            # if group is not None:
            #     line.add_attribute("group", group.split()[0])
            line.remove_attribute("target")
            line.remove_attribute("Target")
            print(line, file=args.out)

    pipe.wait()
    outcode = pipe.returncode
    if outcode != 0:
        raise OSError("rmOutToGff3.pl failed. Returncode: {}".format(outcode))

    temp.close()

    return


if __name__ == "__main__":
    main()
