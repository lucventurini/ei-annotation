#!/usr/bin/env python3


from Mikado.parsers import to_gff
from Mikado.parsers.GFF import GFF3
from Mikado.loci import Gene
from Mikado.utilities.log_utils import create_null_logger
from Mikado.scales.compare import prepare_reference
import argparse
import sys


__doc__ = """Script to parse and remove spurious, ie low-abundance, transcripts"""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-q", "--quantification-tag", required=True, dest="quant",
                        help="""Quantification tag used in the GTF/GFF file (e.g. "FPKM"). Required""")
    parser.add_argument("-m", "--mono", type=float, default=1,
                        help="Minimum quantification threshold for monoexonic transcripts. Default: %(default)s")
    parser.add_argument("-mu", "--multi", type=float, default=0,
                        help="Minimum quantification threshold for multiexonic transcripts. Default: %(default)s")
    parser.add_argument("gf", type=to_gff)
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    args.reference = args.gf
    args.exclude_utr, args.protein_coding = False, False
    genes, _ = prepare_reference(args,
                                 queue_logger=create_null_logger(),  # We do not need the logger in this instance
                                 ref_gff=isinstance(args.reference, GFF3))

    if args.gf.__annot_type__ == "GFF3":
        print("##gff-version\t3", file=args.out)

    for gene in sorted(genes.values()):

        assert isinstance(gene, Gene), type(gene)
        tids = gene.keys()
        to_remove = []
        for tid in tids:
            if args.quant in gene[tid].attributes:
                val = float(gene[tid].attributes[args.quant])

                if gene[tid].monoexonic:
                    if val < args.mono:
                        to_remove.append(tid)
                    else:
                        continue
                elif val < args.multi:
                    to_remove.append(tid)
                else:
                    continue
            else:
                continue

        if len(to_remove) == 0:
            print(gene.format(format_name=args.gf.__annot_type__), file=args.out)
        elif len(to_remove) == len(tids):
            continue
        else:
            [gene.remove(tid) for tid in to_remove]
            print(gene.format(format_name=args.gf.__annot_type__), file=args.out)
        continue

    return


if __name__ == "__main__":
    main()
