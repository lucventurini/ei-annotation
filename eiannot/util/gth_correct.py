import Mikado
import argparse
import sys
from collections import Counter


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("gff3", nargs="?", default=sys.stdin)
    parser.add_argument("out", nargs="?", type=argparse.FileType("wt"), default=sys.stdout)
    args = parser.parse_args()
    
    gCount = Counter()
    current_gene = None

    if args.gff3 is sys.stdin or args.gff3 == "-":
        args.gff3 = sys.stdin
    else:
        args.gff3 = open(args.gff3, "rt")

    for line in args.gff3:
        line = Mikado.parsers.GFF.GffLine(line)
        if line.header:
            print(line, file=args.out)
        elif line.feature == "gene":
            current_gene = line.attributes["Target"].split()[0]
            line.feature = "match"
            suffix = ".{}".format(gCount[current_gene]) if gCount[current_gene] else ''
            gCount[current_gene] += 1
            current_gene = "{current_gene}{suffix}".format(**locals())
            line.id = current_gene
            print(line, file=args.out)
        else:
            line.parent = current_gene
            if line.feature == "exon":
                line.feature = "match_part"
            print(line, file=args.out)


main()
