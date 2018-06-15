import argparse
from ..workflow import AnnotationWorklow
import yaml


__doc__ = """Main entry point for the annotation pipeline."""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('-s', '--samplesheet', required=False,
                        help='Samplesheet of RNA-seq reads to be used for the analysis.')
    parser.add_argument('-g', '--genome', help="Genome FASTA file. Required.", required=True)
    parser.add_argument('-t',
                        '--ref-transcriptome',
                        dest='ref_transcriptome', required=False,
                        help='Optional reference transcriptome to guide RNASeq alignments and predictions.')
    parser.add_argument('-c', '--configuration', required=True)
    parser.add_argument('-o', '--out', default='Snakefile', help='Workflow destination file')
    args = parser.parse_args()

    args.configuration = yaml.load(open(args.configuration))

    workflow = AnnotationWorklow(args.configuration, args.genome, args.samplesheet)
    print(workflow.graph)

    with open(args.out, 'wt') as out:
        print(workflow, file=out)
    return


if __name__ == '__main__':
    main()
