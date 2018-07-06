import argparse
from ..workflow import AnnotationWorklow
from ..rnaseq.alignments.workflow import LongAlignmentsWrapper, ShortAlignmentsWrapper
from ..rnaseq.assemblies.workflow import ShortAssemblerWrapper
from ..rnaseq.mikado.workflow import Mikado
import yaml
import snakemake


__doc__ = """Main entry point for the annotation pipeline."""


def main():

    steps = [
            Mikado.__final_rulename__,
            LongAlignmentsWrapper.__final_rulename__,
            ShortAssemblerWrapper.__final_rulename__,
            ShortAlignmentsWrapper.__final_rulename__
             ]

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('-s', '--samplesheet', required=False,
                        help='Samplesheet of RNA-seq reads to be used for the analysis.')
    parser.add_argument('-g', '--genome', help="Genome FASTA file. Required.", required=True)
    parser.add_argument("-n", "--dryrun", action="store_true", default=False)
    parser.add_argument("--dag", action="store_true", default=False)
    parser.add_argument("--steps", choices=steps, nargs="+")
    parser.add_argument('-t',
                        '--ref-transcriptome',
                        dest='ref_transcriptome', required=False,
                        help='Optional reference transcriptome to guide RNASeq alignments and predictions.')
    parser.add_argument('-c', '--configuration', required=True)
    parser.add_argument('-o', '--out', default='Snakefile', help='Workflow destination file')
    args = parser.parse_args()

    args.configuration = yaml.load(open(args.configuration))

    workflow = AnnotationWorklow(args.configuration, args.genome, args.samplesheet)

    with open(args.out, 'wt') as out:
        print("threads={threads}".format(threads=workflow.threads), file=out, end="\n\n")
        print(workflow, file=out)
    snakemake.snakemake(out.name,
                        dryrun=args.dryrun,
                        printdag=args.dag,
                        until=args.steps)

if __name__ == '__main__':
    main()
