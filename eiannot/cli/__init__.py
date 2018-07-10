import argparse
from ..workflow import AnnotationWorklow
from ..rnaseq.alignments.__init__ import LongAlignmentsWrapper, ShortAlignmentsWrapper
from ..rnaseq.alignments.portcullis import PortcullisWrapper
from ..rnaseq.assemblies.__init__ import ShortAssemblerWrapper
from ..proteins.__init__ import ExonerateProteinWrapper
from ..repeats.__init__ import RepeatMasking
from ..rnaseq.mikado.__init__ import Mikado
import yaml
import snakemake
from itertools import chain
import shutil
import os
import sys
import time


__doc__ = """Main entry point for the annotation pipeline."""


def get_sub_commands(SCHEDULER, prefix, additional, log_folder):
    """Adapted from Daijin."""
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    elif not (os.path.isdir(log_folder) and os.access(log_folder, os.W_OK)):
        raise OSError("Invalid directory specified: {log_folder}".format(**locals()))

    res_cmd = ""
    sub_cmd = ""

    if SCHEDULER.upper() == "LSF":
        sub_cmd = "bsub"
        res_cmd = " ".join([" -R rusage[mem={{cluster.memory}}]span[ptile={{threads}}] -n {{threads}}",
                            "-q {{cluster.queue}} -oo /dev/null",
                            "-J {prefix}_{rule} -oo {log_folder}/{prefix}_{{rule}}_%j.out"]).format(
            **locals())
    elif SCHEDULER.upper() == "PBS":
        sub_cmd = "qsub"
        res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q {cluster.queue}"
    elif SCHEDULER.upper() == "SLURM":
        sub_cmd = "sbatch"
        res_cmd = " ".join([" -N 1 -n 1 -c {{threads}} -p {{params.queue}} --mem={{params.memory}}",
                            "-J {prefix}_{{rule}} -o {log_folder}/{prefix}_{{rule}}_%j.out ",
                            "-e {log_folder}/{prefix}_{{rule}}_%j.err"]).format(**locals())

    res_cmd = "{} {}".format(res_cmd, additional)
    return res_cmd, sub_cmd


def main():

    steps = {
        "mikado": [Mikado.__final_rulename__],
        "align": [
                LongAlignmentsWrapper.__final_rulename__,
                ShortAlignmentsWrapper.__final_rulename__,
                PortcullisWrapper.__final_rulename__
                  ],
        "repeats": [RepeatMasking.__final_rulename__],
        "proteins": [ExonerateProteinWrapper.__final_rulename__],
        "assemble": [ShortAssemblerWrapper.__final_rulename__],
    }

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('-s', '--samplesheet', required=False,
                        help='Samplesheet of RNA-seq reads to be used for the analysis.')
    parser.add_argument('-g', '--genome', help="Genome FASTA file. Required.", required=True)
    parser.add_argument("--steps", choices=steps, nargs="+", default=[])
    parser.add_argument("--logdir", default=os.path.join(".", "eiannot-log"))

    snakeparse = parser.add_argument_group("Snakemake options")
    snakeparse.add_argument("--cores", default=1, type=int)
    snakeparse.add_argument("-j", "--jobs", default=1, type=int)
    snakeparse.add_argument("-n", "--dryrun", action="store_true", default=False)
    parser.add_argument("--dag", action="store_true", default=False)
    parser.add_argument('-t',
                        '--ref-transcriptome',
                        dest='ref_transcriptome', required=False,
                        help='Optional reference transcriptome to guide RNASeq alignments and predictions.')
    parser.add_argument('-c', '--configuration', required=True)
    parser.add_argument('-o', '--out', default='Snakefile', help='Workflow destination file')

    snakeparse.add_argument("--scheduler", choices=["lsf", "slurm", "pbs"], default='')
    snakeparse.add_argument("--prefix", default="eiannot",
                            help="Optional prefix to prepend to job names while using DRMAA.")
    snakeparse.add_argument("--detailed-summary", "-D", action='store_true', default=False,
                            dest="detailed_summary",
                            help="Print detailed summary of all input and output files")
    snakeparse.add_argument("--no_drmaa", "-nd", action='store_true', default=False,
                          help="Use this flag if you wish to run without DRMAA, for example, \
        if running on a HPC and DRMAA is not available, or if running locally on your own machine or server.")
    snakeparse.add_argument("-ad", "--additional-drmaa", default="", type=str,
                          dest="additional_drmaa", help="Additional parameters to be added to the DRMAA command.")
    snakeparse.add_argument("--rerun-incomplete", "--ri", action='store_true', default=False,
                            dest="rerun_incomplete",
                            help="Re-run all jobs the output of which is recognized as incomplete.")
    snakeparse.add_argument("--forcerun", "-R", nargs="+", metavar="TARGET",
                            help="Force the re-execution or creation of the given rules or files. \
                            Use this option if you changed a rule and want to have all its output in your \
                            workflow updated.")
    snakeparse.add_argument("--cleanup-metadata", dest="cleanup_metadata",
                            nargs="+",
                            default=[],
                            help="List of files that are complete even if Daijin has lost \
                            track of them because eg. they were generated by external processes.")
    snakeparse.add_argument("--list", "-l", action='store_true', default=False,
                            help="List resources used in the workflow")
    snakeparse.add_argument("-p", "--printshellcmds", help="Print out the shell commands that will be executed.")
    snakeparse.add_argument("--nolock", action="store_true", default=False,
                            help="Do not lock the working directory. Use with caution!")
    args = parser.parse_args()

    args.configuration = yaml.load(open(args.configuration))

    workflow = AnnotationWorklow(args.configuration, args.genome, args.samplesheet)

    with open(args.out, 'wt') as out:
        print("threads={threads}".format(threads=workflow.threads), file=out, end="\n\n")
        print(workflow, file=out)

    if args.steps:
        execute = list(chain.from_iterable(steps[step] for step in args.steps))
    else:
        execute = []

    res_cmd, sub_cmd = get_sub_commands(args.scheduler, args.prefix, args.additional_drmaa, args.logdir)

    cluster_var = None
    if args.no_drmaa is True and sub_cmd:
        cluster_var = sub_cmd + res_cmd

    drmaa_var = None
    if args.no_drmaa is False and res_cmd:
        try:
            import drmaa
            _ = drmaa.Session()
        except (RuntimeError, ImportError, AttributeError):
            print(
                "WARNING: DRMAA not installed or not configured properly. Switching to local/cluster mode. \
Please use the \"-nd\" flag to run the Ei-Annot pipeline if you do not plan to use DRMAA.",
                file=sys.stderr)
            drmaa_var = None
            args.no_drmaa = True
        else:
            drmaa_var = res_cmd

    snakemake.snakemake(out.name,
                        dryrun=args.dryrun,
                        printdag=args.dag,
                        until=execute,
                        cores=args.cores,
                        nodes=args.jobs,
                        workdir=os.path.abspath("."),
                        # cluster_config=hpc_conf,
                        cluster=cluster_var,
                        drmaa=drmaa_var,
                        printshellcmds=args.printshellcmds,
                        snakemakepath=shutil.which("snakemake"),
                        stats="daijin_tr_" + str(time.ctime()) + ".stats",
                        force_incomplete=args.rerun_incomplete,
                        detailed_summary=args.detailed_summary,
                        list_resources=args.list,
                        latency_wait=60 if args.scheduler else 1,
                        forceall=args.dag,
                        forcerun=args.forcerun,
                        lock=(not args.nolock)
                        )


if __name__ == '__main__':
    main()
