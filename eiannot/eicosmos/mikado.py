from cosmos.api import Cosmos
from .utils import loader, define_output, add_config
import os
import subprocess as sp
from frozendict import frozendict


@loader(programs="mikado")
def mikado_prepare(in_genome, in_conf, out_dir,  load="", procs="", log=None):

    command = """{load} mikado prepare --start-method=spawn --procs={config[threads]} --fasta={in_genome} \
    --json-conf={conf} -od {outdir} 2>&1 > {log}""".format(
        **locals())
    return command


@loader(programs=["mikado"])
def sanitize_blast_database(in_proteins, out_dir, out_put, load=""):
    fastas = " ".join(in_proteins)
    command = """{load} sanitize_blast_db.py --out {out_put} {fastas}""".format(**locals())
    return command


@loader(programs=["mikado"])
def split_mikado_fasta(in_transcripts, out_dir, out_put, chunks, load=""):
    command = "{load} split_fasta.py -m {chunks} {in_transcripts} {out_dir} && touch {out_put}".format(
        **locals())
    return command


@loader(programs=["blast"])
@define_output
def make_protein_blast_db(in_file, log, config=frozendict(), suffix=".pog", load=""):

    command = "{load} makeblastdb -in {in_file} -out {out_file} -dbtype prot -parse_seqids > {log} 2>&1".format(
        **locals())
    return command


@loader(programs=["blast"])
def blastx(in_db, in_pog, in_chunk, out_xml, config=frozendict(), log=""):

    sp.call("mkdir -p {}".format(os.path.dirname(out_xml)), shell=True)
    if os.path.exists(in_chunk):
        command = r"""{load} blastx -num_threads {config["threads"]} -outfmt 5
        -query {in_chunk} -db {in_db} -evalue {config[blastx][evalue]}
        -max_target_seqs {config[blastx][max_target_seqs]} > {out_xml}""".format(**locals())
        return command
    else:
        return """"""


@loader(programs=["diamond"])
@add_config
def diamond_index(in_db, in_dir, out_dmd):
    dbname = os.path.splitext(in_db)[0]
    log = os.path.join(os.path.dirname(in_db), "diamond_index.log")
    command = r"{load} diamond makedb --threads {config[threads]} --in {in_db} --db {dbname} 2> {log} > {log}".format(
        **locals()
    )


@loader(programs=["diamond"])
def diamond(in_db, in_chunk, chunk_id, out_xml, **kwargs):

    config=kwargs["config"]
    blast_dir = ""  # TODO: define this
    log=os.path.join(blast_dir, "logs", "chunk-{chunk_id}.blastx.log")
    command = r"""{kwargs[load]} diamond blastx --threads {config[threads]} --outfmt xml --compress 1
    --out {out_xml} --max-target-seqs {config[blastx][max_target_seqs]} --evalue {config[blastx][evalue]} --db
    {in_db} --salltitles --query {in_chunk} --sensitive > {log} 2> {log}""".format(**locals())
    return command


@loader(programs=["prodigal"])
def prodigal(in_fasta, out_gff, **kwargs):

    """"""


def recipe(*args, **kwargs):

    Cosmos.initdb(args.sqlite)




def main():
    pass


if __name__ == "__main__":
    main()