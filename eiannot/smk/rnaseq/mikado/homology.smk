from eicore.external_process.snakemake_helper import loadPreCmd
import os


rule blast_all:
	input: expand(os.path.join(BLAST_DIR, "xmls", "chunk-{chunk_id}-proteins.xml.gz"), chunk_id=CHUNK_ARRAY)
	output: touch(os.path.join(BLAST_DIR, "blastx.all.done")


rule create_blast_database:
  input: fa=BLASTX_TARGET
  output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
  message: "Creating the BLASTX database"
  params:
    fastas=" ".join(BLASTX_TARGET),
    load=loadPreCmd(config.get("load", dict()).get("mikado", None))
  shell: """{params.load} sanitize_blast_db.py --out {output} {params.fastas}"""

rule split_fa:
	input:
	  tr=rules.mikado_prepare.output.fa
	output: touch(os.path.join(BLAST_DIR, "fastas", "split.done"))
	params:
		outdir=os.path.join(BLAST_DIR, "fastas", "chunk"),
		chunks=config["blastx"]["chunks"],
		load=loadPreCmd(config.get("load", dict()).get("mikado", None))
	threads: 1
	message: "Splitting fasta: {input.tr}"
	shell: "{params.load} split_fasta.py -m {params.chunks} {input.tr} {params.outdir}"


if config.get("mikado", dict()).get("use_diamond", True) is False:
    rule diamond_index:
        input:
            fa=os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
        output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.dmnd")
        params:
            db=os.path.join(BLAST_DIR, "index", "blastdb-proteins"),
            load=loadPreCmd(config.get("load", dict()).get("diamond", None))
        log: os.path.join(BLAST_DIR, "diamond.index.log")
        message: "Making DIAMOND protein database for: {input.fa}"
        threads: THREADS
        shell: "{params.load} diamond makedb --threads THREADS --in {input.fa} --db {params.db} 2> {log} > {log}"

    rule diamond:
        input:
            db=rules.diamond_index.output,
            split=rules.split_fa.output
        output: os.path.join(BLAST_DIR, "xmls", "chunk-{chunk_id}-proteins.xml.gz")
        params:
            load=loadPreCmd(config.get("load", dict()).get("diamond", None)),
            tr=os.path.join(BLAST_DIR, "fastas", "chunk_{chunk_id}.fasta")
        threads: THREADS
        log: os.path.join(BLAST_DIR, "logs", "chunk-{chunk_id}.blastx.log")
        shell: """{params.load}
        if [ -s {params.tr} ]; then
            diamond blastx --threads {threads} --outfmt xml --compress 1 --out {output} --max-target-seqs {BLASTX_MAX_TARGET_SEQS} --evalue {BLASTX_EVALUE} --db {input.db} --salltitles --query {params.tr} --sensitive > {log} 2> {log};
          else
            touch {output};
          fi"""
else:
    rule make_blast:
        input: fa=os.path.join(BLAST_DIR, "index", "blastdb-proteins.fa")
        output: os.path.join(BLAST_DIR, "index", "blastdb-proteins.pog")
        params:
            db=os.path.join(BLAST_DIR, "index", "blastdb-proteins"),
            load=loadPre(config.get("load", dict()).get("blast", None))
        log: os.path.join(BLAST_DIR, "blast.index.log")
        message: "Making BLAST protein database for: {input.fa}"
        shell: "{params.load} makeblastdb -in {input.fa} -out {params.db} -dbtype prot -parse_seqids > {log} 2>&1"

    rule blastx:
        input:
            db=rules.make_blast.output,
            split=rules.split_fa.output
        output: os.path.join(BLAST_DIR, "xmls", "chunk-{chunk_id}-proteins.xml.gz")
        params:
            tr=os.path.join(BLAST_DIR, "fastas", "chunk_{chunk_id}.fasta"),
            db=os.path.join(BLAST_DIR, "/index/blastdb-proteins"),
            load=loadPreCmd(config.get("load", dict()).get("blast", None)),
            uncompressed=os.path.join(BLAST_DIR, "xmls", "chunk-{chunk_id}-proteins.xml"),
        log: os.path.join(BLAST_DIR, "logs", "chunk-{chunk_id}.blastx.log")
        threads: THREADS
        message: "Running BLASTX for mikado transcripts against: {params.tr}"
        shell: "{params.load}
        if [ -s {params.tr} ]; then
          blastx -num_threads {threads} -outfmt 5 -query {params.tr} -db {params.db} -evalue {BLASTX_EVALUE} -max_target_seqs {BLASTX_MAX_TARGET_SEQS} > {params.uncompressed} 2> {log};
        else
          touch {params.uncompressed};
        fi &&
        gzip {params.uncompressed}"
