from eicore.external_process.snakemake_helper import loadPreCmd
import os

if config.get("mikado", dict()).get("use_prodigal", True) is True:
    rule prodigal:
        input: rules.mikado_prepare.output.fa
        output: os.path.join(PROD_DIR, "transcripts.fasta.prodigal.gff3")
        params:
            outdir=PROD_DIR_FULL,
            tr="transcripts.fasta",
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            tr_out="transcripts.fasta.prodigal.gff3",
            load=loadPreCmd(config.get("load", dict()).get("prodigal", None)),
            minprot=config.get("transdecoder", dict()).get("min_protein_len", 300)
        log: os.path.join(PROD_DIR_FULL, "prodigal.log")
        threads: 1
        message: "Running PRODIGAL on Mikado prepared transcripts: {input}"
        shell: "{params.load} mkdir -p {params.outdir} && cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && prodigal -f gff -g 1 -i {params.tr} -o {params.tr_out} > {log} 2>&1"
    orf_out = rules.prodigal.output
else:
  rule transdecoder_lo:
        input: rules.mikado_prepare.output.fa
        output: TDC_DIR+"/transcripts.fasta.transdecoder_dir/longest_orfs.gff3"
        params:
            outdir=TDC_DIR_FULL,
            tr="transcripts.fasta",
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            load=loadPreCmd(config.get("load", dict()).get("transdecoder", None)),
            minprot=config["transdecoder"]["min_protein_len"]
        log: os.path.join(TDC_DIR_FULL, "transdecoder.longorf.log"),
            # ss="-S" if MIKADO_STRAND else ""
        threads: 1
        message: "Running transdecoder longorf on Mikado prepared transcripts: {input}"
        run:
            if config.get("transdecoder", dict()).get("execute", True) is True:
                shell("{params.load} cd {params.outdir} && ln -sf {params.tr_in} {params.tr} && TransDecoder.LongOrfs -m {params.minprot} -t {params.tr} > {log} 2>&1")
            else:
                shell("mkdir -p $(dirname {output}) && touch {output}")

    rule transdecoder_pred:
        input:
            mikado=rules.mikado_prepare.output.fa,
            trans=rules.transdecoder_lo.output
        output: os.path.join(TDC_DIR, "transcripts.fasta.transdecoder.bed")
        params:
            outdir=TDC_DIR_FULL,
            tr_in=os.path.join(MIKADO_DIR_FULL, "mikado_prepared.fasta"),
            lolog=os.path.join(TDC_DIR_FULL, "transdecoder.longorf.log"),
            plog=os.path.join(TDC_DIR_FULL, "transdecoder.predict.log"),
            tr="transcripts.fasta",
            load=loadPreCmd(config.get("load", dict()).get("transdecoder", None))
            # ss="-S" if MIKADO_STRAND else ""
        log: os.path.join(TDC_DIR_FULL, "transdecoder.predict.log")
        threads: THREADS
        message: "Running transdecoder predict on Mikado prepared transcripts: {input}"
        run:
            if config.get("transdecoder", dict()).get("execute", True) is True:
                shell("{params.load} cd {params.outdir} && TransDecoder.Predict -t {params.tr} --cpu {threads} > {log} 2>&1")
            else:
                shell("mkdir -p $(dirname {output}) && touch {output}")
    orf_out = rules.transdecoder_pred.output
