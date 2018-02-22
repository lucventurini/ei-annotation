from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.smk.rnaseq import gmap_intron_lengths

#
#rule asm_trinity_dn:
#  input:
#  output: ASM_DIR+"/trinity-/Trinity-DN.fasta"

# Note: this rule can be executed ONLY for Trinity versions after 2.3.2

rule trinity_dn_all:
  input: expand(os.path.join(ASM_DIR, "output", "trinity_dn_{strand}.gff"), strand=["unstranded", "firststrand", "secondstrand"])
  output: touch(os.path.join(ASM_DIR, "trinity_dn.done")


rule create_trinity_sheet:
  input: 



rule asm_trinity_dn_unstranded:
  input: rules.create_trinity_sheet.output
  output: os.path.join(ASM_DIR, "trinity_unstranded", "Trinity.fasta")
  params:
    load=loadPreCmd(config.get("load", dict()).get("trinity", None)),
    outdir=os.path.join(ASM_DIR, "output"),
    extra=lambda wildcards: config.get("asm_methods", dict()).get("trinity_dn", ""),
    memory=config.get("tgg", dict()).get("max_mem", 5000)
  threads: THREADS
  shell: """{params.load} Trinity {params.extra} --output {params.outdir} --samples_file {input} --full_cleanup --CPU {threads} --max_memory {params.memory}"""


rule asm_map_trinity_dn_unstranded:
	input:
		transcripts=rules.asm_trinity_dn_unstranded.output,
		index=rules.gmap_index.output
	output:
		gff=touch(os.path.join(ASM_DIR, "output", "trinity_dn_unstranded.gff"))
	params:
		load=loadPreCmd(config.get(dict(), "load").get("gmap", None)),
		gff=os.path.join(ASM_DIR, "trinity_unstranded", "trinity_dn_unstranded.gff"),
		link_src=os.path.join("..", "trinity_unstranded", "trinity_dn_unstranded.gff"),
		intron_length=gmap_intron_lengths(config, MAX_INTRON),
		index_dir=os.path.join(ALIGN_DIR, "gmap", "index")
	log: os.path.join(ASM_DIR, "trinitygmap-unstranded.log")
	threads: THREADS
	message: "Mapping trinity transcripts to the genome (run {wildcards.run2}): {input.transcripts}"
	shell: "{params.load} gmap --dir={params.index_dir} --db={NAME} --min-intronlength={MIN_INTRON} {params.intron_length} --format=3 --min-trimmed-coverage={TGG_COVERAGE} --min-identity={TGG_IDENTITY} -n {TGG_NPATHS} -t {THREADS} {input.transcripts} > {params.gff} 2> {log} && ln -sf {params.link_src} {output.gff} && touch -h {output.gff}"

rule asm_trinity_dn_firststrand:
  input: rules.create_trinity_sheet.output
  output: touch(os.path.join(ASM_DIR, "trinity_dn_firststrand", "Trinity.fasta"))
  params:
    load=loadPreCmd(config.get("load", dict()).get("trinity", None)),
    outdir=os.path.join(ASM_DIR, "trinity_dn_firststrand"),
    extra=lambda wildcards: config.get("asm_methods", dict()).get("trinity_dn", ""),
    strand="--SS_lib_type=RF",
    memory=config.get("tgg", dict()).get("max_mem", 5000)
  threads: THREADS
  shell: """{params.load} Trinity {params.extra} {params.strand} --output {params.outdir} --samples_file {input} --full_cleanup --CPU {threads} --max_memory {params.memory}"""


rule asm_map_trinity_dn_firststrand:
	input:
		transcripts=rules.asm_trinity_dn_firststrand.output,
		index=rules.gmap_index.output
	output:
		gff=touch(os.path.join(ASM_DIR, "output", "trinity_dn_firststrand.gff"))
	params:
		load=loadPreCmd(config.get(dict(), "load").get("gmap", None)),
		gff=os.path.join(ASM_DIR, "trinity_firststrand", "trinity_firststrand.gff"),
		link_src=os.path.join("..", "trinity_firststrand", "trinity_firststrand.gff"),
		intron_length=gmap_intron_lengths(config, MAX_INTRON),
		index_dir=os.path.join(ALIGN_DIR, "gmap", "index")
	log: os.path.join(ASM_DIR, "trinitygmap-firststrand.log")
	threads: THREADS
	message: "Mapping trinity transcripts to the genome: {input.transcripts}"
	shell: "{params.load} gmap -z sense_filter --dir={params.index_dir} --db={NAME} --min-intronlength={MIN_INTRON} {params.intron_length} --format=3 --min-trimmed-coverage={TGG_COVERAGE} --min-identity={TGG_IDENTITY} -n {TGG_NPATHS} -t {THREADS} {input.transcripts} > {params.gff} 2> {log} && ln -sf {params.link_src} {output.gff} && touch -h {output.gff}"


rule asm_trinity_dn_secondstrand:
  input: rules.create_trinity_sheet.output
  output: touch(os.path.join(ASM_DIR, "trinity_dn_secondstrand", "Trinity.fasta"))
  params:
    load=loadPreCmd(config.get("load", dict()).get("trinity", None)),
    outdir=os.path.join(ASM_DIR, "trinity_dn_secondstrand"),
    extra=lambda wildcards: config.get("asm_methods", dict()).get("trinity_dn", ""),
    strand="--SS_lib_type=FR",
    memory=config.get("tgg", dict()).get("max_mem", 5000)
  threads: THREADS
  shell: """{params.load} Trinity {params.extra} {params.strand} --output {params.outdir} --samples_file {input} --full_cleanup --CPU {threads} --max_memory {params.memory}"""


rule asm_map_trinity_dn_secondstrand:
	input:
		transcripts=rules.asm_trinity_dn_secondstrand.output,
		index=rules.gmap_index.output
	output:
		gff=touch(os.path.join(ASM_DIR, "output", "trinity_dn_secondstrand.gff"))
	params:
		load=loadPreCmd(config.get(dict(), "load").get("gmap", None)),
		gff=os.path.join(ASM_DIR, "trinity_firststrand", "trinity_secondstrand.gff"),
		link_src=os.path.join("..", "trinity_firststrand", "trinity_secondstrand.gff"),
		intron_length=gmap_intron_lengths(config, MAX_INTRON),
		index_dir=os.path.join(ALIGN_DIR, "gmap", "index")
	log: os.path.join(ASM_DIR, "trinitygmap-secondstrand.log")
	threads: THREADS
	message: "Mapping trinity transcripts to the genome: {input.transcripts}"
	shell: "{params.load} gmap -z sense_filter --dir={params.index_dir} --db={NAME} --min-intronlength={MIN_INTRON} {params.intron_length} --format=3 --min-trimmed-coverage={TGG_COVERAGE} --min-identity={TGG_IDENTITY} -n {TGG_NPATHS} -t {THREADS} {input.transcripts} > {params.gff} 2> {log} && ln -sf {params.link_src} {output.gff} && touch -h {output.gff}"