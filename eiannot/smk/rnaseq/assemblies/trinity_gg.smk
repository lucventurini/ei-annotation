from eicore.external_process.snakemake_helper import loadPreCmd
from eiannot.rnaseq import trinityParameters, extractSample, trinityStrandOption, trinityInput, gmap_intron_lengths
import os


rule asm_trinitygg:
	input:
		bam=os.path.join(ALIGN_DIR, "output", "{alrun}.sorted.bam"),
		align=rules.align_all.output,
		ref=rules.sanitize_reference.output
	output: os.path.join(ASM_DIR, "trinity-{run2,\d+}-{alrun}", "Trinity-GG.fasta")
	params:
		outdir=ASM_DIR+"/trinity-{run2}-{alrun}",
		load=loadPreCmd(config.get("load", dict()).get("trinity", None)),
		extra=lambda wildcards: config.get("asm_methods", dict()).get("trinity", [""]*int(wildcards.run2))[int(wildcards.run2)],
		strand=lambda wildcards: trinityStrandOption(extractSample(wildcards.alrun), SAMPLE_MAP),
		base_parameters=lambda wildcards: trinityParameters(config, extractSample(wildcards.alrun), input.ref, TGG_MAX_MEM)
	log: os.path.join(ASM_DIR, "trinity-{run2}-{alrun}.log")
	threads: THREADS
	message: "Using trinity in genome guided mode to assemble (run {wildcards.run2}): {input.bam}"
	shell: "{params.load} Trinity --seqType=fq {params.strand} --output={params.outdir} {params.extra} --genome_guided_max_intron={MAX_INTRON}  --CPU={threads}  {params.base_parameters}={input.bam} > {log} 2>&1"


rule asm_map_trinitygg:
	input:
		transcripts=rules.asm_trinitygg.output,
		index=rules.gmap_index.output
	output:
		gff=os.path.join(ASM_DIR, "output", "trinity-{run2,\d+}-{alrun}.gff")
	params:
		load=loadPreCmd(config.get("load", dict()).get("gmap", None)),
		gff=os.path.join(ASM_DIR, "trinity-{run2}-{alrun}", "trinity-{run2}-{alrun}.gff"),
		link_src=os.path.join("..", "trinity-{run2}-{alrun}", "/trinity-{run2}-{alrun}.gff"),
		intron_length=gmap_intron_lengths(config, MAX_INTRON)
	log: os.path.join(ASM_DIR, "trinitygmap-{run2}-{alrun}.log")
	threads: THREADS
	message: "Mapping trinity transcripts to the genome (run {wildcards.run2}): {input.transcripts}"
	shell: "{params.load} gmap --dir={ALIGN_DIR}/gmap/index --db={NAME} --min-intronlength={MIN_INTRON} {params.intron_length}  --format=3 --min-trimmed-coverage={TGG_COVERAGE} --min-identity={TGG_IDENTITY} -n {TGG_NPATHS} -t {THREADS} {input.transcripts} > {params.gff} 2> {log} && ln -sf {params.link_src} {output.gff} && touch -h {output.gff}"


rule trinity_all:
	input: expand(os.path.join(ASM_DIR, "output", "trinity-{run2}-{alrun}.gff"), run2=TRINITY_RUNS, alrun=ALIGN_RUNS)
	output: touch(os.path.join(ASM_DIR, "trinity.done"))
