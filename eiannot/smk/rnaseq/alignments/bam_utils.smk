rule bam_sort:
input: bam=ALIGN_DIR+"/output/{align_run}.bam"
output: ALIGN_DIR+"/output/{align_run}.sorted.bam"
params:
    load=loadPre(config, "samtools"),
    temp=ALIGN_DIR+"/sort_{align_run}"
threads: int(THREADS)
message: "Using samtools to sort {input.bam}"
shell: "{params.load} samtools sort -o {output} -O bam -m 1G -T {params.temp} -@ {threads} {input.bam}"


rule bam_index:
input: rules.bam_sort.output
output: ALIGN_DIR+"/output/{align_run}.sorted.bam.bai"
params: load=loadPre(config, "samtools")
threads: 1
message: "Using samtools to index: {input}"
shell: "{params.load} samtools index {input}"


rule bam_stats:
input:
    bam=rules.bam_sort.output,
    idx=rules.bam_index.output
output: ALIGN_DIR+"/output/{align_run}.sorted.bam.stats"
params:
    load=loadPre(config, "samtools"),
    plot_out=ALIGN_DIR+"/output/plots/{align_run}/{align_run}"
threads: 1
message: "Using samtools to collected stats for: {input}"
shell: "{params.load} samtools stats {input.bam} > {output} && plot-bamstats -p {params.plot_out} {output}"
