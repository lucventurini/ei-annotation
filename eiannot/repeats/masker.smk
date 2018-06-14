import os
from eicore.external_process.snakemake_helper import loadPreCmd

def get_input():
  pass


if config["repeats"]["model"] is True:

  rule masker:
    input:
      genome=REF,
      library=os.path.join(outdir, "reference_and_rmodeler.lib")
    output:
      masked=os.path.join(outdir, "genome.fa.masked"),
      rmout=os.path.join(outdir, "genome.fa.out"),
    params:
      load=loadPreCmd(config["load"]["repeatmasker"]),
      dir=outdir
    threads: THREADS
    shell: """{params.load} cd {params.dir} && RepeatMasker -lib {input.library} -pa {threads} -frag 500000 {input.genome}"""

else:

  rule masker:
    input: REF
    output:
      masked=os.path.join(outdir, "genome.fa.masked"),
      rmout=os.path.join(outdir, "genome.fa.out"),
    params:
      load=loadPreCmd(config["load"]["repeatmasker"]),
      species=config["repeats"]["species"],
      dir=outdir
    threads: THREADS
    shell: """{params.load} cd {params.dir} && RepeatMasker -species {params.species} -pa {threads} -frag 500000 {input.genome}"""

