import os
from . import MikadoOp, __modes__
from .serialise import MikadoSerialise


class MikadoPick(MikadoOp):

    def __init__(self, serialise: MikadoSerialise):

        super().__init__(is_long=serialise.is_long)
        self.configuration = serialise.configuration
        self.input = serialise.output
        self.input["gtf"] = serialise.input["gtf"]
        self.input["cfg"] = serialise.input["cfg"]
        self.__serialise_dir = serialise.outdir
        mode = self.mode
        self.output = {"loci": os.path.join(
            self.loci_dir, "mikado-{mode}.loci.gff3").format(**locals()),
                       "link": os.path.join(self.outdir, "mikado.loci.gff3")}
        self.log = os.path.join(self.loci_dir, "mikado-{mode}.pick.err").format(
            **locals()
        )
        self.message = "Running mikado picking stage in mode: {mode}".format(**locals())

    @property
    def mode(self):
        return self.configuration["mikado"]["pick"]["mode"]

    @mode.setter
    def mode(self, mode):
        assert mode in __modes__
        self.__mode = mode

    @property
    def loader(self):
        return ["mikado"]

    @property
    def _rulename(self):
        return "mikado_pick_{mode}".format(mode=self.mode)

    @property
    def cmd(self):
        load = self.load

        cmd = "{load} "
        mode = self.mode
        threads = self.threads
        cmd += "mkdir -p {loci_dir} && mikado pick --source Mikado_{mode} --mode={mode} --procs={threads} "
        input, output = self.input, self.output
        cmd += "--start-method=spawn --json-conf={input[cfg]} "
        loci_out = os.path.basename(self.output["loci"])
        outdir = self.outdir
        loci_dir = self.loci_dir
        log = self.log
        link_dir = os.path.dirname(self.output["link"])
        link_src = os.path.relpath(self.output["loci"],
                                   start=os.path.dirname(self.output["link"]))
        link_dest = os.path.basename(self.output["link"])
        cmd += "-od {loci_dir} --loci_out {loci_out}  -lv INFO -db {input[db]} {input[gtf]} > {log} 2>&1 "
        cmd += " && mkdir -p {link_dir} && cd {link_dir} && ln -s {link_src} {link_dest}"
        cmd = cmd.format(**locals())

        return cmd

    @property
    def loci_dir(self):
        return os.path.join(self.mikado_dir, "pick")

    @property
    def outdir(self):
        return os.path.join(self.mikado_dir, "output")


class IndexMikado(MikadoOp):

    def __init__(self, mikado: MikadoPick):
        super().__init__(is_long=mikado.is_long)
        self.input = mikado.output
        self.output = {"midx": self.input["link"] + ".midx"}
        self.log = os.path.join(os.path.dirname(self.input["link"]), "index_loci.log")
        self.configuration = mikado.configuration

    @property
    def _rulename(self):
        return "index_mikado_loci"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):
        load = self.load
        input, log = self.input, self.log

        cmd = "{load} mikado compare -r {input[link]} -l {log} --index"

        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):

        return os.path.dirname(self.input["loci"])

    @property
    def threads(self):
        return 1


class MikadoStats(MikadoOp):

    def __init__(self, index: IndexMikado):

        super().__init__(is_long=index.is_long)
        self.input = index.output
        self.input.update(index.input)
        self.outdir = index.outdir
        self.configuration = index.configuration
        self.output = {"stats": os.path.splitext(self.input["link"])[0] + ".stats"}
        self.message = "Calculating statistics for Mikado run"
        self.log = self.output["stats"] + ".log"

    @property
    def _rulename(self):
        return "mikado_stats"

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["link"])

    @property
    def cmd(self):

        load = self.load
        input, output, log = self.input, self.output, self.log
        cmd = "{load} mikado util stats {input[link]} {output[stats]} > {log} 2>&1"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1
