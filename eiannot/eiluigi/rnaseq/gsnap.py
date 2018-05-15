import luigi
import luigi.target
import sciluigi
import os
import subprocess
from itertools import product


def create_dir(target: sciluigi.TargetInfo):

    """Quick snippet to create a path given a TargetInfo *instance*"""

    outdir = os.path.dirname(target.path)
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)


class GsnapIndex(sciluigi.SlurmTask):

    aligndir = luigi.Parameter()
    name = luigi.Parameter()
    genome = None

    @property
    def indexdir(self):
        return os.path.join(str(self.aligndir), "gmap", "index")

    @property
    def logfile(self):
        return os.path.join(str(self.aligndir), "gmap.index.log")

    def out_index(self):
        return sciluigi.TargetInfo(self,
                                   os.path.join(self.indexdir,
                                                str(self.name), str(self.name) + ".sachildguide1024")
                                   )

    def run(self):

        cmd = "{load} gmap_build --dir={dir} --db={name} {input} > {log} 2>&1".format(
            load="",
            dir=self.indexdir,
            name=self.name,
            input=self.genome,
            log=self.logfile)
        self.ex(cmd)


class GsnapTask(sciluigi.SlurmTask):

    gmap_loading = luigi.Parameter()  # This could be provided instead by a function
    in_files = None

    def out_bam(self):

        return luigi.target.AtomicLocalFile(
            os.path.join(os.path.join("eiannot"), "rnaseq", "2-align",
                         "gsnap",
                         "{}-{}".format(self.sample, self.run_num),
                         "gsnap.bam"))

    def run(self):

        # os.makedirs(os.path.dirname(self.output()))
        create_dir(self.output())

        # Create the command line
        cmd = """{load} gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {extra} \
        --novelsplicing=1 --localsplicedist={MAX_INTRON} \
        --nthreads={threads} --format=sam \
        --npaths=20 {params.infiles} 2> {log} | \
        samtools view -b -@ {threads} - \
        > {output.bam} && ln -sf {params.link_src} {output.link} && touch -h {output.link}""".format()

        # Now execute the command
        self.ex(cmd)


class Gsnaps(sciluigi.WorkflowTask):

    def output(self):

        return "toucher.txt"

    def run(self):
        with open(self.output(), "wt") as out:
            # Create dummy file
            pass

    def complete(self):
        """
        If the task has any outputs, return ``True`` if all outputs exist.
        Otherwise, return ``False``.

        However, you may freely override this method with custom logic.
        """
        required = luigi.task.flatten(self.requires())
        max_time = max(os.stat(req).st_ctime for req in required)

        outputs = luigi.task.flatten(self.output())
        import warnings
        if len(outputs) == 0:
            warnings.warn(
                "Task %r without outputs has no custom complete() method" % self,
                stacklevel=2
            )
            return False

        return all(map(lambda output: output.exists() and os.stat(output.path).st_ctime > max_time, outputs))


if __name__ == "__main__":

    luigi.build([GsnapTask(sample="A", run_num="0")])
