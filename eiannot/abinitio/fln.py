#!/usr/bin/env python3

from ..abstract import EIWrapper, AtomicOperation
from .abstract import augustus_root_dir, AugustusMethod
from ..rnaseq.mikado.pick import MikadoStats
from ..rnaseq.mikado.downstream import MikadoSequenceExtractor, MikadoSequenceIndexer
from ..rnaseq.mikado.downstream import SelfDiamondP, MikadoDiamondIndex, ConvertMikadoToBed12
from ..rnaseq.mikado import Mikado
from ..preparation import SanitizeProteinBlastDB, BlastxIndex
from ..proteins.abstract import _get_value
from ..preparation import FaidxGenome
import os
import abc
import functools


class FlnWrapper(EIWrapper):

    __final_rulename__ = "fln_done"

    def __init__(self, mikado: Mikado, fai: FaidxGenome):

        super().__init__()

        self.__is_long = mikado.stats.is_long
        self.portcullis = mikado.portcullis
        self.long_alignments = mikado.long_alignments
        self.short_alignments = mikado.short_alignments
        self.faidx = fai

        if mikado.stats:
            # Start extracting the data ...
            self.configuration = mikado.configuration
            if self.dbs:
                self.sanitizer = SanitizeProteinBlastDB(self.configuration, dbs=self.dbs)
                self.blast_index = BlastxIndex(self.sanitizer)
                self.add_edge(self.sanitizer, self.blast_index)
            else:
                self.blast_index = None
            sequence_extractor = MikadoSequenceExtractor(mikado.stats)
            self.add_edge(mikado, sequence_extractor)
            diamond_index = MikadoDiamondIndex(sequence_extractor)
            self.add_edge(sequence_extractor, diamond_index)
            self.self_blast = SelfDiamondP(diamond_index)
            self.add_edge(diamond_index, self.self_blast)

            indexer = MikadoSequenceIndexer(sequence_extractor)
            self.add_edge(sequence_extractor, indexer)

            splitter = SplitMikadoFasta(indexer)
            self.add_edge(indexer, splitter)

            fln_chunks = [FLN(splitter, chunk, self.blast_index) for chunk in range(1, splitter.chunks + 1)]

            self.add_edges_from([(splitter, chunk) for chunk in fln_chunks])
            if self.blast_index:
                self.add_edges_from([(self.blast_index, chunk) for chunk in fln_chunks])

            concat = ConcatenateFLN(fln_chunks)
            self.add_edges_from([(chunk, concat) for chunk in fln_chunks])
            convert_mikado = ConvertMikadoToBed12(mikado.stats)
            self.add_edge(mikado, convert_mikado)
            self.fln_filter = FilterFLN(concat, convert_mikado, self.self_blast)
            self.add_edge(concat, self.fln_filter)
            self.add_edge(convert_mikado, self.fln_filter)
            self.add_edge(self.self_blast, self.fln_filter)
            self.training = FlnCategoryStats(self.fln_filter, "Training", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.training)
            self.testing = FlnCategoryStats(self.fln_filter, "Testing", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.training)
            self.gold = FlnCategoryStats(self.fln_filter, "Gold", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.gold)
            self.silver = FlnCategoryStats(self.fln_filter, "Silver", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.silver)
            self.bronze = FlnCategoryStats(self.fln_filter, "Bronze", is_long=mikado.stats.is_long)
            self.add_edge(self.fln_filter, self.bronze)

            if mikado.stats.is_long:
                self.__final_rulename__ += "_long"

            self.add_final_flag()
        else:
            pass

    @property
    def flank(self):
        return self.configuration.get("mikado", {}).get("pick", {}).get("flank", 200)

    @property
    def training_candidates(self):
        return self.training.input["Training"]

    @property
    def testing_candidates(self):
        return self.testing.input["Testing"]

    @property
    def is_long(self):
        return self.__is_long

    @property
    def flag_name(self):
        return os.path.join(self.fln_filter.outdir, "fln.done")

    @property
    def dbs(self):
        return self.configuration["mikado_homology"]["prot_dbs"]


class FLNOp(AugustusMethod, metaclass=abc.ABCMeta):

    def __init__(self, ancestor, is_long=None):

        super().__init__(ancestor.configuration)
        if is_long is None:
            self.__is_long = ancestor.is_long
        else:
            self.__is_long = is_long
        self.input.update(ancestor.output)

    @property
    @abc.abstractmethod
    def _rulename(self):
        pass

    @property
    def is_long(self):
        return self.__is_long

    @property
    def rulename(self):
        return "{rulename}{long}".format(rulename=self._rulename,
                                         long="_long" if self.is_long else "")

    @property
    def toolname(self):
        return "full_lengther_next"

    @property
    def outdir(self):
        fln = os.path.join("1-FLN{long}".format(long="-long-reads" if self.is_long else ""))
        return os.path.join(self._root_dir, augustus_root_dir, fln, self.__subfolder)


class SplitMikadoFasta(FLNOp):

    def __init__(self, indexer: MikadoSequenceIndexer):
        super().__init__(indexer)
        self.input.update(indexer.input)
        self.log = os.path.join(os.path.dirname(self.outdir), "split.log")
        self.output["split_flag"] = os.path.join(self.outdir, "split.done")
        self.output["chunks"] = ["{outprefix}_{chunk}.fasta".format(
            outprefix=self.outprefix, chunk=str(_).zfill(3)) for _ in range(1, self.chunks + 1)]

    @property
    def chunks(self):
        return self.configuration.get("full_lengther_next", dict()).get("chunks", 1)

    @property
    def _rulename(self):
        return "split_mikado_loci_fasta"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["mikado"]

    @property
    def cmd(self):

        load = self.load
        log = self.log
        input, output = self.input, self.output
        chunks = self.chunks
        outprefix = self.outprefix

        cmd = "{load} split_fasta.py -m {chunks} {input[transcripts]} {outprefix} && touch {output[split_flag]}".format(
            **locals()
        )

        cmd = cmd.format(**locals())
        return cmd

    @property
    def __subfolder(self):
        return os.path.join("Chunks", "Fasta")

    @property
    def outprefix(self):
        return os.path.join(self.outdir, "chunk")

    @property
    def is_small(self):
        return True


@functools.lru_cache(typed=True, maxsize=4)
def get_fln_version(loader):

    return True


class FLN(FLNOp):

    def __init__(self,
                 mikado_split: SplitMikadoFasta,
                 chunk,
                 blastx_index: BlastxIndex=None):

        super().__init__(mikado_split)
        self.__chunk = None
        self.__split = mikado_split
        self.chunk = chunk
        # self.input["fasta_chunk"] = self.fasta_chunk
        self.output = {"dbannotated": os.path.join(self.outdir, "fln_results", "pt_seqs"),
                       "flag": os.path.join(self.outdir, "fln.done")}
        self.touch = False
        self.log = os.path.join(self.outdir, "fln.log")
        if blastx_index is not None:
            self.input["fln_user"] = blastx_index.output["db"]

    @property
    def chunk(self):
        return self.__chunk

    @chunk.setter
    def chunk(self, chunk):

        if not isinstance(chunk, int):
            raise ValueError(chunk)
        if ("{outprefix}_{chunk}.fasta".format(outprefix=self.__split.outprefix, chunk=str(chunk).zfill(3))
                not in self.input["chunks"]):
            raise ValueError("{} not found among inputs!")
        self.__chunk = chunk

    @property
    def user_db(self):
        if self.input.get("fln_user", None):
            prefix = os.path.splitext(self.input.get("fln_user", None))[0]
            loc = os.path.relpath(prefix, start=self.outdir)
            return " -u {loc} ".format(**locals())
        else:
            return " "

    @property
    def fasta_chunk(self):
        return os.path.join(self.__split.outdir, "{outprefix}_{chunk}.fasta".format(
            outprefix=os.path.basename(self.__split.outprefix), chunk=str(self.chunk).zfill(3)))

    @property
    def loader(self):
        return ["full_lengther_next"]

    @property
    def _rulename(self):
        return "fln_chunk_{}".format(self.chunk)

    @property
    def cmd(self):

        taxon = self.taxon
        dbs = self.fln_dbs
        outdir = self.outdir
        log = os.path.basename(self.log)
        fasta_link_src = os.path.relpath(os.path.abspath(self.fasta_chunk),
                                         start=os.path.abspath(self.outdir))
        load = self.load
        threads = self.threads
        output = self.output
        db_verify = os.path.relpath(self.output["dbannotated"],
                                    start=self.outdir)
        flag = os.path.basename(self.output["flag"])

        dbs = "-a {dbs}".format(**locals()) if get_fln_version(self.load) else ""
        user_db = self.user_db

        cmd = "{load} mkdir -p {outdir} && cd {outdir} &&"
        cmd += " if [[ ! -e chunk.fasta ]]; then ln -s {fasta_link_src} chunk.fasta; fi && "
        cmd += "full_lengther_next {user_db} -w {threads} --taxon_group={taxon} {dbs} -f chunk.fasta 2>&1 > {log}"
        cmd += " && [[ -s {db_verify} ]] && touch {flag}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def __subfolder(self):
        return os.path.join("Chunks", "FLN", "Chunk_{}".format(str(self.chunk).zfill(5)))

    @property
    def taxon(self):
        assert self.configuration["programs"]["full_lengther_next"]["taxon"] in ["fungi",
                                                                 "human", "invertebrates",
                                                                 "mammals", "plants", "rodents",
                                                                 "vertebrates"]

        return self.configuration["programs"]["full_lengther_next"]["taxon"]

    @property
    def fln_dbs(self):
        if "dbs" not in self.configuration["programs"]["full_lengther_next"]:
            return "stnp"
        else:
            assert all(item in ("s", "t", "n", "p", "c") for item
                       in self.configuration["programs"]["full_lengther_next"]["dbs"])
            return "".join(self.configuration["programs"]["full_lengther_next"]["dbs"])


class ConcatenateFLN(FLNOp):

    def __init__(self, runs: [FLN]):

        assert runs
        super().__init__(runs[0])
        self.configuration = runs[0].configuration
        self.input = {"dbs": [run.output["dbannotated"] for run in runs]}
        self.output["table"] = os.path.join(self.outdir, "fln_table.txt")

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "concatenate_flns"

    @property
    def loader(self):
        return []

    @property
    def cmd(self):

        input, output = self.input, self.output
        input = " ".join(self.input["dbs"])
        outdir = self.outdir
        cmd = """mkdir -p {outdir} && cat {input} | """
        cmd += """awk 'NR==1 || $1!="Query_id"' > {output[table]}"""
        cmd = cmd.format(**locals())
        return cmd

    @property
    def __subfolder(self):
        return "output"

    @property
    def is_small(self):
        return True


class FilterFLN(FLNOp):

    def __init__(self, concat: ConcatenateFLN,
                 to_bed12: ConvertMikadoToBed12,
                 self_blast: SelfDiamondP):
        super().__init__(concat)
        self.input.update(to_bed12.output)
        self.input.update(self_blast.output)
        # self.input.update(to_bed12.input)
        self.output = {"table": self.outprefix + ".table.txt",
                       # "list": self.outprefix + ".list.txt",
                       }
        for category in ['Training', 'Testing', 'Gold', 'Silver', 'Bronze']:
            self.output[category] = self.outprefix + ".{category}.gff3".format(**locals())

        self.log = os.path.join(self.outdir, "filter_fln.log")

    @property
    def __subfolder(self):
        return "output"

    @property
    def outprefix(self):
        # TODO: change
        return os.path.join(self.outdir, "output")

    @property
    def threads(self):
        return 1

    @property
    def _rulename(self):
        return "filter_fln"

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]

    @property
    def flank(self):
        return self.configuration.get("mikado", {}).get("pick", {}).get("flank", 200)

    @property
    def identity(self):
        return _get_value(self.configuration, None, "identity")

    @property
    def coverage(self):
        return _get_value(self.configuration, None, "coverage")

    @property
    def max_training(self):
        return self.configuration.get("abinitio", {}).get("max_training_models", 2000)

    @property
    def cmd(self):

        load = self.load
        outdir = self.outdir
        input, log = self.input, self.log
        outprefix = self.outprefix
        mikado_loci = os.path.splitext(self.input["bed12"])[0]
        max_intron = self.max_intron
        flank = self.flank
        coverage, identity = self.coverage, self.identity
        max_training = self.max_training

        # parser = argparse.ArgumentParser(__doc__)
        # parser.add_argument("--flank", type=int, default=1000)
        # parser.add_argument("--max-intron", dest="max_intron", type=int, default=10000)
        # parser.add_argument("-cov", "--coverage", type=perc, default=80)
        # parser.add_argument("-id", "--identity", type=perc, default=80)
        # parser.add_argument("--max-training", default=2000, type=int)
        # parser.add_argument("fln")
        # parser.add_argument("mikado")
        # parser.add_argument("blast")
        # parser.add_argument("out_prefix")
        # args = parser.parse_args()

        cmd = "{load} "
        cmd += "mkdir -p {outdir} && filter_fln.py "
        cmd += " --flank {flank} "
        cmd += " --max-intron {max_intron} -cov {coverage} -id {identity} --max-training {max_training} "
        cmd += " {input[table]} {mikado_loci} {input[blast_txt]} {outprefix} > {log} 2> {log}"

        cmd = cmd.format(**locals())

        return cmd

    @property
    def is_small(self):
        return False


class FlnCategoryStats(FLNOp):

    def __init__(self, fln: FilterFLN, category, is_long=False):

        super().__init__(fln)
        self.input = fln.output
        # self.input.update(index.input)
        self.category = category
        self.configuration = fln.configuration
        self.output = {"stats": os.path.splitext(self.input[self.category])[0] + ".stats"}
        self.message = "Calculating statistics for FLN filtering, {self.category} category".format(**locals())
        self.log = self.output["stats"] + ".log"

    @property
    def rulename(self):
        return "stats_fln_{category}{long}".format(category=self.category, long="" if not self.is_long else "_long")

    @property
    def loader(self):
        return ["mikado"]

    @property
    def loci_dir(self):
        return os.path.dirname(self.input["link"])

    @property
    def __subfolder(self):
        return "output"

    @property
    def cmd(self):

        load = self.load
        input, output, log = self.input, self.output, self.log
        inp_file = self.input[self.category]
        cmd = "{load}  if [ -s {inp_file} ]; then "
        cmd += " mikado util stats {inp_file} {output[stats]} > {log} 2>&1; else "
        cmd += " touch {output[stats]}; fi"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def threads(self):
        return 1

    @property
    def is_small(self):
        return True
