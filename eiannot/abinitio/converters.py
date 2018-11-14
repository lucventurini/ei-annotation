from ..abstract import EIWrapper
import abc
# from ..rnaseq.mikado.__init__ import Mikado
from ..rnaseq.alignments.portcullis import PortcullisWrapper
from ..rnaseq.alignments import ShortAlignmentsWrapper
from ..rnaseq.alignments.bam import Bam2BigWig, MergeWigs
from ..repeats.__init__ import RepeatMasking
from ..proteins.__init__ import ProteinWrapper
from ..proteins.abstract import _get_value, FilterAlignments
from .fln import FlnWrapper, FilterFLN
from .abstract import AugustusMethod, augustus_root_dir
import os


__doc__ = """A single file of hints has to be generated.

Portcullis: only "intron"
Mikado:
  - "exon": internal exons
  - "exonpart": beginning/end of the model
  - "intron": the introns present
  - "CDSpart": beginning/end of the model
  - "CDS": internal CDS exons.
Portcullis junctions:
  - "intron"  
Proteins:
  - "intron"
  - "CDS": internal exons
  - "CDSpart": terminal exons
Repeats:
  - "nonexonpart" 

"""


augustus_help = """

              5. USING HINTS (AUGUSTUS+)
              --------------------------
AUGUSTUS can take hints on the gene structure. It currently accepts 16 types of hints:  
start, stop, tss, tts, ass, dss, exonpart, exon, intronpart, intron, CDSpart, CDS, UTRpart, UTR, irpart, nonexonpart.
The hints must be stored in a file in gff format containing one hint per line. 
Example of a hintsfile:

HS04636	mario	exonpart	500	506	.	-	.	source=M
HS04636	mario	exon	966	1017	.	+	0	source=P
HS04636	AGRIPPA	start	966	968	6.3e-239	+	0	group=gb|AAA35803.1;source=P
HS04636	AGRIPPA	dss	2199	2199	1.3e-216	+	.	group=gb|AAA35803.1;source=P
HS04636	mario	stop	7631	7633	.	+	0	source=M
HS08198	AGRIPPA	intron	2000	2000	0	+	.	group=ref|NP_000597.1;source=E
HS08198	AGRIPPA	ass	757	757	1.4e-52	+	.	group=ref|NP_000597.1;source=E

The fields must be separated by a tabulator. In the first column (field) the sequence name is given. 
In this case the hints are together about two sequences. The second field is the name of the program that produced the hint. 
It is ignored here. The third column specifies the type of the hint. The 4th and 5th column specify the begin and end position of the hint.
Positions start at 1. The 6th colum gives a score. The 7th the strand. The 8th the reading frame as defined in the GFF standard. 
The 9th column contains arbitrary extra information but it must contain a string 'source=X' where X
is the source identifier of the hint. Which values for X are possible is specified in the file augustus/config/extrinsic.cfg,
e.g. X=M, E, or P.

AUGUSTUS can follow a hint, i.e. predict a gene structure that is compatible with it, or 
AUGUSTUS can ignore a hint, i.e. predict a gene structure that is not compatible with it. The probability that AUGUSTUS ignores a 
hint is the smaller the more reliable the hints of this type are. 

Run AUGUSTUS using the --hintsfile option. Example:

>augustus --species=human --hintsfile=../examples/hints.gff --extrinsicCfgFile=../config/extrinsic.human.MPEC.cfg ../examples/example.fa 

As an alternative to giving the option --extrinsicCfgFile you can replace augustus/config/extrinsic.cfg with the appropriate file, as this 
file is read by default when the option --extrinsicCfgFile is not given.

Explanation of the file format of the extrinsic.cfg file.

The gff/gtf file containing the hints must contain somewhere in the last
column an entry source=?, where ? is one of the source characters listed in
the line after [SOURCES] above. You can use different sources when you have
hints of different reliability of the same type, e.g. exon hints from ESTs
and exon hints from evolutionary conservation information.

In the [GENERAL] section the entries second column specify a bonus for obeying
a hint and the entry in the third column specify a malus (penalty) for
predicting a feature that is not supported by any hint. The bonus and the
malus is a factor that is multiplied to the posterior probability of gene
structueres. 
Example: 
  CDS     1000  0.7  ....
means that, when AUGUSTUS is searching for the most likely gene structure,
every gene structure that has a CDS exactly as given in a hint gets
a bonus factor of 1000. Also, for every CDS that is not supported the
probability of the gene structure gets a malus of 0.7. Increase the bonus to
make AUGUSTUS obey more hints, decrease the malus to make AUGUSTUS predict few
features that are not supported by hints. The malus helps increasing
specificity, e.g. when the exons predicted by AUGUSTUS are suspicious because
there is no evidence from ESTs, mRNAs, protein databases, sequence
conservation, transMapped expressed sequences.
Setting the malus to 1.0 disables those penalties. Setting the bonus to 1.0
disables the boni. 

      start: translation start (start codon), specifies an interval that contains
             the start codon. The interval can be larger than 3bp, in which case
             every ATG in the interval gets a bonus. The highest bonus is given
             to ATGs in the middle of the interval, the bonus fades off towards the ends.
       stop: translation end  (stop codon), see 'start'
        tss: transcription start site, see 'start'
        tts: transcription termination site, see 'start'
        ass: acceptor (3') splice site, the last intron position, for only approximately known ass an interval can be specified
        dss: donor (5') splice site, the first intron position, for only approximately known dss an interval can be specified
   exonpart: part of an exon in the biological sense. The bonus applies only
             to exons that contain the interval from the hint. Just
             overlapping means no bonus at all. The malus applies to every
             base of an exon. Therefore the malus for an exon is exponential
             in the length of an exon: malus=exonpartmalus^length.
	     Therefore the malus should be close to 1, e.g. 0.99.
       exon: exon in the biological sense. Only exons that exactly match the
             hint get a bonus. Exception: The exons that contain the start
             codon and stop codon. This malus applies to a complete exon
             independent of its length.
 intronpart: introns both between coding and non-coding exons. The bonus
             applies to every intronic base in the interval of the hint.
     intron: An intron gets the bonus if and only if it is exactly as in the hint.
    CDSpart: part of the coding part of an exon. (CDS = coding sequence)
        CDS: coding part of an exon with exact boundaries. For internal exons
             of a multi exon gene this is identical to the biological
             boundaries of the exon. For the first and the last coding exon
             the boundaries are the boundaries of the coding sequence (start, stop).
        UTR: exact boundaries of a UTR exon or the untranslated part of a
             partially coding exon.
    UTRpart: The hint interval must be included in the UTR part of an exon.
     irpart: The bonus applies to every base of the intergenic region. If UTR
             prediction is turned on (--UTR=on) then UTR is considered
             genic. If you choose against the usual meaning the bonus of
             irparts to be much smaller than 1 in the configuration file you
             can force AUGUSTUS to not predict an intergenic region in the
             specified interval. This is useful if you want to tell AUGUSTUS
             that two distant exons belong to the same gene, when AUGUSTUS
             tends to split that gene into smaller genes.
nonexonpart: intergenic region or intron. The bonus applies to very non-exon
             base that overlaps with the interval from the hint. It is
             geometric in the length of that overlap, so choose it close to
             1.0. This is useful as a weak kind of masking, e.g. when it is
             unlikely that a retroposed gene contains a coding region but you
             do not want to completely forbid exons.
  genicpart: everything that is not intergenic region, i.e. intron or exon or UTR if
             applicable. The bonus applies to every genic base that overlaps with the
             interval from the hint. This can be used in particular to make Augustus
             predict one gene between positions a and b if a and b are experimentally
             confirmed to be part of the same gene, e.g. through ESTs from the same clone.
             alias: nonirpart

Any hints of types dss, intron, exon, CDS, UTR that (implicitly) suggest a donor splice
site allow AUGUSTUS to predict a donor splice site that has a GC instead of the much more common GT.
AUGUSTUS does not predict a GC donor splice site unless there is a hint for one.

Starting in column number 4 you can tell AUGUSTUS how to modify the bonus 
depending on the source of the hint and the score of the hint. 
The score of the hints is specified in the 6th column of the hint gff/gtf.
If the score is used at all, the score is not used directly through some
conversion formula but by distinguishing different classes of scores, e.g. low
score, medium score, high score. The format is the following:
First, you specify the source character, then the number of classes (say n), then you
specify the score boundaries that separate the classes (n-1 thresholds) and then you specify
for each score class the multiplicative modifier to the bonus (n factors). 

Examples:

M 1 1e+100
means for the manual hint there is only one score class, the bonus for this
type of hint is multiplied by 10^100. This practically forces AUGUSTUS to obey
all manual hints.

T    2       1.5 1 5e29
For the transMap hints distinguish 2 classes. Those with a score below 1.5 and
with a score above 1.5. The bonus if the lower score hints is unchanged and
the bonus of the higher score hints is multiplied by 5x10^29.

D    8     1.5  2.5  3.5  4.5  5.5  6.5  7.5  0.58  0.4  0.2  2.9  0.87  0.44 0.31  7.3
Use 8 score classes for the DIALIGN hints. DIALIGN hints give a score, a strand and
reading frame information for CDSpart hints. The strand and reading frame are often correct but not
often enough to rely on them. To account for that I generated hints for all
6 combinations of a strand and reading frame and then used 2x2x2=8 different
score classes:
{low score, high score} x {DIALIGN strand, opposite strand} x {DIALIGN reading frame, other reading frame}
This example shows that scores don't have to be monotonous. A higher score
does not have to mean a higher bonus. They are merely a way of classifying the
hints into categories as you wish. In particular, you could get the effect of
having different sources by having just hints of one source and then distinguishing
more scores classes.

Alternative Transcripts / Alternative Splicing (evidence based)
---------------------------------------------------------------

AUGUSTUS can predict alternative splicing or - more general - alternative transcripts that are suggested by evidence given in hints.
The method is very general. But to give an example: If two EST alignments to the same genomic area cannot be explained by a single
transcript then AUGUSTUS can predict a gene with two different splice forms, one splice form
compatible with each of the EST alignments.

Grouping hints:
Each hint can be given a group name, by specifying 'group=goupname;' or 'grp=goupname;' in the last column for the hint
in the gff file. This should be used to group all the hints coming from the alignment of the same
sequence to the genome. For example, if an EST with the name est_xyz aligns to the genome with one gap
suggesting an intron then the hints resulting from that alignment could look like this

HS04636	blat2hints	exonpart	500	599	.	+	.	group=est_xyz; source=E
HS04636	blat2hints	intron		600	700	.	+	.	group=est_xyz; source=E
HS04636	blat2hints	exonpart	701	900	.	+	.	group=est_xyz; source=E

Grouping tells AUGUSTUS that hints belong together. Ideally, all hints of a group are obeyed by a predicted
transcript or the whole group of hints is ignored when making the prediction.

Prioritizing groups:
Hints or hint groups can be given a priority by specifying 'priority=n;' or 'pri=n' in the last
column for the hint in the gff file. For example

HS04636	blat2hints	exonpart	500	599	.	+	.	priority=2; source=E
HS04636	blat2hints	intron		550	650	.	+	.	priority=5; source=mRNA

When two hints or hint groups contradict each other then the hints with the lower priority number
are ignored. This is especially useful if for a genome several sources of hints are available,
where one source should be trusted when in doubt. For example, the rhesus macaque currently has few native ESTs
but human ESTs often also align to rhesus. Giving the hints from native ESTs a higher priority means
that AUGUSTUS uses only them for genes with support by native ESTs and uses the alien EST alignments
when native ESTs alignments are not available for a gene. When the priority is not specified, it is
internally set to -1.

When AUGUSTUS is run with --alternatives-from-evidence=false then all hints are given to AUGUSTUS at
the same time no whether they can be explained with a single transcript per gene. AUGUSTUS will then
choose the most likely transcript variant.

When AUGUSTUS is run with --alternatives-from-evidence=true then AUGUSTUS will predict alternative
transcripts based on the alternatives the hints suggest. This can be any form of alternative
splicing, including nested genes contained in introns of other genes, overlapping genes, alternative
translation starts and variation in UTR.

"""

outdir = os.path.join("3-Hints")


class ConvertToHints(EIWrapper):

    __final_rulename__ = "augustus_conversion_done"

    def __init__(self,
                 run: int,
                 mikado: FlnWrapper,
                 mikado_long: FlnWrapper,
                 repeats: RepeatMasking,
                 proteins: ProteinWrapper):

        super().__init__(configuration=mikado.configuration)

        self.__coverages = {"+": [], "-": [], ".": []}
        self.alignments = mikado.short_alignments
        self.perform_coverages()
        self.__run = None
        self.run = run

        if mikado.fln_filter:
            self.mikado_converter = ConvertMikado(run=self.run, filterer=mikado.fln_filter)
            self.add_node(self.mikado_converter)
        else:
            self.mikado_converter = None
        if mikado_long.fln_filter:
            self.mikado_converter_long = ConvertMikado(run=self.run, filterer=mikado_long.fln_filter)
            self.add_node(self.mikado_converter_long)
        else:
            self.mikado_converter_long = None

        portcullis = ConvertJunctions(run=self.run, junctions=mikado.portcullis)
        self.add_node(portcullis)
        if repeats.execute:
            repeats = ConvertRepeats(run=self.run, repeats=repeats)
            self.add_node(repeats)
        self.proteins = []
        for filterer in proteins.proteins:
            converter = ConvertProteins(run=self.run, proteins=filterer)
            self.proteins.append(converter)

        # Add them to the graph
        [self.add_node(converter) for converter in self.proteins]
        self.add_final_flag()

    def perform_coverages(self):
        for bam in self.alignments.bams:

            if bam.sample.stranded is True:
                plus = Bam2BigWig(bam, strand="+")
                self.add_edge(bam, plus)
                self.__coverages["+"].append(plus)
                minus = Bam2BigWig(bam, strand="-")
                self.add_edge(bam, minus)
                self.__coverages["-"].append(minus)
            else:
                nonstrand = Bam2BigWig(bam, strand=None)
                self.add_edge(bam, nonstrand)
                self.__coverages["."].append(nonstrand)

        # Remove non-stranded if we have stranded
        if self.__coverages["+"] or not self.always_execute_unstranded:
            self.remove_nodes(self.__coverages["."])
            self.__coverages["."] = []

        wig2hints_store = []

        for strand in self.__coverages.keys():
            if not self.__coverages[strand]:
                continue
            merged = MergeWigs(bigwigs=self.__coverages[strand], strand=strand,
                                faidx=self.alignments.prepare_wrapper.fai)
            self.add_edges_from([(coverage, merged) for coverage in self.__coverages[strand]])
            # We will trust the main class to find out the duplications in creating the WIG files
            # and removing them.
            wig2hints = WigToHints(run=self.run, merger=merged)
            self.add_edge(merged, wig2hints)
            wig2hints_store.append(wig2hints)

    @property
    def __always_skip_unstranded(self):
        value = self.configuration.get("abinitio", {}).get("coverage", {}).get(
            "always_skip_unstranded", False
        )

        assert isinstance(value, bool)
        return value

    @property
    def always_skip_unstranded(self):

        return self.__always_skip_unstranded

    @property
    def __always_execute_unstranded(self):

        value = self.configuration.get("abinitio", {}).get("coverage", {}).get(
            "always_execute_unstranded", False
        )

        assert isinstance(value, bool)
        return value

    @property
    def always_execute_unstranded(self):
        return self.__always_execute_unstranded or (not self.always_skip_unstranded)

    @property
    def flag_name(self):
        return os.path.join(self.outdir, "augustus_hints_conversion.done")

    @property
    def outdir(self):
        return os.path.join(self._root_dir, augustus_root_dir,
                            "run-{run}".format(run=self.run),
                            outdir, "output")


class HintMethod(AugustusMethod, metaclass=abc.ABCMeta):

    def __init__(self, run, configuration):

        super().__init__(configuration)
        self.__run = None
        self.run = run

    @property
    def run(self):
        return self.__run

    @run.setter
    def run(self, run):

        assert isinstance(run, int)
        # Check the run is valid
        self.__run = run

    @property
    def __subfolder(self):
        return "3-Hints"

    @property
    def hint_dir(self):

        return os.path.join(self._augustus_dir, "run-{run}".format(run=self.run), self.__subfolder)


class WigToHints(HintMethod):

    __name__ = "wig_to_hints"

    def __init__(self, run, merger: MergeWigs):

        super().__init__(run=run, configuration=merger.configuration)
        self.input = merger.output
        self.__name_prefix = merger.name_prefix
        self.__strand = merger.strand
        self.output["hints"] = os.path.join(self.outdir, self.name_prefix + ".hints.gff")

    @property
    def name_prefix(self):
        return self.__name_prefix

    @property
    def strand(self):
        return self.__strand

    @property
    def outdir(self):
        return os.path.join(self.hint_dir, "coverage")

    @property
    def rulename(self):
        rule = self.__name__
        if self.strand is not None:
            rule += "_" + self.strand
        return rule

    @property
    def loader(self):
        return ["ei-annotation", "augustus"]

    @property
    def width(self):
        return 10

    @property
    def margin(self):
        return 10

    @property
    def threshold(self):
        return 2

    @property
    def minscore(self):
        return 4

    @property
    def prune(self):
        return 0.1

    @property
    def radius(self):
        return 4.5

    @property
    def src(self):
        # TODO change
        return "W"

    @property
    def priority(self):
        # TODO change
        return 3

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        width, margin, threshold = self.width, self.margin, self.threshold
        minscore, prune, radius = self.minscore, self.prune, self.radius
        wig = self.input["wig"]
        out = self.output["hints"]
        outdir = self.outdir
        src = self.src
        source = "coverage"
        if self.strand:
            source += "_" + self.strand
        if self.strand:
            if self.strand == "forward":
                strand = "+"
            else:
                strand = "-"
        else:
            strand = "."
        priority = self.priority

        cmd = "{load} mkdir -p {outdir} && wig2hints.pl --width={width} --margin={margin} --minthresh={threshold} "
        cmd += " --src={src} --strand={strand} --type=exonpart --UCSC={source} --pri={priority}"
        cmd += " --minscore={minscore} --prune={prune} --radius={radius} < {wig} > {out}"

        cmd = cmd.format(**locals())
        return cmd


class MergeHints(HintMethod):

    __doc__ = """This class has the purpose of bringing together all the hints into a single GFF file."""

    def __init__(self, run, hints, with_coverage: bool):

        if len(hints) == 0:
            raise ValueError("No hints provided")

        super().__init__(configuration=hints[0].configuration, run=hints[0].run)
        self.__with_coverage = None
        self.with_coverage = with_coverage
        self.input["hints"] = []
        self.output["hints"] = os.path.join(self.outdir, "hints.gff")

    @property
    def outdir(self):
        return os.path.join(self.hint_dir, "output")

    @property
    def with_coverage(self):
        return self.__with_coverage

    @with_coverage.setter
    def with_coverage(self, flag):
        if flag not in (False, True):
            raise ValueError(flag)
        self.__with_coverage = flag

    @property
    def rulename(self):
        return "merge_and_sort_hints"

    @property
    def loader(self):
        return ["genometools"]

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        out = self.output["hints"]
        gffs = " ".join(self.input["hints"])
        cmd = "{load}"
        cmd += " gt gff3 -tidy -sort -force -o {out} {gffs}"
        cmd = cmd.format(**locals())
        return cmd


class ConvertMikado(HintMethod):

    """Mikado will have to be converted in *four* batches, with *four* different priorities and two different sources:
    - gold: source "M", default priority: 10
    - silver: source "M", default priority: 9
    - bronze: source "E", default priority: 8
    - all: source "E", default priority: 7
    """

    def __init__(self, run, filterer: FilterFLN):

        super().__init__(run=run, configuration=filterer.configuration)
        self.__is_long = filterer.is_long

        self.input["table"] = filterer.output["table"]
        self.input["bed12"] = filterer.input["bed12"]
        self.output["hints"] = os.path.join(self.outdir, "mikado.hints.{long}gff3".format(
            long="long." if self.is_long else ""))

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], outdir, "output")

    @property
    def gold_score(self):
        # TODO CHANGE
        return self.configuration.get("abinitio", {}).get("mikado", {}).get("scores", {}).get("gold", 10)

    @property
    def silver_score(self):
        # TODO CHANGE
        return self.configuration.get("abinitio", {}).get("mikado", {}).get("scores", {}).get("silver", 9)

    @property
    def bronze_score(self):
        # TODO CHANGE
        return self.configuration.get("abinitio", {}).get("mikado", {}).get("scores", {}).get("bronze", 8)

    @property
    def all_score(self):
        # TODO CHANGE
        return self.configuration.get("abinitio", {}).get("mikado", {}).get("scores", {}).get("all", 7)

    @property
    def gold_source(self):
        # TODO CHANGE
        src = self.configuration.get("abinitio", {}).get("mikado", {}).get("sources", {}).get("gold", "M")
        assert src in ("M", "E", "P", "RM", "W")
        return src

    @property
    def silver_source(self):
        # TODO CHANGE
        src = self.configuration.get("abinitio", {}).get("mikado", {}).get("sources", {}).get("silver", "M")
        assert src in ("M", "E", "P", "RM", "W")
        return src

    @property
    def bronze_source(self):
        # TODO CHANGE
        src = self.configuration.get("abinitio", {}).get("mikado", {}).get("sources", {}).get("bronze", "E")
        assert src in ("M", "E", "P", "RM", "W")
        return src

    @property
    def all_source(self):
        # TODO CHANGE
        src = self.configuration.get("abinitio", {}).get("mikado", {}).get("sources", {}).get("all", "E")
        assert src in ("M", "E", "P", "RM", "W")
        return src

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]

    @property
    def cmd(self):

        load = self.load
        input, output = self.input, self.output
        mikado_loci = self.input["bed12"]
        
        cmd = "{load} filtered_fln_to_hints.py "
        gold_score, gold_source = self.gold_score, self.gold_source
        cmd += " -gs {gold_score} --gold-source {gold_source} "
        
        silver_score, silver_source = self.silver_score, self.silver_source
        cmd += " -ss {silver_score} --silver-source {silver_source} "
        
        bronze_score, bronze_source = self.bronze_score, self.bronze_source
        cmd += " -bs {bronze_score} --bronze-source {bronze_source} "

        all_score, all_source = self.all_score, self.all_source
        cmd += " -as {all_score} --all-source {all_source} "

        table = self.input["table"]
        out = self.output["hints"]
        cmd += " {mikado_loci} {table} {out}"

        cmd = cmd.format(**locals())

        return cmd

    @property
    def threads(self):
        return 1

    @property
    def is_long(self):
        return self.__is_long

    @property
    def rulename(self):
        return "prepare_mikado_hints{long}".format(long="_long" if self.is_long else "")


class ConvertRepeats(HintMethod):

    def __init__(self, run: int, repeats: RepeatMasking):

        if not repeats.execute:
            return
        if repeats.execute:
            super().__init__(run=run, configuration=repeats.configuration)
            self.input["table"] = repeats.output["table"]
            self.output["gff3"] = os.path.join(self.outdir, "repeat_hints.gff3")
            self.log = os.path.join(os.path.dirname(self.outdir), "logs", "extract_repeats.log")
            self.logdir = os.path.dirname(self.log)

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], outdir, "output")

    @property
    def loader(self):
        return ["repeatmasker", "ei-annotation"]

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        # TODO CHANGE
        load = self.load
        cmd = "{load} mkdir -p {logdir} && "
        logdir, log = self.logdir, self.log
        input, output = self.input, self.output
        log = self.log
        cmd += "prepare_repeat_hints.py {input.table} {output.gff3} > {log} 2> {log}"

        cmd = cmd.format(**locals())

        return

    @property
    def rulename(self):
        return "convert_repeat_hints"


class ConvertProteins(HintMethod):

    def __init_subclass__(cls, **kwargs):
        pass

    def __init__(self, run: int, proteins: FilterAlignments):

        super().__init__(run=run, configuration=proteins.configuration)
        self.__dbname = proteins.dbname
        self.input = proteins.output
        self.output["hints"] = os.path.join(self.outdir, "{dbname}.hints.gff".format(dbname=self.dbname))

    @property
    def rulename(self):
        return "convert_proteins_to_hints_{dbname}".format(dbname=self.dbname)

    @property
    def loader(self):
        return ["mikado", "ei-annotation"]

    @property
    def outdir(self):
        return os.path.join(self.hint_dir, "output")

    @property
    def dbname(self):
        return self.__dbname

    @property
    def priority(self):
        # TODO implement
        pass

    @property
    def source(self):
        # TODO implement
        pass

    @property
    def threads(self):
        return 1

    @property
    def cmd(self):
        load = self.load
        input, output = self.input, self.output
        priority, source = self.priority, self.source
        cmd = "{load} convert_proteins_to_hints.py -p {priority} -s {source} {input[gff3]} {output[hints]}"
        cmd = cmd.format(**locals())
        return cmd


class ConvertJunctions(HintMethod):

    def __init__(self, run: int, junctions: PortcullisWrapper):

        super().__init__(run=run, configuration=junctions.configuration)
        self.configuration = junctions.configuration
        self.input["tab"] = junctions.junctions_task.output["tab"]
        self.output["gold"] = os.path.join(self.outdir,
                                           "portcullis_junctions.gold.{source}{score}.gff3".format(
                                               source=self.gold_source, score=self.gold_score
                                           ))
        self.output["silver"] = os.path.join(self.outdir,
                                             "portcullis_junctions.silver.{source}{score}.gff3".format(
                                                 source=self.silver_source, score=self.silver_score
                                             ))

    @property
    def rulename(self):
        return "prepare_hints_portcullis"

    @property
    def threads(self):
        return 1

    @property
    def loader(self):
        return ["ei-annotation"]

    @property
    def gold_source(self):
        #TODO implement
        pass

    @property
    def gold_score(self):
        # TODO implement
        pass

    @property
    def silver_source(self):
        # TODO: most probably this will have to go into the configuration
        pass

    @property
    def silver_score(self):
        # TODO implement
        pass

    @property
    def threshold(self):
        # TODO: most probably this will have to go into the configuration
        return 1

    @property
    def cmd(self):

        input, output = self.input, self.output
        outdir = self.outdir
        load = self.load
        outdir = self.outdir
        priority = " ".join([str(self.gold_score), str(self.silver_score)])
        sources = " ".join([str(self.gold_source), str(self.silver_source)])
        threshold = self.threshold
        maxintron = self.max_intron  # TODO: maybe this should be yet another value in the conf?
        prefix = os.path.join(self.outdir, "portcullis_junctions")
        cmd = "{load} mkdir -p {outdir} && filter_portcullis.py -p {priority} -s gold silver -t {threshold}"
        cmd += " --sources {sources} -mi {maxintron} {input[tab]} {prefix}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], outdir, "output")
