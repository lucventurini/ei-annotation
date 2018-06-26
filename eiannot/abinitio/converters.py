from ..abstract import AtomicOperation, EIWrapper
from ..rnaseq.mikado.workflow import Mikado
from ..rnaseq.alignments.portcullis import PortcullisWrapper, PortcullisMerge, PortcullisMergeFailed
from ..repeats.workflow import RepeatMasking
from ..proteins.workflow import ExonerateProteinWrapper
import os


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

outdir = os.path.join( "abinitio", "3-Hints")

class ConvertToHints(EIWrapper):

    def __init__(self):
        pass


class ConvertMikado(AtomicOperation):

    pass

class ConvertRepeats(AtomicOperation):

    pass

class ConvertProteins(AtomicOperation):
    pass


class ConvertJunctions(AtomicOperation):

    def __init__(self, junctions: (PortcullisMerge, PortcullisMergeFailed)):

        super().__init__()
        self.configuration = junctions.configuration
        self.input["gff3"] = junctions.input["gff3"]
        self.is_pass = (isinstance(junctions, PortcullisMerge))
        self.output["gff3"] = os.path.join(self.outdir, "portcullis_junctions_{passed}.gff3".format(
            passed="pass" if self.is_pass else "failed"
        ))

    @property
    def rulename(self):
        return "prepare_hints_junction_{passed}".format(passed="pass" if self.is_pass else "failed")

    @property
    def loader(self):
        return []

    @property
    def priority(self):
        # TODO: most probably this will have to go into the configuration
        if self.is_pass:
            return 6
        else:
            return 4

    @property
    def cmd(self):

        input, output, priority = self.input, self.output, self.priority
        outdir = self.outdir

        cmd = "mkdir -p {outdir} && sed 's/grp/group/; s/;$//; s/$/;pri={priority}/' {input[gff3]} > {output[gff3]}"
        cmd = cmd.format(**locals())
        return cmd

    @property
    def outdir(self):
        return os.path.join(self.configuration["outdir"], outdir, "output")


class GetCoverage(AtomicOperation):
    pass
