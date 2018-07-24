#!/usr/bin/env perl

# Generate GFF  from exonerate gff file


# v0.2 Release Notes
# Added Alias field to the match line to get the protein name



# Gemy George Kaithakottil 10th March 2015
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
my $prog = basename($0);
my $usage = "
Generate GFF from Exonerate protein2genome GFF output file

## Usage:

	perl $prog --in <file.exonerate.out> --fasta_count [INT] --minIdentity 60 --minCoverage 50

## Required options:
  --in=<exonerage_protein2genome.output>  -- Exonerate output file 
                                             See NOTE below to see how exonerate must have ran
  --fasta=n                               -- Input fasta sequence to exonerate
                                             Basically the fasta sequences you have given into exonerate

## Optional:
  --minIdentity=n <INT>                   -- alignments with identity shorter than this are discarded
                                             [default 0 - i.e., No Cutoff]
  --minCoverage=n <INT>                   -- alignments with coverage shorter than this are discarded
                                             [default 0 - i.e., No Cutoff]
  --help                                  -- prints this option menu and quit

## NOTE:
  Exonerate to be ran with options:
  --model protein2genome --showtargetgff yes --showvulgar yes --softmaskquery yes --softmasktarget yes --bestn 10  --minintron 20 --maxintron 50000 --ryo \">%qi\\tlength=%ql\\talnlen=%qal\\tscore=%s\\tpercentage=%pi\\nTarget>%ti\\tlength=%tl\\talnlen=%tal\\n\"
  
 Output to STDOUT

Contact: Gemy.Kaithakottil\@tgac.ac.uk
         gemygk\@gmail.com
";
my $exfile;
my $fasta="";
my $fasta_count = 0;
my $minIdentity = 0;
my $minCoverage = 0;
my $help;
#if ($#ARGV < 2 ) {
#    print "$usage";
#    exit;
#}

GetOptions(
           'in=s'=>\$exfile,
           'fasta:s'=>\$fasta,
           'minIdentity:i'=>\$minIdentity,
           'minCoverage:i'=>\$minCoverage,
           "h|help" => \$help
           );
#if (@ARGV==0) {
#  die $usage;
#}
if ($help) {
    die $usage;
}
unless ($exfile && $fasta ) {die "\n## ERROR:Need input file --in & --fasta options\n$usage\n"};
my %transList=();
my $nextRecord=0;
my $count=1;
my $path=1;
my $transcript;
my $boolean=0;
my @pre_transcript=();


my %full_trans=();
my %full_aligns=();
my %fltr_trans=();
my %fltr_aligns=();
open(XNT, "<$exfile") || die "Couldn't open $exfile\n";
while (<XNT>) {
    next if /^#/;
    next unless /\S/;
	chomp;
   
   if($nextRecord) {
	 # Clearing earlier records
	 #$count=1;
   }
   
   if (/^vulgar:/) {
   		# skip the first record
   		if($boolean) {
   			# The array is in the format
   			# name region
   			# GFF match
   			# GFF match_part
   			# GFF match_part
   			# GFF attributes ... and so on with "name region"

   			# get the alignment region -- first of arrary
   			my ($trans1,$region) = split /\t/, $pre_transcript[0];
   			$region = "aln:" . $region;
   			# get the aln len, len of trans, cov, identity, score - from last of array
   			my ($trans2,$attrib) = split /\t/, $pre_transcript[$#pre_transcript];
   			my $note=$region ."|".$attrib;
   			#print "Alignments attributes:$attrib\n";
   			# the columns are:
   			#alnLen=347|len=430|cov=80.70|identity=76.38|score=1246
   			my @a = split /\|/,$attrib;
   			my ($curr_cov) = $a[2] =~ /cov:(\d+\D\d+)/;
   			my ($curr_id) = $a[3] =~ /id:(\d+\D\d+)/;
   			#print "Current coverage:$curr_cov\tCurrent identity:$curr_id\n";

   			#print "$trans1\t$region|$attrib\n";
   			# Only printing the alignments that are over the cutoff
   			if ($curr_cov >= $minCoverage && $curr_id >= $minIdentity) {
	   			# NOT printing the first and last element of array
	   			for (my $i=1;$i<$#pre_transcript;$i++) {
	   			#for (my $i=0;$i<$#pre_transcript+1;$i++) {
	   				my @f = split /\t/,$pre_transcript[$i];
	   				my $type=$f[2];
	   				if ($type eq "match") {
	   					print "$pre_transcript[$i];Note=$note\n";

	   					# store the filtered transcripts and alignments counts
	   					my @match_line = split /\t/, $pre_transcript[$i];
	   					my ($passFiltr_align) = $match_line[8] =~ /Name=([^;]+)/;
	   					my ($passFiltr_trans) = $match_line[8] =~ /Name=(.*)\.path\d+/;
	   					#print "$match_line[8]\t$passFiltr_align\n";
	   					#print "$match_line[8]\t$passFiltr_trans\n";
						unless (exists $fltr_trans{$passFiltr_trans}) {
							$fltr_trans{$passFiltr_trans} += 1;
						}
						unless (exists $fltr_aligns{$passFiltr_align}) {
							$fltr_aligns{$passFiltr_align} += 1;
						}
	   				} else {
	   					print "$pre_transcript[$i]\n";
	   				}
	   			}
	   		}
   			@pre_transcript=();
   		}
		#print "$_\n";
		$boolean=1;
		my @f = split /\s/;
		my $name=$f[1];
		my $start=$f[2]+1;
		my $end=$f[3];
		#print "$name\t$start-$end\n";
		push (@pre_transcript,"$name\t$start-$end");

		if (exists($transList{$name})) {
			$path++;
		} else {
			$transList{$name} +=1;
			$path=1;
		}
		
		$nextRecord=1;

		# Clearing earlier records
	 	$count=1;
   } 
   if (/exonerate:protein2genome:local/i) {
		#print "$_\n";
		my @f = split /\t/, $_, 9;
		if (@f < 8) { warn "Not gff format"; next }
		my $seqname = $f[0];
		my $source = $f[1];
		my $type = $f[2];
		my $start = $f[3];
		my $end = $f[4];
		my $score = $f[5];
		my $strand = $f[6];
		if ($end < $start) {
			my $tmp = $start;
			$start = $end;
			$end = $tmp;
		}
		if ($type eq "gene") {
			($transcript) = $f[8] =~ /sequence\s*([^;]+)/;
			$transcript =~ s/\s+$//;
			#print "$seqname\t$source\tmatch\t$start\t$end\t$score\t$strand\t.\tID=$transcript;Name=$transcript;Note=$transcript\n";
			push (@pre_transcript,"$seqname\t$source\tmatch\t$start\t$end\t$score\t$strand\t.\tID=$transcript.path$path;Name=$transcript.path$path;Alias=$transcript");

			# store the full transcripts and alignments counts
			unless (exists $full_trans{$transcript}) {
				$full_trans{$transcript} += 1;
			}
			unless (exists $full_aligns{"$transcript.path$path"}) {
				$full_aligns{"$transcript.path$path"} += 1;
			}


		} elsif ($type eq "exon") {
			my ($insertions) = $f[8] =~ /insertions\s*([^;]+)/;
			$insertions =~ s/\s+$//;
			my ($deletions) = $f[8] =~ /deletions\s*([^;]+)/;
			$deletions =~ s/\s+$//;
			#print "$seqname\t$source\tmatch_part\t$start\t$end\t$score\t$strand\t.\tID=$transcript.exon$count;Parent=$transcript;insertions=$insertions;deletions=$deletions\n";
			push (@pre_transcript,"$seqname\t$source\tmatch_part\t$start\t$end\t$score\t$strand\t.\tID=$transcript.path$path.exon$count;Parent=$transcript.path$path;insertions=$insertions;deletions=$deletions");
			$count++;
		}
   }
   if (/^>/) {
		my @f = split /\t/;
		my $name=$f[0];
		$name=~s/>//;
		my ($len)     = $f[1] =~ /length=(\d+)/;
		my ($alnlen)  = $f[2] =~ /alnlen=(\d+)/;
		my ($score)   = $f[3] =~ /score=(\d+)/;
		my ($percent) = $f[4] =~ /percentage=(\d+\D\d+)/;
		my $coverage  = sprintf("%0.2f",( ($alnlen/$len)*100 ) );
		#print "$name\talnLen=$alnlen|length=$len|coverage=$coverage|identity=$percent|score=$score\n";
		push (@pre_transcript,"$name\talnLen:$alnlen|len:$len|cov:$coverage|id:$percent|score:$score");
   }

   if(eof(XNT)) {
		# The array is in the format
		# name region
		# GFF match
		# GFF match_part
		# GFF match_part
		# GFF attributes ... and so on with "name region"

		# get the alignment region -- first of arrary
		my ($trans1,$region) = split /\t/, $pre_transcript[0];
		$region = "aln:" . $region;
		# get the aln len, len of trans, cov, identity, score - from last of array
		my ($trans2,$attrib) = split /\t/, $pre_transcript[$#pre_transcript];
		my $note=$region ."|".$attrib;
		#print "Alignments attributes:$attrib\n";
		# the columns are:
		#alnLen=347|len=430|cov=80.70|identity=76.38|score=1246
		my @a = split /\|/,$attrib;
		my ($curr_cov) = $a[2] =~ /cov:(\d+\D\d+)/;
		my ($curr_id) = $a[3] =~ /id:(\d+\D\d+)/;
		#print "Current coverage:$curr_cov\tCurrent identity:$curr_id\n";

		#print "$trans1\t$region|$attrib\n";
		# Only printing the alignments that are over the cutoff
		if ($curr_cov >= $minCoverage && $curr_id >= $minIdentity) {
			# NOT printing the first and last element of array
			for (my $i=1;$i<$#pre_transcript;$i++) {
			#for (my $i=0;$i<$#pre_transcript+1;$i++) {
				my @f = split /\t/,$pre_transcript[$i];
				my $type=$f[2];
				if (!defined $type) {
					next;
				}

				if ($type eq "match") {
					print "$pre_transcript[$i];Note=$note\n";

					# store the filtered transcripts and alignments counts
					my @match_line = split /\t/, $pre_transcript[$i];
					my ($passFiltr_align) = $match_line[8] =~ /Name=([^;]+)/;
					my ($passFiltr_trans) = $match_line[8] =~ /Name=(.*)\.path\d+/;
					#print "$match_line[8]\t$passFiltr_align\n";
					#print "$match_line[8]\t$passFiltr_trans\n";
				unless (exists $fltr_trans{$passFiltr_trans}) {
					$fltr_trans{$passFiltr_trans} += 1;
				}
				unless (exists $fltr_aligns{$passFiltr_align}) {
					$fltr_aligns{$passFiltr_align} += 1;
				}
				} else {
					print "$pre_transcript[$i]\n";
				}
			}
		}
	}
}
close(XNT);

# Get the count of FASTA sequences
if (not -e "$fasta.fai") {
    `samtools faidx $fasta`
}

$fasta_count = `cat $fasta.fai | wc -l`;

# Print the number of transcripts and alignments
my $input_transcripts = scalar (keys %full_trans);
my $input_coverage = sprintf("%0.2f",( ($input_transcripts/$fasta_count*100 ) ) );
print STDERR "\n## Input Stats:\n";
print STDERR "## Number of transcripts aligned     => $input_transcripts ($input_coverage %)\n";
print STDERR "## Number of transcript alignments   => " . scalar (keys %full_aligns) . "\n";
my $filtered_transcripts = scalar (keys %fltr_trans);
my $filtered_coverage = sprintf("%0.2f",( ($filtered_transcripts/$fasta_count*100 ) ) );
print STDERR "\n## Output Stats (Filtered at coverage >= $minCoverage% and identity >= $minIdentity%):\n";
print STDERR "## Number of transcripts aligned (pass filter)    => $filtered_transcripts ($filtered_coverage %)\n";
print STDERR "## Number of transcript alignments (pass filter)  => " . scalar (keys %fltr_aligns) . "\n\n";
exit;
