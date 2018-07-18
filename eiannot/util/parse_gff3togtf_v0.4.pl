#!/usr/bin/env perl
#
# For given input GFF file convert to GTF
#

# v0.4 Release # Wednesday, 07 September 2016, 11:39AM
# - Add gene_name and transcript_name attributes to all feature types
# - Add exon_id to exon
# - Add protein_id to CDS

# v0.3 Release
# - Added ncRNA feature (Minor update)

# v0.2 Release
# - Added transcript line printing to output (Minor update)
#
#
# AUTHOR: Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk || gemygk@gmail.com)

use strict;
use warnings;
use File::Basename;

my $prog = basename($0);
my $usage = "
	$prog -- For given input GFF file convert to GTF

	Usage: $prog <*.gff3>

	The input GFF3 *MUST* be sorted and can have following types:
		* gene==>[mRNA|ncRNA|transcript]==>[exon|CDS|five_prime_UTR|five_prime_utr|three_prime_UTR|three_prime_utr]
		* match==>match_part

Contact: Gemy.Kaithakottil\@earlham.ac.uk || Gemy.Kaithakottil\@gmail.com
";

if (@ARGV==0) {
  die $usage;
}

my $gff = shift or die $usage;

my $gene="";
my $g_mrna="";
my $exon=0;

my $boolean=1;
my $exon_count=1;
my $cds_count=1;
#print "##gff-version 3\n";
open(ID ,"< $gff") || die "Cannot open $gff:$!\n";
while (<ID>) { # Read lines from file(s) specified on command line. Store in $_.
    next if /^#/; # Remove comments from $_.
    next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
	chomp;

    my @f = split /\t/; # split $_ at tabs separating fields.

    # if gene do not print the line
	if($f[2] eq "gene") {
		# Reset gene counters
		($gene) = $f[8] =~ /ID\s*=\s*([^;]+)/;
	}
	# if match I should print the line
	if($f[2] eq "match" ) {
		# Reset gene counters
		($gene) = $f[8] =~ /ID\s*=\s*([^;]+)/;
		($g_mrna) = $f[8] =~ /ID\s*=\s*([^;]+)/;
		print "$f[0]\t$f[1]\ttranscript\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; gene_name \"$gene\"; transcript_name \"$g_mrna\";\n";
		$exon_count=1;
		$cds_count=1;
	}
	if($f[2] =~ /^(mRNA|ncRNA|transcript)$/ ) {
		($g_mrna) = $f[8] =~ /ID\s*=\s*([^;]+)/;
		print "$f[0]\t$f[1]\ttranscript\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; gene_name \"$gene\"; transcript_name \"$g_mrna\";\n";
		$exon_count=1;
		$cds_count=1;
	}
	if($f[2] eq "match_part" ) {
		($g_mrna) = $f[8] =~ /Parent\s*=\s*([^;]+)/;
		print "$f[0]\t$f[1]\texon\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; exon_number \"$exon_count\"; gene_name \"$gene\"; transcript_name \"$g_mrna\"; exon_id \"$g_mrna.exon_$exon_count\";\n";
		$exon_count++;
	}
	if($f[2] eq "five_prime_utr" || $f[2] eq "five_prime_UTR") {
		print "$f[0]\t$f[1]\t5UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; gene_name \"$gene\"; transcript_name \"$g_mrna\";\n";
	}
	if($f[2] eq "exon") {
		print "$f[0]\t$f[1]\texon\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; exon_number \"$exon_count\"; gene_name \"$gene\"; transcript_name \"$g_mrna\"; exon_id \"$g_mrna.exon_$exon_count\";\n";
		$exon_count++;							# v0.2 update
	}
	if($f[2] eq "CDS") {
		print "$f[0]\t$f[1]\tCDS\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; gene_name \"$gene\"; transcript_name \"$g_mrna\"; protein_id \"$g_mrna.cds_$cds_count\";\n";
		$cds_count++;
	}
	if($f[2] eq "three_prime_utr" || $f[2] eq "three_prime_UTR") {
		print "$f[0]\t$f[1]\t3UTR\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$g_mrna\"; gene_name \"$gene\"; transcript_name \"$g_mrna\";\n";
	}
}
close(ID);