#!/usr/bin/env perl
#
# Script to parse Mikado Gold GFF to TransDecoder GFF

# v0.1 # Saturday, 17 March 2018, 02:47PM
#################################################################################
#################################################################################

# This script is released under the MIT license.

# MIT License

# Copyright (c) 2018 Gemy George Kaithakottil

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

#################################################################################
#################################################################################

#
# AUTHOR: Gemy George Kaithakottil (Gemy.Kaithakottil@gmail.com)

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
my $prog = basename($0);

my $usage = "
	Script to parse Mikado Gold GFF to TransDecoder GFF

	Usage: $prog <source> <genewise.gff>

	where,
	<source> - source field for output GFF3 (column 2)

Contact: Gemy.Kaithakottil\@gmail.com
";

my $source = shift or die $usage;
my $file = shift or die $usage;

open (TEXT_FILE,"< $file") or die "$!";
my %gene_info=();
my $gene_boo=0;
my $cds_count=1;
my $mrna_count=0;
my $exon_count=0;
my $utr5_count=0;
my $utr3_count=0;

while (<TEXT_FILE>) {
    next if /^#/; # Remove comments from $_.
    next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
	chomp;
	s/\r//g; # remove ^M from the line
	my @f = split (/\t/); # split $_ at tabs separating fields.
	my $type = $f[2];
	my $attrib = $f[8];
    next unless ($type =~ /^(mRNA|exon|CDS|five_prime_UTR|three_prime_UTR)$/i ); # only consider main GFF features
	if ($type =~ /^(mRNA)$/i) {
		$mrna_count++;
		$exon_count = 0;
		$utr5_count=0;
		$utr3_count=0;
		my ($id) = $attrib =~ /ID=([^;\n]+)/;
		if ($gene_boo) {
			print "\n";
		}
		$gene_boo=1;

		$attrib = $attrib . ";Note=$id";
		$cds_count=1;
		# print the GFF3 lines
		my $g_line = join "\t",$f[0],$source,"gene",$f[3],$f[4],$f[5],$f[6],$f[7],"ID=$id|g.$mrna_count";
		print $g_line . "\n";
		my $m_line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],"ID=$id|m.$mrna_count;Parent=$id|g.$mrna_count";
		print $m_line . "\n";
	} elsif ($type =~ /^(exon)$/i) {
		$exon_count++;
		my ($id) = $attrib =~ /Parent=([^;\n]+)/;
		# print the GFF3 lines
		my $line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],"ID=$id|m.$mrna_count.exon$exon_count;Parent=$id|m.$mrna_count";
		print $line . "\n";
	} elsif ($type =~ /^(CDS)$/i) {
		my ($id) = $attrib =~ /Parent=([^;\n]+)/;
		# print the GFF3 lines
		my $line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],"ID=cds.$id|m.$mrna_count;Parent=$id|m.$mrna_count";
		print $line . "\n";
	} elsif ($type =~ /^(five_prime_UTR)$/i) {
		$utr5_count++;
		my ($id) = $attrib =~ /Parent=([^;\n]+)/;
		# print the GFF3 lines
		my $line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],"ID=$id|m.$mrna_count.utr5p$utr5_count;Parent=$id|m.$mrna_count";
		print $line . "\n";
	} elsif ($type =~ /^(three_prime_UTR)$/i) {
		$utr3_count++;
		my ($id) = $attrib =~ /Parent=([^;\n]+)/;
		# print the GFF3 lines
		my $line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],"ID=$id|m.$mrna_count.utr3p$utr3_count;Parent=$id|m.$mrna_count";
		print $line . "\n";
	}
	# if ($type =~ /^(cds)$/i) {
	# 	my ($id) = $attrib =~ /ID=([^;\n]+)/;
	# 	my ($name) = $attrib =~ /Name=([^;\n]+)/;
	# 	$type = "match_part";
	# 	if (exists $gene_info{$id}) {
	# 		$attrib = "ID=$id.cds$cds_count;Parent=$id;Name=$name";
	# 		$cds_count++;
	# 	} else {
	# 		die "Error: Transcript's '$id' parent gene not processed before, check the gff3 for the transcript\n";
	# 	}
	# }
	# # print the GFF3 lines
	# my $x_line = join "\t",$f[0],$source,$type,$f[3],$f[4],$f[5],$f[6],$f[7],$attrib;
	# print $x_line . "\n";
}
exit;

