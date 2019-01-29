#!/usr/bin/perl
# this script is part of a series of scripts to prepare genotype files
# by transposing it

use strict; use warnings;

die "usage: snp_raw_to_tsv.pl" unless @ARGV == 1;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

my $linenum = 1;
while (my $line = <IN>)
{
	chomp $line;
	my @stuff = split(" ", $line);
	
	my $outname = "out" . $linenum;
	open (my $out, ">", $outname) or die "Could not open $outname\n";
	
	# format TCGA file name
	# from p09.TCGA-AA-3556-10 to TCGA-AA-3556
	if ($linenum != 1)
	{
		my @tcga = split(/\./, $stuff[1]);
		my @tcga2 = split("-", $tcga[1]);
		my $sample = $tcga2[0] . "-" . $tcga2[1] . "-" . $tcga2[2];
	
		print $out "$sample\n";
	}
	else
	{	
		print $out "$stuff[1]\n";
	}
	
	for (my $i = 6; $i < @stuff; $i++)
	{
		print $out "$stuff[$i]\n";
	}
	close $out;
	$linenum++;
}
close IN;
