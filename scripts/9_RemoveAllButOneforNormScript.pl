#!/usr/bin/env perl -w
#
#output one sequence from the original profile to a new file
#output all assembly seqs from that taxon to same file.

#need to know number of sequnces in profile alignment. 
#It must be the same number in all alignments.

#Usage: 9_RemoveAllButOneforNormScript.pl number_of_sequences_in_ref fasta_file_to_parse output_file_name

$numSeqs = $ARGV[0];
$input=$ARGV[1];
$output=$ARGV[2];


open INFH, "<$input";
open OUT, ">$output";

$count = 0;
while(<INFH>)
{
	if(/^>/)
	{
		$count++;
		$count1 = 0;
		if($count == 1)
		{
			$count1++;
			print OUT $_;
		}
		if($count > $numSeqs)
		{
			$count1++;
			print OUT $_;
		}
	}
	elsif(! /^>/ && $count1 == 1)
	{
		print OUT $_;
	}
}
close INFH;
