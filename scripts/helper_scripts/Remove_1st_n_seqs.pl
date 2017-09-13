#!/usr/bin/env perl -w

#Remove the first n sequeces from a fasta file.

#Usage: 9_RemoveAllButOneforNormScript.pl number_of_sequences_to_remove fasta_file_to_parse output_file_name

$numSeqs = $ARGV[0];
$input=$ARGV[1];
$output=$ARGV[2];


open INFH, "<$input" or die "Can't open $input, $!\n";
open OUT, ">$output" or die "Can't open $output, $!\n";

$count = 0;
while(<INFH>)
{
	if(/^>/)
	{
		$count++;
	}
	if($count > $numSeqs)
	{
		print OUT $_;
	}
}
close INFH;
