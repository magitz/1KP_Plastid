#!/usr/bin/env perl

#This script reads through a bunch of .Seqs files (output from 2b_process.pl)
#It finds any sequences that have BlastX hits in the same region against multiple genes
# 5/12/14 MAG Added some "or die"s on open files for more graceful failure.

$dir=$ARGV[0];

@files = <$dir/*.Seqs>;

if (! @files)  #MAG 05/12/14 added checking for existence of files being added to the array. Exit with -1 if no files match $dir/*.Seqs.
{	
	print "No *.Seqs files to process in $dir, exiting\n";
	exit
}

open OUT, ">$dir/RepeatSeqs_GB.txt", or die;
%filehash = ();
%starthash = ();
%endhash = ();
%badhash = ();

foreach $file (@files)
{
	open FH2, "<$file", or die;  
	while (<FH2>)
	{
		if (/^(\S+)\s+(\d+)\s+(\d+)/)
		{
			$seq = $1;
			$num1 = $2;
			$num2 = $3;
	
			if (! exists $starthash{$seq})
			{
				if ($num1 < $num2)
				{
				$filehash{$seq} = $file;
				$starthash{$seq} = $num1;
				$endhash{$seq} = $num2;
				}
				else 
				{
				$filehash{$seq} = $file;
				$starthash{$seq} = $num2;
				$endhash{$seq} = $num1;
				}   
			}
	
			else
			{
				if (! exists $badhash{$seq})
				{
					$badhash{$seq} = 1;		
					print OUT "$seq\t$filehash{$seq}\t$starthash{$seq}\t$endhash{$seq}\n";
					
					if ($num1 < $num2)
					{
						print OUT "$seq\t$file\t$num1\t$num2\n";
					}
					
					else
					{
						print OUT "$seq\t$file\t$num2\t$num1\n";
					}
				}
	
				else
				{
					if ($num1 < $num2)
					{
						print OUT "$seq\t$file\t$num1\t$num2\n";
					}
					else
					{
						print OUT "$seq\t$file\t$num2\t$num1\n";
					}
				}
			}
		}
	}
}

