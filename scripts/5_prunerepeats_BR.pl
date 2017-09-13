#!/usr/bin/env perl 

#Open up the file that includes a list of the repeat sequences you want to prune from each gene- one sequence per line

#1) make an array of each gene present in file DupContigsToRemove.txt
#2) make another array of all genes that don't need to be altered.
#3) list all gene files in dir

#3) open each file and print to ED seqs folder
#	if the gene file is in the gene array
#	a) open DupContigsToRemove.txt, make a hash of contigs that are listed for that gene 
#	b) go through it and print only seqs NOT in the contig hash
#	if not just print it with a new name
# 5/12/14 MAG Added some "or die"s on open files for more graceful failure.

$dir=$ARGV[0];

open LOG, ">$dir/2out/5_pruneRepeats.log" or die;

#1)
open FH, "<$dir/DupContigsToRemove.txt" or die;
@geneArray = ();
%geneHash = ();
while(<FH>)
{
	if(/^\S+\s+$dir\/(\S+)\.Seqs/)
	{
		$gene = $1;
		if(! exists $geneHash{$gene})
		{
			$geneHash{$gene} = 1;
			push @geneArray, $gene;
		}
	}
}
@geneArray = sort (@geneArray);
close FH;

for $gene (@geneArray)
{print LOG "Gene in array: $gene\n";}

#2
@files= <$dir/*.Seqs>;
`mkdir $dir/EditedSeqs`;

@goodGeneArray = ();
foreach $file (@files)
{
	if($file =~ /^$dir\/(\S+)\.Seqs/)
	{
		$goodGene = $1;
		if (! exists $geneHash{$goodGene})
		{
			push @goodGeneArray, $goodGene;
		}
	}
}


for $gene (@geneArray)
{
	print LOG "Start of $gene\n";
	foreach $file (@files)
	{
		if($file =~/^$dir\/($gene)\.Seqs/)
		{
			print LOG "Changing $gene\.Seqs\n"; 
			$file = $1;
			%contigHash = ();
			print LOG "Gene found: $gene\n";
			open FH, "<$dir/DupContigsToRemove.txt";
			while(<FH>)
			{
				if(/^(\S+)\s+$gene\./)
				{
					#print "$gene found while searching DupContigsToRemove.txt\n";
					$contig = $1;
					if(! exists $contigHash{$contig})
					{
						#print "$contig found while searching $gene in DupContigsToRemove.txt\n";
						$contigHash{$contig} = 1;
					}
				}
			}
			close FH;
			print LOG "opening $file.Seqs\n";
			
			open FH3, "<$dir/$file.Seqs" or die;
			open OUT, ">$dir/EditedSeqs/$file.ED.Seqs" or die;
			while(<FH3>)
			{
				if(/^(\S+)\s+\d+\s+\d+/)
				{
					$contig2 = $1;
					if (! exists $contigHash{$contig2})
					{
						print OUT;
					}
				}
				else
				{
					print OUT;
				}
			}
			close FH3;
			close OUT;
		}
	}
}

for $goodGene (@goodGeneArray)	
{		
	foreach $file (@files)
	{
		if($file =~ /^$dir\/($goodGene)\.Seqs/)
		{
			$file = $1;
			print LOG "Not changing $file\n";
			open FH3, "<$dir/$file.Seqs" or die;
			open OUT, ">$dir/EditedSeqs/$file.ED.Seqs" or die;
			while(<FH3>)
			{print OUT "$_";}
			close FH3;
			close OUT;
		}
	}	
}	