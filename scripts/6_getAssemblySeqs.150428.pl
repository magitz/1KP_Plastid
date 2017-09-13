#!/usr/bin/env perl -w
use File::Basename;

#This takes the output from 2b_process.pl and the assembly files and outputs for 
#each gene, a fasta file with all of the relevant assembly sequences.  These sequences
#should now be in the right orientation and trimmed.

# Steps 3-5 are now in between step 2 and step 6.
# This script has been modified to take most variables as input rather than needing to change them in here.

# 5/12/14 MAG Added some "or die"s on open files for more graceful failure.

$dir=$ARGV[0];	#The name of the directory to process ($my_dir in 2a_process script).
$clade=$ARGV[1]; #The clade we're working on.

@files= <$dir/EditedSeqs/*.Seqs>;

if (! @files)  #MAG 05/12/14 added checking for existence of files being added to the array. Exit if no files match $dir/*.Seqs.
{	
	print "No *.Seqs files to process in $dir/EditedSeqs/, exiting\n";
	exit
}

`mkdir $dir/Pulled`;

foreach $file (@files)
{
#change next line if needed
    if ($file =~ /$dir\/EditedSeqs\/(\S+\.Seqs)/)
    {
	($gene,$directory,$ext) = fileparse($file, qr/\..*/);
	
	$file = $1;
	print "Getting Sequences from $file ...\n";
	open FH2, "<$dir/EditedSeqs/$file" or die;

	#You may want to change the output file name here
	open OUT, ">$dir/Pulled/$file.fasta" or die;

	
	while (<FH2>)
	{
		
	    if (/$gene\.(\S+).blast/)
	    {
	    $sample = $1;
		$assembly = $sample . ".fa";
		print "$assembly\n";
	    }

	    elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/)
	    {
		$seq = $1;
		$start = $2;
		$end = $3;
		#$dstart = $4;
		#$dend = $5;
		if ($start < $end)
		{
		    $orient = "F";
		}
		else {$orient = "R";}
		$len = abs ($end - $start) + 1;
		if ($start < $end)
		{
			$realstart = $start - 1;
		}
		else {$realstart = $end - 1;}
		
		#Change directory below?
		$punt= "$clade/Assemblies\/$assembly";
		open FH3, "<$punt" or print "Can't open Assemply file for reading, $punt: $!\n"; # 2/18/15: MAG: changed from die to print since some assemblies were removed and are no longer there.
		#print "getting $seq from $assembly\n";
		$ok = 1;
		while (<FH3>)
		{
		    if (/^>/)
		    {
				$ok++;
				if ($ok == 1)
				{
					$overlap = substr($dna, $realstart, $len);
					if ($orient eq "F")
					{
					print OUT "$overlap\n";
					}
					else 
					{
					$revoverlap = reverse $overlap;
					$revoverlap =~ tr/ACGT/TGCA/;
					print OUT "$revoverlap\n";
					}
				}

				if (/$seq/)		#Not sure why Brad had commented this out and went with $seq\s+ since that caused problems with the taxa that are Genus_sp. 
										#I think earlier test removed the "." and now it's back and without it, these don't get pulled.
	#			if (/$seq\s+/)  #Brad had this in place and the line above commented out.
				{
					$dna = "";
					$ok = 0;
					print OUT ">$sample\_$seq\n";
				}
		    }
		    elsif (/^(\S+)/ && $ok == 0)
		    {
			$line = $1;
			chomp $line;
			$dna = "$dna" . "$line";
#			print "$realstart\t$len\n$dna\n";
#			$overlap = substr($dna, $realstart, $len);
#			print OUT "$overlap\n";
		    }
		}
		if ($ok == 0)  #Handle cases where last sequence in file was the one we wanted.
		{
			$overlap = substr($dna, $realstart, $len);
			if ($orient eq "F")
			{
			print OUT "$overlap\n";
			}
			else 
			{
			$revoverlap = reverse $overlap;
			$revoverlap =~ tr/ACGT/TGCA/;
			print OUT "$revoverlap\n";
			}
		}
		#close FH3;
	    }
	}
	close FH2;

    }
}
