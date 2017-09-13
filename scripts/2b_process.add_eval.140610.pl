#!/usr/bin/env perl

###########################################################
#
# 2a.process.DATE.pl  e_value_cutoff directory_to_process
#  This script takes the blast output from blasting each sample's assembly to the each gene's reference database
#  and identifies the regions that match.
#
# 
# The outputs are:
#		directory_to_process/refs --Reference sequence for the best blast matches of each gene. Will be used in Genewise steps later.
#		directory_to_process/gene.Seqs with start and stop coordinates for each scaffold blasting to a gene.
######
# 7/16/12
# MAG Modified this to also print file with contig name and name of best blast hit
# This will be used in the genewise step to pull the sequence of the best blast hit to use as the reference in genewise.
#
# 5/21/13
# MAG added a percent identity threshold to help deal with samples that blast to similar genes. 
# I found that often for similar genes (e.g. psaA and psaB), the e-value could be 0 for the wrong gene
# At the same time percent identity was low ~45%. Filter by this. 
#
# 5/23/13
#  MAG added ability to set 2 thresholds for minidentity. Some genes have highly similar ones with cp and also in mtDNA.
#
# 6/28/13
# MAG changed to read list of files into array rather than a separate file.

#Minimum length of nucleotide overlap in query (assembly) sequence to the alignment
$minlength = 100;

#Minimum percentage:  length of overlap in query (assembly) sequence / length of longest sequence in the alignment
$minpercent = 0.25;

#Minimum percent identity to blast targe needed to keep the sample.
$minidentity = 65.0;

@similar_genes=('atpA', 'atpB', 'ndhK', 'ndhC', 'psaA', 'psaB'); #List of genes to use the higher value for minidentity.
																#12/04/13, For some reason despite note above about psaA/B, I hadn't put those in this list!!
$similar_min_identity=80.0;

#
# E-value cutoff to use. Can blast at any value less strignet than what you set here, 
# but can't get results that were excluded in blast
# Eg if blast at 1.0e-5, can be more stringent in pulling results by using  1.0e-10.
#####This is now passed in as a command line option. 
$e_cutoff = $ARGV[0];
$clade= $ARGV[1];
$dir = $ARGV[2];
$ref_path=$ARGV[3]; #5/12/14 To generalize for multiple references, add passing of ref_path. Should not end in /
$blast_out=$ARGV[4]; #5/12/14 added passing of blast output path. Should not end in /

print "Evalue is: $e_cutoff\n";
print "Working in: $dir\n";
@files = <$ref_path/*.faa>;

foreach $file (@files)
{
    if ($file =~ /$ref_path\/(\S+)\.faa/)
    {
		#Getting the length of the longest sequence in the alignment for a gene             
		$longalign = 0;
		$locus = $1;
		#print "$locus\n";
	
		if ($locus ~~ @similar_genes){		#MAG-5/23/13--Allow changing the minidentity based on locus to handle loci that are highly similar amongst themselves and mtDNA copies.
			$locus_minidentity=$similar_min_identity;
		}
		else
		{
			$locus_minidentity=$minidentity;
		}

 
		#Output name
		open OUT, ">$dir/$locus.Seqs" or die "Can't open .Seqs output file: $!";

		#Find longest sequence in alignement
		open FHa, "<$ref_path/$locus.faa" or die "Can't open reference $reaf_path/$locus.faa: $!";
		$len=0;
		$longalign=0;
	
		while (<FHa>)
		{
			if (/^(\S+)/ && ! /^>/)
			{
				$len = $len + length($1); #keep adding the length of each line until next sequence header.
			}
			elsif (/^>/)
			{
				if ($len > $longalign)
				{
					$longalign = $len *3; #These are amino acid lengths, we're going to compare to nt below.
				}
				$len=0; #reset len
			}
		}
		close FHa;
		print "$locus \tLongest= $longalign\n";
		#Find Blast output
		#CHANGE THIS DIRECTORY
		@blast_files=<$blast_out/blasttable/$locus/$locus.*.out>;
		if (! @blast_files)  #MAG 05/12/14 added checking for existence of files being added to the array. Exit if no files match $dir/*.Seqs.
		{	
			print "No blast result files to process in $blast_out/blasttable/$locus/$locus.*.out, exiting\n";
			exit
		}

		
		foreach $blast_file (@blast_files)
		{
			if ($blast_file =~ /$locus\/($locus.\S+).out/)
			{
			$file = $1;  ##MAG: #WOW reusing the file variable in the same loop!!! Careful...I think we're ok, but ???
			#		print "$file\n";
			#Now I'm getting candidate sequences from the assembly
			#CHANGE THIS DIRECTORY
			open FH2, "<$blast_out/blasttable/$locus/$file.out" or die "Can't open $blast_out/blasttable/$locus/$file.out: $!";
			%qstarthash = ();
			%qendhash = ();
			%dstarthash = ();
			%dendhash = ();
			%evaluehash = ();
			@qarray = ();
			$num = 0;
			%ref_hash=();
			while (<FH2>)
			{
				if (/^(\S+)\s+(\S+)\s+(\S+)\s+\d+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/)
				{
					$qseq = $1;
					$ref=$2;
					$pid=$3;
					$qstart = $4;
					$qend = $5;
					$dstart = $6;
					$dend = $7;
					$evalue = $8;
				
					if (($evalue < $e_cutoff) and ( $pid > $locus_minidentity))
					{
						if (! exists $ref_hash{$qseq})
						{
							$ref_hash{$qseq}=$ref;
						}
					
						if (! exists $qstarthash{$qseq})
						{
							$qstarthash{$qseq} = $qstart;
							$qendhash{$qseq} = $qend;
							$dstarthash{$qseq} = $dstart;
							$dendhash{$qseq} = $dend;
							$evaluehash{$qseq}= $evalue;							
							push @qarray, $qseq;
							$num = 0;
						}
					
						else
						{
							if (($qstart < $qend) && ($qend > $qendhash{$qseq}))
							{
								$qendhash{$qseq} = $qend;
				#				$dendhash{$qseq} = $dend;
							}
							elsif (($qstart > $qend) && ($qend < $qendhash{$qseq}))
							{
								$qendhash{$qseq} = $qend;
				#				$dendhash{$qseq} = $dend;
							}
					
							if (($qstart < $qend) && ($qstart < $qstarthash{$qseq}))
							{
								$qstarthash{$qseq} = $qstart;
				#				$dendhash{$qseq} = $dstart; 
							}
							elsif (($qstart > $qend) && ($qstart > $qstarthash{$qseq}))
							{
								$qstarthash{$qseq} = $qstart;
				#				$dendhash{$qseq} = $dstart; 
							}
							
							if ($evalue <$evaluehash{$qseq})
							{
								$evaluehash{$qseq}=$evalue;
							}
						}
					}
				}
			}
			close FH2;
		
			open REF_OUT, ">$dir/refs/$locus/ref.$file.out";
		
			while(($k, $v) = each (%ref_hash) ) 
			{
				print REF_OUT "$k\t$v\n";
			}
		  
			print OUT "\n$file\n";
		
			for $seq(@qarray)
			{
				$querylength = abs ($qendhash{$seq} - $qstarthash{$seq});
				if (($querylength >= $minlength) && ($querylength >= ($minpercent*$longalign) ))
				{
					
					print OUT "$seq\t$qstarthash{$seq}\t$qendhash{$seq}\t$dstarthash{$seq}\t$dendhash{$seq}\t$evaluehash{$seq}\n";
				}
			}
			}
	#	    close FH;
		}
    }
}
