#!/usr/bin/perl

# Input is the 6_getAssemblySeqs output:
#	For each gene, there is a file with all of the fragments for all of the species.
#	Fasta headers will look like: CODE_scaffold-name or Genus_species_scaffold-name
#	 put codes in each file in an array
# 	 for each item in array loop through file and spit out those seqs.to individual files.
# 	 these file should be labeled by code and gene

#####
##
# 7/16/12
# MAG--modified to add pulling the best blast hit using the file made in step 2.
# use this to run genewise
###
###
$dir=$ARGV[0];
$clade=$ARGV[1];
$AA_path=$ARGV[2];  #5/12/14 To generalize for multiple references, add passing of AA_path. Will not have trailing /.

########
# Turn debuging on and off
# Set to 2 for lots of logging, 1 for normal debugging, 0 for no debugging
$debug=1;
########
if($debug > 0)
{	print "Running with degug log\n";
	open DEBUG, ">$dir/debug.log" or die "Can't open debugging log file: $!\n";
}


#Get $TMPDIR from environment and check it's set--should be set if run through scheduler.
$tmpdir = ''; # Getting core dumps using local scratch disk, so skip for now $ENV{'TMPDIR'};
if (length($tmpdir) == 0)
{
	system "mkdir $dir/temp"; #If no TMPDIR set, make directory and set TMPDIR.
	print "No TMPDIR set in enviroment, unsing $dir/temp instead\n";
	$tmpdir="$dir/temp";
}
else 
{
	print "Using temp dir from PBS: $tmpdir \n";
}

#make a list of all filenames (each gene)
@files = <$dir/Pulled_Good/*seqs.cleaned.fasta> ;


system "mkdir $dir/Genewise";

foreach $file (@files)
{      
	if (! @files)  #MAG 05/12/14 added checking for existence of files being added to the array. Exit with -1 if no files match $dir/*.Seqs.
	{	
		print "No *seqs.cleaned.fasta files to process in $dir/Pulled_Good/, exiting\n";
		exit
	}

     
	if($debug > 1){print DEBUG "Making hash for: $file\n";}

	$gene = `basename $file .seqs.cleaned.fasta`;
	chomp($gene);
	
	if($debug > 1){print DEBUG "Gene is: $gene\n";}

	$ok = 0;		
	@taxon_array = ();
	%taxon_hash = ();
	
	# 1) open it
	open FH2, "<$file" or die "Can't open file: $!";
	
	# 2) create an array of unique taxon codes
	while(<FH2>)
	{
		if(/^>(\S+)_\S+/)	#04/30/15: Make non-greedy to get up to the last _ (Use the whole name not just the first bit before _).
		{
			$taxon_code = $1;
			if (exists $taxon_hash{$taxon_code})
			{
				$taxon_hash{$1}++;

			}
			else
			{
				push @taxon_array, $taxon_code;
				$taxon_hash{$1} = 1;
				if($debug > 1){print DEBUG "\tAdding $1 new value is: $taxon_hash{$1}\n";}
			}
		}
	}	
	close FH2;

#	Go back through the files, and for each sequence, run genewise, parse the results into a file that
#     has the longest cdna portion for all of the contigs for each taxon for each gene in one file.		
	
	#make folders for each gene
	system "mkdir $dir/Genewise/$gene";
	
	for $taxon_code(@taxon_array)
	{
		open FH2, "<$file" or die "Can't open file: $!"; 
		if($debug > 0){print DEBUG "\n\n\nWorking on taxon:$taxon_code \n";}	
		####Make a hash with the names of contigs and names of their best blast hits.
		#Open the file that has the contig names and the name of the best blast hit for each.
		# append the gene and taxon specific stuff to the $ref_path variable
		$ref_path = "$dir/refs/";

		$ref_path=$ref_path . $gene . "\/ref." . $gene . "." . $taxon_code . ".blast.out";
		
		@ref_file=glob($ref_path);
		if($debug > 1){print DEBUG "Ref path: $ref_file[0]\n";}

		$num_refs=@ref_file;
		if($debug > 0){print DEBUG "There are $num_refs matches for ref_file: $ref_path\n";}
		if($num_refs > 1 )
		{
			print "WARNING !!!!!!!!!!!!!!!!\n";
			print "More than 1 ref file matched for search with $ref_path\n";
		}
		open REF, "<$ref_file[0]" or die "Can't open reference file $ref_file[0]: $!\n";
		%ref_hash=();
		#read that file into a hash
		while (<REF>)
		{
			if($debug > 1){print DEBUG "Getting reference from $ref_file[0]\n";}
			chomp $_; #remove line break
			@line_bits=split(/\t/,$_);
			
			$ref_hash{ $line_bits[0] } = $line_bits[1];
			if($debug > 1){print DEBUG "\tAdding $ref_hash{$line_bits[0]} to ref_hash\n";}

		}
		close REF;
		$orig_seq_num=0; #counter for number of sequences for a taxon at a gene
		$trimmed_seq_num=0; #counter for number of sequences after genewise--should be the same.

		#Go through the file again, this time finding the ref seq and running genewise.
		while (<FH2>)
		{
			#If line is fasta description:
			if (/^>$taxon_code/)
			{
				$orig_seq_num++;
				$header=$_;
				# Find the reference sequence for the contig in question and make a 
				# temp file with that reference AA sequence.
				while ( ($key, $value) = each(%ref_hash) ) {
					# The original contig names used to make the ref_hash are part of 
					# the name the contig has now. But the name has had some extra 
					# stuff added to it, so this seemed like an inefficient, but
					# workable way to find the match.
						if($debug > 1){print DEBUG "Looking for header: $header with key: $key\n";}

					 if ($header=~/$key/)
					 {
						if($debug > 1){print DEBUG "\tKey is part of header...getting ref\n";}
					
					$AA_file=$AA_path . "/" . $gene . ".faa";
						open REF, "<$AA_file" or die "Can't open AA genefile $AA_file: $!\n";

					 if($debug > 1){print DEBUG "Reading $AA_file, looking for $ref_hash{$key}\n";}
					while (<REF>)
						{
							if (/$ref_hash{$key}/) #Found header we're looking for. ## 12/24/14: MAG: removed ^> from start of regex as names now have orders added to them so they weren't matching from the start.
							{
								###GET that sequence
								if($debug > 1){print DEBUG "Found header with $ref_hash{$key}\n";}
								$inseq=1;
								$ref_seq="";
								while ($inseq == 1)
								{
									$line=<REF>; #get the next line
									if($line=~ /^>/) {  #Keep adding lines until next fasta header
										$inseq=2;
									}
									elsif (eof(REF)) { #found end of file--this was up above as an ||, but didn't work where there was no line break at end of file.
										chomp($line);
										$ref_seq=$ref_seq . $line;  #still want to get that last line added to the ref_seq.
										$inseq=2;
									}
									else{
										chomp($line);
										$ref_seq=$ref_seq . $line;
									}
								}
										
								if($debug > 1){print DEBUG "sequence is $ref_seq\n";}
								open REF_TEMP, ">$tmpdir/temp.ref.fa" or die "Can't open REF_TEMP:$tmpdir/temp.ref.fa: $!\n";
								print REF_TEMP ">" . $ref_hash{$key} . "\n";
								print REF_TEMP $ref_seq;
							}
						}
					}
				}
				#Make a temp fasta file with the nt contig sequence being looked at.		
				open SEQ_TEMP, ">$tmpdir/temp.seq.fa" or die "Can't open SEQ_TEMP:$!\n";
				print SEQ_TEMP ">$taxon_code\n";  # 5/20/13: Had been using $header, but the long names caused segmentation faults. # 6/27/13--want to get beck to the full name for recovering the scaffold later--See printing PEPOUT below.
				if($debug > 1){print DEBUG "Sequence header: $header\n";}
				#get the next line in FH2, which has the sequence ###7/28/13 MAG: Fix this so that we aren't assuming single line sequences.
				$inseq=1;			
				$tempseq_seq="";
				while ($inseq == 1)
				{
					$line=<FH2>; #get the next line
					if($line=~ /^>/) {  #Keep adding lines until next fasta header
						$inseq=2;
					}
					elsif (eof(FH2)) { #found end of file--this was up above as an ||, but didn't work where there was no line break at end of file.
						chomp($line);
						$tempseq_seq=$tempseq_seq . $line;  #still want to get that last line added to the ref_seq.
						$inseq=2;
					}
					else{
						chomp($line);
						$tempseq_seq=$tempseq_seq . $line;
					}
				}
								
				
				chomp($tempseq_seq);
				if($debug > 1){print DEBUG "Sequence: $tempseq_seq\n";}
				print SEQ_TEMP $tempseq_seq ;
				
				if($debug > 0){print DEBUG "Running genewise for $gene with $taxon_code\n";}
				#run genewise using the temp files
				$genewise=`genewise $tmpdir/temp.ref.fa $tmpdir/temp.seq.fa -cdna -sum > $dir/Genewise/$gene/$gene.genewise.out`; 
			
				#Parse the genewise output
				open PEPOUT, ">>$dir/Genewise/$gene/$taxon_code.genewise.out";
				open GW, "<$dir/Genewise/$gene/$gene.genewise.out";

				
				$current_seq="";
				$longest=0;
				$in_seqs=0;
				$num_frags=0;
				$long_name;
				#Parse genewise output to find the longest and print that to the final file.
				while (<GW>)
				{
					if(/\/\//)
					{}
					elsif (/^>/)
					{
						if($debug > 0){print DEBUG "Found genewise result header\n";}
						$num_frags++;
						$in_seqs=1; #now past all the header information and into the sequence output
						if (length($current_seq) > $longest)
						{
							if($debug > 0){print DEBUG "New longest seq: $current_name of length" . length($current_seq) . " \n";}

							$long_seq=$current_seq;
							$long_name=$current_name;
						}
						
						$current_name=$_;
						$current_seq="";
					}
					
					elsif ($in_seqs== 1)
					{
						chomp($_);
						$current_seq=$current_seq . $_;
						if($debug > 0){print DEBUG "getting genewise results sequences: result $num_frags\n";}
					}
				}
				close GW;
				#handle the last sequence in the file.
				if (length($current_seq) > $longest)
				{
					$long_seq=$current_seq;
					$long_name=$current_name;
				}
				if($debug > 0)
				{	
					print DEBUG "For $taxon_code, $long_name had a total of $num_frags\n";
				}
				print PEPOUT $header; ##Was $long_name;--6/27/13 changed to header to keep full name in output for blasting later.
				print PEPOUT $long_seq . "\n";
				$trimmed_seq_num++;

				if ($num_frags lt 1)
									{
											chomp($gene);
											chomp($taxon_code);
											chomp($header);
											$header=~s/\>//;
											#Genewise failed for some reason. Move tmp files to 2error folder to help diagnose issue.
											print "Genewise failed for $header at $gene with $taxon_code, with error: $genewise\n";
											`mv $tmpdir/temp.seq.fa $dir/2error/$gene.$taxon_code.$header.temp.seq.fa`;
											`mv $tmpdir/temp.ref.fa $dir/2error/$gene.$taxon_code.$header.temp.ref.fa`;

									}
				else
				{
					`rm $tmpdir/temp.seq.fa`;
					`rm $tmpdir/temp.ref.fa`;
					`rm $dir/Genewise/$gene/$gene.genewise.out`;
				}
			}
		}
		if ($orig_seq_num != $trimmed_seq_num)
		{
			print "WARNING LOST SOME SEQUENCES IN GENEWISE STEP FOR $gene with $taxon_code\n";
			print "Started with $orig_seq_num, but only got $trimmed_seq_num in the end!!!!\n";		
		}
	}
	
}


