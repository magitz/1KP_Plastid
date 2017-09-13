#/usr/bin/env perl 

#Script will open the RepeatSeqs_GB.txt file which contains all contigs that blasted to 
# more than one gene. It will then pick out those regions that are overlapping
# and print them to a file so they can be removed with a subsequent script 
# Any non-overlapping regions will be saved.
# 5/12/14 MAG Added some "or die"s on open files for more graceful failure.


$dir=$ARGV[0];

#input
open FH, "<$dir/RepeatSeqs_GB.txt" or die;

#output
open LOG, ">$dir/2out/log.txt" or die;
open OUT, ">$dir/DupContigsToRemove.txt" or die;

%contigHash = ();
@contigArray = ();
$totalSeqCount = 0;

#read RepeatSeqs_GB.txt and put the unique contig names in an array
while (<FH>)
{
	if(/^(\S+)\t\S+\t\d+\t\d+\s+/)
	{
		$contig = $1;
		if(! exists $contigHash{$contig})
		{
			push @contigArray, $contig;
			$contigHash{$contig} = 1; 
		}
		elsif(exists $contigHash{$contig})
		{
			$contigHash{$contig}++;
		}
	}
}
close FH;

# for each contig name determine if any of the positive blast regions on that contig
# are overlapping and print out the name of any regions that overlap with any others to a
# file that will be used to exclude those contigs/regions from the appropriate gene files
foreach $item (@contigArray)
{
	@RemoveArray = ();
	@beginArray = ();
	@endArray = ();
	%beginHash = ();
	%endHash = ();
	
	#make 2 hashes for the contig. The first is the beginning of each gene in that contig
	#the second is the end of each gene in the contig. 
	#In each hash Key is gene name. Value is start/end position.
	open FH, "<$dir/RepeatSeqs_GB.txt" or die;
	while (<FH>)
	{
		if(/^$item\t(\S+)\t(\d+)\t(\d+)\s+/)
		{
			$gene = $1;
			$begin = $2;
			$end = $3;

			if(! exists $beginHash{$gene})
			{
				$beginHash{$gene} = $begin;
				$endHash{$gene} = $end;
			}
		}
	}
	close FH;


	#sort beginHash by starting position of each gene
	#populate %positionHash with position of gene in %beginHash as key and gene name as value
	#for each gene in %beginHash push its start position to beginArray
	#and push the corresponding end value of that gene from end hash to end array
	# now for each gene in order, the start value ans end value are in the same order in two 
	#different arrays.

	sub hashValueAscendingNum 
	{
	   $beginHash{$a} <=> $beginHash{$b};
	}
	#print "$item\n";
	%positionHash = ();
	$posNum = 0;
	foreach $key (sort hashValueAscendingNum  (keys (%beginHash)))
	{
	#	print "$beginHash{$key}\t\t $key\n";
		$positionHash{$posNum} = $key;
		push @beginArray, $beginHash{$key};
		push @endArray, $endHash{$key};
		$posNum++;
	}
	
	
	# open a hash to store genes on contigs you want to remove and push unique genes to array
	# start with the end value of the first gene in the endArray and compare it to all other
	# genes beginning values to see if it overlaps with any genes. IF it does print both
	# names so you can remove them.move to the second genes end value and repeat process. Etc. 
	%generemoveHash = ();
	$numItemsInArrays = scalar (@endArray);
	print LOG "Working on $item\n";
	print  LOG "Number items in Array $numItemsInArrays\n";
	$start = 1;
	for (0 .. $numItemsInArrays)
	{
		$posEnd = $_; # starts at 0 (first pos in Array)
		print LOG "End of first seq in position $posEnd is at ";
		$end = $endArray[$posEnd]; #end equals 600 etc.
		print LOG "$end\n";
		for($start .. ($numItemsInArrays-1))
		{
			print LOG "Start is equal to $start\n";
			$posBegin = $_; #starts at the second item in the end array.
			print LOG "Start of second seq in position $start is at $beginArray[$posBegin]\n";
			print LOG "IF $end is > $beginArray[$posBegin]\n";
			if( $end > $beginArray[$posBegin])
			{		
 				if(! exists $generemoveHash{$positionHash{$posBegin}})
 				{
 					push @RemoveArray, $positionHash{$posBegin};
 					print LOG "Remove this region: $positionHash{$posBegin}\n";
 					$generemoveHash{$positionHash{$posBegin}} = 1;
 				}
 				if (! exists $generemoveHash{$positionHash{$posEnd}})
 				{
 					push @RemoveArray, $positionHash{$posEnd};
 					print LOG "Remove this region: $positionHash{$posEnd}\n";
 					$generemoveHash{$positionHash{$posEnd}} = 1;
 				}
 			}
		}
		$start++;
	}
	
	#print out each contigname and each gene to remove from that contig.
	foreach $position (@RemoveArray) 
	{
		print OUT "$item\t";
		print OUT "$position\n";
	}
	print LOG "\n";
}
