#!/usr/bin/perl
# =====================================================
# 1). This scripts creates a scaffold out of the mafft 
# alignments of othorgroup contigs (reference not included).
# It assumes that the reference sequence is always the last
# in the alignment - mafft fasta alignment format.
# 2). Computes the coverage of the reference by the 
# constituent contigs excluding the "N" and "-" characters.
# => coverage = bases covered / reference length
# 3). Computes the divergence based on mismatches in the 
# overlapping regions of the contigs
# 4). Gaps opened within contigs (excluding overhangs) are
# not padded with "Ns", but deleted. Gaps opened between 
# adjacent contigs are padded "Ns".
# => divergence = mismtaches / total overlap length   
# INPUT = Parent directory containing ortho directories
# OUTPUT = ortho scaffolds with coverage and divergence
# on the definition line 
# 
# Eric Wafula
# 08/107/10
# =====================================================

#USAGE: 10_create_mafft_scaffold_MG.fixed.mvref.pl <SUMMARY_FILE> <infile> <outfile> <path for temp files> <file of divergent scaffolds>
# run in directory labeld 1_Results above
# first put all gene folders in a folder called "genes" and make a dir called "output" for output

use strict;

# ============= Check for input parameters =========
if (!$ARGV[2]) {
    print "USAGE: 10_create_mafft_scaffold_MG.fixed.mvref.pl <SUMMARY_FILE> <infile> <outfile> <path for temp files> <file of divergent scaffolds>\n";
    exit(1);
}

my $summary_file= $ARGV[0]; #  summary output file
my $input= $ARGV[1];
my $outfile=$ARGV[2];
my $temp_path=$ARGV[3];
my $divergent=$ARGV[4];

#print "Input file is: $input\n";
#print "Outputfile is: $outfile\n";
#print "Temp files going to: $temp_path\n";
#print "Sammary information in: $summary_file\n";


#print "START - "; system "date"; #Based on prior runs, this takes about 1 second to run, we don't need all of this junk in the logs.

# ====== open mafft alignment file and get the alignements ======
open (IN, "<$input") or die "Couldn't open $input file, $!";
my @align = ();
while (my $line = <IN>){
	if ($line =~ /^>/){
		START:
		my $seq = "";
		while (<IN>){
			if ($_ =~ /^>/){last;}
			chomp($_);
			$seq .= $_;
		} # end while for inner <IN>
		if ($seq ne ""){
		$seq =~ s/\s+//g;
		push (@align, $seq);

		goto START;
		}
	}
} # end while for outer <IN>
close IN;

my $tempref = $align[0]; #MAG --move ref to last
shift(@align); #MAG-- remove the first (reference) sequence
push(@align, $tempref);


# ====== create an alignment base position matrix ========
open (MATRIX, ">$temp_path/temp") or die "Couldn't open temp file, $!";
foreach my $i (0..$#align-1){  
	my @elements = split(//,$align[$i]);
	foreach my $j (0..$#elements){
		if($j == $#elements){ print MATRIX "$elements[$j]\n";}
		else {print MATRIX "$elements[$j] ";}
	}
}
close MATRIX;
# ====== Transpose the alignment base position matrix =======
open (TMATRIX, ">$temp_path/temp2") or die "Couldn't open temp2 file, $!";
open (MATRIX, "<$temp_path/temp") or die "Couldn't open temp file, $!";
my $unequal=0;
$_=<MATRIX>;
$_ =~ s/\r?\n//;
my @out_rows = split (/\s+/, $_);
my $num_out_rows = $#out_rows+1;
while(<MATRIX>) {
	$_ =~ s/\r?\n//;
	my @F = split (/\s+/, $_);
	foreach my $k (0 .. $#F) {
		$out_rows[$k] .= " $F[$k]";
	}
	if ($num_out_rows != $#F+1) {
		$unequal=1;
		warn "\nWARNING! Rows in input had different numbers of columns\n" if $unequal; 
		warn "\nTransposed table: result has $. columns and $num_out_rows rows\n\n"; 
		exit(1);    
	}
}
foreach my $row (@out_rows) {
	print TMATRIX "$row\n"; 
}
close TMATRIX;
close MATRIX;
# ====== compute consesus sequence out alignment ======
open (TMATRIX, "<$temp_path/temp2") or die "Couldn't open temp2 file, $!";
my $consensus = "";
my $mismatches = 0;
my $overlap_len = 0;

while (<TMATRIX>){
	chomp($_);
	$_ =~ s/\s//g;
	my %hash2;
	my $count = 0;
	my $tracker = -1;
	my $tie = 0;
	my ($x, $c, $g, $t, $e) = (0, 0, 0, 0, 0);
	while ($_ =~ /a/ig){$x++;}
	while ($_ =~ /c/ig){$c++;}
	while ($_ =~ /g/ig){$g++;}
	while ($_ =~ /t/ig){$t++;}
	while ($_ =~ /[^acgt]/ig){$e++;}
	if (($x+$c+$g+$t) > 1){$overlap_len++;} 
	if ($e == $#align){$consensus .= "N"; next;}
	else {
		my %hash = ("a", $x, "c", $c, "g", $g, "t", $t);
		foreach my $keys (sort {$hash{$b} <=> $hash{$a}} keys %hash){
			if ($tracker == -1){
				$count++;
				$tracker = $hash{$keys};
				$hash2{$keys}=$count;
			}
			else{
				if ($hash{$keys} == $tracker){
				$hash2{$keys}=$count;
				$tie++;
				if ($hash{$keys} > 0 ) {$mismatches++;} #MAG new method
				}
				else{
					$count++;
					$count += $tie;
					$tie = 0;
					$hash2{$keys}=$count;
					$tracker = $hash{$keys};
					#$mismatches += $hash{$keys}; #original method 
					if ($hash{$keys} > 0 ) {$mismatches++;} #MAG new method
				}
			}
		}
		if ($hash2{'a'} < $hash2{'c'} and $hash2{'a'} < $hash2{'g'} and $hash2{'a'} < $hash2{'t'}){$consensus .= 'a';}
		if ($hash2{'c'} < $hash2{'a'} and $hash2{'c'} < $hash2{'g'} and $hash2{'c'} < $hash2{'t'}){$consensus .= 'c';}
		if ($hash2{'g'} < $hash2{'a'} and $hash2{'g'} < $hash2{'c'} and $hash2{'g'} < $hash2{'t'}){$consensus .= 'g';}
		if ($hash2{'t'} < $hash2{'a'} and $hash2{'t'} < $hash2{'c'} and $hash2{'t'} < $hash2{'g'}){$consensus .= 't';}
		if ($hash2{'a'} == $hash2{'c'} and $hash2{'a'} < $hash2{'g'} and $hash2{'a'} < $hash2{'t'}){$consensus .= 'M';}
		if ($hash2{'a'} == $hash2{'g'} and $hash2{'a'} < $hash2{'c'} and $hash2{'a'} < $hash2{'t'}){$consensus .= 'R';}
		if ($hash2{'a'} == $hash2{'t'} and $hash2{'a'} < $hash2{'c'} and $hash2{'a'} < $hash2{'g'}){$consensus .= 'W';}
		if ($hash2{'c'} == $hash2{'g'} and $hash2{'c'} < $hash2{'a'} and $hash2{'c'} < $hash2{'t'}){$consensus .= 'S';}
		if ($hash2{'c'} == $hash2{'t'} and $hash2{'c'} < $hash2{'a'} and $hash2{'c'} < $hash2{'g'}){$consensus .= 'Y';}
		if ($hash2{'g'} == $hash2{'t'} and $hash2{'g'} < $hash2{'a'} and $hash2{'g'} < $hash2{'c'}){$consensus .= 'K';}
		if ($hash2{'a'} == $hash2{'c'} and $hash2{'a'} == $hash2{'g'} and $hash2{'a'} < $hash2{'t'}){$consensus .= 'V';}
		if ($hash2{'a'} == $hash2{'c'} and $hash2{'a'} == $hash2{'t'} and $hash2{'a'} < $hash2{'g'}){$consensus .= 'H';}
		if ($hash2{'a'} == $hash2{'g'} and $hash2{'a'} == $hash2{'t'} and $hash2{'a'} < $hash2{'c'}){$consensus .= 'D';}
		if ($hash2{'c'} == $hash2{'g'} and $hash2{'c'} == $hash2{'t'} and $hash2{'c'} < $hash2{'a'}){$consensus .= 'B';}
		if ($hash2{'a'} == $hash2{'c'} and $hash2{'a'} == $hash2{'g'} and $hash2{'a'} == $hash2{'t'}){$consensus .= 'X';}
	}
} # end if while TMATRIX
close TMATRIX;

# ====== compute scaffold reference coverage and divergence  ======
# compute reference coverage by scaffold excluding Ns and dashes ("-")
my $coverage = 0;

my @scfld = split (//, $consensus);
my @ref = split (//, $align[$#align]); 
if ($#scfld == $#ref){
	foreach my $l (0..$#scfld){if($scfld[$l] ne "N" and $ref[$l] ne "-"){$coverage++;}} 
}
else {print "SCAFFOLD AND REFERENCE ALIGNMENTS LENGTHS NOT PROPORTIONAL!!!..\n $#scfld \n $#ref\n"; exit(1);}
my $reference = $align[$#align];
#		$reference =~ s/-//g;
#		my $ref_len = length($reference);
my $ref_len = 1;
$coverage = $coverage/$ref_len;
$coverage = sprintf("%.3f", $coverage);
# compute scaffold divergence - based on mismatches in overapping regions of contigs
my $divergence;
if ($overlap_len == 0){$divergence = 0;}
else {$divergence = $mismatches/$overlap_len;}
$divergence = sprintf("%.3f", $divergence);
# get scaffolds with diverenge threshold below 5%
open (FILE, ">>$summary_file") or die "Couldn't open $summary_file, $!";
open (DIVERGE, ">>$divergent") or die "Coudn't open $divergent, $!";

if ($divergence < 0.05){
		if ($overlap_len == 0){print FILE "$input\t$coverage\tNO\tYES\n";}
	else {print FILE "$input\t$coverage\t$divergence\tYES\n";}
	
	# # ============ open output fasta file ============
	open (OUT, ">>$outfile") or die "Couldn't open $outfile file, $!";

	
	print OUT ">$input"."Coverage="."$coverage Divergence "."$divergence\n";
	# remove within contig gaps
	my %gaps;
	foreach my $n (0..$#align-1){
		my $str = "$align[$n]";
		my @array = split (//, $str);
		for (my $p = 0; $p < $#array-1; $p++){if($array[$p] eq "-"){$array[$p] = "N";}else{last;}}
		for (my $q = $#array; $q > -1; $q--){if($array[$q] eq "-"){$array[$q] = "N";}else{last;}}		 
		foreach my $m (0..$#array){
			if ($array[$m] eq "-"){$gaps{$m}=$m;}
		} 
	}
	foreach my $r (keys %gaps){
		my %overhang;
		foreach my $t (0..$#align-1){
			my $str2 = "$align[$t]";
			my @array2 = split (//, $str2);
			foreach my $s (0..$#array2){if($array2[$s] ne "-"){$overhang{$s}=$s;}}
		}
		if (!$overhang{$r}){delete ($scfld[$gaps{$r}]);}
	}
	$consensus = join('', @scfld);
	# trim scaffold leading and trailing Ns
#			$consensus =~ s/^N+//;
	$consensus =~ s/N+$//;
	$consensus =~ s/X/N/g;
	$consensus = uc ($consensus);
	my $len += length($consensus);
	$consensus =~ s/.{80}(?=.)/$&\n/g;
	print OUT "$consensus\n";
}
else{
	print FILE "$input\t$coverage\t$divergence\tNO\n";
	print   "Too divergent ($divergence), consensus not made for $input\n";
	print DIVERGE "$divergence \t $input\n";

}
close FILE;
close OUT;
system ("rm $temp_path/temp*");
#print "STOP - "; system "date";
exit(0)
