#!/usr/bin/perl
# =====================================================
#  Takes input of aligned contigs and makes consensus
#  Uses BioPerl consensus_iupac method code stolen from:
#  https://github.com/josephhughes/Sequence-manipulation/blob/master/Consensus.pl
#
# Matt Gitzendanner
# 7/10/13
# =====================================================

#USAGE: USAGE: compute_consensus_from_mafft.pl <infile> <outfile>

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;



# ============= Check for input parameters =========
if (!$ARGV[1]) {
    print "USAGE: compute_consensus_from_mafft.pl <infile> <outfile>\n";
    exit(1);
}

my $input= $ARGV[0];
my $outfile=$ARGV[1];

#print "Input file is: $input\n";
#print "Outputfile is: $outfile\n";

my $filename=$1."_consensus" if $input=~/(.+)\..+/;
my $in = Bio::AlignIO->new(-file => "$input" ,
                         -format => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');
my $ali = $in->next_aln();
my $iupac_consensus = $ali->consensus_iupac(); # dna/rna alignments only
my $seq = Bio::Seq->new(-seq => "$iupac_consensus",
                        -display_id => $filename);

$out->write_seq($seq);