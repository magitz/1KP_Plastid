#!/bin/sh

#Submits the blast jobs for a clade:
# sh scripts/1a.do_blast.140814.sh Saxifragales

#Get clade passed on commandline
clade=$1

##################################################################
# Variables that likely need to be changed for each specific case
#
# Set reference gene dataset location here needs tailing / --It's used in ${REF}*${REF_ENDING} in 1b.
REF=mtDNA.original/
# Set ending name for reference files, eg. .fna, .fasta, etc.
REF_ENDING=".faa"
# Set ouptut path, where you want the results to go
OUT_PATH=$clade/mtDNA_blastx/
#
#Where are the assemblies located within the clade directory?
assembly_path="Assemblies"
##################################################################

for i in $clade/$assembly_path/*.fa
do 
	queryname=`basename $i -SOAPdenovo-Trans-assembly.fa`
	echo Submitting $queryname with: qsub scripts/1b_gene_blast_AA.qsub -v CLADE=$clade,QUERY=$i,REF=$REF,REF_ENDING=$REF_ENDING,OUT_PATH=$OUT_PATH -N $queryname.$clade
	qsub scripts/1b_gene_blast_AA.qsub -v CLADE=$clade,QUERY=$i,REF=$REF,REF_ENDING=$REF_ENDING,OUT_PATH=$OUT_PATH -N $queryname.$clade
	sleep 0.5
done