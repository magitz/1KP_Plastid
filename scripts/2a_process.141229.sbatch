#!/bin/bash
#SBATCH --account=ufhpc
#SBATCH --qos=ufhpc
#SBATCH --mail-user=magitz@ufl.edu
#SBATCH --mail-type=FAIL,END
##SBATCH --output=$CLADE/Parse_out.%j.o
##SBATCH --error=$CLADE/Parse_error.%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=64:00:00



############################################################
# qsub scripts/2a_process.qsub -v StartStep=step_to_start_at,CLADE=clade_folder_to_process,blast_out=folder with blast results,[my_dir=set my_dir from commandline, my_final_dir=set my_final_dir from commandline]

#This worked on HPG2. Can't get output/error names set
#sbatch --export=CLADE=Caryophyllales,blast_out=Caryophyllales/Blastx/,mydir=Caryophyllales/cpDNA_June2016/Parse.e25.160623,StartStep=9 --output=Caryophyllales/Parse_out.%j.o scripts/2a_process.141229.sbatch 


# This script takes the blast output from step 1 (blast of all 1kp scaffolds to each gene reference file)
# It does all of the post blast processing using several scripts and applications
#
################
#
# 6/27/13- MAG
# * Added a second blast step before combining pulled scaffold fragments. We found that there were often odd fragments
# 	that were being included. All fragments have been selected based on blast similarity to the gene of interest, using
#	a blast dataset that consists of just the reference taxa from the clade of interest.
#	Now we reblast each fragment to the reference dataset for the gene for Brads whole dataset.
#	If the best blast hit is outside the clade of interest, throw out that fragment.
#
# * Added StartStep variable to allow restarting partial analyses at desired step.
#
# * Also redid paths/variable so that scripts can be run from a single copy in base_dir/scripts and analyze each clade in
#	base_dir/clade. Before this, I had copies in each directory that needed updating each time a change was made.
#
# 07/16/13: MAG
# * Change method of calculating consensus of scaffolds to get single sequence for each sample.
#
# 07/19/13: MAG
# *  Modified script blast_to_full_dataset.v2.py, adding output of summary stats for good and bad blast hits.
#
# 05/08/14: MAG
# * Modified to handle multiple independent sets of reference data (cpDNA, mtDNA, etc.).
#
# 12/29/14: MAG
# * Got rid of the prerequisite of having a list of 4-letter codes for all the samples in the clade. This list is now
#   automatically generated based on the samples in assembly directory.
#
# 1/5/15: MAG
# * Added option to stop after step 7. This allows rest of steps to be done more in parallel for faster processing.

###########################################################

###########################################################
# Set some manually configured values here
###########################################################

#These values are good for the first time through, but best to use next section for reruns.
#Make sure to change BOTH of the values below. This saves a bunch of parsing, so just change the 2 to correspond.
evalue=e25		#used for directory naming
full_evalue=1.0e-25 #used for parsing.

# Path to full set of references used in blast step at Step 7
# For the first round, this should be your set of genes taken from GenBank for example. For subsequent rounds, it will
# include newly found genes from your data.

# Path to the reference dataset to be used. Should have .faa files here and .fna files in nucleotide/
# For cpDNA, the nucleotide files should have also been converted to blast databases.
#full_dataset_ref_path=/scratch/lfs/magitz/1kp/Summer2013/Brads_dataset/wOrders
full_dataset_ref_path=/ufrc/ufhpc/magitz/research/1kp/Summer2013/cpDNA_round1_problems_removed/wBrads  #For second round, has first round results and Brads.
#full_dataset_ref_path=/scratch/lfs/magitz/1kp/Summer2013/mtDNA/no_gap

#Path to sub-alignements for aligning fragments. For cpDNA, this is per clade alignments. Used in steps 7-9.
sub_alignments=$CLADE/AA_RefGeneFiles

# Name for output directory: will be $CLADE/out_name/
out_name=cpDNA_June2016

#Concatenated file of all the scaffold stats:
stats_file=All.scaffold-stats.141229.txt

#Assembly path
Assemblies=$CLADE/Assemblies

#Stop after step 7?
Step7Stop=0
###########################################################

#Print the date and set a date stamp
date
date_stamp=`date "+%y%m%d"`
#Create output directory.
echo "Results of this run will be in $CLADE/$out_name"
mkdir $CLADE/$out_name

###########################################################
# Handle required command line options here
###########################################################

#Need to pass clade name to work on
if [[ -z $CLADE ]]
then
echo "You must specify the clade name with CLADE= option"
exit -1
fi

#Need to pass the location of the blast output to work on (from step 1b: BLAST_PATH="mtDNA/Blast/blasttable/")
if [[ -z $blast_out ]]
then
echo "You must specify the location of blast results with blast_out= option"
exit -1
fi
###########################################################


###########################################################
# Handle optional command line options
###########################################################

#if not passed in qsub command, which is useful for reprocessing partial results, set my_dir and my_final_dir
if [[ -z $my_dir ]]
then
my_dir=$CLADE/$out_name/Parse.$evalue.$date_stamp
fi

if [[ -z $my_final_dir ]]
then
my_final_dir=$CLADE/$out_name/Final_scaffolds.$evalue.$date_stamp
fi
###########################################################




#Load needed modules for things called directly from here.
module load perl
module load mafft
module load python

# ncbi_blast and wise modules are loaded in step 7 if run.


###########################################################
# Step 2: Create directories, run 2b_process.140512.pl
###########################################################
if [[ $StartStep -le 2 ]]
then
#create output directories
mkdir $my_dir
mkdir $my_dir/2out
mkdir $my_dir/2error
mkdir $my_dir/refs

for i in $full_dataset_ref_path/*.faa
do
subject_base=`basename $i .faa`
mkdir $my_dir/refs/$subject_base
done

echo Directories made, starting step 2b
date
perl scripts/2b_process.140512.pl $full_evalue $CLADE $my_dir $full_dataset_ref_path $blast_out
if [[ $? -ne 0 ]]					#Exit if step 2 returns error code.
then
echo "Step 2, 2b_process.140512.pl, failed (code=$?). Stopping."
exit -1
fi

fi


###########################################################
# Step 3 (and 4, 5): Find duplicates and prune them
###########################################################

if [[ $StartStep -le 3 ]]
then
echo Starting steps 3-5
date
perl scripts/3_finddups.pl $my_dir
if [[ $? -ne 0 ]]					#Exit if step 3 returns error code.
then
echo "Step 3, 3_finddups.pl, failed (code=$?). Stopping."
exit -1
fi

if [[ ! -s $my_dir/RepeatSeqs_GB.txt ]] #Exit if output of step 3 is zero length
then
echo "Step 3, 3_finddups.pl, resulted in zero length file at $my_dir/RepeatSeqs_GB.txt. Stopping."
exit -1
fi


echo Starting step 4
date
perl scripts/4_ListDuplicates2Remove.pl $my_dir #Step 4
if [[ $? -ne 0 ]]					#Exit if step 4 returns error code.
then
echo "Step 4, 4_ListDuplicates2Remove.pl, failed (code=$?). Stopping."
exit -1
fi

echo Starting step 5
date
perl scripts/5_prunerepeats_BR.pl $my_dir #Step 5
if [[ $? -ne 0 ]]					#Exit if step 5 returns error code.
then
echo "Step 5, 5_prunerepeats_BR.pl, failed (code=$?). Stopping."
exit -1
fi

fi


###########################################################
# Step 6: Get the sequences from matched scaffolds
###########################################################

if [[ $StartStep -le 6 ]]
then
echo Starting step 6
date
perl scripts/6_getAssemblySeqs.pl $my_dir $CLADE
if [[ $? -ne 0 ]]					#Exit if step 6 returns error code.
then
echo "Step 6, 6_getAssemblySeqs.pl, failed (code=$?). Stopping."
exit -1
fi
fi


###########################################################
# Step 7: 6/27/13 Add blast_to_full_dataset.py
# Step 7 blasts pulled data to full dataset and parses to
# remove fragments with distant blast hits. Then runs Genewise on keepers.
###########################################################

if [[ $StartStep -le 7 ]]
then
echo Starting step 7
date
#Need ncbi_blast for blast script.
module load ncbi_blast

#Set the velues of $full_dataset_ref_path and $stats_file above in manually configured options section.
python scripts/blast_to_full_dataset.v2.py -d $my_dir -c $CLADE -b $full_dataset_ref_path/nucleotide -s $stats_file

if [[ $? -ne 0 ]]					#Exit if step 7 returns error code.
then
echo "Step 7, blast_to_full_dataset.v2.py, failed (code=$?). Stopping."
exit -1
fi

#Need GeneWise for step 7.
module load wise
perl scripts/7_SeparateGenewise.new.pl $my_dir $CLADE $full_dataset_ref_path
if [[ $? -ne 0 ]]					#Exit if step 7 returns error code.
then
echo "Step 7, 7_SeparateGenewise.new.pl, failed (code=$?). Stopping."
exit -1
fi

fi

if [[ Step7Stop -eq 0 ]]
then
#Prep directories for last steps.
mkdir $my_dir/mafft
mkdir $my_final_dir


#There is no step 8....
for i in $sub_alignments/*.faa
do
GENE=`basename $i .faa`

###########################################################
# Step 9: Run mafft to align pulled data to references
# Then compute the consensus to get one sequence for each gene for a taxon.
###########################################################
if [[ $StartStep -le 9 ]]
then
if [[ -e $my_final_dir/$GENE ]] # we've already processed this gene
then
echo "Skipping $gene, already done."
else
echo  -n "Working on gene $GENE: " #Print gene to log without newline (-n)

mkdir $my_dir/mafft/$GENE
mkdir $my_final_dir/$GENE

if [[ ! -s $sub_alignments/nucleotide/$GENE.fna ]] #Warn if reference sequence file zero length
then
echo "No reference sequences found at: $sub_alignments/nucleotide/$GENE.fna. Stopping."
#exit -1
fi

#Count how many sequences are in the reference allignment for the gene--used in step 9.
ref_seqs=`grep -c \>  $sub_alignments/nucleotide/$GENE.fna`

for i in $my_dir/Genewise/$GENE/*.out # For each taxon file in the $GENE directory, run mafft to get alignement
do
if [ -e "$i" ]  # if there's nothing there, i will equal "$my_dir/Genewise/$GENE/*.out" and useless files are made.
then
code_genename=`basename $i .out`

echo -n "." #add some dots to same line as GENE above for each mafft run.
#Uses the nucleotide alignment file generated by prep_partitions_for_blast.sh
mafft-linsi --thread $SLURM_CPUS_PER_TASK --add $i  $sub_alignments/nucleotide/$GENE.fna > $my_dir/mafft/$GENE/$code_genename.mafft.out

#####
# MAG: 07/16/13: This used a script that was written by Eric Wafula to make the consensus.
# I found one bug in the script, and thought it was working. However there still seem to be
# some odd results. For example the divergence is calculated over only the overlapping portions
# so if 2 scaffolds overlap by 3bp and differ, they are too divergent and won't be combined.
# I also found lots of problems with Claudia's datasets where there were dozens of overalapping
# scaffolds, this calculated >100% divergence. Something is not right here, so let's use BioPerl.
##Run step 10, which takes the alignments of scaffolds to the single reference and creates the consensus sequence for the taxon.
##USAGE: 10_create_mafft_scaffold_MG.fixed.mvref.pl <SUMMARY_FILE> <infile> <outfile> <path for temp files> <file of divergent scaffolds>
#echo "Command is:	perl scripts/10_create_mafft_scaffold_MG.fixed.mvref.pl $my_final_dir/scaffold_summary.out $my_dir/mafft/$GENE/$code_genename.mafft.singleref.out $my_final_dir/$GENE/$code_genename.$GENE.fa $TMPDIR $my_final_dir/Divegent_scaffold_summary.txt "
#
#perl scripts/10_create_mafft_scaffold_MG.fixed.mvref.pl $my_final_dir/scaffold_summary.out $my_dir/mafft/$GENE/$code_genename.mafft.singleref.out $my_final_dir/$GENE/$code_genename.$GENE.fa $TMPDIR $my_final_dir/Divegent_scaffold_summary.txt

# Modify step 9 to remove all refercences, leaving only the aligned sample data: add 1 to $ref_seqs.
#Run step 9 to remove all but one of the reference sequences to prep for step 10.
# Nothing called Step 10 anymore.
perl scripts/Remove_1st_n_seqs.pl $ref_seqs $my_dir/mafft/$GENE/$code_genename.mafft.out $my_dir/mafft/$GENE/$code_genename.mafft.noref.out
if [[ $? -ne 0 ]]					#Exit if step 9 returns error code.
then
echo "Step 9, Remove_1st_n_seqs.pl, failed (code=$?). Stopping."
exit -1
fi

perl scripts/compute_consensus_from_mafft.pl $my_dir/mafft/$GENE/$code_genename.mafft.noref.out $my_final_dir/$GENE/$code_genename.$GENE.fa
if [[ $? -ne 0 ]]					#Exit if step 9 returns error code.
then
echo "Step 9, compute_consensus_from_mafft.pl failed (code=$?) while working on $GENE for $code_genename. Stopping."
exit -1
fi

else
echo "No data for $GENE, removing directories"
rm -r $my_dir/mafft/$GENE
rm -r $my_final_dir/$GENE
fi
done
echo "Done with $GENE."
fi
fi

###########################################################
# Step 11: Finalizing all the output.
# Rarely a need to start here, so just keeping it in the same loop at stop 9 for convenience.
###########################################################

if [[ $StartStep -le 11 ]]
then
#Create the clade_code_list based on the samples in the Assemblies folder.
#This is a list of all the 4-letter codes of samples in this clade.
#Modify path of assemblies if they are not in $CLADE/Assemblies.
ls $Assemblies | cut -f1 -d"-" > $CLADE/clade_code_list.$date_stamp
#Link this to the clade_code_list file
cd $CLADE
rm clade_code_list
ln -s clade_code_list.$date_stamp clade_code_list
cd ..

#For each gene, run the finalize script.
echo Running finalize on gene $GENE
python scripts/finalize.py -i 1kp_taxon_list_full.150101.txt  -g $my_final_dir/$GENE -f $CLADE/clade_code_list
if [[ $? -ne 0 ]]					#Exit if step 9 returns error code.
then
echo "Step 11, finalize.py , failed (code=$?). Stopping."
exit -1
fi

mafft-linsi --thread $SLURM_CPUS_PER_TASK --add $my_final_dir/results/$GENE.fna  $sub_alignments/nucleotide/$GENE.fna > $my_final_dir/results/$GENE.mafft.fna
fi

done
fi

echo Done.


date
