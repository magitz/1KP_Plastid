#!/bin/bash

##Runs steps 9-11 for the 2a_process scripts in parallel with each gene in each clade as a separate job.

#ref_path=$1		#The reference file path.--No use CLADE/AA_RefGeneFiles
result_dir=$1 	#Not including CLADE: cpDNA_Dec2014_results

my_parse="Parse.e25.150219"
my_final="Final_scaffolds.e25.150227"

date_stamp=`date "+%y%m%d"`

#clade_list="Algae Asparagales_Liliales Basal_Angios Basal_Eudicots Campanulids Caryophyllales Conifers Cornales_Ericales Fabids Lamiids Malvids Monocots Non_seed Saxifragales"
clade_list="Saxifragales"


echo "Processing all the genes in $ref_path"

for clade in $clade_list
do
	ref_path=$clade/AA_RefGeneFiles
	
	my_dir=$clade/$result_dir/$my_parse
	my_final_dir=$clade/$result_dir/$my_final
	
	for i in $ref_path/*.faa
	do
		gene=`basename $i .faa`
		num_results=`ls $my_dir/Genewise/$gene/ | wc -l`
		if [ $num_results -gt 0 ]  #Are there results to parse? 
		then
			echo "Submitting for $clade and $gene with: qsub scripts/9b_runLastStepsInParallel.sh -N $clade.$gene -v CLADE=$clade,GENE=$gene,ref_path=$ref_path,my_dir=${my_dir},my_final_dir=${my_final_dir}"	

			qsub scripts/9b_runLastStepsInParallel.sh -N $clade.$gene -v CLADE=$clade,GENE=$gene,ref_path=$ref_path,my_dir=${my_dir},my_final_dir=${my_final_dir}		
			sleep 0.1
		else #skip if not
			echo "No results for $gene in $clade, skipping"
			echo "No results for $gene in $clade, skipping" >> skipped_genes.$date_stamp.log
		fi
	done
done
