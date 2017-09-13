#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import subprocess
import os
import fnmatch
import shutil

#  This script goes through the result directories (FinalReslts.e25.YYMMDD) and combine from all clades 

#Parse commandline options.

parser = argparse.ArgumentParser()
parser.add_argument("-e", help="evalue", default="e25")
parser.add_argument("-o", help="Output folder")
parser.add_argument("-r", help="result path within each clade folder")

args = parser.parse_args()

evalue = args.e
out_dir = args.o
results = args.r

clade_dirs= [ 'Algae',
				'Asparagales_Liliales',
				'Basal_Angios',
				'Basal_Eudicots',
				'Campanulids',
				'Caryophyllales',
				'Conifers',
				'Cornales_Ericales',
				'Fabids',
				'Lamiids',
				'Malvids',
				'Monocots',
				'Non_seed',
				'Saxifragales' ]

for clade in clade_dirs :
	print "Working on clade: %s" %(clade)
	
	final_glob= clade + '/' + results + '/Final_scaffolds.' + evalue + '*' 	#Glob all of the Final_scaffold directories.
	final_dirs=glob.glob(final_glob)
	
	if len(final_dirs) > 0:
		#Sort and take the last one, should be the latest.
		final_dirs.sort()
		result_dir= final_dirs[-1] + "/results/"
	
		result_glob=result_dir + '*.fna'
	
		summary_file=result_dir + 'Final.lengths.summary.txt'		#Copy the summaries to dest dir.
		dest_summary= out_dir + '/' + clade + '_Final.lengths.summary.txt'
		shutil.copyfile(summary_file, dest_summary)
	
		for gene in glob.glob(result_glob) :
			if not fnmatch.fnmatch(gene, '*.mafft.fna') :   #Skip the mafft alignments which also in in .fna
				try:
					GENE=open(gene, 'r')
				except IOError:
					print "Can't open file", gene
			
			
				gene_name= os.path.basename(gene)
				gene_name= re.sub('.fna', '', gene_name)
			
				out_file= out_dir + '/' + gene_name + '.fna'  #Handle the nucleotide files (*.fna)
			
				try:
					OUT=open(out_file, 'a')
				except IOError:
					print "Can't open file", out_file
			
				for record in SeqIO.parse(GENE, "fasta") :
					record.id = '>' + clade + '_' + record.id # Add clade name to the record id.

					seq_sub=re.sub('^N+', '', str(record.seq)) #remove any Ns at the start of the sequence
					record.seq=seq_sub
			
					if len(record.seq) > 0 :	#Only write out sequences that are non-zero in length.
						OUT.write(record.id)	#Write the record to the gene file
						OUT.write ('\n')
						OUT.write(record.seq)	
						OUT.write ('\n')
			

				aa_file=re.sub('.fna', '.faa', gene)	#Handle the amino acid files (*.faa)
			
				try:
					AA_FILE=open(aa_file, 'r')
				except IOError:
					print "Can't open file", aa_file
			
	
				out_file= out_dir + '/' + gene_name + '.faa'
			
				try:
					OUT=open(out_file, 'a')
				except IOError:
					print "Can't open file", out_file
			
				for record in SeqIO.parse(AA_FILE, "fasta") :
					record.id = '>' + clade + '_' + record.id # Add clade name to the record id.

					seq_sub=re.sub('^X+', '', str(record.seq)) #remove any Xs at the start of the sequence
					record.seq=seq_sub
			
					if len(record.seq) > 0 :	#Only write out sequences that are non-zero in length.
						OUT.write(record.id)	#Write the record to the gene file
						OUT.write ('\n')
						OUT.write(record.seq)	
						OUT.write ('\n')
	else:
		print "No data in %s, skipping." %(final_glob)
		