#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import os

# This script goes through the output of 10_create_mafft_scaffold_MG.fixed.mvref.pl and cleans things up.
#	Combines files into one file per gene
#	Aligns the new sequences to the profile alignment
#	Creates a table that has the number of non-N bases for each gene for each taxon

# It also attempts to recover the longest sequence from samples that had something pulled, but it was lost along the way
# 	Sequences tend to be lost when there are multiple scaffolds that overlap, but are too divergent. In most caseso

#Parse commandline options.

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="1kp list of taxa file (all 1,400+ taxa)")
parser.add_argument("-f", help="Full list of species that were processed for this group--4 letter codes as used in prep steps.")
parser.add_argument("-g", help="Gene directory to parse")
parser.add_argument("-b", help="Backup directory with seqeunces to get if main set fails for a taxon--Parse_e##/cpPulled/")

args = parser.parse_args()

kp_list_file = args.i
full_taxon_set=args.f
gene_dir=args.g
cpPulled_dir=args.b

KP_LIST = open (kp_list_file, 'r')

taxon_hash={}

for Line in KP_LIST :
	Line = Line.strip('\n')
	tax_code=Line[:4]
	taxon=Line
	taxon_hash[tax_code] = taxon

gene_bits=re.split("/",gene_dir) #eg. Final_scaffolds/accD/

try:
	os.makedirs(gene_bits[0] + "/results/")
except OSError as exception:
	pass

gene_out_name= gene_bits[0] + "/results/" + gene_bits[1] + ".fna"
gene_out_aa_name=gene_bits[0] + "/results/" + gene_bits[1] + ".faa"

GENE_OUT=open(gene_out_name, 'w')
GENE_AA_OUT= open(gene_out_aa_name, 'w')
SUM_FILE=open(gene_bits[0] + '/results/Final.lengths.summary.txt', 'a')


sum_ntlen_hash={}
sum_aalen_hash={}
sum_aastoplen_hash={}

ALL_TAXA=open(full_taxon_set, 'r') #Initialize the summary hashes with all taxa so that all will be reported in the summary even if there wasn't data for that taxon.
for Line in ALL_TAXA:
	Line = Line.strip('\n')
	if len(Line) == 4 :
		sum_ntlen_hash[taxon_hash[Line]]=0
		sum_aalen_hash[taxon_hash[Line]]=0
		sum_aastoplen_hash[taxon_hash[Line]]=0
	elif Line != '':
		print "Error: Line in Full taxon list (-f) is %s" %(Line)
		
#print a header line in the summary file.
SUM_FILE.write ('Gene\tStatistic\t')
for key in sorted(sum_ntlen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(key))
SUM_FILE.write ('\n')
ALL_TAXA.cose()

ALL_TAXA=open(full_taxon_set, 'r') # Go back through the list of all the taxa that we should have data for.	
for LINE in ALL_TAXA :
	Line = Line.strip('\n')
	if len(Line) == 4 :
		file = gene_dir + Line + "*"
		try:
			GENE_IN= open(file, 'rU')
			#We got some data for this taxon/gene pair, so parse it.				
			for record in SeqIO.parse(GENE_IN, "fasta") :
				identifier= record.id		#Get the name as it is e.g: >Parse/mafft/accD/ADHK.genewise.mafft.singleref.outCoverage=339.000 Divergence 0.000
				id_bits=re.split("/",identifier) #Split by /s
				code=id_bits[3][:4] #First 4 chars of id_bits[3] is 4-letter code.
		
				#Write fasta file to GENE_OUT with format of Genus_spp_4letter_gene
				out_name=taxon_hash[code]+ "_"+ code
		
				record.id=out_name
				record.description='' #Clear the description field
				SeqIO.write(record, GENE_OUT, "fasta")
				
				#Generate translation and look at length to 1st stop codon
				protein=record.seq.translate(to_stop=True)
				sum_aastoplen_hash[taxon_hash[code]]=len(protein)
		
				#Then get full length translation and output to file
				protein=record.seq.translate()
				sum_aalen_hash[taxon_hash[code]]=len(protein)
		
				#write out amino acid file 		
				out_name=">"+out_name
				GENE_AA_OUT.write (out_name)
				GENE_AA_OUT.write ("\n")
				GENE_AA_OUT.write (str(protein))
				GENE_AA_OUT.write ("\n")
		
				#Get some summary stats for the sequences.
				no_Ns=str(record.seq)
				no_Ns=no_Ns.replace("N","")
		
				sum_ntlen_hash[taxon_hash[code]]=len(no_Ns)
		except:
		#file doesn't exist, which means we have no data for this species. Find out why.
			try:
				get_cpPulled=cpPulled_dir + gene_bits[1] + "*"
				GENE_IN=open(get_cpPulled, 'r')
				for record in SeqIO.parse(GENE_IN, "fasta") :
					identifier= record.id		#Get the name as it is e.g: >HQRJ_scaffold-HQRJ-2109559-Loropetalum_chinense
					code=identifier[:4] #First 4 chars is 4-letter code.
					if code == Line : #We have a sequence of the current taxon/gene combo
						

			
#Print summary information to file.		
#  Gene	Statistic	Taxon_1	Taxon_2	Taxon_3	etc...

SUM_FILE.write ('%s \t Length nt \t' %(gene_bits[1]))
for key in sorted(sum_ntlen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_ntlen_hash[key]))
	
SUM_FILE.write ('\n%s \t Length aa \t' %(gene_bits[1]))
for key in sorted(sum_aalen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_aalen_hash[key]))
	
SUM_FILE.write ('\n%s \t Length aa to 1st stop \t' %(gene_bits[1]))
for key in sorted(sum_aastoplen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_aastoplen_hash[key]))
	
SUM_FILE.write ('\n')


		

		
		 
		