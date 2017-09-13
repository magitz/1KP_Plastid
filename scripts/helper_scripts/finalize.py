#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import os, errno

# This script goes through the output of 10_create_mafft_scaffold_MG.fixed.mvref.pl and cleans things up.
#	Combines files into one file per gene
#	Aligns the new sequences to the profile alignment
#	Creates a table that has the number of non-N bases for each gene for each taxon


#Parse commandline options.

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="1kp list of taxa file (all 1,400+ taxa)")
parser.add_argument("-f", help="Full list of species that were processed for this group--4 letter codes as used in prep steps.")
parser.add_argument("-g", help="Gene directory to parse")
args = parser.parse_args()

kp_list_file = args.i
full_taxon_set=args.f
gene_dir=args.g

try:
	KP_LIST = open (kp_list_file, 'r')
except IOError:
	print "Can't open file:", kp_list_file

taxon_hash={}

for Line in KP_LIST :
	Line = Line.strip('\n')
	tax_code=Line[:4]
	taxon=Line[5:] + "_" + tax_code
	taxon_hash[tax_code] = taxon


gene_name= os.path.basename(gene_dir) #get gene name from last part of path.

result_path= os.path.dirname(gene_dir) + "/results/"  #Get the directory one up from the passed directory and add /results/

gene_dir = gene_dir + "/*" ##E.g.: Saxifragales/Final_scaffolds.e25.130628/ycf4
file_list=glob.glob(gene_dir)

try:
	os.makedirs(result_path)
except OSError, err:
    # Reraise the error unless it's about an already existing directory 
    if err.errno != errno.EEXIST or not os.path.isdir(result_path): 
        raise

gene_out_name= result_path + gene_name + ".fna"
gene_out_aa_name=result_path  + gene_name + ".faa"

try:
	GENE_OUT=open(gene_out_name, 'w')
except IOError:
	print "Can't open file:", gene_out_name

try:
	GENE_AA_OUT= open(gene_out_aa_name, 'w')
except IOError:
	print "Can't open file:", gene_out_aa_name
	
try:
	SUM_FILE=open(result_path + 'Final.lengths.summary.txt', 'a')
except IOError:
	print "Can't open file:", result_path, "Final.lengths.summary.txt"

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
	
	
for file in file_list :
	print file
	GENE_IN= open(file, 'rU')
	for record in SeqIO.parse(GENE_IN, "fasta") :
		identifier= record.id		#Get the name as it is e.g: >Saxifragales/Parse/mafft/accD/ADHK.genewise.mafft.singleref.outCoverage=339.000 Divergence 0.000
		#id_bits=re.split("/",identifier) #Split by /s
		#code=id_bits[4][:4] #First 4 chars of id_bits[3] is 4-letter code.
		
		code=os.path.basename(record.id)[:4]  #First 4 chars of last part of name is 4-letter code.
		print record.id
		print code
		
		#Write fasta file to GENE_OUT with format of Genus_spp_4letter_gene
		record.id=taxon_hash[code]
		record.description='' #Clear the description field
		SeqIO.write(record, GENE_OUT, "fasta")
				
		#Generate translation and look at length to 1st stop codon
		#MAG 8/6/13 add the ungap to remove gaps as codons that were all gaps were causing problems. Won't handle out of frame gaps, but it's a start.
		protein=record.seq.ungap("-").translate(to_stop=True)
		sum_aastoplen_hash[taxon_hash[code]]=len(protein)
		
		#Then get full length translation and output to file
		#MAG 8/6/13 add the ungap to remove gaps as codons that were all gaps were causing problems. Won't handle out of frame gaps, but it's a start.
		protein=record.seq.ungap("-").translate()
		sum_aalen_hash[taxon_hash[code]]=len(protein)
		
		#write out amino acid file 		
		out_name=">"+taxon_hash[code]
		GENE_AA_OUT.write (out_name)
		GENE_AA_OUT.write ("\n")
		GENE_AA_OUT.write (str(protein))
		GENE_AA_OUT.write ("\n")
		
		#Get some summary stats for the sequences.
		no_Ns=str(record.seq)
		no_Ns=no_Ns.replace("N","")
		
		sum_ntlen_hash[taxon_hash[code]]=len(no_Ns)

#Print summary information to file.		
#  Gene	Statistic	Taxon_1	Taxon_2	Taxon_3	etc...

SUM_FILE.write ('%s \t Length nt \t' %(gene_name))
for key in sorted(sum_ntlen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_ntlen_hash[key]))
	
SUM_FILE.write ('\n%s \t Length aa \t' %(gene_name))
for key in sorted(sum_aalen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_aalen_hash[key]))
	
SUM_FILE.write ('\n%s \t Length aa to 1st stop \t' %(gene_name))
for key in sorted(sum_aastoplen_hash.iterkeys()):
	SUM_FILE.write ('%s \t' %(sum_aastoplen_hash[key]))
	
SUM_FILE.write ('\n')


		

		
		 
		