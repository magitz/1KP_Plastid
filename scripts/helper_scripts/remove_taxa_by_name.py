#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

#This only works on fasta files. Don't try to do this on phylip as each sequence is written one by one and you get phylip header for each.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file")
parser.add_argument("-r", help="list of names to remove")
parser.add_argument("-o", help="output file")

args = parser.parse_args()

infile = args.i
remove_list = args.r
outfile = args.o

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	REMOVE=open(remove_list, 'r')
except IOError:
	print "Can't open file", remove_list

try:
	OUT=open(outfile, 'w')
except IOError:
		print "Can't open file", outfile
	
taxa_to_remove = []

for Line in REMOVE :
	Line = Line.strip('\n')
	taxa_to_remove.append(Line)	#Add the whole line, which should be the name of the taxon to remove
	
for record in SeqIO.parse(IN, "fasta") :
	
	if record.id in taxa_to_remove :
		print "Removing %s from %s" %(record.id, infile)
	else :
		OUT.write(record.format("fasta"))	
	
