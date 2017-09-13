#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re


parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file")
parser.add_argument("-r", help="list of 4-letter codes to remove")
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
	OUT=open(outfile, 'a')
except IOError:
		print "Can't open file", outfile
	
taxa_to_remove = []

for Line in REMOVE :
	Line = Line.strip('\n')
	taxa_to_remove.append(Line[:4])	#1st 4 characters of the line should be 1KP 4-letter code.
	
for record in SeqIO.parse(IN, "fasta") :
	tax_code=record.id[-4:]
	if tax_code in taxa_to_remove :
		print "Removing %s from %s" %(record.id, infile)
	else :
		OUT.write(record.format("fasta"))	
	
