#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

#sorts two fasta files in the same (alphabetical) order to run through pal2nal.

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input aa file")
parser.add_argument("-n", help="input nt file")
parser.add_argument("-o", help="output aa file")
parser.add_argument("-f", help="output nt file")

args = parser.parse_args()

infile = args.i
outfile = args.o

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	NT=open(args.n, 'r')
except IOError:
	print "Can't open file", args.n

try:
	OUTAA=open(outfile, 'w')
except IOError:
		print "Can't open file", outfile

try:
	OUTNT=open(args.f, 'w')
except IOError:
		print "Can't open file", args.f


aa_seqs= SeqIO.to_dict(SeqIO.parse(IN, "fasta"))
nt_seqs= SeqIO.to_dict(SeqIO.parse(NT, "fasta"))

for key in sorted(aa_seqs):

	#Some AA seq names have _1, _2, _3, from the frame used in the translation. That's not in the nt name, so need to remove that.
	frame=["_1", "_2", "_3"]
	if aa_seqs[key].id[-2:] in frame:
		aa_seqs[key].id=aa_seqs[key].id[:-2]
		

	SeqIO.write(aa_seqs[key], OUTAA, "fasta")
	SeqIO.write(nt_seqs[aa_seqs[key].id], OUTNT, "fasta") #Odd way to do this, but we are changing the .id, not the dictionary key above, so need to keep key as it is for aa_seqs. Use the new .id in aa_seqs to get correct key for nt_seqs
