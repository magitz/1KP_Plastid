#!/usr/bin/env python

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re


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

KeepRecordList=[]

for record in SeqIO.parse(IN, "fasta") :
	if "*" in record.seq or "J" in record.seq:			
		print ("Skipping %s, has *s or Js in translated AA sequence." %(record.id))
		
	else:
		SeqIO.write(record, OUTAA, "fasta")
		
		#Some AA seq names have _1, _2, _3, from the frame used in the translation. That's not in the nt name, so need to remove that.
		frame=["_1", "_2", "_3"]
		if record.id[-2:] in frame:
			record.id=record.id[:-2]

		KeepRecordList.append(record.id)

for record in SeqIO.parse(NT, "fasta") :
	if record.id in KeepRecordList:
		SeqIO.write(record, OUTNT, "fasta")
