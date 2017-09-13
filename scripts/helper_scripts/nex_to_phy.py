#!/usr/bin/env python

import argparse
from Bio import AlignIO

# A simple converter from Nexus to Phylip format using BioPython.
# Matt Gitzendanner
# University of Florida

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file")
parser.add_argument("-o", help="output file")

args = parser.parse_args()

infile = args.i
outfile = args.o

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	OUT=open(outfile, 'a')
except IOError:
		print "Can't open file", outfile
		
alignment = AlignIO.read(IN, "nexus")
AlignIO.write([alignment], OUT, "phylip-relaxed")
