#!/usr/bin/env python

import argparse
from Bio import SeqIO

#
# Matt Gitzendanner
# University of Florida
#
# 1kp_trim_taxa.py -i file.fa -o file.phy
#
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input fasta file")
parser.add_argument("-o", help="output fasta file")
parser.add_argument("-l", help="List of 4-letter codes to remove")
parser.add_argument("-m", help="Remove individuals samples where there is a combined sample? Default=1", default=1)

args = parser.parse_args()

InfileName = args.i
OutfileName = args.o
ListFileName = args.l
PruneMultipleTissue= args.m

def GetSpeciesName(RecordName):
	DashSplit=Record.id.split("-")
	DashSplitBits=DashSplit[0].split("_")
	UnderscoreList=["Asparagales", "Basal", "Cornales", "Non"]
	if Record_bits[0] in UnderscoreList: #Handle Basal_Angios, Basal_Eudicots and others with _ in name of clade			
		SpeciesBits=DashSplitBits[2:]
	else:
		SpeciesBits=DashSplitBits[1:]
	
	SpeciesName="_".join(SpeciesBits)
	return(SpeciesName)
	
try:
	InFile=open(InfileName, 'r')
except IOError:
	print "Can't open input fasta file: ", InfileName

try:
	OutFile=open(OutfileName, 'w')
except IOError:
		print "Can't open output fasta file", OutfileName


try:
	ListFile=open(ListFileName, 'r')
except IOError:
	print "Can't open list of 4-letter codes: ", ListFileName



# Build the list of taxa to exclude based on list file.
ExcludeList=[]

for Line in ListFile:
	Line=Line.strip('\n')
	ExcludeList.append(Line)


# Build list of taxa in the file:

CombinedSpeciedDict={} # Key=Genus_species, value=CODE :4-letter code of the combined sample.

for Record in SeqIO.parse(InFile, "fasta"):
	Record_bits = Record.id.split("_")		#split the record ID by _ 
	if "combined" in Record_bits:				#check if "combined" is part of the name.
		Code=Record_bits[-1] #last element will be 4-letter code
		
		if Code not in ExcludeList:
			Species=GetSpeciesName(Record.id)
		
			CombinedSpeciedDict[Species]=Code
		
InFile.close() #Close and then reopen, starting at top.

try:
	InFile=open(InfileName, 'r')
except IOError:
	print "Can't open input fasta file: ", InfileName


for Record in SeqIO.parse(InFile, "fasta"):
	Record_bits = Record.id.split("_")		#split the record ID by _ 
	Code=Record_bits[-1] #last element will be 4-letter code
	Species=GetSpeciesName(Record.id)
	if Code not in ExcludeList:
		if Species in CombinedSpeciedDict:
			if Code == CombinedSpeciedDict[Species]:
				SeqIO.write(Record, OutFile, "fasta")
			else:
				if PruneMultipleTissue == 1:
					print "Excluding single tissues of combined sample: \t%s\t%s" %(Species, Code)
				else:
					SeqIO.write(Record, OutFile, "fasta")
		else:
			SeqIO.write(Record, OutFile, "fasta")
	else:
		print "Excluding specified sample: \t%s\t%s" %(Species, Code)
			
			
	


	
	