###############################################################################
# create_fastas_given_seqIDs
# Copyright (c) 2016, Joshua J Hamilton, Robin R Rohwer, and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This script takes the sequence IDs specified in the input file and finds
# them in the original query fasta file. It then creates a new fasta file 
# containing just the desired sequences.
#
# The format for typing into the terminal is:
# $ python fetch_fastas_with_seqIDs.py idFile fastaFile outputFile
# where 	fetch_fastas_with_seqIDs.py		path to this script
#		idFile                           the file with \n separated seqIDs
#		fastaFile                        the fasta file containing all the seqIDs
#		outputFile                       the new fasta file the script generates

#%%#############################################################################
### Import packages
################################################################################
from Bio import SeqIO
import os
import sys

#%%#############################################################################
### Read input arguments from command line into variable names
################################################################################
idFile = sys.argv[1]
fastaFile = sys.argv[2]
outputFile = sys.argv[3]

#%%#############################################################################
### Process the list of all seqIDs in the input file and retrieve the resulting
### sequences from the fasta file. Create a new fasta file containing just
### those sequences.
################################################################################

# Delete the output file if it already exists so that you append to a blank file
if os.path.isfile(outputFile) :
	os.remove(outputFile)
	print("Existing file with your chosen output file name was deleted.")


# Create hash of all SeqIDs in the idFile
pidentIDs = {}
with open(idFile) as ID:
    for line in ID:
        pidentIDs[line.strip()] = None
            
# Create hash of all SeqIDs and sequences in query fasta file
allIDs = {}
fileHandle = open(fastaFile)
for record in SeqIO.parse(fileHandle, "fasta") :
    allIDs[record.id] = record
fileHandle.close()
            
# Generate the output file
SeqIDsFile = open(outputFile,"w")
for key in pidentIDs :
    SeqIO.write(allIDs[key], SeqIDsFile, "fasta")
SeqIDsFile.close()



