"""
RRR 8-25-15

This script takes the sequence IDs selected from the blast output 
and finds them in the original query fasta file.  
It then creates a new fasta file containing just the desired sequences.

The format for typing into the terminal is:

$ python fetch_fastas_with_seqIDs.py idfile fastafile outputfile

where 	fetch_fastas_with_seqIDs.py is the path to this script
		idfile  is the file with \n separated seqIDs
		fastafile is the fasta file containing all the seqIDs
		outputfile is the new fasta file the script generates
"""

import os	# a library for manipulating files on your computer
import sys	# a library for reading arguments from the command line

# Read input arguments from command line into variable names
idfile = sys.argv[1]
fastafile = sys.argv[2]
outputfile = sys.argv[3]

# delete the output file if it already exists so that you append to a blank file
if os.path.isfile(outputfile) :
	os.remove(outputfile)
	print("Existing file with your chosen output file name was deleted.")


# Create hash of all SeqIDs in the idfile
pidentIDs = {}
with open(idfile) as ID:
    for line in ID:
            pidentIDs[line.strip()] = None
            
# Create hash of all SeqIDs and sequences in query fasta file
allIDs = {}
with open(fastafile) as fasta:
    for line in fasta:
        if str.startswith(line, '>') :
            key = line.strip()[1:]
        else :
            val = line.strip()
            allIDs[key] = val
            
## Generate the output file
SeqIDsFile = open(idfile,"r")

with open(outputfile, "a") as ResultFile :
    for key in pidentIDs :
        ResultFile.write('>'+key+'\n')
        ResultFile.write(allIDs[key]+'\n')
SeqIDsFile.close()

