###############################################################################
# removeBadLineages.py
# Copyright (c) 2018, Joshua J Hamilton, Robin R Rohwer, and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Remove any sequence not assigned to a freshwater-specific lineage.
# This script requires four arguments, read from the command line:
# inputFasta = Input fasta file, from which to remove lineages: i.e., FreshTrain25Jan2018.fasta
# inputTaxon = Input taxonomy file, from which to remove lineages: i.e., FreshTrain25Jan2018.taxonomy
# outputFasta = Output fasta file, i.e., custom.fasta
# outputTaxon = Output fasta file, i.e., custom.fasta
# From the input files, this script removes any sequences not assigned to a
# freshwate-specific lineage, and write the remaining sequences to the output
# files
#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import sys

#%%#############################################################################
### User-defined files and other declarations
################################################################################

inputFasta = sys.argv[1]
inputTaxon = sys.argv[2]
outputFasta = sys.argv[3]
outputTaxon = sys.argv[4]

# Valid lineages. We discard any sequence not belonging to this lineage.
lineageList = ['acI','acIII','acIV','acSTL','acTH1','acTH2','acV',
               'alfI','alfI-A','alfI-B','alfII','alfIII','alfIV','alfV','alfVI','alfVII','alfVIII',
               'bacI','bacII','bacIII','bacIV','bacV','bacVI',
               'betI','betII','betIII','betIV','betV','betVII',
               'gamI','gamII','gamIII','gamIV','gamV',
               'LD19','Luna1','Luna3','verI-A','verI-B',
]

# List of sequence IDs that don't have a valid lineage. Will get populated
# when correcting the taxonomy file.
badLineageList = []

#%%#############################################################################
### Read in the taxonomy file and remove any sequences that don't have a valid
### lineage. Record the seequence IDs.
################################################################################

inTaxonHandle = open(inputTaxon, 'r')
outTaxonHandle = open(outputTaxon, 'w')

for curLine in inTaxonHandle.readlines():
    curLine = curLine.strip()    

    # Split along a tab to separate seqID from taxonomy
    # Grab the lineage from the taxonomy information
    [seqID, seqTaxon] = curLine.split('\t')
    seqTaxon = seqTaxon.split(';')
    lineage = seqTaxon[4]
    
    # If the lineage is valid, write the new taxonomy file
    if lineage in lineageList:
        newTaxonString = seqID+'\t'+';'.join(seqTaxon)
        outTaxonHandle.write(newTaxonString+'\n')    

    # Otherwise, record the seqID
    else:
        badLineageList.append(seqID)

outTaxonHandle.close()

#%%#############################################################################
### Read in the sequence file and drop all sequences that don't have a valid
### lineage
################################################################################

inFastaHandle = open(inputFasta, 'r')
outFastaHandle = open(outputFasta, 'w')

for curRecord in SeqIO.parse(inFastaHandle, 'fasta') :

    # If the lineage is valid, write the new taxonomy file
    if curRecord.name not in badLineageList:        
        outFastaHandle.write('>'+curRecord.id+'\n')
        outFastaHandle.write(str(curRecord.seq)+'\n')
    
outFastaHandle.close()