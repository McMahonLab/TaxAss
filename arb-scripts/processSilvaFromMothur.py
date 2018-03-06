###############################################################################
# processSilvaFromMothur.py
# Copyright (c) 2018, Joshua J Hamilton, Robin R Rohwer, and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Before running this script, download the mothur-compatible SILVA files 
# and reformat for use with TaxAss. As of 2018-01-25, these files can be 
# accessed from:
#   https://www.mothur.org/wiki/Silva_reference_files
# You will want to download the "Full length sequences and taxonomy references."

# This script requires four arguments, read from the command line:
# inputFasta = Input fasta file, to make compliant with TaxAss
# inputTaxon = Input taxonomy file, to make compliant with TaxAss
# outputFasta = Output fasta file, i.e., custom.fasta
# outputTaxon = Output fasta file, i.e., custom.txonomy
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import re
import sys

#%%#############################################################################
### User-defined files
################################################################################

inputFasta = sys.argv[1]
inputTaxon = sys.argv[2]
outputFasta = sys.argv[3]
outputTaxon = sys.argv[4]

#%%#############################################################################
### Read in the alignment file and reformat
################################################################################
inFastaHandle = open(inputFasta, 'r')
outFastaHandle = open(outputFasta, 'w')

for record in SeqIO.parse(inFastaHandle, 'fasta') :

    # Remove record description so fasta gets printed without it
    record.id = re.sub('\.', '', record.id)

    # Write the results to file
    record.seq = re.sub('\.', '', str(record.seq))
    record.seq = re.sub('-', '', str(record.seq))
    
    outFastaHandle.write('>'+record.id+'\n')
    outFastaHandle.write(record.seq+'\n')
    
outFastaHandle.close()
        
#%%#############################################################################
### Read in the taxonomy file and reformat
################################################################################
inTaxonHandle = open(inputTaxon, 'r')
outTaxonHandle = open(outputTaxon, 'w')

for curLine in inTaxonHandle.readlines():
    curLine = curLine.strip()    

    # 1st part is the seqID, remove the '.' for consistency
    seqID = curLine.split('\t')[0]
    seqID = re.sub('\.', '', seqID)

    taxonString = curLine.split('\t')[1]
    taxonArray = taxonString.split(';')
    
    if len(taxonArray) == 7:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+taxonArray[2]+';o__'+taxonArray[3]+';f__'+taxonArray[4]+ \
        ';g__'+taxonArray[5]+';s__'+taxonArray[6]+';'
    elif len(taxonArray) == 6:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+taxonArray[2]+';o__'+taxonArray[3]+';f__'+taxonArray[4]+ \
        ';g__'+taxonArray[5]+';s__'+';'
    elif len(taxonArray) == 5:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+taxonArray[2]+';o__'+taxonArray[3]+';f__'+taxonArray[4]+ \
        ';g__'+';s__'+';'
    elif len(taxonArray) == 4:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+taxonArray[2]+';o__'+taxonArray[3]+';f__'+';g__'+';s__'+';'
    elif len(taxonArray) == 3:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+taxonArray[2]+';o__'+';f__'+';g__'+';s__'+';'
    elif len(taxonArray) == 2:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+taxonArray[1]+ \
        ';c__'+';o__'+';f__'+';g__'+';s__'+';'
    elif len(taxonArray) == 1:
        newTaxonString = seqID+'\t'+'k__'+taxonArray[0]+';p__'+';c__'+';o__'+';f__'+';g__'+';s__'+';'

    outTaxonHandle.write(newTaxonString+'\n')
        
outTaxonHandle.close()