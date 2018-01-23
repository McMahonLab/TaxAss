###############################################################################
# convertFreshTrain.py
# Copyright (c) 2018, Joshua J Hamilton, Robin R Rohwer, and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Update the taxonomic strucutre of the FreshTrain reference database to be 
# compatabile with the mothur-compatible SILVA reference files.
#
# Tested on FreshTrain18Aug2016 and the mothur-compatible reference files from
# SILVA v132.
#
# Performs the following operations:
# Removes any sequence not assigned to a freshwater-specific lineage. We were
# too lazy to map the taxonomy of these individual sequences to SILVA. For the
# most part they aren't freshwater-specific, anyway.
#
# Update the taxonomy of any freshwater-specific lineages, as follows:
# SILVA Classification	FreshTrain Classification	Applies to Lineage
# c__Bacteroidia	c__[Saprospirae]	bacI, bac IV
# o__Chitinophagales	o__[Saprospirales]	bacI, bac IV
# c__Bacteroidia	c__Flavobacteriia	bacII, bacV
# c__Bacteroidia	c__Cytophagia	bacIII
# c__Bacteroidia	c__Sphingobacteriia	bacVI
# c__Gammaproteobacteria	c__Betaproteobacteria	betI, betII, betIII, betIV, betV, betVII
# o__Betaproteobacteriales	o__Burkholderiales	betI, betII, betIII, bet VII
# o__Betaproteobacteriales	o__Methylophilales	betIV
# o__Betaproteobacteriales	o__undefined	betV
# c__Verrucomicrobiae	c__[Methylacidiphilae]	LD19
# c__Verrucomicrobiae	c__[Spartobacteria]	verI-A, verI-B
# o__Chthoniobacterales	o__[Chthoniobacterales]	verI-A, verI-B
# o__Frankiales	o__Actinomycetales	acI, acSTL, acTH1
# o__Micrococcales	o__Actinomycetales	acIII
# o__Microtrichales	o__Acidimicrobiales	acIV
# o__Corynebacteriales	o__Actinomycetales	acTH2
# o__uncultured	o__Acidimicrobiales	acV
# o__SAR11_clade	o__Rickettsiales	alfV
# o__Acetobacterales	o__Rhodospirillales	alfVIII
# o__Micrococcales	o__Actinomycetales	Luna1, Luna3
# Note: SILVA represents the betaproteobacteria as an order within the
# gammaproteobacteria, so all 'bet' lineages get classified as gammaprotoes
# at the class level.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import os

#%%#############################################################################
### User-defined files and other declarations
################################################################################

inputFasta = 'custom.fasta'
inputTaxon = 'custom.taxonomy'
outputFasta = 'temp.fasta'
outputTaxon = 'temp.taxonomy'

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

# Dictionary of class-level disagreements. For each lineage, will update the
# class with the appropriate SILVA class-level taxonomic assignment
newClassDict = {
               'bacI': 'c__Bacteroidia',
               'bacII': 'c__Bacteroidia',
               'bacIII': 'c__Bacteroidia',
               'bacIV': 'c__Bacteroidia',
               'bacV': 'c__Bacteroidia',
               'bacVI': 'c__Bacteroidia',
               'betI': 'c__Gammaproteobacteria',
               'betII': 'c__Gammaproteobacteria',
               'betIII': 'c__Gammaproteobacteria',
               'betIV': 'c__Gammaproteobacteria',
               'betV': 'c__Gammaproteobacteria',
               'betVII': 'c__Gammaproteobacteria',
               'LD19': 'c__Verrucomicrobiae',
               'verI-A': 'c__Verrucomicrobiae',
               'verI-B': 'c__Verrucomicrobiae'
        }

# Dictionary of order-level disagreements. For each lineage, will update the
# order with the appropriate SILVA order-level taxonomic assignment
newOrderDict = {
        'acI': 'o__Frankiales',
        'acIII': 'o__Micrococcales',
        'acIV': 'o__Microtrichales',
        'acSTL': 'o__Frankiales',
        'acTH1': 'o__Frankiales',
        'acTH2': 'o__Corynebacteriales',
        'acV': 'o__uncultured',
        'alfV': 'o__SAR11_clade',
        'alfVIII': 'o__Acetobacterales',
        'bacI': 'o__Chitinophagales',
        'bacIV': 'o__Chitinophagales',
        'betI': 'o__Betaproteobacteriales',
        'betII': 'o__Betaproteobacteriales',
        'betIII': 'o__Betaproteobacteriales',
        'betIV': 'o__Betaproteobacteriales',
        'betV': 'o__Betaproteobacteriales',
        'betVII': 'o__Betaproteobacteriales',
        'Luna1': 'o__Micrococcales',
        'Luna3': 'o__Micrococcales',
        'verI-A': 'o__Chthoniobacterales',
        'verI-B': 'o__Chthoniobacterales'
        }


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
    seqID = curLine.split('\t')[0]
    seqTaxon = curLine.split(';')
    lineage = seqTaxon[4]
    
    # Check if the class needs to be updated and update id
    if lineage in newClassDict.keys():
        seqTaxon[2] = newClassDict[lineage]
        
    # Check if the class needs to be updated and update id
    if lineage in newOrderDict.keys():
        seqTaxon[3] = newOrderDict[lineage]
    
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

#%%#############################################################################
### Clean up your files
################################################################################

os.remove(inputFasta)
os.remove(inputTaxon)
os.move(outputFasta, inputFasta)
os.move(outputTaxon, inputTaxon)