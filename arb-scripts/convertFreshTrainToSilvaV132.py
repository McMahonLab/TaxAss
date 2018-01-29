###############################################################################
# convertFreshTrainToSilvaV132.py
# Copyright (c) 2018, Joshua J Hamilton, Robin R Rohwer, and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Update the taxonomic strucutre of the FreshTrain reference database to be 
# compatabile with the mothur-compatible SILVA V132 reference files.
#
# Tested on FreshTrain25Jan2018 and the mothur-compatible reference files from
# SILVA v132.
#
# Updates the taxonomy of any freshwater-specific lineages, as follows:
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
# Note: SILVA V132 represents the betaproteobacteria as an order within the
# gammaproteobacteria, so all 'bet' lineages get classified as gammaproteobacteria
# at the class level.

# This script requires two arguments, read from the command line:
# inputTaxon = Input taxonomy file, to make compliant with FreshTrain
# outputTaxon = Output fasta file, i.e., custom.txonomy
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import sys

#%%#############################################################################
### User-defined files and other declarations
################################################################################

#inputTaxon = sys.argv[1]
#outputTaxon = sys.argv[2]

inputTaxon = '/Users/joshamilton/Documents/Research/2015-TagTaxonomyWorkflow/2018-01-18-SILVA-FreshTrain/FreshTrain-files/FreshTrain25Jan2018NoBadLineage/FreshTrain25Jan2018NoBadLineage.taxonomy'
outputTaxon = 'custom.taxonomy'

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
### Read in the taxonomy file and update the class and lineage as needed.
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
    
    # Check if the class needs to be updated and update id
    if lineage in newClassDict.keys():
        seqTaxon[2] = newClassDict[lineage]
        
    # Check if the class needs to be updated and update id
    if lineage in newOrderDict.keys():
        seqTaxon[3] = newOrderDict[lineage]
     
    newTaxonString = seqID+'\t'+';'.join(seqTaxon)
    outTaxonHandle.write(newTaxonString+'\n')    
    

outTaxonHandle.close()