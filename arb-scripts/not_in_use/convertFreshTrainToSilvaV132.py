# 2018-4-30 Robin Rohwer and Josh Hamilton

# Update the coarse-level taxonomic nomenclature of the FreshTrain reference database to be 
# compatabile with SILVA v132. Changes identified using Database_Improvement_Workflow.Rmd directions.

# Syntax:
# python convertFreshTrainToSilvaV128.py oldFreshTrain.taxonomy newFreshTrain.taxonomy

# oldFreshTrain.taxonomy is the gg-formatted mothur-formatted taxonomy file and
# newFreshTrain.taxonomy is the silva 128-formatted mothur-formatted output file


#%%#############################################################################
### Import packages
################################################################################

import sys

#%%#############################################################################
### User-defined files and other declarations
################################################################################

inputTaxon = sys.argv[1]
outputTaxon = sys.argv[2]

# inputTaxon = '/Users/joshamilton/Documents/Research/2015-TagTaxonomyWorkflow/2018-01-18-SILVA-FreshTrain/FreshTrain-files/FreshTrain25Jan2018NoBadLineage/FreshTrain25Jan2018NoBadLineage.taxonomy'
# outputTaxon = 'custom.taxonomy'

# Dictionary of class-level disagreements. For each lineage, will update the
# class with the appropriate SILVA class-level taxonomic assignment
newClassDict = {
               'bacI': 'Bacteroidia',
               'bacII': 'Bacteroidia',
               'bacIII': 'Bacteroidia',
               'bacIV': 'Bacteroidia',
               'bacV': 'Bacteroidia',
               'bacVI': 'Bacteroidia',
               'betI': 'Gammaproteobacteria',
               'betII': 'Gammaproteobacteria',
               'betIII': 'Gammaproteobacteria',
               'betIV': 'Gammaproteobacteria',
               'betV': 'Gammaproteobacteria',
               'betVII': 'Gammaproteobacteria',
               'LD19': 'Verrucomicrobiae',
               'verI-A': 'Verrucomicrobiae',
               'verI-B': 'Verrucomicrobiae'
        }

# Dictionary of order-level disagreements. For each lineage, will update the
# order with the appropriate SILVA order-level taxonomic assignment
newOrderDict = {
        'acI': 'Frankiales',
        'acIII': 'Micrococcales',
        'acIV': 'Microtrichales',
        'acSTL': 'Frankiales',
        'acTH1': 'Frankiales',
        'acTH2': 'Corynebacteriales',
        'alfV': 'SAR11_clade',
        'alfVIII': 'Acetobacterales',
        'bacI': 'Chitinophagales',
        'bacIV': 'Chitinophagales',
        'betI': 'Betaproteobacteriales',
        'betII': 'Betaproteobacteriales',
        'betIII': 'Betaproteobacteriales',
        'betIV': 'Betaproteobacteriales',
        'betV': 'Betaproteobacteriales',
        'betVII': 'Betaproteobacteriales',
        'Luna1': 'Micrococcales',
        'Luna3': 'Micrococcales',
        'verI-A': 'Chthoniobacterales',
        'verI-B': 'Chthoniobacterales'
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