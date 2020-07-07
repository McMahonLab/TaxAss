# 2018-4-30 Robin Rohwer and Josh Hamilton

# Update the coarse-level taxonomic nomenclature of the FreshTrain reference database to be 
# compatabile with SILVA v128. Changes identified using Database_Improvement_Workflow.Rmd directions.

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

# Dictionary of class-level disagreements. For each lineage, will update the
# class with the appropriate SILVA class-level taxonomic assignment
newClassDict = {
               'bacI': 'Sphingobacteriia',
               'bacIV': 'Sphingobacteriia',
               'verI-A': 'Spartobacteria',
               'verI-B': 'Spartobacteria',
               'LD19': 'Verrucomicrobia_Incertae_Sedis'
        }

# Dictionary of order-level disagreements. For each lineage, will update the
# order with the appropriate SILVA order-level taxonomic assignment
newOrderDict = {
        'acI': 'Frankiales',
        'acIII': 'Micrococcales',
        'acSTL': 'Frankiales',
        'acTH1': 'Frankiales',
        'acTH2': 'Corynebacteriales',
        'alfV': 'SAR11_clade',
        'bacI': 'Sphingobacteriales',
        'bacIV': 'Sphingobacteriales',
        'betV': 'Nitrosomonadales',
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