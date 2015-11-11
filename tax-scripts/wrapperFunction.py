###############################################################################
# wrapperFunction.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to fun the 16S taxonomic assignment pipeline on the
# iTag data from the Lake Mendota time-series
################################################################################

#%%#############################################################################
### Import packages
################################################################################
import glob
import operator
import os
import re
import shutil
import subprocess

#%%#############################################################################
### User-defined files and folder structure
################################################################################
# Define data folders
otuCountsFolder = 'rawData'
processedDataFolder = 'processedData'

# Define user-specified variables
pident = 90
numSeqs = 5
queryFasta = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/MendotaTags/rawData/otus.fasta'
customFasta = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/16STaxAss-GitHubRepo/databases/FWonly_7Sept2015.fasta'
customTaxonomy = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/16STaxAss-GitHubRepo/databases/FWonly_7Sept2015.taxonomy'
generalFasta = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/16STaxAss-GitHubRepo/databases/gg_13_5.fasta'
generalTaxonomy = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/16STaxAss-GitHubRepo/databases/gg_13_5.taxonomy'

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
pathToScripts = '/Users/joshamilton/Documents/Research/2015d-TagTaxonomyWorkflow/16STaxAss-GitHubRepo/tax-scripts'

# Names for output files generated during script operation

# BLASTing
customDB = re.search('(.+)\.fasta', customFasta).group(1)
customDBstring = re.search('.+/(.+)\.fasta', customFasta).group(1)
generalDBstring = re.search('.+/(.+)\.fasta', generalFasta).group(1)
otuFileName = re.search('.+/(.+)\.fasta', queryFasta).group(1)
queryBlastResults = processedDataFolder+'/'+otuFileName+'.custom.blast'
formattedBlastResults = processedDataFolder+'/'+otuFileName+'.custom.table'

# Finding queries to be classifed by curated and general databases
queryAbovePident = processedDataFolder+'/'+otuFileName+'.above.'+str(pident)
hitsAbovePident = processedDataFolder+'/'+otuFileName+'.above.'+str(pident)+'.hitsAbove.csv'
otusAbovePident = processedDataFolder+'/'+otuFileName+'.above.'+str(pident)+'.fasta'
queryBelowPident = processedDataFolder+'/'+otuFileName+'.below.'+str(pident)
hitsBelowPident = processedDataFolder+'/'+otuFileName+'.below.'+str(pident)+'.hitsBelow.csv'
otusBelowPident = processedDataFolder+'/'+otuFileName+'.below.'+str(pident)+'.fasta'
queryMissing = processedDataFolder+'/'+otuFileName+'.missing.'+str(pident)

#%%#############################################################################
### Step 0 - Data Formatting
################################################################################

# See dataProcessing.py file
        
#%%#############################################################################
### Step 1 - Make BLAST Database
################################################################################

subprocess.call(['makeblastdb',
    '-dbtype', 'nucl',
    '-in', customFasta,
    '-input_type', 'fasta',
    '-parse_seqids',
    '-out', customDB])

#%%#############################################################################
### Step 2 - Run BLAST
################################################################################

subprocess.call(['blastn',
                 '-query', queryFasta,
                 '-task', 'megablast',
                 '-db', customDB,
                 '-out', queryBlastResults,
                 '-outfmt', '11',
                 '-max_target_seqs', str(numSeqs)])
                 
#%%#############################################################################
### Step 3 - Reformat BLAST Database
################################################################################
    
subprocess.call(['blast_formatter',
                 '-archive', queryBlastResults,
                 '-outfmt', '6 qseqid pident length qlen qstart qend',
                 '-out', formattedBlastResults])

#%%#############################################################################
### Step 4 - Filter BLAST Results
################################################################################

subprocess.call(['Rscript', pathToScripts+'/find_seqIDs_with_pident.R',
                 formattedBlastResults,
                 queryAbovePident, hitsAbovePident, str(pident), 'TRUE'])
                 
subprocess.call(['Rscript', pathToScripts+'/find_seqIDs_with_pident.R',
                 formattedBlastResults,
                 queryBelowPident, hitsBelowPident, str(pident), 'FALSE'])
                 
#%%#############################################################################
### Step 4.5 - Construct plots
################################################################################
                 
subprocess.call(['Rscript', pathToScripts+'/plot_blast_hit_stats.R',
                 hitsAbovePident, str(pident)])
                 
#%%#############################################################################
### Step 5 - Recover sequence IDs with no BLAST hit
################################################################################
subprocess.call(['python', pathToScripts+'/fetch_seqIDs_blast_removed.py',
                 queryFasta, formattedBlastResults, queryBelowPident])
                 
#%%#############################################################################
### Step 6 - Filter BLAST Results
################################################################################
                 
subprocess.call(['python', pathToScripts+'/fetch_fastas_with_seqIDs.py',
                 queryAbovePident, queryFasta, otusAbovePident])

subprocess.call(['python', pathToScripts+'/fetch_fastas_with_seqIDs.py',
                 queryBelowPident, queryFasta, otusBelowPident])
                 
#%%#############################################################################
### Step 7 - Run mothur
################################################################################

subprocess.call(['mothur',
                 '#classify.seqs(fasta='+otusAbovePident+', template='
                 +customFasta+', taxonomy='+customTaxonomy
                 +', method=wang, probs=T, processors=2)'])


subprocess.call(['mothur',
                 '#classify.seqs(fasta='+otusBelowPident+', template='
                 +generalFasta+', taxonomy='+generalTaxonomy
                 +', method=wang, probs=T, processors=2)'])
                 
#%%#############################################################################
### Step 8 - Concatenate taxonomy files
################################################################################

fileList = [processedDataFolder+'/'+otuFileName+'.above.'+str(pident)+'.'+customDBstring+'.wang.taxonomy', processedDataFolder+'/'+otuFileName+'.below.'+str(pident)+'.'+generalDBstring+'.wang.taxonomy']
with open(processedDataFolder+'/'+otuFileName+'.taxonomy', 'w') as outfile:
    for fname in fileList:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
                
#%%#############################################################################
### Step 9 - Run mothur again
################################################################################

subprocess.call(['mothur',
                 '#classify.seqs(fasta='+otuCountsFolder+'/'+otuFileName+'.fasta, template='
                 +generalFasta+', taxonomy='+generalTaxonomy
                 +', method=wang, probs=T, processors=2)'])

# Move output files to the proper directory

for file in glob.glob(otuCountsFolder+'/*gg*'):                                                                                                               
    shutil.move(file, processedDataFolder)
    
#%%#############################################################################
### Step 10 - Replace all tabs with semicolons in taxonomy files
################################################################################
    
fileList = [processedDataFolder+'/'+otuFileName+'.taxonomy', processedDataFolder+'/'+otuFileName+'.'+generalDBstring+'.wang.taxonomy', ]

for fname in fileList:
# Read in the file
    tempFile = None
    with open(fname, 'r') as file :
        tempFile = file.read()

# Replace the target string
    tempFile = tempFile.replace('\t', ';')

# Write the file out again
    with open(fname, 'w') as file:
        file.write(tempFile)