#!/bin/bash

# This is a sourceable bash script that runs step 16 of the workflow.  
# Clean-up step: intermediate files will be deleted, and all other files 
# will be sorted into 4 folders: analysis, data, databases, scripts.
# RRR 2/11/16

# Make sure you're really done before deleting stuff.
# This leaves everything you'd need to re-analyze from the beginning,
# But removes things you'd need to re-analyze from the middle.

# 16
rm custom.db.* custom.8mer custom.custom* custom.tree* general.8mer general.general* general.tree* *wang* mothur.*.logfile otus.custom.blast* ids* otus.below*.fasta otus.above*.fasta otus.[0-9][0-9].taxonomy otus.[0-9][0-9][0-9].taxonomy otus.general.taxonomy otus.custom.taxonomy otus.custom.[0-9]* custom.general* *pvalues total* final*names
mkdir scripts ; mv *.py *.R *.sh *.md *.Rmd *.html scripts
mkdir analysis ; mv conflicts* *.txt plots analysis
mkdir data ; mv otus* *.rds data
mkdir databases ; mv *.taxonomy *.fasta databases

printf 'All Done! \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0