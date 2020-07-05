#!/bin/bash

# RRR 2/11/16
# This is a sourceable bash script that runs steps 15-15.5 of the workflow.  
# NOTE that the pident chosen here generate the final otu taxonomy table.
# You will need to have run steps 1-14 with multiple pidents to make the plots in 15.5a

# command line syntax:
# ./RunStep_15.sh otus custom general 98 80 80 2

# ---- Receive input from terminal arguments --------------------------------------------------------------------------------------------

otus=$1
custom=$2
general=$3
pident=($4)
fwbootstrap=$5
ggbootstrap=$6
processors=$7

printf "Running TaxAss steps 15, 15.5a, and 15.5b.\n"
printf "\notu filenames: $otus.fasta and $otus.abund\n"
printf "custom database filenames: $custom.fasta and $custom.taxonomy\n"
printf "general database filenames: $general.fasta and $general.taxonomy\n"
printf "\nCreating the final taxonomy table with pident: $pident\n"
printf "\nClassification bootstrap confidence cutoffs are $fwbootstrap %% for the custom classification and $ggbootstrap %% for the general classification.\n"
printf "Using $processors processors.\n"

# ---- Make final taxonomy table --------------------------------------------------------------------------------------------------------
# 15
Rscript find_classification_disagreements.R ${otus}.$pident.taxonomy ${otus}.${general}.taxonomy ids.above.$pident conflicts_$pident $pident $fwbootstrap $ggbootstrap final &&

# ---- Plot how awesome your TaxAss taxonomy is -----------------------------------------------------------------------------------------
# 15.5.a
Rscript plot_classification_improvement.R final.taxonomy.pvalues final.general.pvalues total.reads.per.seqID.csv plots final.taxonomy.names final.general.names ids.above.$pident &&
# 15.5.b
mothur "#classify.seqs(fasta=${otus}.fasta, template=${custom}.fasta, taxonomy=${custom}.taxonomy, method=wang, probs=T, processors=$processors, cutoff=0)" &&
cat ${otus}.${custom}.wang.taxonomy > ${otus}.${custom}.taxonomy &&
sed 's/[[:blank:]]/\;/' <${otus}.${custom}.taxonomy >${otus}.${custom}.taxonomy.reformatted &&
mv ${otus}.${custom}.taxonomy.reformatted ${otus}.${custom}.taxonomy &&
mkdir conflicts_forcing
Rscript find_classification_disagreements.R ${otus}.${custom}.taxonomy ${otus}.$pident.$fwbootstrap.$ggbootstrap.taxonomy ids.above.$pident conflicts_forcing NA $fwbootstrap $ggbootstrap forcing &&
Rscript plot_classification_disagreements.R ${otus}.abund plots conflicts_forcing ${otus}.${custom}.$fwbootstrap.taxonomy ${otus}.$pident.$fwbootstrap.$ggbootstrap.taxonomy &&

# ---- done (except step 16 cleanup) ----------------------------------------------------------------------------------------------------
printf 'Steps 15, 15.5.a., and 15.5.b. have finished running.  When you are finished with all your analysis you can tidy up your working directory with step 16. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'

exit 0