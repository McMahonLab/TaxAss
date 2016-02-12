#!/bin/bash

# This is a sourceable bash script that runs steps 15-16 of the workflow.  
# NOTE that the pident chosen here is the final taxonomy one, you may want to change it.
# RRR 2/11/16

# Choose final pident.

pident=("98")

# 15
Rscript find_classification_disagreements.R otus.$pident.taxonomy otus.general.taxonomy ids.above.$pident conflicts_$pident $pident 85 70 final
# 16
Rscript plot_classification_improvement.R final.taxonomy.pvalues final.general.pvalues total.reads.per.seqID.csv plots
mothur "#classify.seqs(fasta=otus.fasta, template=custom.fasta, taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
cat otus.custom.wang.taxonomy > otus.custom.taxonomy
sed 's/[[:blank:]]/\;/' <otus.custom.taxonomy >otus.custom.taxonomy.reformatted
mv otus.custom.taxonomy.reformatted otus.custom.taxonomy
mkdir conflicts_forcing
Rscript find_classification_disagreements.R otus.custom.taxonomy otus.$pident.85.70.taxonomy ids.above.$pident conflicts_forcing NA 85 70 forcing
Rscript plot_classification_disagreements.R otus.abund plots NA conflicts_forcing otus.custom.85.taxonomy NA NA

printf 'Steps 15-16 have finished running.  When you\'re finished with all your analysis you can tidy up your working directory with step 17. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0