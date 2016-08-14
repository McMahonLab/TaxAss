#!/bin/bash

# This is a sourceable bash script that runs steps 1-13 of the workflow.  
# In addition to avoiding problems that happened with pasting blocks of code,
# This lets you run multiple pidents in paralelle which makes it faster.
# The actual commands in this script are identical to the workflow commands
# that you type directly into the terminal command line.
# The only additional code is an array of pidents and a loop to run them all.
# RRR 1/12/16

# Do this after checking for bugs with RunSteps_1-13_single_pident.sh

# Choose pidents to test over.

pident=("100" "99" "98" "97" "96" "95")


# Next Run steps 5, 7-10, & 12-13 with different pident cutoffs
# Define a function called runagain since you repeat this part many times in paralelle

runagain () {
   # 5
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.$1 $1 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.$1 $1 FALSE
   # 7 b
   cat ids.below.$1 ids.missing > ids.below.$1.all
   # 8
   python create_fastas_given_seqIDs.py ids.above.$1 otus.fasta otus.above.$1.fasta
   python create_fastas_given_seqIDs.py ids.below.$1.all otus.fasta otus.below.$1.fasta
   # 9
   mothur "#classify.seqs(fasta=otus.above.$1.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.$1.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   # 10
   cat otus.above.$1.custom.wang.taxonomy otus.below.$1.general.wang.taxonomy > otus.$1.taxonomy
   # 12 a,b
   sed 's/[[:blank:]]/\;/' <otus.$1.taxonomy >otus.$1.taxonomy.reformatted
   mv otus.$1.taxonomy.reformatted otus.$1.taxonomy
   # 13
   mkdir conflicts_$1
   Rscript find_classification_disagreements.R otus.$1.taxonomy otus.general.taxonomy ids.above.$1 conflicts_$1 $1 85 70
}

# Using function runagain run the additional pidents in paralelle

length=${#pident[*]}

for p in ${pident[*]:1:$length} 
do 
	runagain $p &
done 
wait

# the & lets the levels of the loop run in parallel
# the wait makes sure all the loops finish before the script exits

# Next run step 14- plotting everything together to choose final pident


printf 'Steps 1-13 have finished running for the rest of the pidents. Next run step 14 to choose your cutoff. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0
