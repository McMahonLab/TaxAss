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

# First Run steps 1-11 to generate databases and folders exactly following workflow
# Note: still gotta do the reformatting on your own (step 0)


# Next Run steps 4-9 and 11-12 with different pident cutoffs
# Define a function called runagain since you repeat this part many times in paralelle

runagain () {
   # 5
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.$1 $1 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.$1 $1 FALSE
   # 6 b
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified $1 plots
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

# Now run step 14- plotting everything together to choose final pident
# First generate the arguments for the command call:

args_string=""
for p in ${pident[*]}
do
   args_string+=" conflicts_$p ids.above.$p $p"
done

# 14
Rscript plot_classification_disagreements.R otus.abund plots regular NA $args_string

printf 'Steps 1-14 have finished running.  Now analysze the plots from step 14 to choose your final pident and generate your final taxonomy file in step 15.  Optionally you can compare to how your taxonomy would have been in step 16. At the end tidy up your working directory with step 17. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0
