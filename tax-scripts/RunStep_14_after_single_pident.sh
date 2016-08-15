#!/bin/bash

# This is a sourceable bash script that runs steps 1-13 of the workflow.  
# In addition to avoiding problems that happened with pasting blocks of code,
# This lets you run multiple pidents in paralelle which makes it faster.
# The actual commands in this script are identical to the workflow commands
# that you type directly into the terminal command line.
# The only additional code is an array of pidents and a loop to run them all.
# RRR 1/12/16

# Do this after checking for bugs with RunSteps_1-13_single_pident.sh 
# and then running 1-13 in paralelle with RunSteps_1-13_after_single_pident.sh

# Choose pidents to test over.

pident=("100" "99" "98" "97" "96" "95")


# Now run step 14- plotting everything together to choose final pident
# First generate the arguments for the command call:

args_string=""
for p in ${pident[*]}
do
   args_string+=" conflicts_$p ids.above.$p $p"
done

# 14
Rscript plot_classification_disagreements.R otus.abund plots regular NA $args_string

printf 'Step 14 has finished running.  Now analysze the plots from step 14 to choose your final pident, then generate your final taxonomy file in step 15.  Optionally you can look at how much better your taxonomy is now in step 15.5. At the end tidy up your working directory and delete large intermediate files with step 16. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0
