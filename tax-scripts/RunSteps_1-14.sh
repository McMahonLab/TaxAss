#!/bin/bash

# This is a sourceable bash script that runs steps 1-13 of the workflow.  
# In addition to avoiding problems that happened with pasting blocks of code,
# This lets you run multiple pidents in paralelle which makes it faster.
# The actual commands in this script are identical to the workflow commands
# that you type directly into the terminal command line.
# The only additional code is an array of pidents and a loop to run them all.
# RRR 1/12/16

# Choose pidents to test over.

pident=("100" "99" "98" "97" "96" "95")
fwbootstrap=("80")
ggbootstrap=("80")

# First Run steps 1-13 to generate databases and folders exactly following workflow
# Note: still gotta do the reformatting manually (step 0)

# 1
makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db &&
# 2
blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast -outfmt 11 -max_target_seqs 5 &&
# 3
blast_formatter -archive otus.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table &&
# 4
Rscript calc_full_length_pident.R otus.custom.blast.table otus.custom.blast.table.modified &&
#5
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.${pident[0]} ${pident[0]} TRUE &&
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.${pident[0]} ${pident[0]} FALSE &&
# 6
mkdir plots &&
RScript plot_blast_hit_stats.R otus.custom.blast.table.modified ${pident[0]} plots &&
# 7
python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table.modified ids.missing &&
cat ids.below.${pident[0]} ids.missing > ids.below.${pident[0]}.all &&
# 8
python create_fastas_given_seqIDs.py ids.above.${pident[0]} otus.fasta otus.above.${pident[0]}.fasta &&
python create_fastas_given_seqIDs.py ids.below.${pident[0]}.all otus.fasta otus.below.${pident[0]}.fasta &&
# 9
mothur "#classify.seqs(fasta=otus.above.${pident[0]}.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2, cutoff=0)" &&
mothur "#classify.seqs(fasta=otus.below.${pident[0]}.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2, cutoff=0)" &&
# 10
cat otus.above.${pident[0]}.custom.wang.taxonomy otus.below.${pident[0]}.general.wang.taxonomy > otus.${pident[0]}.taxonomy &&
# 11
mothur "#classify.seqs(fasta=otus.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2, cutoff=0)" &&
cat otus.general.wang.taxonomy > otus.general.taxonomy &&
# optional step 11.5 database comparison not run here.
# 12
sed 's/[[:blank:]]/\;/' <otus.${pident[0]}.taxonomy >otus.${pident[0]}.taxonomy.reformatted &&
mv otus.${pident[0]}.taxonomy.reformatted otus.${pident[0]}.taxonomy &&
sed 's/[[:blank:]]/\;/' <otus.general.taxonomy >otus.general.taxonomy.reformatted &&
mv otus.general.taxonomy.reformatted otus.general.taxonomy &&
# 13
mkdir conflicts_${pident[0]} &&
Rscript find_classification_disagreements.R otus.${pident[0]}.taxonomy otus.general.taxonomy ids.above.${pident[0]} conflicts_${pident[0]} ${pident[0]} $fwbootstrap $ggbootstrap &&

# Next Run steps 4-9 and 11-12 with different pident cutoffs
# Define a function called runagain since you repeat this part many times in paralelle

# ** should I put && at the ends of the lines here or would that mess up the paralellization in the runagain loop?
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
   mothur "#classify.seqs(fasta=otus.above.$1.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2, cutoff=0)"
   mothur "#classify.seqs(fasta=otus.below.$1.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2, cutoff=0)"
   # 10
   cat otus.above.$1.custom.wang.taxonomy otus.below.$1.general.wang.taxonomy > otus.$1.taxonomy
   # 12 a,b
   sed 's/[[:blank:]]/\;/' <otus.$1.taxonomy >otus.$1.taxonomy.reformatted
   mv otus.$1.taxonomy.reformatted otus.$1.taxonomy
   # 13
   mkdir conflicts_$1
   Rscript find_classification_disagreements.R otus.$1.taxonomy otus.general.taxonomy ids.above.$1 conflicts_$1 $1 $fwbootstrap $ggbootstrap
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
# the && makes it only run the next step if the previous step exited successfully

# Now run step 14- plotting everything together to choose final pident
# First generate the arguments for the command call:

args_string=""
for p in ${pident[*]}
do
   args_string+=" conflicts_$p ids.above.$p $p"
done

# 14
Rscript plot_classification_disagreements.R otus.abund plots regular NA NA $args_string &&

printf 'Steps 1-14 have finished running.  Now analysze the plots from step 14 to choose your final pident and generate your final taxonomy file in step 15.  Optionally you can then check out your classification improvement in step 15.5. At the end delete intermediate files with step 16. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0
