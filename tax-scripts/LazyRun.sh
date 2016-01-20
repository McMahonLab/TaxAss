#!/bin/bash

# This is a sourceable bash script that runs steps 1-13 of the workflow.  
# In addition to avoiding problems that happened with pasting blocks of code,
# This lets you run multiple pidents in paralelle which makes it faster.
# The actual commands in this script are identical to the workflow commands
# that you type directly into the terminal command line.
# The only additional code is an array of pidents and a loop to run them repeatedly.
# RRR 1/12/16

# Choose pidents to test over.

pident=("100" "98" "96" "94")

# First Run steps 1-12 to generate databases and folders exactly following workflow
# Note: still gotta do the reformatting on your own (step 0)

makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db
blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast -outfmt 11 -max_target_seqs 5
blast_formatter -archive otus.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table
Rscript calc_full_length_pident.R otus.custom.blast.table otus.custom.blast.table.modified
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.${pident[0]} ${pident[0]} TRUE 
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.${pident[0]} ${pident[0]} FALSE
mkdir plots
RScript plot_blast_hit_stats.R otus.custom.blast.table.modified ${pident[0]} plots
python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
cat ids.below.${pident[0]} ids.missing > ids.below.${pident[0]}.all
python create_fastas_given_seqIDs.py ids.above.${pident[0]} otus.fasta otus.above.${pident[0]}.fasta
python create_fastas_given_seqIDs.py ids.below.${pident[0]}.all otus.fasta otus.below.${pident[0]}.fasta
mothur "#classify.seqs(fasta=otus.above.${pident[0]}.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
mothur "#classify.seqs(fasta=otus.below.${pident[0]}.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
cat otus.above.${pident[0]}.custom.wang.taxonomy otus.below.${pident[0]}.general.wang.taxonomy > otus.${pident[0]}.taxonomy
mothur "#classify.seqs(fasta=otus.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
cat otus.general.wang.taxonomy > otus.general.taxonomy
mothur "#classify.seqs(fasta=custom.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
mv custom.general.wang.taxonomy custom.general.taxonomy
sed 's/[[:blank:]]/\;/' <otus.${pident[0]}.taxonomy >otus.${pident[0]}.taxonomy.reformatted
mv otus.${pident[0]}.taxonomy.reformatted otus.${pident[0]}.taxonomy
sed 's/[[:blank:]]/\;/' <otus.general.taxonomy >otus.general.taxonomy.reformatted
mv otus.general.taxonomy.reformatted otus.general.taxonomy
sed 's/[[:blank:]]/\;/' <custom.general.taxonomy >custom.general.taxonomy.reformatted
mv custom.general.taxonomy.reformatted custom.general.taxonomy
sed 's/[[:blank:]]/\;/' <custom.taxonomy >custom.custom.taxonomy
mkdir conflicts_${pident[0]}
Rscript find_classification_disagreements.R otus.${pident[0]}.taxonomy otus.general.taxonomy ids.above.${pident[0]} conflicts_${pident[0]} ${pident[0]} 70
mkdir conflicts_database
Rscript find_classification_disagreements.R custom.custom.taxonomy custom.general.taxonomy NA conflicts_database NA 70 database 

# Next Run steps 4-9, the first half of 11, and 12 with different pident cutoffs
# Define a function called runagain since you repeat this part many times

runagain () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.$1 $1 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.$1 $1 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified $1 plots
   cat ids.below.$1 ids.missing > ids.below.$1.all
   python create_fastas_given_seqIDs.py ids.above.$1 otus.fasta otus.above.$1.fasta
   python create_fastas_given_seqIDs.py ids.below.$1.all otus.fasta otus.below.$1.fasta
   mothur "#classify.seqs(fasta=otus.above.$1.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.$1.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.$1.custom.wang.taxonomy otus.below.$1.general.wang.taxonomy > otus.$1.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.$1.taxonomy >otus.$1.taxonomy.reformatted
   mv otus.$1.taxonomy.reformatted otus.$1.taxonomy
   mkdir conflicts_$1
   Rscript find_classification_disagreements.R otus.$1.taxonomy otus.general.taxonomy ids.above.$1 conflicts_$1 $1 70
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

# Now run step 13- plotting everything together to choose final pident

args_string=""
for p in ${pident[*]}
do
   args_string+=" conflicts_$p ids.above.$p $p"
done

Rscript plot_classification_disagreements.R otus.abund plots conflicts_database $args_string

printf 'Steps 1-13 have finished running.  Now analysze the plots from step 13 to choose your final pident and generate your final taxonomy file in step 14.  Then tidy up your working directory with step 15. \n \a'
sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; sleep .1; printf '\a'; 

exit 0
