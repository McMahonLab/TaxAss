#!/bin/bash

# First Run steps 1-12 to generate databases and folders exactly following workflow
# Note: still gotta do the reformatting on your own (step 0)

makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db
blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast -outfmt 11 -max_target_seqs 5
blast_formatter -archive otus.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.100 100 TRUE 
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.100 100 FALSE
mkdir plots
RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 100 plots
python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
cat ids.below.100 ids.missing > ids.below.100.all
python create_fastas_given_seqIDs.py ids.above.100 otus.fasta otus.above.100.fasta
python create_fastas_given_seqIDs.py ids.below.100.all otus.fasta otus.below.100.fasta
mothur "#classify.seqs(fasta=otus.above.100.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
mothur "#classify.seqs(fasta=otus.below.100.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
cat otus.above.100.custom.wang.taxonomy otus.below.100.general.wang.taxonomy > otus.100.taxonomy
mothur "#classify.seqs(fasta=otus.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
cat otus.general.wang.taxonomy > otus.general.taxonomy
sed 's/[[:blank:]]/\;/' <otus.100.taxonomy >otus.100.taxonomy.reformatted
mv otus.100.taxonomy.reformatted otus.100.taxonomy
sed 's/[[:blank:]]/\;/' <otus.general.taxonomy >otus.general.taxonomy.reformatted
mv otus.general.taxonomy.reformatted otus.general.taxonomy
mkdir conflicts_100
Rscript find_classification_disagreements.R otus.100.taxonomy otus.general.taxonomy ids.above.100 conflicts_100 100 70

# Next Run steps 4-9, the first half of 11, and 12 with different pident cutoffs
# Define a function called pident since you repeat this part many times

pident () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.$1 $1 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.$1 $1 FALSE
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

# complete the next steps for all of the pidents you include below

for i in 99 98 97 96 95 94 93 92 91 90 
do 
	pident $i &
done 

wait

# the & lets the levels of the loop run in parallel
# the wait makes sure all the loops finish before the script exits

# Now run step 13- plotting everything together to choose final pident
# for all the pidents you're comparing



Rscript plot_classification_disagreements.R otus.abund plots conflicts_100 100 conflicts_99 99 conflicts_98 98 conflicts_97 97 conflicts_96 96 conflicts_95 95 conflicts_94 94 conflicts_93 93 conflicts_92 92 conflicts_91 91 conflicts_90 90 

