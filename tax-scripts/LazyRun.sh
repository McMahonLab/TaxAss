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
# Define functions for these repeated parts so that each block can run in paralelle

Pident99 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.99 99 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.99 99 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 99 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.99 ids.missing > ids.below.99.all
   python create_fastas_given_seqIDs.py ids.above.99 otus.fasta otus.above.99.fasta
   python create_fastas_given_seqIDs.py ids.below.99.all otus.fasta otus.below.99.fasta
   mothur "#classify.seqs(fasta=otus.above.99.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.99.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.99.custom.wang.taxonomy otus.below.99.general.wang.taxonomy > otus.99.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.99.taxonomy >otus.99.taxonomy.reformatted
   mv otus.99.taxonomy.reformatted otus.99.taxonomy
   mkdir conflicts_99
   Rscript find_classification_disagreements.R otus.99.taxonomy otus.general.taxonomy ids.above.99 conflicts_99 99 70
}

Pident98 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.98 98 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.98 98 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 98 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.98 ids.missing > ids.below.98.all
   python create_fastas_given_seqIDs.py ids.above.98 otus.fasta otus.above.98.fasta
   python create_fastas_given_seqIDs.py ids.below.98.all otus.fasta otus.below.98.fasta
   mothur "#classify.seqs(fasta=otus.above.98.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.98.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.98.custom.wang.taxonomy otus.below.98.general.wang.taxonomy > otus.98.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.98.taxonomy >otus.98.taxonomy.reformatted
   mv otus.98.taxonomy.reformatted otus.98.taxonomy
   mkdir conflicts_98
   Rscript find_classification_disagreements.R otus.98.taxonomy otus.general.taxonomy ids.above.98 conflicts_98 98 70
}

Pident97 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.97 97 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.97 97 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 97 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.97 ids.missing > ids.below.97.all
   python create_fastas_given_seqIDs.py ids.above.97 otus.fasta otus.above.97.fasta
   python create_fastas_given_seqIDs.py ids.below.97.all otus.fasta otus.below.97.fasta
   mothur "#classify.seqs(fasta=otus.above.97.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.97.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.97.custom.wang.taxonomy otus.below.97.general.wang.taxonomy > otus.97.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.97.taxonomy >otus.97.taxonomy.reformatted
   mv otus.97.taxonomy.reformatted otus.97.taxonomy
   mkdir conflicts_97
   Rscript find_classification_disagreements.R otus.97.taxonomy otus.general.taxonomy ids.above.97 conflicts_97 97 70
}

Pident96 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.96 96 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.96 96 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 96 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.96 ids.missing > ids.below.96.all
   python create_fastas_given_seqIDs.py ids.above.96 otus.fasta otus.above.96.fasta
   python create_fastas_given_seqIDs.py ids.below.96.all otus.fasta otus.below.96.fasta
   mothur "#classify.seqs(fasta=otus.above.96.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.96.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.96.custom.wang.taxonomy otus.below.96.general.wang.taxonomy > otus.96.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.96.taxonomy >otus.96.taxonomy.reformatted
   mv otus.96.taxonomy.reformatted otus.96.taxonomy
   mkdir conflicts_96
   Rscript find_classification_disagreements.R otus.96.taxonomy otus.general.taxonomy ids.above.96 conflicts_96 96 70
}

Pident95 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.95 95 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.95 95 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 95 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.95 ids.missing > ids.below.95.all
   python create_fastas_given_seqIDs.py ids.above.95 otus.fasta otus.above.95.fasta
   python create_fastas_given_seqIDs.py ids.below.95.all otus.fasta otus.below.95.fasta
   mothur "#classify.seqs(fasta=otus.above.95.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.95.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.95.custom.wang.taxonomy otus.below.95.general.wang.taxonomy > otus.95.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.95.taxonomy >otus.95.taxonomy.reformatted
   mv otus.95.taxonomy.reformatted otus.95.taxonomy
   mkdir conflicts_95
   Rscript find_classification_disagreements.R otus.95.taxonomy otus.general.taxonomy ids.above.95 conflicts_95 95 70
}

Pident94 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.94 94 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.94 94 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 94 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.94 ids.missing > ids.below.94.all
   python create_fastas_given_seqIDs.py ids.above.94 otus.fasta otus.above.94.fasta
   python create_fastas_given_seqIDs.py ids.below.94.all otus.fasta otus.below.94.fasta
   mothur "#classify.seqs(fasta=otus.above.94.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.94.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.94.custom.wang.taxonomy otus.below.94.general.wang.taxonomy > otus.94.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.94.taxonomy >otus.94.taxonomy.reformatted
   mv otus.94.taxonomy.reformatted otus.94.taxonomy
   mkdir conflicts_94
   Rscript find_classification_disagreements.R otus.94.taxonomy otus.general.taxonomy ids.above.94 conflicts_94 94 70
}

Pident93 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.93 93 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.93 93 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 93 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.93 ids.missing > ids.below.93.all
   python create_fastas_given_seqIDs.py ids.above.93 otus.fasta otus.above.93.fasta
   python create_fastas_given_seqIDs.py ids.below.93.all otus.fasta otus.below.93.fasta
   mothur "#classify.seqs(fasta=otus.above.93.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.93.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.93.custom.wang.taxonomy otus.below.93.general.wang.taxonomy > otus.93.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.93.taxonomy >otus.93.taxonomy.reformatted
   mv otus.93.taxonomy.reformatted otus.93.taxonomy
   mkdir conflicts_93
   Rscript find_classification_disagreements.R otus.93.taxonomy otus.general.taxonomy ids.above.93 conflicts_93 93 70
}

Pident92 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.92 92 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.92 92 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 92 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.92 ids.missing > ids.below.92.all
   python create_fastas_given_seqIDs.py ids.above.92 otus.fasta otus.above.92.fasta
   python create_fastas_given_seqIDs.py ids.below.92.all otus.fasta otus.below.92.fasta
   mothur "#classify.seqs(fasta=otus.above.92.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.92.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.92.custom.wang.taxonomy otus.below.92.general.wang.taxonomy > otus.92.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.92.taxonomy >otus.92.taxonomy.reformatted
   mv otus.92.taxonomy.reformatted otus.92.taxonomy
   mkdir conflicts_92
   Rscript find_classification_disagreements.R otus.92.taxonomy otus.general.taxonomy ids.above.92 conflicts_92 92 70
}

Pident91 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.91 91 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.91 91 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 91 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.91 ids.missing > ids.below.91.all
   python create_fastas_given_seqIDs.py ids.above.91 otus.fasta otus.above.91.fasta
   python create_fastas_given_seqIDs.py ids.below.91.all otus.fasta otus.below.91.fasta
   mothur "#classify.seqs(fasta=otus.above.91.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.91.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.91.custom.wang.taxonomy otus.below.91.general.wang.taxonomy > otus.91.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.91.taxonomy >otus.91.taxonomy.reformatted
   mv otus.91.taxonomy.reformatted otus.91.taxonomy
   mkdir conflicts_91
   Rscript find_classification_disagreements.R otus.91.taxonomy otus.general.taxonomy ids.above.91 conflicts_91 91 70
}

Pident90 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.90 90 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.90 90 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 90 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.90 ids.missing > ids.below.90.all
   python create_fastas_given_seqIDs.py ids.above.90 otus.fasta otus.above.90.fasta
   python create_fastas_given_seqIDs.py ids.below.90.all otus.fasta otus.below.90.fasta
   mothur "#classify.seqs(fasta=otus.above.90.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.90.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.90.custom.wang.taxonomy otus.below.90.general.wang.taxonomy > otus.90.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.90.taxonomy >otus.90.taxonomy.reformatted
   mv otus.90.taxonomy.reformatted otus.90.taxonomy
   mkdir conflicts_90
   Rscript find_classification_disagreements.R otus.90.taxonomy otus.general.taxonomy ids.above.90 conflicts_90 90 70
}

Pident88 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.88 88 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.88 88 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 88 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.88 ids.missing > ids.below.88.all
   python create_fastas_given_seqIDs.py ids.above.88 otus.fasta otus.above.88.fasta
   python create_fastas_given_seqIDs.py ids.below.88.all otus.fasta otus.below.88.fasta
   mothur "#classify.seqs(fasta=otus.above.88.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.88.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.88.custom.wang.taxonomy otus.below.88.general.wang.taxonomy > otus.88.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.88.taxonomy >otus.88.taxonomy.reformatted
   mv otus.88.taxonomy.reformatted otus.88.taxonomy
   mkdir conflicts_88
   Rscript find_classification_disagreements.R otus.88.taxonomy otus.general.taxonomy ids.above.88 conflicts_88 88 70
}

Pident86 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.86 86 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.86 86 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 86 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.86 ids.missing > ids.below.86.all
   python create_fastas_given_seqIDs.py ids.above.86 otus.fasta otus.above.86.fasta
   python create_fastas_given_seqIDs.py ids.below.86.all otus.fasta otus.below.86.fasta
   mothur "#classify.seqs(fasta=otus.above.86.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.86.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.86.custom.wang.taxonomy otus.below.86.general.wang.taxonomy > otus.86.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.86.taxonomy >otus.86.taxonomy.reformatted
   mv otus.86.taxonomy.reformatted otus.86.taxonomy
   mkdir conflicts_86
   Rscript find_classification_disagreements.R otus.86.taxonomy otus.general.taxonomy ids.above.86 conflicts_86 86 70
}

Pident84 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.84 84 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.84 84 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 84 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.84 ids.missing > ids.below.84.all
   python create_fastas_given_seqIDs.py ids.above.84 otus.fasta otus.above.84.fasta
   python create_fastas_given_seqIDs.py ids.below.84.all otus.fasta otus.below.84.fasta
   mothur "#classify.seqs(fasta=otus.above.84.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.84.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.84.custom.wang.taxonomy otus.below.84.general.wang.taxonomy > otus.84.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.84.taxonomy >otus.84.taxonomy.reformatted
   mv otus.84.taxonomy.reformatted otus.84.taxonomy
   mkdir conflicts_84
   Rscript find_classification_disagreements.R otus.84.taxonomy otus.general.taxonomy ids.above.84 conflicts_84 84 70
}

Pident82 () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.82 82 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.82 82 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified 82 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
   cat ids.below.82 ids.missing > ids.below.82.all
   python create_fastas_given_seqIDs.py ids.above.82 otus.fasta otus.above.82.fasta
   python create_fastas_given_seqIDs.py ids.below.82.all otus.fasta otus.below.82.fasta
   mothur "#classify.seqs(fasta=otus.above.82.fasta, template=custom.fasta,  taxonomy=custom.taxonomy, method=wang, probs=T, processors=2)"
   mothur "#classify.seqs(fasta=otus.below.82.fasta, template=general.fasta, taxonomy=general.taxonomy, method=wang, probs=T, processors=2)"
   cat otus.above.82.custom.wang.taxonomy otus.below.82.general.wang.taxonomy > otus.82.taxonomy
   sed 's/[[:blank:]]/\;/' <otus.82.taxonomy >otus.82.taxonomy.reformatted
   mv otus.82.taxonomy.reformatted otus.82.taxonomy
   mkdir conflicts_82
   Rscript find_classification_disagreements.R otus.82.taxonomy otus.general.taxonomy ids.above.82 conflicts_82 82 70
}

pident () {
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.above.$1 $1 TRUE 
   Rscript filter_seqIDs_by_pident.R otus.custom.blast.table otus.custom.blast.table.modified ids.below.$1 $1 FALSE
   RScript plot_blast_hit_stats.R otus.custom.blast.table.modified $1 plots
   python find_seqIDs_blast_removed.py otus.fasta otus.custom.blast.table ids.missing
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

# Now Run all the functions (other pidents) you want to in paralelle:
# just comment out any pidents you don't want to do this time

pident 80

wait

Rscript plot_classification_disagreements.R otus.abund plots conflicts_100 100 conflicts_99 99 conflicts_98 98 conflicts_97 97 conflicts_96 96 conflicts_95 95 conflicts_94 94 conflicts_93 93 conflicts_92 92 conflicts_91 91 conflicts_90 90 conflicts_88 88 conflicts_86 86 conflicts_84 84 conflicts_82 82 conflicts_80 80
