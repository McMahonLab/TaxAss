#!/bin/bash

# This performs the minimum amount of processing to get a final taxonomy file.
# That means that you do not try different percent identity cutoffs to choose the best one.
# That might make sense for you if you have already made a similarity choice, for example by
# choosing a cutoff to cluster OTUs. Then just have pident match that cutoff.
# Note: this also skips the BLAST check (step 6). You could go back and just do that one.
# Note: still run step 16 to tidy up.

# Choose pident.

pident=("98")
fwbootstrap=("80")
ggbootstrap=("80")

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
# 6 skip: blast stats plot
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
# 11 skip: general-only classification
# 12
sed 's/[[:blank:]]/\;/' <otus.${pident[0]}.taxonomy >otus.${pident[0]}.taxonomy.reformatted &&
mv otus.${pident[0]}.taxonomy.reformatted otus.${pident[0]}.taxonomy &&
# 13 skip: coarse taxonomy conflict checking
# 14 skip: compare and choose a pident
# 15
Rscript find_classification_disagreements.R otus.$pident.taxonomy quickie ids.above.$pident conflicts_$pident $pident $fwbootstrap $ggbootstrap final &&


printf 'All done. If this finished without errors, RunStep_16.sh to delete intermediate files and organize the rest into convenient folders. \n \a'
sleep .2; printf '\a'; sleep .2; printf '\a'; sleep .2; printf '\a'; sleep .2; printf '\a'; sleep .2; printf '\a'; 

exit 0
