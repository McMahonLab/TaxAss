#!/bin/bash

# RRR 1/17/18
# This is to address reviewer 2's comments
# This script si based on the surv fresh trim_ref_align_to_primer_region.sh script
# pair it with TaxAss_full-length_16S.R

# This trims a fasta file of full-length 16S sequences to the primers in the oligos file.
# It uses mothur and file formats match the mothur program requirements

# syntax:
# ./trim_full-length_to_primer_region.sh primers.oligos aligned-fasta-without-filetype name-for-trimmed.fasta

# define input files:
oligosfile="$1" # the primers
refalign="$2"   # the aligned sequences to-be-trimmed. file name without .fasta extension. ex. mara.fasta would be entered as argument "mara"
# an ecoli full length 16S sequence named exactly ecoli_16S.fasta
createdfile="$3" # a name for the trimmed reference file, like mara_v4.fasta

mothur "#pcr.seqs(fasta=ecoli_16S.fasta, oligos=$oligosfile, keepprimer=F, keepdots=T)" &&

mothur "#align.seqs(fasta=ecoli_16S.pcr.fasta, reference=$refalign.fasta)" &&

mothur "#summary.seqs(fasta=ecoli_16S.pcr.align)" &&

startnum=`awk '{print $2}' ecoli_16S.pcr.summary | sed '2q;d'`
# for mothur v1.39.5 correct for a bug in the starting position by subtracting 1
correctstart=$((startnum-1))
endnum=`awk '{print $3}' ecoli_16S.pcr.summary | sed '2q;d'`
echo $startnum $correctstart $endnum

# note this is using the corrected start num for mother 1.39.5
mothur "#pcr.seqs(fasta=$refalign.fasta, start=$correctstart, end=$endnum, keepdots=F, nomatch=reject, pdiffs=0, rdiffs=0)" &&

# note: double check that there are still 289 fasta sequences, since non-matching ones are tossed

rm mothur.*.logfile
rm *.8mer
rm ecoli_16S.pcr.align ecoli_16S.pcr.align.report ecoli_16S.pcr.fasta ecoli_16S.pcr.summary

sed 's/-//g' <$refalign.pcr.fasta >$createdfile
rm $refalign.pcr.fasta