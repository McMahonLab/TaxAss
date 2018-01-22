#!/bin/bash

# RRR 1/17/18
# This is to address reviewer 2's comments
# This script si based on the surv fresh trim_ref_align_to_primer_region.sh script
# pair it with TaxAss_full-length_16S.R

# This trims a fasta file of full-length 16S sequences to the primers in the oligos file.
# It uses mothur and file formats match the mothur program requirements

# syntax:
# ./trim_full-length_to_primer_region.sh full-length.fasta primers.oligos ecoli.fasta

mothur "#pcr.seqs(fasta=ecoli_16S.fasta, oligos=$oligosfile, keepprimer=T, keepdots=T)" &&
