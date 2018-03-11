#!/bin/bash

# RRR Mar 11, 2018

# this un-aligns the mothur-formatted silva file, and removed long taxonomies from seq names.
# For silva.nr_v132.align it goes from 10.68 GB to 318.14 MB. Dang!! save some space!
# but FYI it does take a few minutes(20-30) to run b/c the aligned file is so big.

# also have to remove the dot from the seqID names in the taxonomy file, since all dots 
# are removed from the fasta file that keeps the seqIDs matching btwn files.

# the first argument is the mothur-formatted silva .align file downloaded from mothur website
# the second argument is a name for the new fasta file this script creates.
# the third argument is the mothur-formatted silva .tax file downlaoded form mothur website
# the fourth argument is a name for the new taxonomy file this script creates


# Example Syntax:
# ./un-align_silva.sh silva.nr_v132.align general.fasta silva.nr_v132.tax general.taxonomy

silva="$1"
smallsilva="$2"
tax="$3"
periodlesstax="$4"

# sed 's/\.//g' <$silva >nodots
# printf "removed dots \n"
# sed 's/-//g' <nodots >nodashes
# rm nodots
# printf "removed dashes \n"
# sed 's/	.*$//g' <nodashes >$smallsilva
# rm nodashes
# printf "removed taxonomies after seqIDs \n \a fasta done \n"

sed 's/\.//' <$tax >$periodlesstax
printf "removed dots from seqID names in taxonomy file too \n \a now actually done \n"

exit