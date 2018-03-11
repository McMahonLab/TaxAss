#!/bin/bash

# RRR Mar 11, 2018

# this un-aligns the mothur-formatted silva file, and removed long taxonomies from seq names.
# For silva.nr_v132.align it goes from 10.68 GB to 318.14 MB. Dang!! save some space!
# but FYI it does take a few minutes(20-30) to run.
# the first argument is the mothur-formatted silva file downloaded from mothur website
# the second argument is a name for the new file this script creates.

# Example Syntax:
# ./un-align_silva.sh silva.nr_v132.align silva_132_unaligned.fasta

silva="$1"
smallsilva="$2"

sed 's/\.//g' <$silva >nodots
printf "removed dots \n"
sed 's/-//g' <nodots >nodashes
rm nodots
printf "removed dashes \n"
sed 's/	.*$//g' <nodashes >$smallsilva
rm nodashes
printf "removed taxonomies after seqIDs \n \a done \n \a"

exit