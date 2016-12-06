#!/bin/bash

# RRR 12-6-16
# This reformats the greengenes download taxonomy file into a mothur (and TaxAss) compatible format:
# This is explained in step 0 of TaxAss directions

ggfile=$"$1"

sed 's/ //g' <$ggfile >NoSpaces
sed 's/$/;/' <NoSpaces >EndLineSemicolons
mv EndLineSemicolons $ggfile
rm NoSpaces

exit

# notes to self b/c I suck at bash:
# $1 is a built in bash variable that is the first argument following the script name when sourced.
