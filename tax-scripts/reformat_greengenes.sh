#!/bin/bash

# RRR 12-6-16
# This reformats the greengenes download taxonomy file into a mothur (and TaxAss) compatible format:
# It also renames the reformatted file to general.taxonomy and deletes the original
# This is explained in step 0 of TaxAss directions

ggfile=$"$1"

sed 's/ //g' <$ggfile >NoSpaces
rm $ggfile
sed 's/$/;/' <NoSpaces >EndLineSemicolons
rm NoSpaces
mv EndLineSemicolons general.taxonomy


exit

# notes to self b/c I suck at bash:
# $1 is a built in bash variable that is the first argument following the script name when sourced.
