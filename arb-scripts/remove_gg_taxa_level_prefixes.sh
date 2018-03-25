#!/bin/bash

# RRR Mar 25, 2018
# Script to remove the p__ headers from taxa names in FreshTrain so that they are compatible with SILVA

withprefix="$1"
reformatted="$2"

sed 's/k__//g' <$withprefix | \
sed 's/p__//g' | \
sed 's/c__//g' | \
sed 's/o__//g' >$reformatted

exit


