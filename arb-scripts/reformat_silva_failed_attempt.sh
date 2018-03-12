#!/bin/bash

# RRR 12-6-16
# This reformats the silva download fasta file into a blast-compatible fasta file (shorter simple IDs)
# and a mothur (and TaxAss) compatible taxonomy file:
# This is explained in step 0 of TaxAss directions

silvfile=$"$1"

# create a separte taxonomy file (same name as input with .taxonomy extension)
grep ">" <$silvfile >silva.taxonomy
sed 's/>//g' <silva.taxonomy >nocarrots
sed 's/ /	/' <nocarrots >firstspacetab
sed 's/.*Eukaryota.*//' <firstspacetab >noeuks
grep -v "^$" <noeuks >onlyendmessedupnow
# and you see that this is hopeless fucked b/c silva provides a bunch of other shit at the
# end of (only some) of the rows with the same delimination as btwn the ranks. damnit.
# so that means you have to either learn ARB or use the outdated mothur/qiime versions.
# fuck you silva, I'm not learning arb. You did this on purpose, didn't you? So we'd
# all have to use your stupid baby ARB. 


# create the fasta with only IDs as names (same name as input)
sed 's/ .*//g' <$silvfile >silva.fasta

exit

# notes to self b/c I suck at bash:
# $1 is a built in bash variable that is the first argument following the script name when sourced.
# in regex ^$ signifies an empty line 