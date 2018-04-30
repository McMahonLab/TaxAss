# RRR 2018-4-30
# This script removes references from the FreshTrain that we manually looked at first.
# These references end up conflicted in SILVA because the FreshTrain gives them the same lineage
# name while silva has 2 different order names. To fix the problem, we remove half the lineage 
# from the freshtrain so that it only has one order name. This is a simple approach based on
# that fact that there are very few references for the lineages in question so they are not
# well enough flushed out to justify keeping them and changing them in the FreshTrain.
# I used the Database_Improvement workflow and classified FreshTrain against silva. Then I
# manually looked at the unique conflicts in excel and then in the taxonomy files to identify
# these sequences to remove.

# ---- define paths ----

userprefs <- commandArgs(trailingOnly = TRUE)

start.fasta.file <- userprefs[1]
start.taxonomy.file <- userprefs[2]
end.fasta.file <- userprefs[3]
end.taxonomy.file <- userprefs[4]
silva.version <- userprefs[5]

cat("\n\ncrap re-comment file paths\n\n")
start.fasta.file <- "../FreshTrain-files/FreshTrain25Jan2018Greengenes13_5/FreshTrain25Jan2018Greengenes13_5.fasta"
start.taxonomy.file <- "../FreshTrain-files/FreshTrain25Jan2018Greengenes13_5/FreshTrain25Jan2018Greengenes13_5.taxonomy"
end.fasta.file <- "~/Desktop/test.fasta"
end.taxonomy.file <- "~/Desktop/test.taxonomy"
silva.version <- "v132"


# ---- go ----

