# Table 1 BLAST must be recalculated example:
# RRR 8-19-16 ----

# draft 3- this is now Table 2. did not recalculate, just copied from poster.

# find an example of 5 blast hits 
# that get re-calculated to something else 
# and the best hit is not the first one.

# ---- define file paths ----

blast.file.path <- "../../take_mendota_clust_2/otus.custom.blast.table"
blast.file.path <- "../../poster_mend_unclust/otus.custom.blast.table"

# !!! also source all the functions in calc_full_length_pident.R

# ---- dig around ----

blast <- import.BLAST.data(File = blast.file.path)
str(blast)
blast <- format.BLAST.data(BlastTable = blast)
str(blast)
blast <- calc.full.pIDs(BlastTable = blast)
str(blast)

blast.unique.ids <- choose.best.hit(BlastTable = blast, OutputFile = "~/Desktop/throwmeout.blast")
str(blast.unique.ids)

index.not.1 <- which(blast.unique.ids$hit.num.best.ids != 1)
blast.not.1 <- blast.unique.ids[index.not.1, ]
str(blast.not.1)

index.also.good <- which(blast.not.1$true.pids >= 98)
blast.also.good <- blast.not.1[index.also.good, , drop = FALSE]
str(blast.also.good)
# narrowed that down!!!

ex.seqIDs <- NULL
for (s in 1:nrow(blast.also.good)){
  ex.seqIDs <- c(ex.seqIDs, blast.also.good$qseqid)
}
str(ex.seqIDs)

index.full.table <- list(NULL)
for (n in 1:length(ex.seqIDs)){
  index.full.table[[n]] <- which(blast$qseqid == ex.seqIDs[n])
}
str(index.full.table)

ex.tables <- list(NULL)
for (s in 1:length(ex.seqIDs)){
  ex.tables[[s]] <- blast[index.full.table[[s]], ]
}

for (s in 1:length(ex.tables)){
  names(ex.tables)[s] <- ex.tables[[s]][1,1]
  ex.tables[[s]] <- ex.tables[[s]][ ,-c(1,5)]
}
ex.tables

# there they are! 

# for poster:

file.name <- paste("~/Desktop/TaxonomyTrainingSets/BLASTing/take_mendota_clust_2/plots/blast_recalc_example_seqID_", names(ex.tables)[1], ".csv", sep = "")
write.csv(x = ex.tables[[1]], file = file.name, quote = FALSE, row.names = FALSE)

# in the fig_5 plot script- start with th cyano example seqIDs that were 100%
cyano.ex.seqids
ex.seqIDs <- cyano.ex.seqids
# source the other scropts functions
# run through the above stuff
# the cyanos not that interesting b/c all teh hits the same.





