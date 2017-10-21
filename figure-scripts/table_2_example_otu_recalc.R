# Table 1 BLAST must be recalculated example:
# RRR 8-19-16 ----

# find an example of 5 blast hits 
# that get re-calculated to something else 
# and the best hit is not the first one.

# ---- define file paths ----

all.hits.blast.file.path <- "~/Desktop/TaxAss-BatchFiles-go/Mendota/Quickie-TaxAss-Mendota/otus.custom.blast.table"

blast.file.path <- "~/Desktop/TaxAss-BatchFiles-go/Mendota/Quickie-TaxAss-Mendota/otus.custom.blast.table.modified"

# ---- define functions ----

# (copied from calc_full_length_pident.R)
import.BLAST.data <- function(File){
  blast.file.path <- File
  blast <- read.table(file = blast.file.path, sep = "\t", stringsAsFactors = F, colClasses = "character")
  colnames(blast) <- c("qseqid","pident","length","qlen","qstart","qend")
  return(blast)
}

# (copied from calc_full_length_pident.R)
format.BLAST.data <- function(BlastTable){
  blast <- BlastTable
  
  # format numbers as numberic instead of integer
  blast[,2:6] <- apply(blast[,2:6], 2, as.numeric)
  
  # add a query alignment length column 
  q.align <- blast$qend - blast$qstart + 1   # add 1 b/c start is 1 instead of 0.  for 4 aligned pairs, start=1, end=4, length=4-1+1=4
  blast <- cbind(blast,q.align)              
  
  # remove the now unnecessary qend/qstart columns
  blast <- blast[,-(5:6)]
  
  # are most of the alignments full length?
  #   cat("\nalignment length \tmean:", mean(blast$length), "\t\tmin:", min(blast$length), "\t\tmax:",max(blast$length),
  #       "\nthe query length \tmean:", mean(blast$qlen), "\t\tmin:", min(blast$qlen), "\t\tmax:",max(blast$qlen),"\n")
  
  return(blast)
}

# (copied from calc_full_length_pident.R)
calc.full.pIDs <- function(BlastTable){
  blast <- BlastTable
  
  # define function to apply to each row
  pid.fun <- function(NumericBlastMatrixRow){
    b <- NumericBlastMatrixRow
    pid <- (b[1] * b[2]) / (b[2] - b[4] + b[3]) # (pident * length) / (length - q.align + qlen)
    return(pid)
  }
  
  # create a numeric blast matrix to use the apply funciton
  b <- blast[,-1]
  b <- as.matrix(b)
  
  # find the "true" pid for the full length sequence of each blast hit
  true.pids <- apply(b, 1, pid.fun)
  
  # add this to the blast table
  blast <- cbind(blast, true.pids)
  
  return(blast)
}

# (copied from plot_blast_hit_stats.R)
import.BLAST.hits.table <- function(FilePath){
  blast.file <- FilePath
  blast <- read.table(file = blast.file, header = F, sep = "", colClasses = "character")
  colnames(blast) <- c("qseqid", "pident", "length", "qlen", "q.align", "true.pids", "hit.num")
  blast[ ,2:7] <- apply(blast[ ,2:7], 2, as.numeric)
  return(blast)
}

find.when.lower.hit.better <- function(full.table){
  index.not.1 <- which(full.table$hit.num != 1)
  blast.not.1 <- full.table[index.not.1, ]
  cat(length(index.not.1), "total (", length(index.not.1) / nrow(full.table) * 100, "% ) of blast results\nhad better recalculated percent identity when a lower hit was used.")
  return(blast.not.1)
}

select.only.high.initial.pidents <- function(blast.table, pident){
  index <- which(blast.table$pident >= pident)
  initial.good <- blast.table[index, ]
  cat(length(index), "or (", length(index) / nrow(blast.table) * 100, "% ) of the", nrow(blast.table), "blast results where 1st hit was not best\nhad an initial blast pident >", pident )
  return(initial.good)
}

select.only.low.initial.pidents <- function(blast.table, pident){
  index <- which(blast.table$pident < pident)
  initial.bad <- blast.table[index, ]
  cat(length(index), "or (", length(index) / nrow(blast.table) * 100, "% ) of the", nrow(blast.table), "blast results where 1st hit was not best\nhad an initial blast pident <", pident )
  return(initial.bad)
}

select.only.high.recalc.pidents <- function(blast.table, pident){
  index <- which(blast.table$true.pids >= pident)
  recalc.good <- blast.table[index, ]
  cat(length(index), "or (", length(index) / nrow(blast.table) * 100, "% ) of the", nrow(blast.table), "blast results where 1st hit was not best\nhad a recalculated pident >", pident )
  return(recalc.good)
}

select.only.low.recalc.pidents <- function(blast.table, pident){
  index <- which(blast.table$true.pids < pident)
  recalc.bad <- blast.table[index, ]
  cat(length(index), "or (", length(index) / nrow(blast.table) * 100, "% ) of the", nrow(blast.table), "blast results where 1st hit was not best\nhad a recalculated pident <", pident )
  return(recalc.bad)
}

pull.out.1st.hit.only <- function(full.table){
  index.1 <- which(full.table$hit.num == 1)
  blast.1 <- full.table[index.1, ]
  cat("Pulled out only the 1st hits from the all the blast results, regardless of whether they were best.")
  return(blast.not.1)
}

find.shared.seqIDs <- function(table.1, table.2){
  set.1 <- table.1$qseqid
  set.2 <- table.2$qseqid
  combined <- c(set.1, set.2)
  index <- duplicated(combined) # duplicated marks only the SECOND instance as true (why order you combined matters)
  index <- index[-(1:length(set.1))]
  table.3 <- table.2[index, ]
  return(table.3)
}

# ---- use functions ----

# look only at the best recalculated hits:

blast <- import.BLAST.hits.table(FilePath = blast.file.path)

blast.not.1 <- find.when.lower.hit.better(full.table = blast)

good.before <- select.only.high.initial.pidents(blast.table = blast.not.1, pident = 98)

bad.before <- select.only.low.initial.pidents(blast.table = blast.not.1, pident = 98)

good.after <- select.only.high.recalc.pidents(blast.table = blast.not.1, pident = 98)

bad.after <- select.only.low.recalc.pidents(blast.table = blast.not.1, pident = 98)

# look at all the hits:

blast.all <- import.BLAST.data(File = all.hits.blast.file.path)
blast.all <- format.BLAST.data(BlastTable = blast.all)
blast.all <- calc.full.pIDs(BlastTable = blast.all)

hit.1 <- pull.out.1st.hit.only(full.table = blast.all)

hit.1.bad <- select.only.low.recalc.pidents(blast.table = hit.1, pident = 98)

# find examples where: 
#   1st hit was not best hit
#   best hit recalc made pident cutoff
#   1st hit recalc didn't make pident cutoff

# There are none...
length(unique(hit.1.bad$qseqid)) # 
length(unique(good.after$qseqid))
length(unique(hit.1.bad$qseqid)) + length(unique(good.after$qseqid)) - length(unique(c(hit.1.bad$qseqid, good.after$qseqid)))
# seqIDs in hit.1.bad + seqIDs in good after are same total number if combine them and make unique- there is no overlap!

# yup still none shared
find.shared.seqIDs(table.1 = hit.1.bad, table.2 = good.after)

# look yet another way, there really is no overlap
good.after.seqids <- good.after[ ,1, drop = FALSE]
merge(x = good.after, y = hit.1.bad, by = "qseqid")


# ---- paper table ----

# RIP

# ---- poster table ----

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





