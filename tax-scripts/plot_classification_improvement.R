# RRR 1/26/16



taxonomy.pvalues.path <- "../../take9c/final.taxonomy.pvalues"
reads.table.path <- "../../take9c/total.reads.per.seqID"
gg.pvalues.path <- "../../take9c/final.general.pvalues"


#####
# Define functions to import data
#####

import.pvalues <- function(FilePath){
  taxonomy.pvalues.path <- FilePath
  pvals <- read.table(file = taxonomy.pvalues.path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  pvals[,2:8] <- apply(X = pvals[,2:8], MARGIN = 2, FUN = as.numeric)
  pvals[,1] <- as.character(pvals[,1])
  return(pvals)
}

import.read.counts <- function(FilePath){
  reads.table.path <- FilePath
  reads <- read.table(file = reads.table.path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  reads[,2] <- as.numeric(reads[,2])
  reads[,1] <- as.character(reads[,1])
  return(reads)
}


#####
# Define functions to manipulate the data
#####

order.by.seqID.and.combine.tables <- function(ReadsTable, PvaluesTable){
  reads <- ReadsTable
  pvals <- PvaluesTable
  
  # first check there are no duplicated names- can't have any ties in the order!
  if (length(pvals$seqID) != length(unique(pvals$seqID)) | length(reads$seqID) != length(unique(reads$seqID)) | length(pvals$seqID) != length(reads$seqID)){
    cat("uh oh major problem, the number of seqIDs is different in your two tables!")
  }
  
  # this orders numbers as characters, so 1, 10, 100.  oh well, only way to know it will work with character names for seqIDs
  index.reads <- order(reads[ ,1])
  index.pvals <- order(pvals[ ,1])
  reads <- reads[index.reads, ]
  pvals <- pvals[index.pvals, ]
  
  # check that they are in the same order now
  if (all.equal(reads[,1], pvals[,1])){
    cat("orders match- good.")
  }else{
    cat("crap something's wrong, they didn't end up in the same order.  must fix.")
  }
  
  pvals.with.reads <- cbind(reads, pvals[,2:8])
  
  return(pvals.with.reads)
}

convert.to.name.presence.absence <- function(PvaluesTable){
  pvals <- PvaluesTable
  
  pvals[ ,3:9] <- pvals[ ,3:9] > 0
  
  return(pvals)
}

convert.to.reads.presence.absence <- function(TrueFalseTable){
  pvals <- TrueFalseTable
  
  for (c in 3:9){
    pvals[ ,c] <- pvals[ ,c] * pvals[ ,2]
  }
  
  return(pvals)
}


#####
# Define functions to plot the data
#####

plot.num.classified


#####
# Use functions!
#####

fw.pvals <- import.pvalues(FilePath = taxonomy.pvalues.path)
gg.pvals <- import.pvalues(FilePath = gg.pvalues.path)

reads <- import.read.counts(FilePath = reads.table.path)

fw.pvals <- order.by.seqID.and.combine.tables(ReadsTable = reads, PvaluesTable = fw.pvals)
gg.pvals <- order.by.seqID.and.combine.tables(ReadsTable = reads, PvaluesTable = gg.pvals)

otus.named.fw <- convert.to.name.presence.absence(PvaluesTable = fw.pvals)
otus.named.gg <- convert.to.name.presence.absence(PvaluesTable = gg.pvals)

reads.named.fw <- convert.to.reads.presence.absence(TrueFalseTable = otus.named.fw)
reads.named.gg <- convert.to.reads.presence.absence(TrueFalseTable = otus.named.gg)







