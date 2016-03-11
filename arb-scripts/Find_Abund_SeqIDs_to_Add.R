# RRR 3-1-16
# This is a script to find high read seqIDs that are missing from FW database


#####
# Define File Paths and Input Variables
#####

# input
seqID.reads.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/total.reads.per.seqID.csv"
ave.perc.abund.cutoff <- .1
ids.classified.by.fw.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/ids.above.99"
taxonomy.table.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/otus.99.taxonomy"
blast.pidents.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/otus.custom.blast.table.modified"

# output
priority.taxa.file.path <- "~/Dropbox/Trina/3-3-16 priority seqIDs to add to ARB/mendota_priority_seqID_taxonomies.csv"
priority.seqIDs.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/mendota_priority_seqIDs_for_step7"

#####
# Define functions for import/export
#####

import.seqID.reads <- function(FilePath){
  seqID.reads.file.path <- FilePath
  
  seqID.reads <- read.csv(file = seqID.reads.file.path, stringsAsFactors = FALSE)
  seqID.reads[ ,1] <- as.character(seqID.reads[ ,1])
  seqID.reads[ ,2] <- as.numeric(seqID.reads[ ,2])
  
  return(seqID.reads)
}

import.fw.ids <- function(FilePath){
  ids.above.path <- FilePath
  
  # scan is only able to import numbers, throws error if values are character
  fw.ids <- read.table(file = ids.above.path, stringsAsFactors = FALSE)
  fw.ids <- fw.ids[ ,1]
  fw.ids <- as.character(fw.ids)
  
  return(fw.ids)
}

import.taxonomy.assignments <- function(FilePath){
  taxonomy.table.file.path <- FilePath
  tax.table <- read.table(file = taxonomy.table.file.path, header = FALSE, sep = ";", stringsAsFactors = FALSE)
  tax.table <- tax.table[ ,-9]
  tax.table[ ,1] <- as.character(tax.table[ ,1])
  return(tax.table)
}

import.blast.pidents <- function(FilePath){
  blast.pidents.file.path <- FilePath
  
  blast.stats <- read.table(file = blast.pidents.file.path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  seqIDs.pidents <- blast.stats[ ,c(1,6)]
  colnames(seqIDs.pidents) <- c("seqID","full.pident")
  
  return(seqIDs.pidents)
}

export.priority.taxa.for.Trina <- function(PriorityTaxaTable, FilePath){
  priority.taxa.file.path <- FilePath
  priority.taxa <- PriorityTaxaTable
  
  write.csv(x = priority.taxa, file = priority.taxa.file.path, row.names = FALSE, quote = FALSE)
  cat("\nWrote file ", priority.taxa.file.path, 
      "\nThe format is comma delimited columns: seqID, average % abundance, corrected blast pident, taxonomy assignments",
      "\nNote that average % abundance is in percents- 0.1 means 0.1 % not 10 %",
      "\nNote that taxonomy names include all assignments- only trust p-values that are at least 70 %",
      "\nNote that pident values of NA are for sequences with e-values too low for plast to return a hit")
}

export.priority.seqIDs.to.grab.fastas <- function(PrioritySeqIDsTable, FilePath){
  priority.seqIDs.file.path <- FilePath
  priority.seqIDs <- PrioritySeqIDsTable
  
  write.table(x = priority.seqIDs[ ,1], file = priority.seqIDs.file.path, row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat("\nWrote file ", priority.seqIDs.file.path, 
      "\nThis is a list of the seqIDs with high average (over all samples) relative abundance that are not similar to anything currently in the FW database.",
      "\nThe format is accepted by step 7 of the taxonomy assignment workflow, which creates a fasta file for a list of seqIDs\n")
}

#####
# Define functions for analysis
#####

find.abund.seqIDs <- function(TotReadsTable, PercAbundCutoff){
  seqID.reads <- TotReadsTable
  cutoff <- PercAbundCutoff
  
  # Change to average percent abundance (averaged over all samples)
  tot.reads <- sum(seqID.reads$reads)
  seqID.reads[ ,2] <- seqID.reads[ ,2] / tot.reads * 100
  colnames(seqID.reads)[2] <- "AvePercAbund"
  
  # Remove IDs below cutoff, removing here makes tables smaller so next step faster
  index <- order(seqID.reads[ ,2], decreasing = TRUE)
  seqID.reads <- seqID.reads[index, ]
  
  index <- which(seqID.reads[ ,2] >= cutoff)
  seqID.reads <- seqID.reads[index, ]
  
  return(seqID.reads)
}

remove.seqIDs.in.FW.already <- function(TotReadsTable, FWseqIDs){
  abund.seqIDs <- TotReadsTable
  fw.ids <- FWseqIDs
  
  # Using counters f and a: f = each fw ID, a = each abund read
  index <- NULL
  for (f in 1:length(fw.ids)){
    for (a in 1:nrow(abund.seqIDs)){
      if (fw.ids[f] == abund.seqIDs[a,1]){
        index <- c(index,a)
      }
    }
  }
  abund.seqIDs <- abund.seqIDs[-index, ]
  return(abund.seqIDs)
}

find.pidents.of.priority.seqIDs <- function(SeqIDsTable, BlastTable){
  priority.seqIDs <- SeqIDsTable
  seqIDs.pidents <- BlastTable
  
  index <- NULL
  for(s in 1:nrow(priority.seqIDs)){
    temp.index <- which(seqIDs.pidents[ ,1] == priority.seqIDs[s,1])
    index <- c(index, temp.index)
  }
  priority.pidents <- seqIDs.pidents[index, ]
  
  return(priority.pidents)
}

find.tax.assignments.of.priority.seqIDs <- function(SeqIDsTable, TaxonomyTable){
  priority.seqIDs <- SeqIDsTable
  tax.table <- TaxonomyTable
  
  index <- NULL
  for (s in 1:nrow(priority.seqIDs)){
    temp.index <- which(tax.table[ ,1] == priority.seqIDs[s,1])
    index <- c(index, temp.index)
  }
  tax.table <- tax.table[index, ]
  colnames(tax.table) <- c("seqID","kingdom","phylum","class","order","family","genus","species")
  return(tax.table)
}

combine.taxonomy.abundance.pidents <- function(TaxonomyTable, AbundanceTable, BlastTable){
  priority.seqIDs <- AbundanceTable
  priority.taxa <- TaxonomyTable
  priority.pidents <- BlastTable
  
  # the order of seqIDs must match to combine tables
  index.1 <- order(priority.seqIDs[ ,1])
  index.2 <- order(priority.taxa[ ,1])
  priority.seqIDs <- priority.seqIDs[index.1, ]
  priority.taxa <- priority.taxa[index.2, ]
  
  # the blast results do not include all seqIDs, add NA for missing pidents
  all.priority.pidents <- data.frame(seqIDs = "placeholder", pidents = 0, stringsAsFactors = FALSE)
  for (s in 1:nrow(priority.seqIDs)){
    index <- which(priority.seqIDs[s,1] == priority.pidents[ ,1])
    if(length(index) > 0){
      all.priority.pidents[s, ] <- priority.pidents[index, ]
    }else{
      all.priority.pidents[s, ] <- c(priority.seqIDs[s,1], NA)
    }
  }
  
  cat("\ndo taxonomy and read seqIDs match: ", all.equal(priority.seqIDs[ ,1], priority.taxa[ ,1]))
  cat("\ndo taxonomy and pident seqIDs match: ", all.equal(all.priority.pidents[ ,1], priority.taxa[ ,1]))
  
  # add tables together so Trina can see taxonomy, abundance, and pident all together.
  priority.taxa <- cbind(priority.seqIDs, all.priority.pidents[ ,-1, drop = FALSE], priority.taxa[ ,-1])
  index <- order(priority.taxa[ ,2], decreasing = TRUE)
  priority.taxa <- priority.taxa[index, ]
  
  return(priority.taxa)
}

#####
# Use functions
#####

# import data
seqID.reads <- import.seqID.reads(FilePath = seqID.reads.file.path)

fw.ids <- import.fw.ids(FilePath = ids.classified.by.fw.file.path)

tax.table <- import.taxonomy.assignments(FilePath = taxonomy.table.file.path)

seqIDs.pidents <- import.blast.pidents(FilePath = blast.pidents.file.path)


# process data
abund.seqIDs <- find.abund.seqIDs(TotReadsTable = seqID.reads, PercAbundCutoff = ave.perc.abund.cutoff)

priority.seqIDs <- remove.seqIDs.in.FW.already(TotReadsTable = abund.seqIDs, FWseqIDs = fw.ids)

priority.pidents <- find.pidents.of.priority.seqIDs(SeqIDsTable = priority.seqIDs, BlastTable = seqIDs.pidents)

priority.taxa <- find.tax.assignments.of.priority.seqIDs(SeqIDsTable = priority.seqIDs, TaxonomyTable = tax.table)

priority.taxa <- combine.taxonomy.abundance.pidents(TaxonomyTable = priority.taxa, AbundanceTable = priority.seqIDs, BlastTable = priority.pidents)


#export data
export.priority.taxa.for.Trina(PriorityTaxaTable = priority.taxa, FilePath = priority.taxa.file.path)

export.priority.seqIDs.to.grab.fastas(PrioritySeqIDsTable = priority.seqIDs, FilePath = priority.seqIDs.file.path)

