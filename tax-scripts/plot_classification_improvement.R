# RRR 1/26/16

# syntax for command line:

# $ Rscript plot_classification_improvement.R final.taxonomy.pvalues final.general.pvalues total.reads.per.seqID plots

#####
# Accept arguments from the command line:
#####

userprefs <- commandArgs(trailingOnly = TRUE)

taxonomy.pvalues.path <- userprefs[1]
gg.pvalues.path <- userprefs[2]
reads.table.path <- userprefs[3]
path.to.plots.folder <- userprefs[4]

# taxonomy.pvalues.path <- "../../take9c/final.taxonomy.pvalues"
# gg.pvalues.path <- "../../take9c/final.general.pvalues"
# reads.table.path <- "../../take9c/total.reads.per.seqID"
# path.to.plots.folder <- "../../take9c/plots"

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

plot.num.classified <- function(GGTable, FWTable, Reads = TRUE, FolderPath){
  ggpvals <- GGTable
  fwpvals <- FWTable
  plots.path <- FolderPath
  
  if (Reads == TRUE){
    normalizer <- sum(fwpvals[ ,2])
    title.word <- "Reads"
  }else{
    normalizer <- nrow(ggpvals)
    title.word <- "OTUs"
  }
  
  gg.classified <- colSums(ggpvals[ ,3:9]) / normalizer * 100
  fw.classified <- colSums(fwpvals[ ,3:9]) / normalizer * 100
  
  class.table <- rbind(gg.classified,fw.classified)
  
  # Save plot as .png file
  png(filename = paste(plots.path, "/Improvement_In_Classification_by_", title.word, ".png", sep = ""), 
      width = 7, height = 5, units = "in", res = 100)
  
  barplot(class.table, beside = TRUE, axes = FALSE, col = c("plum4", "orange2"), main = paste("Improvement in", title.word, "Classified"), ylab = paste("Percent total", title.word, "classified"), xlab = "Classifications Assigned at Each Taxonomic Level")
  axis(2, at = seq.int(from = 0, to = 100, by = 20), labels = seq.int(from = 0, to = 100, by = 20), xpd = T)
  legend(x = "topright", legend = c("General", "Workflow"), fill = c("plum4", "orange2"), bty = "n")
  
  unnecessary.message <- dev.off()
}


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

plot.num.classified(FWTable = otus.named.gg, GGTable = otus.named.gg, Reads = FALSE, FolderPath = path.to.plots.folder)
plot.num.classified(FWTable = reads.named.fw, GGTable = reads.named.gg, Reads = TRUE, FolderPath = path.to.plots.folder)





