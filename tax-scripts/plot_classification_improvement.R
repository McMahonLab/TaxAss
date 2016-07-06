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

# taxonomy.pvalues.path <- "../../take13/final.taxonomy.pvalues"
# gg.pvalues.path <- "../../take13/final.general.pvalues"
# reads.table.path <- "../../take13/total.reads.per.seqID.csv"
# path.to.plots.folder <- "../../take13/plots"

#####
# Define functions to import data
#####

import.pvalues <- function(FilePath){
  taxonomy.pvalues.path <- FilePath
  pvals <- read.table(file = taxonomy.pvalues.path, header = TRUE, sep = ",", colClasses = "character")
  pvals[,2:8] <- apply(X = pvals[,2:8], MARGIN = 2, FUN = as.numeric)
  return(pvals)
}

import.read.counts <- function(FilePath){
  reads.table.path <- FilePath
  reads <- read.table(file = reads.table.path, header = TRUE, sep = ",", colClasses = "character")
  reads[,2] <- as.numeric(reads[,2])
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

plot.num.classified <- function(GGTable, FWTable, Reads = TRUE, FolderPath, Tribe = FALSE){
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
  
  # Don't include tribe on plot if you know tribe-level assignments are inaccurate
  if (Tribe == FALSE){
    class.table <- class.table[ ,-7]
  }
  
  
  # Save plot as .png file
  png(filename = paste(plots.path, "/Improvement_In_Classification_by_", title.word, ".png", sep = ""), 
      width = 7, height = 5, units = "in", res = 100)
  
  barplot(class.table, beside = TRUE, axes = FALSE, col = c("plum4", "orange2"), main = paste("Improvement in", title.word, "Classified"), ylab = paste("Percent total", title.word, "classified"), xlab = "Classifications Assigned at Each Taxonomic Level")
  axis(2, at = seq.int(from = 0, to = 100, by = 20), labels = seq.int(from = 0, to = 100, by = 20), xpd = T)
  legend(x = "topright", legend = c("General", "Workflow"), fill = c("plum4", "orange2"), bty = "n")
  
  unnecessary.message <- dev.off()
}

fancy.barplot <- function(BesideData, StackedData, BarSpacing){
  y <- BesideData
  z <- StackedData
  bar.space.y <- BarSpacing
  bar.width <- 1  # spacing and axis limits will determine what width 1 looks like, no need to change
  
  # find all values from the basic "beside" plot:
  num.sections <- ncol(y)
  bar.spots <- barplot(y, col = rainbow(2, s = .3), beside = TRUE, width = bar.width, space = bar.space.y, plot = FALSE)
  tot.x <- max(bar.spots) + .5 * bar.width + bar.space.y[2]
  max.y <- max(y)
  
  # calculate bar spacing
  bar.space.beside <- c(bar.space.y[2], rep(bar.space.y[1] + bar.space.y[2] + bar.width, length.out <- num.sections - 1))
  bar.space.stacked <- bar.space.y[2] + bar.width + bar.space.y[1]
  
  # make plots
  barplot(y[1, ], col = rainbow(1, s = .3), border = "black", beside = FALSE, width = bar.width, space = bar.space.beside, xlim = c(0, tot.x), ylim = c(0, max.y))
  barplot(z, add = TRUE, col = rainbow(3, s = .3), border = "black", beside = FALSE, width = bar.width, space = bar.space.stacked, xlim = c(0, tot.x), ylim = c(0, max.y))
  
  # Uncomment this to make the basic "beside" plot to check against- sets up empty bars and you watch them fill correctly:
  # bar.spots <- barplot(y, col = rainbow(2, s = 0), beside = TRUE, width = bar.width, space = bar.space.y, xlim = c(0, tot.x), ylim = c(0, max.y))
  # barplot(y[1, ], add = TRUE, col = rainbow(1, s = .3), border = NA, beside = FALSE, width = bar.width, space = bar.space.beside, xlim = c(0, tot.x), ylim = c(0, max.y))
  # barplot(z, add = TRUE, col = rainbow(3, s = .3), border = NA, beside = FALSE, width = bar.width, space = bar.space.stacked, xlim = c(0, tot.x), ylim = c(0, max.y))
}
# fancy.barplot(BesideData = y, StackedData = z, BarSpacing = c(.1,1))

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





