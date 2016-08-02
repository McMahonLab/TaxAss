# ---- RRR 1/26/16 ----

# syntax for command line:

# $ Rscript plot_classification_improvement.R final.taxonomy.pvalues final.general.pvalues total.reads.per.seqID plots


# ---- Accept arguments from the command line: ----

# userprefs <- commandArgs(trailingOnly = TRUE)
# 
# taxonomy.pvalues.path <- userprefs[1]
# gg.pvalues.path <- userprefs[2]
# reads.table.path <- userprefs[3]
# path.to.plots.folder <- userprefs[4]

taxonomy.pvalues.path <- "../../take18playwith/final.taxonomy.pvalues"
gg.pvalues.path <- "../../take18playwith/final.general.pvalues"
reads.table.path <- "../../take18playwith/total.reads.per.seqID.csv"
path.to.plots.folder <- "../../take18playwith/plots"
taxonomy.names.path <- "../../take18playwith/final.taxonomy.names"
gg.names.path <- "../../take18playwith/final.general.names"


# ---- Define functions to import and format data ----

import.pvalues <- function(FilePath){
  pvalues.path <- FilePath
  pvals <- read.table(file = pvalues.path, header = TRUE, sep = ",", colClasses = "character")
  pvals[,2:8] <- apply(X = pvals[,2:8], MARGIN = 2, FUN = as.numeric)
  return(pvals)
}

import.read.counts <- function(FilePath){
  reads.table.path <- FilePath
  reads <- read.table(file = reads.table.path, header = TRUE, sep = ",", colClasses = "character")
  reads[,2] <- as.numeric(reads[,2])
  return(reads)
}

import.and.order.names <- function(FilePath){
  names.path <- FilePath
  tax.names <- read.table(file = names.path, header = TRUE, sep = ",", colClasses = "character")
  
  index <- order(tax.names[ ,1])
  tax.names <- tax.names[index, ]
  
  return(tax.names)
}

check.orders.match <- function(Vector1, Vector2){
  # first check there are no duplicated names- can't have any ties in the order!
  if (length(Vector1) != length(unique(Vector1)) | length(Vector2) != length(unique(Vector2)) | length(Vector1) != length(Vector2)){
    cat("uh oh major problem, the number of seqIDs is different in your two tables!")
  }
  # second check that they are in the same order
  if (all.equal(Vector1, Vector2)){
    cat("orders match- good.")
  }else{
    cat("crap something's wrong, they didn't end up in the same order.  must fix.")
  }
}

order.by.seqID.and.combine.tables <- function(ReadsTable, PvaluesTable){
  reads <- ReadsTable
  pvals <- PvaluesTable
  
  # this orders numbers as characters, so 1, 10, 100.  oh well, only way to know it will work with character names for seqIDs
  index.reads <- order(reads[ ,1])
  index.pvals <- order(pvals[ ,1])
  reads <- reads[index.reads, ]
  pvals <- pvals[index.pvals, ]
  
  check.orders.match(Vector1 = pvals$seqID, Vector2 = reads$seqID)
  
  pvals.with.reads <- cbind(reads, pvals[,2:8])
  
  return(pvals.with.reads)
}



# ---- Define functions to manipulate the data ----

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

convert.to.unchanged.true.false <- function(FWnames, GGnames, TrueFalseTable){
  fw.names <- FWnames
  gg.names <- GGnames
  otus.named.fw <- TrueFalseTable
  check.orders.match(Vector1 = fw.names$seqID, Vector2 = gg.names$seqID)
  check.orders.match(Vector1 = fw.names$seqID, Vector2 = otus.named.fw$seqID)
  
  unchanged <- otus.named.fw[ ,1:2]
  # create a T/F table- T if the names match, F if they don't
  unchanged <- cbind(unchanged, fw.names[ ,-1] == gg.names[ ,-1])
  # adjust T/F table so also F if it was unclassified
  unchanged[ ,-(1:2)] <- unchanged[ ,-(1:2)] * otus.named.fw[ ,-(1:2)]
  
  return(unchanged)
}

convert.to.renamed.true.false <- function(FWnames, GGnames, FW_TrueTalseTable, GG_TrueFalseTable){
  fw.names <- FWnames
  gg.names <- GGnames
  otus.named.fw <- FW_TrueTalseTable
  otus.named.gg <- GG_TrueFalseTable
  
  renamed <- otus.named.fw[ ,1:2]
  # create a T/F table- T if the names DON'T match, F if they match
  renamed <- cbind(renamed, fw.names[ ,-1] != gg.names[ ,-1])
  # adjust T/F table so also F if it was unclassified by either
  renamed[ ,-(1:2)] <- renamed[ ,-(1:2)] * otus.named.fw[ ,-(1:2)]
  renamed[ ,-(1:2)] <- renamed[ ,-(1:2)] * otus.named.gg[ ,-(1:2)]
  
  return(renamed)
}

convert.to.newly.named.true.false <- function(FWnames, GGnames, FW_TrueTalseTable, GG_TrueFalseTable){
  fw.names <- FWnames
  gg.names <- GGnames
  otus.named.fw <- FW_TrueTalseTable
  otus.named.gg <- GG_TrueFalseTable
  
  # create table of those NOT named by gg only- T if it was unclassified by GG
  otus.NOT.named.gg <- otus.named.gg[ ,1:2]
  otus.NOT.named.gg <- cbind(otus.NOT.named.gg, otus.named.gg[ ,-(1:2)] == FALSE)
  # the table of all given a name by FW workflow exists- T if it was classified by FW
  newnamed <- otus.named.fw[ ,1:2]
  # now adjust so T if NOT named by GG but NOW named by FW
  newnamed <- cbind(newnamed, otus.named.fw[ ,-(1:2)] * otus.NOT.named.gg[ ,-(1:2)])
  
  return(newnamed)
}

get.beside.data <- function(GGTable, FWTable, Reads = TRUE){
  ggpvals <- GGTable
  fwpvals <- FWTable
  
  if (Reads == TRUE){
    normalizer <- sum(fwpvals[ ,2])
  }else{
    normalizer <- nrow(ggpvals)
  }
  
  gg.classified <- colSums(ggpvals[ ,3:9]) / normalizer * 100
  fw.classified <- colSums(fwpvals[ ,3:9]) / normalizer * 100
  
  class.table <- rbind(gg.classified,fw.classified)
  
  return(class.table)
}

get.stacked.data <- function(Same, Better, New, Reads = TRUE){
  unchanged <- Same
  renamed <- Better
  newly.named <- New
  
  if (Reads == TRUE){
    normalizer <- sum(unchanged[ ,2]) # could be any of the 3 here
  }else{
    normalizer <- nrow(unchanged)
  }
  
  tot.unchanged <- colSums(unchanged[ ,-(1:2)]) / normalizer * 100
  tot.renamed <- colSums(renamed[ ,-(1:2)]) / normalizer * 100
  tot.newly.named <- colSums(newly.named[ ,-(1:2)]) / normalizer * 100
  
  class.table <- rbind(tot.unchanged, tot.renamed, tot.newly.named)
  return(class.table)
}

check.numbers.add.up <- function(StackedData, BesideData){
  stacked <- StackedData
  beside <- BesideData
  
  checker <- all.equal(colSums(stacked), beside[2, ])
  cat("check sums of stacked match workflow bar of beside: ", checker)
}

# Define functions to plot the data ----

# original "besides"-only plot
plot.num.classified <- function(GGTable, FWTable, Reads = TRUE, FolderPath, Tribe = FALSE, Kingdom = FALSE){
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
  
  # Don't include kingdom on plot since they really should all be the same kingdom and non-bacteria should be filtered out of the data
  if (Kingdom == FALSE){
    class.table <- class.table[ ,-1]
  }
  
  # Save plot as .png file
  png(filename = paste(plots.path, "/Improvement_In_Classification_by_", title.word, ".png", sep = ""), 
      width = 7, height = 5, units = "in", res = 100)
  
  barplot(class.table, beside = TRUE, axes = FALSE, col = c("plum4", "orange2"), main = paste("Improvement in", title.word, "Classified"), ylab = paste("Percent total", title.word, "classified"), xlab = "Classifications Assigned at Each Taxonomic Level")
  axis(2, at = seq.int(from = 0, to = 100, by = 20), labels = seq.int(from = 0, to = 100, by = 20), xpd = T)
  legend(x = "topright", legend = c("General", "Workflow"), fill = c("plum4", "orange2"), bty = "n")
  
  unnecessary.message <- dev.off()
}

# the fancy barplot has both the total classified by general and workflow as bars beside each other,
# but the total classified by the workflow is further split into a stacked chart of the type of classification
fancy.barplot <- function(BesideData, StackedData, BarSpacing){
  y <- BesideData
  z <- StackedData
  bar.space.y <- BarSpacing
  bar.width <- 1  # spacing and axis limits will determine what width 1 looks like, no need to change
  col.y <- "grey"
  col.z <- c("grey", "orange", "red")
  
  # find all values from the basic "beside" plot:
  num.sections <- ncol(y)
  bar.spots <- barplot(y, beside = TRUE, width = bar.width, space = bar.space.y, plot = FALSE)
  tot.x <- max(bar.spots) + .5 * bar.width + bar.space.y[2]
  max.y <- max(y)
  empty.labels <- rep(x = "", times = num.sections)
  
  # calculate bar and label spacing
  bar.space.beside <- c(bar.space.y[2], rep(bar.space.y[1] + bar.space.y[2] + bar.width, length.out = num.sections - 1))
  bar.space.stacked <- bar.space.y[2] + bar.width + bar.space.y[1]
  loc.labels <- bar.spots[seq(from = 1, to = length(bar.spots), by = 2)] + .5 * bar.width + .5 * bar.space.y[1]
  
  # make plots
  barplot(y[1, ], col = col.y, border = "black", beside = FALSE, width = bar.width, space = bar.space.beside, xlim = c(0, tot.x), ylim = c(0, 100), names.arg = empty.labels)
  barplot(z, add = TRUE, col = col.z, border = "black", beside = FALSE, width = bar.width, space = bar.space.stacked, xlim = c(0, tot.x), ylim = c(0, 100), names.arg = empty.labels, axes = FALSE)           
  mtext(text = colnames(z), side = 1, line = 1, at = loc.labels)
}

# these are .csv files of the data used to make the plots.
export.summary.table <- function(Summary, FilePath, PlotType, DataType){
  file.name <- paste(FilePath, "/WorkflowImprovement-", PlotType, "Data-", DataType, ".csv", sep = "")
  write.csv(x = Summary, file = file.name, quote = FALSE)
}


# ---- Use functions! ----

# get data to plot the overall improvement (beside) part of the barplot
fw.pvals <- import.pvalues(FilePath = taxonomy.pvalues.path)
gg.pvals <- import.pvalues(FilePath = gg.pvalues.path)

reads <- import.read.counts(FilePath = reads.table.path)

fw.pvals <- order.by.seqID.and.combine.tables(ReadsTable = reads, PvaluesTable = fw.pvals)
gg.pvals <- order.by.seqID.and.combine.tables(ReadsTable = reads, PvaluesTable = gg.pvals)

otus.named.fw <- convert.to.name.presence.absence(PvaluesTable = fw.pvals)
otus.named.gg <- convert.to.name.presence.absence(PvaluesTable = gg.pvals)

reads.named.fw <- convert.to.reads.presence.absence(TrueFalseTable = otus.named.fw)
reads.named.gg <- convert.to.reads.presence.absence(TrueFalseTable = otus.named.gg)

beside.data.otus <- get.beside.data(FWTable = otus.named.fw, GGTable = otus.named.gg, Reads = FALSE)
beside.data.reads <- get.beside.data(FWTable = reads.named.fw, GGTable = reads.named.gg, Reads = TRUE)

# get data to plot the type of improvement (stacked) part of the barplot
fw.names <- import.and.order.names(FilePath = taxonomy.names.path)
gg.names <- import.and.order.names(FilePath = gg.names.path)
check.orders.match(Vector1 = fw.names$seqID, Vector2 = gg.names$seqID)

otus.unchanged <- convert.to.unchanged.true.false(FWnames = fw.names, GGnames = gg.names, TrueFalseTable = otus.named.fw)
reads.unchanged <- convert.to.reads.presence.absence(TrueFalseTable = otus.unchanged)

otus.reclassified <- convert.to.renamed.true.false(FWnames = fw.names, GGnames = gg.names, FW_TrueTalseTable = otus.named.fw, GG_TrueFalseTable = otus.named.gg)
reads.reclassified <- convert.to.reads.presence.absence(TrueFalseTable = otus.reclassified)

otus.newly.named <- convert.to.newly.named.true.false(FWnames = fw.names, GGnames = gg.names, FW_TrueTalseTable = otus.named.fw, GG_TrueFalseTable = otus.named.gg)
reads.newly.named <- convert.to.reads.presence.absence(TrueFalseTable = otus.newly.named)

stacked.data.otus <- get.stacked.data(Same = otus.unchanged, Better = otus.reclassified, New = otus.newly.named, Reads = FALSE)
stacked.data.reads <- get.stacked.data(Same = reads.unchanged, Better = reads.reclassified, New = reads.newly.named, Reads = TRUE)

check.numbers.add.up(StackedData = stacked.data.otus, BesideData = beside.data.otus)
check.numbers.add.up(StackedData = stacked.data.reads, BesideData = beside.data.reads)

# Generate the Plot!

fancy.barplot(BesideData = beside.data.otus[ ,-c(1,7)], StackedData = stacked.data.otus[ ,-c(1,7)], BarSpacing = c(0,1))
fancy.barplot(BesideData = beside.data.reads[ ,-c(1,7)], StackedData = stacked.data.reads[ ,-c(1,7)], BarSpacing = c(0,1))

# Export the Data!

export.summary.table(Summary = beside.data.otus, FilePath = path.to.plots.folder, PlotType = "Beside", DataType = "OTUs")
export.summary.table(Summary = beside.data.reads, FilePath = path.to.plots.folder, PlotType = "Beside", DataType = "Reads")
export.summary.table(Summary = stacked.data.otus, FilePath = path.to.plots.folder, PlotType = "Stacked", DataType = "OTUs")
export.summary.table(Summary = stacked.data.reads, FilePath = path.to.plots.folder, PlotType = "Stacked", DataType = "Reads")


# Check the fancy plot by looking at the simple component plots

# par(mfrow = c(2,2))
# barplot(beside.data.reads, beside = T, main = "reads", ylim = c(0,100))
# barplot(beside.data.otus, beside = T, main = "otus", ylim = c(0,100))
# barplot(stacked.data.reads[ ,4:7], legend =T, ylim = c(0,100), border = NA, main = "reads")
# barplot(stacked.data.otus[ ,4:7], legend =T, ylim = c(0,100), border = NA, main = "otus")


# Legacy plot- the original simple "beside-only" barplot showing improvement.

# plot.num.classified(FWTable = otus.named.gg, GGTable = otus.named.gg, Reads = FALSE, FolderPath = path.to.plots.folder)
# plot.num.classified(FWTable = reads.named.fw, GGTable = reads.named.gg, Reads = TRUE, FolderPath = path.to.plots.folder)





