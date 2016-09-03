# RRR 8-15-16 ----

# fig 3 demonstrates that the workflow is better than using only an 
# ecosystem-specific database by showing the extent of forcing.

# fig 3a shows how forcing can mess up your fine-level taxa analysis
# fig 3b shows how forcing can skew coarse-level taxa analyses

# How am I defining forcing?
# 1. classify everything with custom database
# 2. look at the freshwater classifications of everything that the workflow put in greengenes
# 3. ignore everything that's unclassified in freshwater classifications
# 4. find all the classifications that are different from the workflow classifications-
#    these are the ones I call "forced" into the wrong classification.

# for fig 3a: baseline is freshwater rank abundance
# 1. make a rank abundance plot of the top lineages as defined by the custom-only taxonomy
#    (this leaves out "unclassified" from the taxa names)
# 2.color the bars red for the OTUs that were not classified that way by the workflow

# for fig 3b: baseline is workflow rank abundance
# 1. make a rank abund of the top phyla in freshwater-only classifications
# 2. make a rank abund of the top phyla in workflow classifications
# 3. compare how different freshwater only is from workflow


# ---- file paths ----

mendota.unclust.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_mend_unclust/plots/ForcedTaxonomyGroups/"
mendota.unclust.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_mend_unclust/plots/FinalTaxonomyGroups/"

bogs.epi.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_epi/plots/ForcedTaxonomyGroups/"
bogs.epi.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_epi/plots/FinalTaxonomyGroups/"

bogs.hypo.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_hypo/plots/ForcedTaxonomyGroups/"
bogs.hypo.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_hypo/plots/FinalTaxonomyGroups/"

michigan.hiseq.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_michigan_hiseq/plots/ForcedTaxonomyGroups/"
michigan.hiseq.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_michigan_hiseq/plots/FinalTaxonomyGroups/"

danube.10.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_danube_10/plots/ForcedTaxonomyGroups/"
danube.10.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_danube_10/plots/FinalTaxonomyGroups/"

# choose one:
file.path.forced <- danube.10.forced
file.path.final <- danube.10.final
#


# ---- define functions ----

import.grouped.folder <- function(FolderPath){
  grouped <- list(kingdom=NULL, phylum=NULL, class=NULL, order=NULL, lineage=NULL, clade=NULL, tribe=NULL)
  files <- list.files(FolderPath)
  for (t in 1:7){
    grouped[[t]] <- read.csv(file = paste(FolderPath, files[t], sep = ""), colClasses = "character")
    grouped[[t]][ ,t+1] <- as.numeric(grouped[[t]][ ,t+1])
  }
  return(grouped)
}

# these are copied from the plot_classification_disagreements.R script

# Find the most abundant taxa at each taxonomy level - and decide if you want to include unclassifieds in your rank abund calcs
find.top.taxa.by.total.reads <- function(TaxonomyList, NumberTopTaxa = "all", RemoveUnclass){
  grouped.taxa <- TaxonomyList
  
  # trim or don't trim the total number of results
  num.taxa <- NULL
  if (NumberTopTaxa == "all"){
    for (t in 1:length(grouped.taxa)){
      num.taxa[t] <- nrow(grouped.taxa[[t]])
    }
  }else{
    num.taxa <- rep.int(x = NumberTopTaxa, times = length(grouped.taxa))
  }
  
  # arrange highest to lowest total reads
  grouped.taxa.ord <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  for (t in 1:7){
    index <- order(grouped.taxa[[t]][ ,(t + 1)], decreasing = TRUE)
    grouped.taxa.ord[[t]] <- grouped.taxa[[t]][index, ]
  }
  
  if (RemoveUnclass == TRUE){
    # remove unclassified taxa b/c those really can't be compared on the same taxa level, and likely wouldn't be included in a "top taxa" analysis anyway.
    not.unclassifieds <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
    for (t in 1:7){
      index <- grep(x = grouped.taxa.ord[[t]][ ,t], pattern =  "unclassified.*", value = FALSE )
      # this is necessary because you can not use -0 as an index
      if (length(index) != 0){
        not.unclassifieds[[t]] <- grouped.taxa.ord[[t]][-index, ]
      }else{
        not.unclassifieds[[t]] <- grouped.taxa.ord[[t]]
      }
    }
  }else{
    # name is misleading now but this is easier than chaning all the names:
    not.unclassifieds <- grouped.taxa.ord
  }
  
  # look just at the top 20 levels
  grouped.taxa.top <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  for (t in 1:7){
    if (nrow(not.unclassifieds[[t]]) < num.taxa[t]){
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]]
    }else{
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]][1:num.taxa[t], ]
    }
  }
  return(grouped.taxa.top)
}

remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

find.forcing.diffs <- function(TopFinalList, AllForcedList){
  top.final <- TopFinalList
  all.forced <- AllForcedList
  
  for (t in 1:length(top.final)){
    fw.reads <- NULL
    difference <- NULL
    for (n in 1:nrow(top.final[[t]])){
      index <- which(all.forced[[t]][ ,t] == top.final[[t]][n,t])
      if (length(index) > 0){
        temp.reads <- all.forced[[t]][index,(t+1)]
      }else{
        temp.reads <- 0
      }
      fw.reads <- c(fw.reads, temp.reads)
      difference <- c(difference, temp.reads - top.final[[t]][n,(t+1)])
    }
    top.final[[t]] <- cbind(top.final[[t]], fw.reads, difference)
  }
  
  for (t in 1:length(top.final)){
    x <- top.final[[t]]
    grey.bars <- x$reads
    red.bars <- x$difference
    blue.bars <- abs(x$difference)
    
    for (r in 1:nrow(x)){
      if (x$difference[r] > 0){
        blue.bars[r] <- 0
      }else if (x$difference[r] < 0){
        red.bars[r] <- 0
        grey.bars[r] <- grey.bars[r] - blue.bars[r]
      }
    }
    
    top.final[[t]] <- cbind(top.final[[t]], grey.bars, red.bars, blue.bars)
  }
  return(top.final)
}

# this removes taxa that are low abundance based on their max of red, grey, and blue heights. (so narrows in on interesting bars while maintaining origional ranks)
filter.out.low.abund <- function(TaxaList, CutoffVector){
  for (t in 1:7){
    max.heights <- NULL
    for (r in 1:nrow(TaxaList[[t]])){
      max.heights[r] <- max(TaxaList[[t]][r,(t + 4):(t + 6)])
    }
    index <- which(max.heights < CutoffVector[t])
    if (length(index) > 0){
      TaxaList[[t]] <- TaxaList[[t]][-index, ]
    }
  }
  return(TaxaList)
}

plot.forcing.diffs <- function(TopTaxaList, NumBars, PlottingLevels = 1:length(TopTaxaList)){
  top.taxa <- TopTaxaList
  
  # pull out and format just the data for the plot
  stacked.data <- list(NULL) 
  for (t in PlottingLevels){
    num.bars <- NumBars
    if (nrow(top.taxa[[t]]) < num.bars){
      num.bars <- nrow(top.taxa[[t]])
    }
    stacked.data[[t]] <- top.taxa[[t]][1:num.bars, c(t + 4, t + 5, t + 6)]
    row.names(stacked.data[[t]]) <- top.taxa[[t]][1:num.bars ,t]
    stacked.data[[t]] <- as.matrix(stacked.data[[t]])
    stacked.data[[t]] <- t(stacked.data[[t]])
    names(stacked.data)[t] <- names(top.taxa)[t]
  }
  
  # export plots and data!
  for (t in PlottingLevels){
    
    # # long names go off the plot
    # taxa.names <- sub(pattern = ".*__", replacement = "", x = colnames(stacked.data[[t]]))
    # taxa.names <- substr(x = taxa.names, start = 1, stop = 20)
    # 
    # # make the y axis not have crazy big numbers on it, put them in rounded percents
    # max.bar <- max(stacked.data[[t]][1, ] + stacked.data[[t]][2, ])
    # y.axis.ticks <- c(0, max.bar * (1/4), max.bar * (1/2), max.bar * (3/4), max.bar)
    # y.axis.labels <- round(x = y.axis.ticks / tot.reads * 100, digits = 0)
    
    par(mar = c(10,5,5,2))
    barplot(height = stacked.data[[t]], beside = FALSE, col = c("grey","red","blue"), main = names(stacked.data)[t], las = 2, ylab = "Relative Abundance (% reads)", border = NA)
    legend(x = "topright", legend = c("Gained from forcing", "Lost from forcing"), fill = c("red", "blue"), border = FALSE, bty = "n", inset = .05)
    
    
  }
}



# ---- use functions ----

grouped.forced.taxa <- import.grouped.folder(FolderPath = file.path.forced)

grouped.final.taxa <- import.grouped.folder(FolderPath = file.path.final)

top.final.taxa <- find.top.taxa.by.total.reads(TaxonomyList = grouped.final.taxa, NumberTopTaxa = "all", RemoveUnclass = FALSE) # all here to export all data, numbars (workflow rank abund) is determined in plot function.

top.final.taxa <- find.forcing.diffs(TopFinalList = top.final.taxa, AllForcedList = grouped.forced.taxa)

top.final.taxa <- filter.out.low.abund(TaxaList = top.final.taxa, CutoffVector = c(0, .5, .5, .5, .5, .5, .5))

plot.forcing.diffs(TopTaxaList = top.final.taxa, NumBars = 800, PlottingLevels = c(2,5))





# ---- POSTER ----
# choose mendota

mendota.unclust.forced <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_mend_unclust/plots/ForcedTaxonomyGroups/"
mendota.unclust.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_mend_unclust/plots/FinalTaxonomyGroups/"
file.path.forced <- mendota.unclust.forced
file.path.final <- mendota.unclust.final

top.taxa <- top.final.taxa
stacked.data <- list(NULL) 
for (t in c(2,5)){
  stacked.data[[t]] <- top.taxa[[t]][ ,c(t + 4, t + 5, t + 6)]
  row.names(stacked.data[[t]]) <- top.taxa[[t]][ ,t]
  stacked.data[[t]] <- as.matrix(stacked.data[[t]])
  stacked.data[[t]] <- t(stacked.data[[t]])
  names(stacked.data)[t] <- names(top.taxa)[t]
}
phylum <- stacked.data$phylum
lineage <- stacked.data$lineage

phylum.file <- "~/Dropbox/Trina/8-20-16_ISME16_figures/blue_red_forcing_phylum.png"
lineage.file <- "~/Dropbox/Trina/8-20-16_ISME16_figures/blue_red_forcing_lineage.png"

phylum.title <- expression(bold("Top Phyla in Lake Mendota"))
y.label <- expression(bold("Relative Abundance (% reads)"))
empty.names <- rep(x = "", times = ncol(phylum))
phylum.names <- sub(pattern = ".*__", replacement = "", x = colnames(phylum))
y.ax.phy <- c(0,10,20,30,40)
empty.names.lin <- rep(x = "", times = ncol(lineage))
lineage.names <- sub(pattern = ".*__", replacement = "", x = colnames(lineage))
lineage.title <- expression(bold("Top Families and Lineages in Lake Mendota"))
y.ax.lin <- c(0,5,10,15,20, 25)
legend.text <- expression(bold("Removed Inaccuracy"), bold("Restored Diversity"),bold("No Change"))


png(filename = phylum.file, width = 6.86, height = 5.92, units = "in", res = 300)
par(mar = c(7,5,2,0))
loc.labels <- barplot(height = phylum, beside = FALSE, col = c("grey","red","blue"), main = "", las = 2, ylab = "", border = NA, names.arg = empty.names, axes = F)
text(x = loc.labels, y = -2, labels = phylum.names, adj = 1, cex = 1.5, xpd = T, srt = 40)
mtext(text = phylum.title, side = 3, line = 0, at = 1.9, adj = 0, cex = 2)
axis(side = 2, at = y.ax.phy, labels = F, tick = T, line = 0, lwd = 3, xpd = T)
mtext(side = 2, at = y.ax.phy - c(0,0,0,0,.5), text = y.ax.phy, line = .5, cex = 1.5)
mtext(text = y.label, side = 2, line = 3, cex = 2, at = 40, adj = 1)
dev.off()




png(filename = lineage.file, width = 9.53, height = 5.92, units = "in", res = 300)
par(mar = c(6.5,1,1,0))
loc.labels <- barplot(height = lineage, beside = FALSE, col = c("grey","red","blue"), main = "", las = 2, ylab = "", border = NA, names.arg = empty.names.lin, axes = F)

text(x = loc.labels, y = -.5, labels = lineage.names, adj = 1, cex = 1.1, xpd = T, srt = 40)

mtext(text = lineage.title, side = 3, line = -1, at = 4.3, adj = 0, cex = 2)

axis(side = 2, at = y.ax.lin, labels = F, tick = T, line = -1, lwd = 3, xpd = T)
mtext(side = 2, at = y.ax.lin, text = y.ax.lin, line = -.5, cex = 1.5)

mtext(text = legend.text, side = 3, line = c(-6,-8,-10), at = 24.7, adj = 0, cex = 2, col = c("red","blue","grey"))

dev.off()
