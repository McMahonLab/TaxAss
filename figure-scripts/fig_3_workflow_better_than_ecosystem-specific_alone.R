# RRR ----

# fig 3 demonstrates that the workflow is better than using only an 
# ecosystem-specific database by showing the extent of forcing.

# fig 3a shows how increased unclassifieds decrease alpha diversity
# fig 3b shows how forcing disproportionally affects your fine-level taxa analysis
# The two panels are the same plot at different taxa levels.

# How I define forcing in fig 3: baseline is TaxAss (ie "correct") rank abundance
# 1. classify dataset with TaxAss
# 2. classify dataset with only FreshTrain
# 3. classifications that disagree were "forced" when using FreshTrain alone
#     - loss of diversity of forced into "unclassified"
#     - added inaccuracy if forced into confident classification

# TaxAss generates the data fed into this plot in optional step _____.
# It has already been filtered to only include the ___ most abundant OTUs.

# ---- file paths ----

# note: the exported csv's are already filtered to be only the top taxa. cutoff = .5 % rel abundance (max bar any color) I think.

mendota.forcing.folder <- "~/Desktop/TaxonomyTrainingSets/BLASTing/ME_GG/analysis/plots/Forcing_csvs/"

danube.forcing.folder <- "~/Desktop/TaxonomyTrainingSets/BLASTing/danube/plots/"

mouse.forcing.folder <- "~/Desktop/TaxonomyTrainingSets/BLASTing/mouse2/analysis/plots/"

# choose one:
forcing.folder <- mendota.forcing.folder

forcing.files <- paste(forcing.folder, list.files(path = forcing.folder), sep = "")

# ---- define functions ----

import.forcing.files <- function(FilePaths){
  forcing <- list(NULL)
  for(t in 1:length(FilePaths)){
    forcing[[t]] <- read.csv(file = FilePaths[t], header = TRUE, colClasses = "character")
    forcing[[t]][ ,-(1:t)] <- apply(X = as.matrix(forcing[[t]][ ,-(1:t)]), MARGIN = 2, FUN = as.numeric)
  }
  names(forcing) <- c("kingdom", "phylum", "class", "order", "lineage", "clade", "tribe")
  return(forcing)
}

# modified very slightly from plot_classification_disagreement's plot.forcing.diffs:
# diffs are: split into three funcitons, not exporting the plot automatically

extract.plot.data <- function(TopTaxaList, NumBars = 800, PlottingLevels = 1:length(TopTaxaList)){
  top.taxa <- TopTaxaList
  stacked.data <- list(NULL) 
  
  # plot specified number of bars or fewer if that many don't exist
  for (t in PlottingLevels){
    num.bars <- NumBars
    if (nrow(top.taxa[[t]]) < num.bars){
      num.bars <- nrow(top.taxa[[t]])
    }
    
    # Pull out only red/blue/grey bar heights and the lowest-level taxa name
    stacked.data[[t]] <- top.taxa[[t]][1:num.bars, c(t + 4, t + 5, t + 6)]
    row.names(stacked.data[[t]]) <- top.taxa[[t]][1:num.bars,t]
    stacked.data[[t]] <- as.matrix(stacked.data[[t]])
    stacked.data[[t]] <- t(stacked.data[[t]])
    names(stacked.data)[t] <- names(top.taxa)[t]
  }
  return(stacked.data)
}

shorten.taxa.names <- function(PlotData, PlottingLevels = 1:length(PlotData), MaxLettersBarLabels = 80){
  stacked.data <- PlotData

  for (t in PlottingLevels){
    # remove p__ from beginning of taxa names
    taxa.names <- sub(pattern = ".*__", replacement = "", x = colnames(stacked.data[[t]]))
    
    # shorten names to a max number of characters- looks inconsistent b/c not monospaced font
    max.length <- MaxLettersBarLabels
    for (n in 1:length(taxa.names)){
      if (nchar(taxa.names[n], type = "chars") > max.length){
        taxa.names[n] <- substr(x = taxa.names[n], start = 1, stop = max.length - 1)
        taxa.names[n] <- paste(taxa.names[n], "-", sep = "")
      }
    }
    
    # replace old names in plotting data
    colnames(stacked.data[[t]]) <- taxa.names
  }
  return(stacked.data)
}

plot.forcing.diffs <- function(PlotData, FolderPath = NULL, PlottingLevels = 1:length(PlotData)){  
  stacked.data <- PlotData
  for (t in PlottingLevels){
    # # adjust y-axis to go to the max value
    # max.bar <- max(stacked.data[[t]][1, ] + stacked.data[[t]][2, ]) # grey bar + red bar
    # y.axis.ticks <- c(0, max.bar * (1/4), max.bar * (1/2), max.bar * (3/4), max.bar)
    # y.axis.labels <- round(x = y.axis.ticks, digits = 0)
    
    # export plot if folder specified
    if (!is.null(FolderPath)){
      plot.name <- paste("fig_3__", t, "_", names(stacked.data)[t], "_forcing.png", sep = "")
      plot.name <- paste(FolderPath, "/", plot.name, sep = "")
      png(filename = plot.name, width = 7, height = 5, units = "in", res = 100)
    }
    
    # the plot
    par(mar = c(10,5,5,2))
    barplot(height = stacked.data[[t]], beside = FALSE, col = c("grey","red","blue"), main = names(stacked.data)[t], las = 2, ylab = "Relative Abundance (% reads)", border = NA)
    legend(x = "topright", legend = c("Prevented Inaccuracy", "Maintained Diversity"), fill = c("red", "blue"), border = FALSE, bty = "n", inset = .05)
    
    # finish exporting plot if folder specified
    if (!is.null(FolderPath)){
      unnecessary.message <- dev.off()
      cat("made plot: ", plot.name, "\n")
    }
  }
}


# ---- use functions for quick look at all taxa levels ----

forcing.data <- import.forcing.files(FilePaths = forcing.files)

forcing.data <- extract.plot.data(TopTaxaList = forcing.data)

forcing.data <- shorten.taxa.names(PlotData = forcing.data, MaxLettersBarLabels = 30)

plot.forcing.diffs(PlotData = forcing.data)

# ---- PAPER ----

plot.forcing.diffs(PlotData = forcing.data, FolderPath = "~/Dropbox/PhD/Write It/draft 3/draft_3_figure_files/")

# ---- POSTER ----
# note: ran scripts from plot_classification_disagreements.R to re-generate data to plot,
# but now plot_classification_disagreements.R is updated to export a nice csv of the data
# so I deleted those no-longer-used scripts from this plotting script.

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
