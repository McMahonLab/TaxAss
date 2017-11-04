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

# TaxAss generates the data fed into this plot in optional step 15.5.b.
# It has already been filtered to only include bars with a max height (grey + blue + red) > .5 % abund
# The paper figure further caps the total # bars to 20 (32 > .5% at lineage)

# ---- file paths ----

# note: the exported csv's are already filtered to be only the top taxa. cutoff = .5 % rel abundance (max bar any color) I think.

mendota.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/Mendota/TaxAss-Mendota/analysis/plots/step_15_5b_Improvement_over_custom-only/"

danube.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/Danube/TaxAss-Danube/analysis/plots/step_15_5b_Improvement_over_custom-only/"

mouse.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/MouseGut/TaxAss-MouseGut/analysis/plots/step_15_5b_Improvement_over_custom-only/"

michigan.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/Michigan/TaxAss-Michigan/analysis/plots/step_15_5b_Improvement_over_custom-only/"

bog.epi.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/TroutBogEpi/TaxAss-TroutBogEpi/analysis/plots/step_15_5b_Improvement_over_custom-only/"

bog.hypo.forcing.folder <- "~/Desktop/TaxAss-BatchFiles-go/TroutBogHypo/TaxAss-TroutBogHypo/analysis/plots/step_15_5b_Improvement_over_custom-only/"

# choose one:
forcing.folder <- mendota.forcing.folder




# ---- define functions ----

get.csv.filenames <- function(Folder){
  csv.files <- list.files(path = Folder)
  index <- grep(pattern = "*.csv", x = csv.files)
  forcing.files <- paste(Folder, csv.files[index], sep = "")
  forcing.files <- forcing.files[1:7] # don't include additional summary files
  return(forcing.files)
}

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

forcing.files <- get.csv.filenames(Folder = forcing.folder)

forcing.data <- import.forcing.files(FilePaths = forcing.files)

forcing.data <- extract.plot.data(TopTaxaList = forcing.data)

forcing.data <- shorten.taxa.names(PlotData = forcing.data, MaxLettersBarLabels = 30)

plot.forcing.diffs(PlotData = forcing.data)

# ---- PAPER ----

save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Figure_3.pdf"
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .27)) # bottom, left, top, right

# ---- 3a ----
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40)
y.tick.labels <- c(0,"", 10, "", 20,"",30,"",40)
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 40

# stays same for both panels ----
bar.labels <- colnames(panel.data)
empty.labels <- rep("", length(bar.labels))
bar.colors <- c("grey", "red", "blue")
loc.labels <- barplot(height = panel.data, beside = F, col = bar.colors, las = 2, border = NA, names.arg = empty.labels, axes = F, ylim = c(0,y.max))

# X axis----
text(x = loc.labels - .1, y = -1, labels = bar.labels, srt = -30, xpd = NA, cex = .95, adj = 0)
# Y axis
axis(side = 2, at = y.axis.ticks, labels = F, line = -.2, xpd = T, tck = -.02)
mtext(text = y.tick.labels, side = 2, at = y.axis.ticks, line = .2, las = 1, cex = .7)
# Title
mtext(text = plot.title, side = 3, line = .8, cex = 1, at = .5, adj = 0)
# Y label
mtext(text = y.axis.label, side = 2, line = 2.1, cex = 1)
# 
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)

# ---- 3b ----
panel.data <- forcing.data$lineage
panel.data <- panel.data[ ,1:20]
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30)
y.tick.labels <- c(0,"",10,"",20,"",30)

# stays same for both panels ----
bar.labels <- colnames(panel.data)
empty.labels <- rep("", length(bar.labels))
bar.colors <- c("grey", "red", "blue")
loc.labels <- barplot(height = panel.data, beside = F, col = bar.colors, las = 2, border = NA, names.arg = empty.labels, axes = F)

# X axis----
lab.y <- c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1) # put lowercase ones higher
text(x = loc.labels - .1, y = -.9, labels = bar.labels, srt = -30, xpd = NA, cex = .9, adj = 0)
# Y axis
axis(side = 2, at = y.axis.ticks, labels = F, line = -.75, xpd = T, tck = -.02)
mtext(text = y.tick.labels, side = 2, at = y.axis.ticks, line = -.35, las = 1, cex = .7)
# Title
mtext(text = plot.title, side = 3, line = .8, cex = 1, at = 3, adj = 0)
# ----
legend.text <- c("Prevented Inaccuracy", "Maintained Richness")
text(x = 17, y = c(20.8,17.8), labels = legend.text, adj = 0, xpd = NA, cex = 1.1)
rect(xleft = 15.5, xright = 16.5, ybottom = 20, ytop = 22, col = bar.colors[2], border = F, xpd = NA)
rect(xleft = 15.5, xright = 16.5, ybottom = 17, ytop = 19, col = bar.colors[3], border = F, xpd = NA)

# box(which = "inner", col=adjustcolor("red", alpha.f = .5), lwd = 3)
# box(which = "outer", col=adjustcolor("blue", alpha.f = .5), lwd = 3)
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)
#----
dev.off()


# ---- Paper Supplemental Figures ----

# Michigan ----
forcing.folder <- michigan.forcing.folder
# Run "quick look" section to load data
save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_3_Michigan.pdf"
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .4)) # bottom, left, top, right
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30)
y.tick.labels <- c(0,"", 10, "", 20,"",30)
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 30
# Run paper figure line 167 to 186
panel.data <- forcing.data$lineage
panel.data <- panel.data[ ,1:20]
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35)
y.tick.labels <- c(0,"",10,"",20,"",30,"")
# Run paper figure lines 195 to 219
mtext(text = expression(bold("Lake Michigan")), side = 3, line = -5, at = 18, cex = 2)
dev.off()

# Danube ----
forcing.folder <- danube.forcing.folder
save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_3_Danube.pdf"
# Run "quick look" section to load data
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .2)) # bottom, left, top, right
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30)
y.tick.labels <- c(0,"", 10, "", 20,"",30)
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 30
# Run paper figure line 167 to 186
panel.data <- forcing.data$lineage
panel.data <- panel.data[ ,1:20]
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40)
y.tick.labels <- c(0,"",10,"",20,"",30,"",40)
# Run paper figure lines 195 to 219
mtext(text = expression(bold("Danube River")), side = 3, line = -5, at = 18, cex = 2)
dev.off()

# Bog Epi ----
forcing.folder <- bog.epi.forcing.folder
save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_3_Bog_Epi.pdf"
# Run "quick look" section to load data
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .6)) # bottom, left, top, right
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40,45,50)
y.tick.labels <- c(0,"", 10, "", 20,"",30,"",40,"",50)
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 50
# Run paper figure line 167 to 186
panel.data <- forcing.data$lineage
panel.data <- panel.data[ ,1:20]
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40,45)
y.tick.labels <- c(0,"",10,"",20,"",30,"",40,"")
# Run paper figure lines 195 to 219
mtext(text = expression(bold("Bog Epilimnion")), side = 3, line = -5, at = 18, cex = 2)
dev.off()

# Bog Hypo ----
forcing.folder <- bog.epi.forcing.folder
save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_3_Bog_Hypo.pdf"
# Run "quick look" section to load data
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .6)) # bottom, left, top, right
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40,45,50)
y.tick.labels <- c(0,"", 10, "", 20,"",30,"",40,"",50)
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 50
# Run paper figure line 167 to 186
panel.data <- forcing.data$lineage
panel.data <- panel.data[ ,1:20]
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,5,10,15,20,25,30,35,40,45)
y.tick.labels <- c(0,"",10,"",20,"",30,"",40,"")
# Run paper figure lines 195 to 219
mtext(text = expression(bold("Bog Hypolimnion")), side = 3, line = -5, at = 18, cex = 2)
dev.off()

# Mouse Gut ----
forcing.folder <- mouse.forcing.folder
save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_3_Mouse_Gut.pdf"
# Run "quick look" section to load data
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,2), nrow = 1))
par(mai = c(.55, .2, .24, 0), omi = c(0, .27, .05, .6)) # bottom, left, top, right
panel.data <- forcing.data$phylum
plot.title <- "Phylum Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[3, ]) # grey bar + BLUE bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,10,20,30,40,50,60,70)
y.tick.labels <- c(0,"",20,"",40,"",60,"")
y.axis.label <- "Relative Abundance (% reads)"
y.max <- 70
# Run paper figure line 167 to 186
panel.data <- forcing.data$lineage
plot.title <- "Family/Lineage Rank Abundance"
max.bar <- max(panel.data[1, ] + panel.data[2, ]) # grey bar + red bar
max.bar # choose axis max based on this
y.axis.ticks <- c(0,10,20,30,40,50,60,70,80,90,100)
y.tick.labels <- c(0,"",20,"",40,"",60,"",80,"",100)
# Run paper figure lines 195 to 202
text(x = loc.labels - .1, y = -1.1, labels = bar.labels, srt = -30, xpd = NA, cex = .9, adj = 0)
# Run paper figure lines 204 to 219
# just paste the stupid legend on later.
mtext(text = expression(bold("Mouse Gut")), side = 3, line = -5, at = 9, cex = 2)
dev.off()

#
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
