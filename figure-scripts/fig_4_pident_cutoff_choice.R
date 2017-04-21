# RRR 8-15-16 ----

# fig 4 is the plot showing how you can identify a good percent identity cutoff

# fig 4a shows the percent identity cutoff that results in the maximum reads classified at different taxa levels
# fig 4b shows the diminishing returns of lowering the cutoff as fewer reads at added to the FW classification

# Supp Table 1 is the sanity check- make sure there aren't any phyla or classes being forced

# ---- Define File Paths ----

file.path.otu.summs <- "../../poster/poster_mend_unclust/plots/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "../../poster/poster_mend_unclust/plots/conflict_summary_by_percent_OTUs.csv"
file.path.read.perc.summs <- "../../poster/poster_mend_unclust/plots/conflict_summary_by_percent_reads.csv"

output.folder.supp <- "~/Desktop/test/supp"
output.folder.fig4 <- "~/Desktop/test/fig4"


# ---- Define Functions ----

import.summary <- function(FilePath){
  sumry <- read.csv(file = FilePath, header = FALSE, colClasses = "character")
  sumry[1,1] <- "pident"
  sumry[ ,-1] <- apply(X = sumry[ ,-1], MARGIN = 2, FUN = as.numeric)
  return(sumry)
}

make.supplemental.table.1 <- function(ConflictsSum, FolderPath = NULL){
  otus <- ConflictsSum[1:4, ]
  # only make a file if a folder is specified
  if (!is.null(FolderPath)){
    file.name <- paste(FolderPath, "Supplemental_Table_1.csv", sep = "")
    write.csv(x = otus, file = file.name, quote = FALSE, row.names = FALSE)
    cat("made: ", file.name)
  }
  return(otus)
}

trim.to.perc.classified <- function(ConflictSum){
  perc.class <- ConflictSum[-c(2:6), ]
  future.colnames <- perc.class[ ,1]
  perc.class <- perc.class[ ,-1]
  perc.class <- apply(perc.class, 2, as.numeric)
  perc.class <- t(perc.class)
  colnames(perc.class) <- future.colnames
  return(perc.class)
}

plot.perc.classified <- function(PercClass, Cutoff, FilePath = NULL){
  pidents <- PercClass[ ,1]
  perc.class <- PercClass[ ,2]
  
  # only save plot if folder specified
  if (!is.null(FilePath)){
    png(filename = FilePath, width = 7, height = 5, units = "in", res = 100)
  }
  
  # Set up an empty plot
  plot.title <- "Percent of Data Classified by Ecosystem-Specific Database"
  x.label <- "percent identity"
  y.label <- paste("Total Reads (%)")
  plot(x = pidents, y = perc.class, type = "n", main = plot.title, cex.main = 1, xlab = x.label, ylab = y.label)
  
  # Fill plot with beautiful data
  lines(x = pidents, y = perc.class, col = "grey", lwd = 3)
  points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = 1)
  abline(v = Cutoff, col = "red", xpd = F)
  
  # stop exporting plot if you started to
  if (!is.null(FilePath)){
    unnecessary.message <- dev.off()
    cat("made plot: ", FilePath)
  }
}

plot.max.classified.........

# ---- Use Functions for quick looks ----

otus <- import.summary(FilePath = file.path.otu.summs)

make.supplemental.table.1(ConflictsSum = otus)

reads <- import.summary(FilePath = file.path.read.perc.summs)

reads.class <- trim.to.perc.classified(ConflictSum = reads)

plot.perc.classified(PercClass = reads.class, Cutoff = 98)


# ---- PAPER ----

draft.3.data <- "../../ME_GG/analysis/plots/conflict_summary_by_percent_reads.csv"
draft.3.fig.4b <- "~/Dropbox/PhD/Write It/draft 3/draft_3_figure_files/fig_4b_percent_fw.png"

reads <- import.summary(FilePath = draft.3.data)
reads.class <- trim.to.perc.classified(ConflictSum = reads)
plot.perc.classified(PercClass = reads.class, Cutoff = 99, FilePath = draft.3.fig.4b)


# ---- ISME16 POSTER ----

file.path.read.perc.summs <- "../../poster_mend_unclust/plots/conflict_summary_by_percent_reads.csv"
output.folder.fig4 <- "~/Dropbox/Trina/8-20-16_ISME16_figures/pident_choice_mendota_unclust.png"

# define functions ----

import.summary <- function(FilePath){
  sumry <- read.csv(file = FilePath, header = FALSE, colClasses = "character")
  sumry[1,1] <- "pident"
  sumry[ ,-1] <- apply(X = sumry[ ,-1], MARGIN = 2, FUN = as.numeric)
  return(sumry)
}

trim.to.perc.classified <- function(ConflictSum){
  perc.class <- ConflictSum[-c(2:6), ]
  future.colnames <- perc.class[ ,1]
  perc.class <- perc.class[ ,-1]
  perc.class <- apply(perc.class, 2, as.numeric)
  perc.class <- t(perc.class)
  colnames(perc.class) <- future.colnames
  return(perc.class)
}

# use functions ----

reads <- import.summary(FilePath = file.path.read.perc.summs)

reads.class <- trim.to.perc.classified(ConflictSum = reads)

PercClass = reads.class
FolderPath = output.folder.fig4

pidents <- PercClass[ ,1]
perc.class <- PercClass[ ,2]

pidents <- pidents[1:6]
perc.class <- perc.class[1:6]

plot.title <- c(expression(bold("Percent Classified by")), expression(bold("the FreshTrain")))
x.label <- expression(bold("Percent Identity Cutoff"))
y.label <- expression(bold("Reads (%)"))
y.min <- min(perc.class)
y.max <- max(perc.class)
x.min <- min(pidents)
x.max <- max(pidents)


png(filename = FolderPath, width = 2.8, height = 3, units = "in", res = 300)
par(mar = c(2.3,2.65,1.8,.51))
plot(x = pidents, y = perc.class, type = "n", ann = F, axes = F, ylim = c(y.min, y.max))
polygon(x = c(pidents[2], pidents[2:4], pidents[4]), y = c(y.min, perc.class[2:4], y.min), xpd=T, border = NA, col = adjustcolor(col = "red", alpha.f = .2))
lines(x = pidents, y = perc.class, col = "grey", lwd = 5)
points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = .75)
lines(x = c(98,98), y = c(y.min, perc.class[3]), col = "red", lwd = 3)
lines(x = c(x.min, 98), y = c(perc.class[3], perc.class[3]), col = "blue", lwd = 3)
axis(side = 1, at = pidents, labels = F, tick = TRUE, line = -.5, lwd = 3, xpd = T)
mtext(text = pidents, side = 1, line = .2, at = pidents - .1, cex = 1.5)
mtext(text = x.label, side = 1, line = 1.5, cex = 1.5, at = 93.7, adj = 0)
axis(side = 2, at = c(50,60,70,80), labels = F, tick = T, line = -.5, lwd = 3, xpd =T)
mtext(text = c(50,60,70,80), at = c(50,60,70,80) + c(.75,0,0,-.75), side = 2, line = .1, cex = 1.5)
mtext(text = y.label, side = 2, line = 1.2, cex = 1.5)
mtext(text = plot.title, side = 3, at = c(94, 95.7), adj = 0, line = c(.4,-.6), cex = 1.5, col = "red")
dev.off()









