# RRR 8-15-16 ----

# fig 4 is the plot showing how you can choose the percent identity cutoff to use
#       based on the diminishing returns of percent database classified by ecosystem-specific

# fig 4a shows the choice for mendota unclustered reads (this looks smooth)
# fig 4b shows the same unclustered plot for other ecosystems (some more dissimilar ones are rougher)
# fig 4c shows the impact of clustering reads on the plot shape (they get rough where the cluster percent is)
# Supp Table 1 is the sanity check- make sure there aren't any phyla or classes being forced

# ---- Define File Paths ----

file.path.otu.summs <- "../../take_mendota_uniq/plots/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "../../take_mendota_uniq/plots/conflict_summary_by_percent_OTUs.csv"
file.path.read.perc.summs <- "../../take_mendota_uniq/plots/conflict_summary_by_percent_reads.csv"
output.folder.supp <- "~/Desktop/figures_8-16-16/SupTable1/mendota_unique_"
output.folder.fig4 <- "~/Desktop/figures_8-16-16/fig4c/mendota_100_"

file.path.otu.summs <- "../../take_bogs_deblur/plots/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "../../take_bogs_deblur/plots/conflict_summary_by_percent_OTUs.csv"
file.path.read.perc.summs <- "../../take_bogs_deblur/plots/conflict_summary_by_percent_reads.csv"
output.folder.supp <- "~/Desktop/figures_8-16-16/SupTable1/bogs_"
output.folder.fig4 <- "~/Desktop/figures_8-16-16/fig4b/bogs_"

file.path.otu.summs <- "../../take_danube_10/plots/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "../../take_danube_10/plots/conflict_summary_by_percent_OTUs.csv"
file.path.read.perc.summs <- "../../take_danube_10/plots/conflict_summary_by_percent_reads.csv"
output.folder.supp <- "~/Desktop/figures_8-16-16/SupTable1/danube_"
output.folder.fig4 <- "~/Desktop/figures_8-16-16/fig4b/danube_"

file.path.otu.summs <- "../../take_mendota_clust/plots/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "../../take_mendota_clust/plots/conflict_summary_by_percent_OTUs.csv"
file.path.read.perc.summs <- "../../take_mendota_clust/plots/conflict_summary_by_percent_reads.csv"
output.folder.supp <- "~/Desktop/figures_8-16-16/SupTable1/mendota_clust_"
output.folder.fig4 <- "~/Desktop/figures_8-16-16/fig4c/mendota_98_"

# ---- Define Functions ----

import.summary <- function(FilePath){
  sumry <- read.csv(file = FilePath, header = FALSE, colClasses = "character")
  sumry[1,1] <- "pident"
  sumry[ ,-1] <- apply(X = sumry[ ,-1], MARGIN = 2, FUN = as.numeric)
  return(sumry)
}

make.supplemental.table.1 <- function(ConflictsSum, FolderPath){
  otus <- ConflictsSum[1:4, ]
  file.name <- paste(FolderPath, "Supplemental_Table_1.csv", sep = "")
  write.csv(x = otus, file = file.name, quote = FALSE, row.names = FALSE)
  cat("made: ", file.name)
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

plot.perc.classified <- function(PercClass, Cutoff, FolderPath){
  pidents <- PercClass[ ,1]
  perc.class <- PercClass[ ,2]
  
  plot.name <- paste(FolderPath, ".png", sep = "")
  png(filename = plot.name, width = 7, height = 5, units = "in", res = 100)
  
  # Set up and empty plot
  plot.title <- "Percent of Data Classified by Ecosystem-Specific Database"
  x.label <- "percent identity"
  y.label <- paste("Total Reads (%)")
  plot(x = pidents, y = perc.class, type = "n", main = plot.title, cex.main = 1, xlab = x.label, ylab = y.label)
  
  # Fill plot with beautiful data
  lines(x = pidents, y = perc.class, col = "grey", lwd = 1.5)
  points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = .5)
  abline(v = Cutoff, col = "red", xpd = F)
  
  dev.off()
  cat("made plot: ", plot.name)
}

# ---- Use Functions ----

otus <- import.summary(FilePath = file.path.otu.summs)

make.supplemental.table.1(ConflictsSum = otus, FolderPath = output.folder.supp)

reads <- import.summary(FilePath = file.path.read.perc.summs)

reads.class <- trim.to.perc.classified(ConflictSum = reads)

plot.perc.classified(PercClass = reads.class, Cutoff = 98, FolderPath = output.folder.fig4)



