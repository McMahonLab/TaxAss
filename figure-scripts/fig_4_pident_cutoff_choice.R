# RRR 8-15-16 ----

# fig 4 is the plot showing how you can identify a good percent identity cutoff

# fig 4a shows the percent identity cutoff that results in the maximum reads classified at different taxa levels
# fig 4b shows the diminishing returns of lowering the cutoff as fewer reads at added to the FW classification

# Supp Table 1 is the sanity check- make sure there aren't any phyla or classes being forced

# ---- Define File Paths ----

# file.path.otu.summs <- "../../poster/poster_mend_unclust/plots/conflict_summary_by_OTUs.csv"
# file.path.otu.perc.summs <- "../../poster/poster_mend_unclust/plots/conflict_summary_by_percent_OTUs.csv"

file.path.read.perc.summs <- "../../ME_GG/analysis/plots/conflict_summary_by_percent_reads.csv"
file.path.reads.class <- "../../ME_GG/analysis/plots/Percent_Reads_Classified_by_Pident.csv"

# output.folder.supp <- "~/Desktop/test/supp"
# output.folder.fig4 <- "~/Desktop/test/fig4"


# ---- Define Functions ----

import.conflict.summary <- function(FilePath){
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

import.classified.summary <- function(FilePath){
  sumry <- read.csv(file = FilePath, header = TRUE, colClasses = "character")
  sumry[ ,1] <- substr(x = sumry[ ,1], start = 8, stop = 10)
  colnames(sumry)[1] <- "pident"
  sumry <- apply(X = sumry, MARGIN = 2, FUN = as.numeric)
  return(sumry)
}

plot.total.classified <- function(SummaryMatrix, PidentValues, FilePath = NULL){
  pidents <- PidentValues
  sum.named <- SummaryMatrix
  
  line.col <- rainbow(n = ncol(sum.named), v = .4)
  x.lim <- c(min(pidents), max(pidents))
  y.lim <- c(40,100)
  y.label <- "Percent Classified (Reads)"
  x.label <- "Percent Identity Cutoff"
  plot.title <- expression(bold("Percent of Dataset Classified"))
  taxa.levels <- sub(pattern = ".fw", replacement = "", x = colnames(sum.named))
  
  if(!is.null(FilePath)){
    png(filename = FilePath, width = 10, height = 5, units = "in", res = 100)
  }
  
  par(mfrow = c(1,ncol(sum.named)), omi = c(.4,.3,.3,.1), mai = c(.2,.3,.3,0))
  for (t in 1:ncol(sum.named)){
    ass <- sum.named[ ,t]
    # basic plot
    plot(x = pidents, y = ass, col = line.col[t], type = "l", ann = F, lwd = 3, axes = F)
    mtext(text = taxa.levels[t], side = 3, line = .5, outer = F, cex = 1.2, col = line.col[t])
    
    # vertical max line
    index <- which(ass == max(ass))
    max.names <- pidents[index]
    abline(v = max.names, col = adjustcolor(col = line.col[t], alpha.f = .3), lwd = 3)
    
    # x axis labels
    x.lab.cols <- rep("black", times = length(pidents))
    x.lab.cols[index] <- line.col[t]
    x.lab.cex <- rep(.7, times = length(pidents))
    x.lab.cex[index] <- 2
    x.lab.line <- rep(.5, times = length(pidents))
    x.lab.line[index] <- 1.5
    empty.x.labels <- rep("", times = length(pidents))
    axis(side = 1, at = pidents, labels = empty.x.labels)
    mtext(text = pidents, side = 1, line = x.lab.line, at = pidents, col = x.lab.cols, cex = x.lab.cex)
    
    # y axis labels
    span <- max(ass) - min(ass)
    y.ax <- c(min(ass), min(ass) + (span * 1/3), min(ass) + (span * 2/3), max(ass))
    y.ax.lab <- round(x = y.ax, digits = 0)
    empty.y.labels <- rep("", times = length(y.ax))
    axis(side = 2, at = y.ax, labels = empty.y.labels)
    mtext(text = y.ax.lab, side = 2, line = .7, at = y.ax)
  } 
  mtext(text = plot.title, side = 3, line = .5, outer = T, cex = 1.2)
  mtext(text = x.label, side = 1, line = 1.5, outer = T, cex = 1.2)
  mtext(text = y.label, side = 2, line = .5, outer = T, cex = 1.2)
  
  if(!is.null(FilePath)){
    unnecessary.message <- dev.off()
    cat("made plot: ", FilePath, "\n")
  }
}


# ---- Use Functions for quick looks ----

otus <- import.conflict.summary(FilePath = file.path.otu.summs)
make.supplemental.table.1(ConflictsSum = otus)

reads <- import.conflict.summary(FilePath = file.path.read.perc.summs)
reads.fw.class <- trim.to.perc.classified(ConflictSum = reads)
plot.perc.classified(PercClass = reads.fw.class, Cutoff = 98)

reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
pident.values <- reads.tot.class[ ,1]
reads.tot.class.plot <- reads.tot.class[ ,-c(1:3,8)]
plot.total.classified(SummaryMatrix =reads.tot.class.plot, PidentValues = pident.values)

# ---- PAPER ----
save.to <- "~/Dropbox/PhD/Write It/draft 4/fig_4.pdf"
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
layout(mat = matrix(c(1,2,3,4,5,5), nrow = 1))
par(omi = c(0,.4,0,0)) # bottom, left, top, right
# 4a ----
par(mai = c(.3,0,.4,.2)) # bottom, left, top, right
# ----
pidents <- pident.values
sum.named <- reads.tot.class.plot

line.col <- "grey" 
pointer.col <- adjustcolor(col = "red", alpha.f = .3)
x.lim <- c(min(pidents), max(pidents))
y.label <- "Reads Classified (%)"
x.label <- "Percent Identity Cutoff"
taxa.levels <- c("Class","Order","Family/\nLineage","Genus/\nClade")
big.title <- "Total Classifications"

# Class ----
ass <- sum.named[ ,1]
plot.title <- taxa.levels[1]
min(ass)
max(ass)
y.lim <- c(90,95)
y.ticks <- c(90,91,92,93,94,95)
y.tick.labs <- c("",91,"",93,"",95)
repeat.these <- function(){
  plot(x = pidents, y = ass, col = line.col, type = "n", ann = F, lwd = 3, axes = F, ylim = y.lim)
  points(x = pidents, y = ass, col = line.col, pch = 19)
  lines(x = pidents, y = ass, col = line.col, lwd = 3)
  # vertical max line
  index <- which(ass == max(ass))
  max.names <- pidents[index]
  lines(x = c(pidents[index],pidents[index]),y = c(y.lim[1], ass[index]), col = pointer.col, lwd = 3)
  # x axis 
  x.lab.text <- c(100,expression(bold("99")),98,97,96,95)
  axis(side = 1, at = pidents, labels = F, tck = -.035, line = -.5)
  mtext(text = x.lab.text, side = 1, line = -.3, at = pidents, col = "black", cex = .6)
  # Y Axis
  axis(side = 2, at = y.ticks, labels = F, tck = -.03, line = .2)
  mtext(text = y.tick.labs, side = 2, at = y.ticks, las = 2, line = .5, cex = .6)
  # plot title
  mtext(text = plot.title, side = 3, line = 1, outer = F, cex = .8, padj = 1)
}
repeat.these()

# ----
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)

# Order ----
ass <- sum.named[ ,2]
plot.title <- taxa.levels[2]
min(ass)
max(ass)
y.lim <- c(85,90)
y.ticks <- c(85,86,87,88,89,90)
y.tick.labs <- c("",86,"",88,"",90)
repeat.these()
# Family ----
ass <- sum.named[ ,3]
plot.title <- taxa.levels[3]
min(ass)
max(ass)
y.lim <- c(75,85)
y.ticks <- c(75,77,79,81,83,85)
y.tick.labs <- c("",77,"",81,"",85)
repeat.these()
# Genus ----
ass <- sum.named[ ,4]
plot.title <- taxa.levels[4]
min(ass)
max(ass)
y.lim <- c(57.5,70)
y.ticks <- c(57.5,60,62.5,65,67.5,70)
y.tick.labs <- c("",60,"",65,"",70)
repeat.these()

# ----
mtext(text = x.label, side = 1, line = -1.1, outer = T, at = .3, cex = .8)
mtext(text = y.label, side = 2, line = 1.9, outer = T, cex = .8)
mtext(text = big.title, side = 3, line = -1.3, outer = T, at = .3)

# box(which = "inner", col=adjustcolor("red", alpha.f = .5), lwd = 3)
# box(which = "outer", col=adjustcolor("blue", alpha.f = .5), lwd = 3)
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)

# 4b ----
par(mai = c(.3,.6,.27,.05)) # bottom, left, top, right
# ----
perc.class <- reads.fw.class[ ,2]
chosen.cutoff <- 99
plot.title <- "FreshTrain Classifications"
x.label <- "Percent Identity Cutoff"
y.label <- paste("Reads in FreshTrain Group (%)")
min(perc.class)
max(perc.class)
y.lim <- c(45,80)
y.ticks <- c(45,50,55,60,65,70,75,80)
y.tick.labs <- c("",50,"",60,"",70,"",80)

plot(x = pidents, y = perc.class, type = "n", axes = F, ann = F, ylim = y.lim)
lines(x = pidents, y = perc.class, col = "grey", lwd = 3)
points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = 1)
# red line
index <- which(pidents == chosen.cutoff)
cutoff.perc <- perc.class[index]
lines(x = c(chosen.cutoff,chosen.cutoff),y = c(y.lim[1], cutoff.perc), col = pointer.col, lwd = 3)
# ----
# X Axis
boldeable.text <- c(100,expression(bold("99")),98, 97, 96, 95)
axis(side = 1, at = pidents, labels = F, tck = -.03, line = -.5)
mtext(text = boldeable.text, side = 1, line = -.2, at = pidents, cex = .7)
# Y Axis
axis(side = 2, at = y.ticks, labels = F, tck = -.03, line = 0, tck = -.03)
mtext(text = y.tick.labs, side = 2, at = y.ticks, las = 2, line = .5, cex = .7)
# Labels
mtext(text = plot.title, side = 3, line = 1.6, at = 100.2, cex = 1, adj = 1, padj = 1)
mtext(text = x.label, side = 1, line = 1.1, at = 95, adj = 0, cex = .8)
mtext(text = y.label, side = 2, line = 2, at = 48, adj = 0, cex = .8)
#----
dev.off()

#
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









