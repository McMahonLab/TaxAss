# RRR 8-15-16 ----

# fig 4 is the plot showing how you can identify a good percent identity cutoff
# it shows the percent identity cutoff that results in the maximum reads classified at different taxa levels

# fig 4b was cut from the paper- shows the diminishing returns of lowering the cutoff as fewer reads are added to the FW classification

# Supp Table 1 is the sanity check- make sure there aren't any phyla or classes being forced

# ---- Define File Paths ----

file.path.otu.summs <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_14_Choose_pident_cutoff/conflict_summary_by_OTUs.csv"
file.path.otu.perc.summs <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_14_Choose_pident_cutoff/conflict_summary_by_percent_OTUs.csv"

file.path.read.perc.summs <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_14_Choose_pident_cutoff/conflict_summary_by_percent_reads.csv"
file.path.reads.class <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"

output.folder.supp <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/"
output.folder.fig4 <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/"


# ---- Define Functions ----

import.conflict.summary <- function(FilePath){
  sumry <- read.csv(file = FilePath, header = FALSE, colClasses = "character")
  sumry[1,1] <- "pident"
  sumry[ ,-1] <- apply(X = sumry[ ,-1], MARGIN = 2, FUN = as.numeric)
  return(sumry)
}

make.supplemental.table.1 <- function(ConflictsSum, FolderPath = NULL){
  otus <- ConflictsSum[1:4, ]
  index <- order(otus[1,-1])
  otus <- otus[ ,c(1,index+1)]
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
    
    # make y axis at least 5 percentage points range so plots less misleading
    y.max <- max(ass)
    y.min <- min(ass)
    y.range <- y.max - y.min
    if (y.range < 5){
      half.range <- y.range / 2
      extend.by <- 2.5 - half.range
      y.min <- y.min - extend.by
      y.max <- y.max + extend.by
    }
    
    # basic plot
    plot(x = pidents, y = ass, col = line.col[t], type = "l", ann = F, lwd = 3, axes = F, ylim = c(y.min, y.max))
    mtext(text = taxa.levels[t], side = 3, line = .5, outer = F, cex = 1.2, col = line.col[t])
    
    # vertical max line
    index <- which(ass == max(ass))
    max.names <- pidents[index]
    lines(x = c(max.names, max.names), y = c(0, max(ass)), col = adjustcolor(col = line.col[t], alpha.f = .3), lwd = 3)
    
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
    span <- y.max - y.min
    y.ax <- c(y.min, y.min + (span * 1/3), y.min + (span * 2/3), y.max)
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
plot.perc.classified(PercClass = reads.fw.class, Cutoff = 99)

reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
pident.values <- reads.tot.class[ ,1]
reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
plot.total.classified(SummaryMatrix =reads.tot.class.plot, PidentValues = pident.values)

# ---- PAPER resubmit version ----

mendota <- import.classified.summary(FilePath = "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv")
michigan <- import.classified.summary(FilePath = "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Michigan/TaxAss-Michigan/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv")
danube <- import.classified.summary(FilePath = "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Danube/TaxAss-Danube/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv")
bogepi <- import.classified.summary(FilePath = "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/TroutBogEpi/TaxAss-TroutBogEpi/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv")
boghypo <- import.classified.summary(FilePath = "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/TroutBogHypo/TaxAss-TroutBogHypo/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv")

pidents <- mendota[ ,1]
ecosystem.list <- list("Lake Mendota" = mendota[ ,-(1:2)], "Lake Michigan" = michigan[ ,-(1:2)], "Danube River" = danube[ ,-(1:2)], "Bog Epilimnion" = bogepi[ ,-(1:2)], "Bog Hypolimnion" = boghypo[ ,-(1:2)])

temp.color.options <- adjustcolor(col = c(rainbow(n = 20, v = .7), "burlywood4"), alpha.f = .8)
ecosystem.colors.lines <- c(temp.color.options[8], temp.color.options[12], temp.color.options[16], temp.color.options[1], temp.color.options[21])
temp.color.options <- adjustcolor(col = c(rainbow(n = 20, v = .7), "burlywood4"), alpha.f = 1)
ecosystem.colors.labels <- c(temp.color.options[8], temp.color.options[12], temp.color.options[16], temp.color.options[1], temp.color.options[21])

# # look at colors
# plot(1:21, col = temp.color.options, cex = 5, pch = 19)
# text(1:21, as.character(1:21))
# plot(1:5, pch = 19, col = ecosystem.colors.lines, cex = 20, xpd = T, axes = F)
# plot(1:5, pch = 19, col = ecosystem.colors.labels, cex = 20, xpd = T, axes = F)
# # check them on this: http://www.color-blindness.com/coblis-color-blindness-simulator/

x.lim <- c(min(pidents), max(pidents))
y.label <- "Reads Classified (%)"
x.label <- "Percent Identity Cutoff"
taxa.levels <- c("Phylum","Class","Order","Family/\nLineage","Genus/\nClade","Species/\nTribe")
num.taxa <- 6

max.eco.val <- as.vector(x = rep(0, times = num.taxa))
for (t in (7 - num.taxa):6){
  for (e in 1:length(ecosystem.list)){
    temp.max <- max(ecosystem.list[[e]][3,t])
    if(temp.max > max.eco.val[t]){
      max.eco.val[t] <- temp.max
    }
  }
}

y.axis.ranges <- matrix(data = c(85,100,85,100,80,95,65,85,45,65,14,45),nrow = 2)
colnames(y.axis.ranges) <- c("p","c","o","l","c","t")
y.axis.ticks <- list("p" = c(85,90,95,100), "c" = c(85,90,95,100), "o" = c(80,85,90,95), "lf" = c(65,70,75,80,85), "cg" = c(45,50,55,60,65), "ts" = c(15,20,25,30,35,40,45))

# ---- Fig 4 ----
save.to <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/pident_all.pdf"
pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")

par(mfrow = c(1,num.taxa), omi = c(.05,.18,.1,.1), mai = c(.3,.3,.4,0)) # bottom, left, top, right

for (t in (7 - num.taxa):6){
  plot(x = pidents, y = pidents, type = "n", ann = F, lwd = 3, axes = F, ylim = y.axis.ranges[ ,t])

  axis(side = 1, at = pidents, labels = F, tck = -.035, line = -.5)
  mtext(text = pidents, side = 1, line = -.2, at = pidents, cex = .6, xpd = T)
  
  axis(side = 2, at = y.axis.ticks[[t]], labels = F, tck = -.03, line = .2)
  mtext(text = y.axis.ticks[[t]], side = 2, at = y.axis.ticks[[t]], las = 2, line = .65, cex = .7, xpd = NA)
  
  mtext(text = taxa.levels[t], side = 3, line = 1, outer = F, cex = .8, padj = 1)
  
  
  lines(x = c(98,98), y = c(y.axis.ranges[1,t], max.eco.val[t]), col = adjustcolor("black", alpha.f = .1), lwd = 3)
  
  for (e in 1:length(ecosystem.list)){
    points(x = pidents, y = ecosystem.list[[e]][ ,t], col = ecosystem.colors.lines[e], pch = 19)
    lines(x = pidents, y = ecosystem.list[[e]][ ,t], col = ecosystem.colors.lines[e], lwd = 3)
  }
}


# Phylum ----
ass <- sum.named[ ,1]
plot.title <- taxa.levels[1]
min(ass)
max(ass)
y.lim <- c(85,100)
y.ticks <- c(85,90,95,100)
y.tick.labs <- c(85,90,95,100)
repeat.these <- function(){
  plot(x = pidents, y = ass, col = line.col, type = "n", ann = F, lwd = 3, axes = F, ylim = y.lim)
  points(x = pidents, y = ass, col = line.col, pch = 19)
  lines(x = pidents, y = ass, col = line.col, lwd = 3)
  # vertical max line
  index <- which(ass == max(ass))
  max.names <- pidents[index]
  lines(x = c(pidents[index],pidents[index]),y = c(y.lim[1], ass[index]), col = pointer.col, lwd = 3)
  # x axis
  x.lab.text <- c(100,99,98,97,96,95) # expression(bold("98")) if wanted to bold the chosen pident
  axis(side = 1, at = pidents, labels = F, tck = -.035, line = -.5)
  mtext(text = x.lab.text, side = 1, line = -.2, at = pidents, col = "black", cex = .6)
  # Y Axis
  axis(side = 2, at = y.ticks, labels = F, tck = -.03, line = .2)
  mtext(text = y.tick.labs, side = 2, at = y.ticks, las = 2, line = .65, cex = .7)
  # plot title
  mtext(text = plot.title, side = 3, line = 1, outer = F, cex = .8, padj = 1)
}
repeat.these()

# Class ----
ass <- sum.named[ ,2]
plot.title <- taxa.levels[2]
min(ass)
max(ass)
y.lim <- c(85,100)
y.ticks <- c(85,90,95,100)
y.tick.labs <- c(85,90,95,100)
repeat.these()

# ----
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)

# Order ----
ass <- sum.named[ ,3]
plot.title <- taxa.levels[3]
min(ass)
max(ass)
y.lim <- c(80,95)
y.ticks <- c(80,85,90,95)
y.tick.labs <- c(80,85,90,95)
repeat.these()
# Family ----
ass <- sum.named[ ,4]
plot.title <- taxa.levels[4]
min(ass)
max(ass)
y.lim <- c(65,85)
y.ticks <- c(65,70,75,80,85)
y.tick.labs <- c(65,70,75,80,85)
repeat.these()
# Genus ----
ass <- sum.named[ ,5]
plot.title <- taxa.levels[5]
min(ass)
max(ass)
y.lim <- c(45,65)
y.ticks <- c(45,50,55,60,65)
y.tick.labs <- c(45,50,55,60,65)
repeat.these()
# Species ----
ass <- sum.named[ ,6]
plot.title <- taxa.levels[6]
min(ass)
max(ass)
y.lim <- c(15,45)
y.ticks <- c(15,20,25,30,35,40,45)
y.tick.labs <- c(15,20,25,30,35,40,45)
repeat.these()


# ----
mtext(text = x.label, side = 1, line = -.8, outer = T, cex = .8)
mtext(text = y.label, side = 2, line = 0, outer = T, cex = .8)
mtext(text = big.title, side = 3, line = -.7, outer = T, at = .15, adj = 0)

# box(which = "inner", col=adjustcolor("red", alpha.f = .5), lwd = 3)
# box(which = "outer", col=adjustcolor("blue", alpha.f = .5), lwd = 3)
# box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)

#----
dev.off()


# # ---- PAPER first version ----
# 
# # ---- Fig 4 ----
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Figure_4.pdf"
# pdf(file = save.to, width = 6.875, height = 3, family = "Helvetica", title = "TaxAss Fig 2", colormodel = "srgb")
# # layout(mat = matrix(c(1,2,3,4), nrow = 1))
# par(mfrow = c(1,6), omi = c(.05,.18,.1,.1), mai = c(.3,.3,.4,0)) # bottom, left, top, right
# 
# # ----
# pidents <- pident.values
# sum.named <- reads.tot.class.plot
# 
# line.col <- "purple4" 
# pointer.col <- adjustcolor(col = "purple4", alpha.f = .1) 
# x.lim <- c(min(pidents), max(pidents))
# y.label <- "Reads Classified (%)"
# x.label <- "Percent Identity Cutoff"
# taxa.levels <- c("Phylum","Class","Order","Family/\nLineage","Genus/\nClade","Species/Tribe")
# big.title <- "Percent Identity Cutoff Where Classifications Are Maximized"
# 
# # Phylum ----  
# ass <- sum.named[ ,1]
# plot.title <- taxa.levels[1]
# min(ass)
# max(ass)
# y.lim <- c(90,100)
# y.ticks <- c(90,92.5,95,97.5,100)
# y.tick.labs <- c("90","","95","","100")
# repeat.these <- function(){
#   plot(x = pidents, y = ass, col = line.col, type = "n", ann = F, lwd = 3, axes = F, ylim = y.lim)
#   points(x = pidents, y = ass, col = line.col, pch = 19)
#   lines(x = pidents, y = ass, col = line.col, lwd = 3)
#   # vertical max line
#   index <- which(ass == max(ass))
#   max.names <- pidents[index]
#   lines(x = c(pidents[index],pidents[index]),y = c(y.lim[1], ass[index]), col = pointer.col, lwd = 3)
#   # x axis 
#   x.lab.text <- c(100,99,98,97,96,95) # expression(bold("98")) if wanted to bold the chosen pident
#   axis(side = 1, at = pidents, labels = F, tck = -.035, line = -.5)
#   mtext(text = x.lab.text, side = 1, line = -.2, at = pidents, col = "black", cex = .6)
#   # Y Axis
#   axis(side = 2, at = y.ticks, labels = F, tck = -.03, line = .2)
#   mtext(text = y.tick.labs, side = 2, at = y.ticks, las = 2, line = .65, cex = .7)
#   # plot title
#   mtext(text = plot.title, side = 3, line = 1, outer = F, cex = .8, padj = 1)
# }
# repeat.these()
# 
# # Class ----
# ass <- sum.named[ ,2]
# plot.title <- taxa.levels[2]
# min(ass)
# max(ass)
# y.lim <- c(85,95)
# y.ticks <- c(85, 87.5, 90, 92.5, 95)
# y.tick.labs <- c("85","","90","","95")
# repeat.these()
# 
# # ----
# # box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# # box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)
# 
# # Order ----
# ass <- sum.named[ ,3]
# plot.title <- taxa.levels[3]
# min(ass)
# max(ass)
# y.lim <- c(80,90)
# y.ticks <- c(80,82.5,85,87.5,90)
# y.tick.labs <- c("80","","85","","90")
# repeat.these()
# # Family ----
# ass <- sum.named[ ,4]
# plot.title <- taxa.levels[4]
# min(ass)
# max(ass)
# y.lim <- c(70,80)
# y.ticks <- c(70,72.5,75,77.5,80)
# y.tick.labs <- c("70","","75","","80")
# repeat.these()
# # Genus ----
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(55,65)
# y.ticks <- c(55,57.5,60,62.5,65)
# y.tick.labs <- c("55","",60,"",65)
# repeat.these()
# # Species ----
# ass <- sum.named[ ,6]
# plot.title <- taxa.levels[6]
# min(ass)
# max(ass)
# y.lim <- c(35,45)
# y.ticks <- c(35,37.5,40,42.5,45)
# y.tick.labs <- c("35","",40,"",45)
# repeat.these()
# 
# 
# # ----
# mtext(text = x.label, side = 1, line = -.8, outer = T, cex = .8)
# mtext(text = y.label, side = 2, line = 0, outer = T, cex = .8)
# mtext(text = big.title, side = 3, line = -.7, outer = T, at = .15, adj = 0)
# 
# # box(which = "inner", col=adjustcolor("red", alpha.f = .5), lwd = 3)
# # box(which = "outer", col=adjustcolor("blue", alpha.f = .5), lwd = 3)
# # box(which = "plot", col=adjustcolor("purple", alpha.f = .5), lwd = 3)
# # box(which = "figure", col=adjustcolor("orange", alpha.f = .5), lwd = 3)
# 
# #----
# dev.off()
# 
# #
# 
# # ---- Supplemental Figure 4 ----
# 
# # these figures source and load their own data, then jump btwn following figure lines and different y axes.
# 
# # Michigan ----
# file.path.reads.class <- "~/Desktop/TaxAss-BatchFiles-go/Michigan/TaxAss-Michigan/analysis/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_4_Michigan.pdf"
# reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
# pident.values <- reads.tot.class[ ,1]
# reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
# # Run Figure Script lines 176-243
# # Family 
# ass <- sum.named[ ,4]
# plot.title <- taxa.levels[4]
# min(ass)
# max(ass)
# y.lim <- c(65,75)
# y.ticks <- c(65,67.5,70,72.5,75)
# y.tick.labs <- c("65","","70","","75")
# repeat.these()
# # Genus 
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(50,60)
# y.ticks <- c(50,52.5,55,57.5,60)
# y.tick.labs <- c("50","",55,"",60)
# repeat.these()
# # Run Figure Script lines 261-274
# big.title <- "Lake Michigan"
# mtext(text = big.title, side = 3, line = -.7, outer = T)
# dev.off()
# 
# # Danube ----
# file.path.reads.class <- "~/Desktop/TaxAss-BatchFiles-go/Danube/TaxAss-Danube/analysis/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_4_Danube.pdf"
# reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
# pident.values <- reads.tot.class[ ,1]
# reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
# # Run Figure Script lines 176-220
# # Class
# ass <- sum.named[ ,2]
# plot.title <- taxa.levels[2]
# min(ass)
# max(ass)
# y.lim <- c(90,100)
# y.ticks <- c(90,92.5,95,97.5,100)
# y.tick.labs <- c("90","","95","","100")
# repeat.these()
# # Order
# ass <- sum.named[ ,3]
# plot.title <- taxa.levels[3]
# min(ass)
# max(ass)
# y.lim <- c(85,95)
# y.ticks <- c(85, 87.5, 90, 92.5, 95)
# y.tick.labs <- c("85","","90","","95")
# repeat.these()
# # Family
# ass <- sum.named[ ,4]
# plot.title <- taxa.levels[4]
# min(ass)
# max(ass)
# y.lim <- c(65,75)
# y.ticks <- c(65,67.5,70,72.5,75)
# y.tick.labs <- c("65","","70","","75")
# repeat.these()
# # Genus 
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(35,55)
# y.ticks <- c(35,37.5,40,42.5,45,47.5,50,52.5,55)
# y.tick.labs <- c("35","",40,"",45,"",50,"",55)
# repeat.these()
# # Species
# ass <- sum.named[ ,6]
# plot.title <- taxa.levels[6]
# min(ass)
# max(ass)
# y.lim <- c(15,40)
# y.ticks <- c(15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40)
# y.tick.labs <- c("15","",20,"",25,"",30,"",35,"",40)
# repeat.these()
# # Run Figure Script lines 273-274
# big.title <- "Danube River"
# mtext(text = big.title, side = 3, line = -.7, outer = T)
# dev.off()
# 
# # Bog Epi ----
# file.path.reads.class <- "~/Desktop/TaxAss-BatchFiles-go/TroutBogEpi/TaxAss-TroutBogEpi/analysis/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_4_Bog_Epi.pdf"
# reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
# pident.values <- reads.tot.class[ ,1]
# reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
# # Run Figure Script lines 176-234
# # Order
# ass <- sum.named[ ,3]
# plot.title <- taxa.levels[3]
# min(ass)
# max(ass)
# y.lim <- c(75,85)
# y.ticks <- c(75,77.5,80,82.5,85)
# y.tick.labs <- c("75","","80","","85")
# repeat.these()
# # Family
# ass <- sum.named[ ,4]
# plot.title <- taxa.levels[4]
# min(ass)
# max(ass)
# y.lim <- c(60,70)
# y.ticks <- c(60,62.5,65,67.5,70)
# y.tick.labs <- c("60","","65","","70")
# repeat.these()
# # Genus
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(40,50)
# y.ticks <- c(40,42.5,45,47.5,50)
# y.tick.labs <- c("40","","45","","50")
# repeat.these()
# # Species
# ass <- sum.named[ ,6]
# plot.title <- taxa.levels[6]
# min(ass)
# max(ass)
# y.lim <- c(20,35)
# y.ticks <- c(20,22.5,25,27.5,30,32.5,35)
# y.tick.labs <- c("20","",25,"",30,"",35)
# repeat.these()
# # Run Figure Script lines 273-274
# big.title <- "Trout Bog Epilimnion"
# mtext(text = big.title, side = 3, line = -.7, outer = T)
# dev.off()
# 
# # Bog Hypo ----
# file.path.reads.class <- "~/Desktop/TaxAss-BatchFiles-go/TroutBogHypo/TaxAss-TroutBogHypo/analysis/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_4_Bog_Hypo.pdf"
# reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
# pident.values <- reads.tot.class[ ,1]
# reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
# # Run Figure Script lines 176-252
# # Genus
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(40,50)
# y.ticks <- c(40,42.5,45,47.5,50)
# y.tick.labs <- c("40","",45,"",50)
# repeat.these()
# # Species
# ass <- sum.named[ ,6]
# plot.title <- taxa.levels[6]
# min(ass)
# max(ass)
# y.lim <- c(20,30)
# y.ticks <- c(20,22.5,25,27.5,30)
# y.tick.labs <- c("20","",25,"",30)
# repeat.these()
# # Run Figure Script lines 273-274
# big.title <- "Trout Bog Hypolimnion"
# mtext(text = big.title, side = 3, line = -.7, outer = T)
# dev.off()
# 
# # Mouse Gut ----
# file.path.reads.class <- "~/Desktop/TaxAss-BatchFiles-go/MouseGut/TaxAss-MouseGut/analysis/plots/step_14_Choose_pident_cutoff/Percent_Reads_Classified_by_Pident.csv"
# save.to <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Figure_4_Mouse_Gut.pdf"
# reads.tot.class <- import.classified.summary(FilePath = file.path.reads.class)
# pident.values <- reads.tot.class[ ,1]
# reads.tot.class.plot <- reads.tot.class[ ,-(1:2)]  #[ ,-c(1:3,8)]
# # Run Figure Script lines 176-220
# # Class
# ass <- sum.named[ ,2]
# plot.title <- taxa.levels[2]
# min(ass)
# max(ass)
# y.lim <- c(90,100)
# y.ticks <- c(90,92.5,95,97.5,100)
# y.tick.labs <- c("90","","95","","100")
# repeat.these()
# # Order
# ass <- sum.named[ ,3]
# plot.title <- taxa.levels[3]
# min(ass)
# max(ass)
# y.lim <- c(90,100)
# y.ticks <- c(90,92.5,95,97.5,100)
# y.tick.labs <- c("90","","95","","100")
# repeat.these()
# # Family
# ass <- sum.named[ ,4]
# plot.title <- taxa.levels[4]
# min(ass)
# max(ass)
# y.lim <- c(80,90)
# y.ticks <- c(80,82.5,85,87.5,90)
# y.tick.labs <- c("80","","85","","90")
# repeat.these()
# # Genus
# ass <- sum.named[ ,5]
# plot.title <- taxa.levels[5]
# min(ass)
# max(ass)
# y.lim <- c(15,25)
# y.ticks <- c(15,17.5,20,22.5,25)
# y.tick.labs <- c("15","",20,"",25)
# repeat.these()
# # Species
# ass <- sum.named[ ,6]
# plot.title <- taxa.levels[6]
# min(ass)
# max(ass)
# y.lim <- c(0,10)
# y.ticks <- c(0,2.5,5,7.5,10)
# y.tick.labs <- c("0","",5,"",10)
# repeat.these()
# # Run Figure Script lines 273-274
# big.title <- "Mouse Gut"
# mtext(text = big.title, side = 3, line = -.7, outer = T)
# dev.off()
# 
# # 
# # ---- old version of Supp Table 1 ----
# save.to <- "~/Dropbox/PhD/Write It/draft 3/draft_3_figure_files/sup_table_1_forcing_check.csv"
# supp.table.1 <- make.supplemental.table.1(ConflictsSum = otus, FolderPath = save.to)
# supp.table.1
# # actually now this is only at the chosen taxa level.
# # Directions in their own figure-script
# 
# # ---- deceased panel 4b (killed by coauthors. apparently it's incomprehensible) ----
# 
# par(mai = c(.3,.6,.27,.05)) # bottom, left, top, right
# # ----
# perc.class <- reads.fw.class[ ,2]
# chosen.cutoff <- 99
# plot.title <- "FreshTrain Classifications"
# x.label <- "Percent Identity Cutoff"
# y.label <- paste("Reads in FreshTrain Group (%)")
# min(perc.class)
# max(perc.class)
# y.lim <- c(45,80)
# y.ticks <- c(45,50,55,60,65,70,75,80)
# y.tick.labs <- c("",50,"",60,"",70,"",80)
# 
# plot(x = pidents, y = perc.class, type = "n", axes = F, ann = F, ylim = y.lim)
# lines(x = pidents, y = perc.class, col = "grey", lwd = 3)
# points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = 1)
# # red line
# index <- which(pidents == chosen.cutoff)
# cutoff.perc <- perc.class[index]
# lines(x = c(chosen.cutoff,chosen.cutoff),y = c(y.lim[1], cutoff.perc), col = pointer.col, lwd = 3)
# # ----
# # X Axis
# boldeable.text <- c(100,expression(bold("99")),98, 97, 96, 95)
# axis(side = 1, at = pidents, labels = F, tck = -.03, line = -.5)
# mtext(text = boldeable.text, side = 1, line = -.2, at = pidents, cex = .7)
# # Y Axis
# axis(side = 2, at = y.ticks, labels = F, tck = -.03, line = 0, tck = -.03)
# mtext(text = y.tick.labs, side = 2, at = y.ticks, las = 2, line = .5, cex = .7)
# # Labels
# mtext(text = plot.title, side = 3, line = 1.6, at = 100.2, cex = 1, adj = 1, padj = 1)
# mtext(text = x.label, side = 1, line = 1.1, at = 95, adj = 0, cex = .8)
# mtext(text = y.label, side = 2, line = 2, at = 48, adj = 0, cex = .8)
# 
# #
# # ---- ISME16 POSTER ----
# 
# file.path.read.perc.summs <- "../../poster_mend_unclust/plots/conflict_summary_by_percent_reads.csv"
# output.folder.fig4 <- "~/Dropbox/Trina/8-20-16_ISME16_figures/pident_choice_mendota_unclust.png"
# 
# # define functions for poster ----
# 
# import.summary <- function(FilePath){
#   sumry <- read.csv(file = FilePath, header = FALSE, colClasses = "character")
#   sumry[1,1] <- "pident"
#   sumry[ ,-1] <- apply(X = sumry[ ,-1], MARGIN = 2, FUN = as.numeric)
#   return(sumry)
# }
# 
# trim.to.perc.classified <- function(ConflictSum){
#   perc.class <- ConflictSum[-c(2:6), ]
#   future.colnames <- perc.class[ ,1]
#   perc.class <- perc.class[ ,-1]
#   perc.class <- apply(perc.class, 2, as.numeric)
#   perc.class <- t(perc.class)
#   colnames(perc.class) <- future.colnames
#   return(perc.class)
# }
# 
# # use functions for poster ----
# 
# reads <- import.summary(FilePath = file.path.read.perc.summs)
# 
# reads.class <- trim.to.perc.classified(ConflictSum = reads)
# 
# PercClass = reads.class
# FolderPath = output.folder.fig4
# 
# pidents <- PercClass[ ,1]
# perc.class <- PercClass[ ,2]
# 
# pidents <- pidents[1:6]
# perc.class <- perc.class[1:6]
# 
# plot.title <- c(expression(bold("Percent Classified by")), expression(bold("the FreshTrain")))
# x.label <- expression(bold("Percent Identity Cutoff"))
# y.label <- expression(bold("Reads (%)"))
# y.min <- min(perc.class)
# y.max <- max(perc.class)
# x.min <- min(pidents)
# x.max <- max(pidents)
# 
# 
# png(filename = FolderPath, width = 2.8, height = 3, units = "in", res = 300)
# par(mar = c(2.3,2.65,1.8,.51))
# plot(x = pidents, y = perc.class, type = "n", ann = F, axes = F, ylim = c(y.min, y.max))
# polygon(x = c(pidents[2], pidents[2:4], pidents[4]), y = c(y.min, perc.class[2:4], y.min), xpd=T, border = NA, col = adjustcolor(col = "red", alpha.f = .2))
# lines(x = pidents, y = perc.class, col = "grey", lwd = 5)
# points(x = pidents, y = perc.class, col = "grey", pch = 19, cex = .75)
# lines(x = c(98,98), y = c(y.min, perc.class[3]), col = "red", lwd = 3)
# lines(x = c(x.min, 98), y = c(perc.class[3], perc.class[3]), col = "blue", lwd = 3)
# axis(side = 1, at = pidents, labels = F, tick = TRUE, line = -.5, lwd = 3, xpd = T)
# mtext(text = pidents, side = 1, line = .2, at = pidents - .1, cex = 1.5)
# mtext(text = x.label, side = 1, line = 1.5, cex = 1.5, at = 93.7, adj = 0)
# axis(side = 2, at = c(50,60,70,80), labels = F, tick = T, line = -.5, lwd = 3, xpd =T)
# mtext(text = c(50,60,70,80), at = c(50,60,70,80) + c(.75,0,0,-.75), side = 2, line = .1, cex = 1.5)
# mtext(text = y.label, side = 2, line = 1.2, cex = 1.5)
# mtext(text = plot.title, side = 3, at = c(94, 95.7), adj = 0, line = c(.4,-.6), cex = 1.5, col = "red")
# dev.off()









