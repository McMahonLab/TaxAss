# Marathonas Validation Figure
# RRR

# ---- paths ----
file.path.summary.table.v4 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results/summary_table.rds"
file.path.summary.table.v3v4 <- "../arb-scripts/Marathonas_test/v3v4_mara/v3v4_results/summary_table.rds"
file.path.summary.table.v4v5 <- "../arb-scripts/Marathonas_test/v4v5_mara/v4v5_results/summary_table.rds"


# file.path.created.plot <- "~/Dropbox/PhD/TaxAss_manuscript/draft 7/re-submission_figures/mara_validation-referencelabel.pdf"
# file.path.created.table <- "~/Dropbox/PhD/Write\ It/draft\ 7/res-submission_figures/mara_validation.csv"

# ---- functions ----

make.stacked.bar.for.manuscript <- function(){ # laxy calls from global env
  # evenly spaced bar labels to match table
  # colored differently ("under" is correct, yellow is db error, red is known error)
  
  label.cex <- 1.2
  axis.cex <- 1.1
  title.cex <- 1.2
  excel.cex <- 1.7
  
  col.correct <- adjustcolor(col = "darkgreen")
  col.under <- adjustcolor(col = rainbow(n = 20, v = .7)[4])
  col.wrong <- adjustcolor(col = "darkred")
  col.vector <- c(col.correct, col.correct, col.correct, col.correct, col.under, col.wrong, col.wrong, col.wrong)
  
  line.loc <- cumsum(results[ ,7])
  line.loc <- line.loc - (.5 * results[ ,7])
  label.loc <- seq(from = 20, to = 275, along.with = line.loc)
  
  max.val <- sum(results[ ,1])
  
  par(mar = c(2,4,2.5,16))
  bar.loc <- barplot(results[ ,-(1:4)], col = col.vector, axisnames = F, axes = F)
  
  mtext(text = colnames(results)[-(1:4)], side = 1, line = .5, at = bar.loc - c(.4,-.05,-.15), cex = label.cex)
  
  axis(side = 2, at = c(0, max.val), labels = F)
  tic.loc <- axis(side = 2, labels = F)
  label.scooch <- c(tic.loc + c(4, rep.int(x = 0, times = length(tic.loc) - 1)) , max.val + 15) # 0 up, 285 up
  mtext(text = c(tic.loc,max.val), side = 2, line = .5, outer = F, at = label.scooch, cex = axis.cex)
  mtext(text = "Number of Sequences", side = 2, line = 2.4, cex = label.cex)
  
  text(x = rep.int(x = 4.2, times = nrow(results)), y = label.loc, labels = row.names(results), xpd = NA, adj = 0, cex = excel.cex, col = col.vector)
  for (n in 1:nrow(results)){
    lines(x = c(3.65,4.15), y = c(line.loc[n], label.loc[n]), xpd = T, col = col.vector[n], lwd = 1.8)
  }
  
  return(label.loc)
} 

make.example.table <- function(){ # all lazy calls to global env
  
  title.cex <- 1.2
  label.cex <- 1.2
  table.cex <- 1.4
  
  half.spacing <- (word.heights[2] - word.heights[1]) * 1/2
  box.tops <- word.heights + half.spacing
  box.bottoms <- word.heights - half.spacing
  box.lefts <- 9.2
  box.middles <- 12.4
  box.rights <- 15.6
  
  col.correct <- adjustcolor(col = "darkgreen", alpha.f = .3)
  col.under <- adjustcolor(col = rainbow(n = 20, v = .7)[4], alpha.f = .5)
  col.wrong <- adjustcolor(col = "darkred", alpha.f = .3)
  col.vector <- c(col.correct, col.correct, col.correct, col.correct, col.under, col.wrong, col.wrong, col.wrong)
  
  # these are actual examples I found in the results, not just made-up examples
  left.text.vect <- c("acI-A1", "bacI-unclassified", "n/a", "acI-A", "acI-B1", "n/a", "betI-A","acI-C")
  right.text.vect <- c("acI-A1","bacI-unclassified", "Microcystaceae","acI", "Actinobacteria", "LD19","betI-B","acI-C1")
  
  # left boxes
  for (b in 1:length(word.heights)){
    rect(xleft = box.lefts, xright = box.middles, ybottom = box.bottoms[b], ytop = box.tops[b], xpd = NA, col = col.vector[b])
    text(x = box.lefts + .2, y = word.heights[b], labels = left.text.vect[b], xpd = NA, adj = 0, cex = table.cex)
  }
  # right boxes
  for (b in 1:length(word.heights)){
    rect(xleft = box.middles, xright = box.rights, ybottom = box.bottoms[b], ytop = box.tops[b], xpd = NA, col = col.vector[b])
    text(x = box.middles + .2, y = word.heights[b], labels = right.text.vect[b], xpd = NA, adj = 0, cex = table.cex)
  }
  
  middle.left <- box.lefts + (box.middles - box.lefts) * (1/2)
  middle.right <- box.middles + (box.rights - box.middles) * (1/2)
  mtext(text = c("Reference","TaxAss"), side = 1, line = .5, at = c(middle.left, middle.right), cex = label.cex)
  mtext(text = "Examples of Each Category:", side = 3, line = 1, at = box.middles, cex = label.cex)
  
}

make.csv.table.for.manuscript <- function(v4,v5,v3){
  x <- cbind(v4[ ,5:7], v5[ ,5:7], v3[ ,5:7])
  x <- x[nrow(x):1, ] # flip to match plot order
  x <- rbind(colnames(x), x)
  x <- rbind(rep(c("Simulated V4 Tags", "Simulated V4-V5 Tags", "Simulated V3-V4 Tags"), each = 3), x)
  return(x)
}

# ---- go ----

results <- readRDS(file = file.path.summary.table.v4)

pdf(file = file.path.created.plot, width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
word.heights <- make.stacked.bar.for.manuscript()
make.example.table()
dev.off()

v4 <- readRDS(file = file.path.summary.table.v4)
v5 <- readRDS(file = file.path.summary.table.v4v5)
v3 <- readRDS(file = file.path.summary.table.v3v4)
supp.table <- make.csv.table.for.manuscript(v4 = v4, v5 = v5, v3 = v3)
write.csv(x = supp.table, file = file.path.created.table)
# in excel:
# delete top row, merge primer labels into 1 cell each, add cell outlines
# change font to helvetica, center text, color background to match figure with eyedropper
# save as pdf in excel
# open powerpoint and insert the pdf, crop off white space
# resize the image using the format picture pane- lock aspect ratio and set width to 6.87
# right click and re-save the image as a pdf



