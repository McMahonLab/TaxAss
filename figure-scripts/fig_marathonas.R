# Marathonas Validation Figure
# RRR

# ---- paths ----
file.path.summary.table <- "../arb-scripts/Marathonas_test/v4_mara/v4_results//summary_table.rds"

file.path.created.plot <- "~/Desktop/mara.pdf"

# ---- functions ----

make.stacked.bar.for.manuscript <- function(){ # laxy calls from global env
  # evenly spaced bar labels to match table
  # colored differently ("under" is correct, yellow is db error, red is known error)
  
  label.cex <- 1.2
  axis.cex <- 1.1
  title.cex <- 1.2
  excel.cex <- 1.7
  
  col.correct <- adjustcolor(col = "darkgreen")
  col.under <- adjustcolor(col = rainbow(n = 20, v = .8)[4])
  col.wrong <- adjustcolor(col = "darkred")
  col.vector <- c(col.correct, col.correct, col.correct, col.correct, col.under, col.under, col.wrong, col.wrong)
  
  line.loc <- cumsum(results[ ,7])
  line.loc <- line.loc - (.5 * results[ ,7])
  label.loc <- seq(from = 20, to = 275, along.with = line.loc)
  
  max.val <- sum(results[ ,1])
  
  par(mar = c(2,4,2.5,16))
  bar.loc <- barplot(results[ ,-(1:4)], col = col.vector, axisnames = F, axes = F)
  
  mtext(text = colnames(results)[-(1:4)], side = 1, line = .5, at = bar.loc - c(.4,0,-.1), cex = label.cex)
  
  axis(side = 2, at = c(0, max.val), labels = F)
  tic.loc <- axis(side = 2, labels = F)
  label.scooch <- c(tic.loc + c(4, rep.int(x = 0, times = length(tic.loc) - 1)) , max.val + 15) # 0 up, 285 up
  mtext(text = c(tic.loc,max.val), side = 2, line = .5, outer = F, at = label.scooch, cex = axis.cex)
  mtext(text = "Number of Sequences", side = 2, line = 2.3, cex = label.cex)
  
  text(x = rep.int(x = 4.2, times = nrow(results)), y = label.loc, labels = row.names(results), xpd = NA, adj = 0, cex = excel.cex, col = col.vector)
  for (n in 1:nrow(results)){
    lines(x = c(3.62,4.15), y = c(line.loc[n], label.loc[n]), xpd = T)
  }
  
  # mtext(text = "TaxAss Performance on Simulated Tags", side = 3, line = 3.2, cex = title.cex)
  
  return(label.loc)
} 

make.example.table <- function(){ # all lazy calls to global env
  
  title.cex <- 1.2
  label.cex <- 1.2
  
  half.spacing <- (word.heights[2] - word.heights[1]) * 1/2
  box.tops <- word.heights + half.spacing
  box.bottoms <- word.heights - half.spacing
  box.lefts <- 9.2
  box.middles <- 12.4
  box.rights <- 15.6
  
  # left boxes
  for (b in 1:length(word.heights)){
    rect(xleft = box.lefts, xright = box.middles, ybottom = box.bottoms[b], ytop = box.tops[b], xpd = NA)
  }
  # right boxes
  for (b in 1:length(word.heights)){
    rect(xleft = box.middles, xright = box.rights, ybottom = box.bottoms[b], ytop = box.tops[b], xpd = NA)
  }
  
  middle.left <- box.lefts + (box.middles - box.lefts) * (1/2)
  middle.right <- box.middles + (box.rights - box.middles) * (1/2)
  mtext(text = c("ARB","TaxAss"), side = 1, line = .5, at = c(middle.left, middle.right), cex = label.cex)
  mtext(text = "Examples of Each Category:", side = 3, line = 1, at = box.middles, cex = label.cex)
  
}


# ---- go ----

results <- readRDS(file = file.path.summary.table)

pdf(file = file.path.created.plot, width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
word.heights <- make.stacked.bar.for.manuscript()
make.example.table()
dev.off()
