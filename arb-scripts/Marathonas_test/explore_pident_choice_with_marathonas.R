# 2018-3-8 RRR

# Look at how pident choice changes database assignment errors

# Re-Run step 15 on the v4_mara or whichever folder using different pidents.

# Run through main R script to get summary tables or error types exported

# Those tables are the input into this script.

# Do number of seqs incorrectly-assigned to GG/FT change?

# ---- File Paths ----

file.96 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_96/summary_table.csv"
file.97 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_97/summary_table.csv"
file.98 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_98/summary_table.csv"
file.99 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_99/summary_table.csv"
file.100 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_100/summary_table.csv"

output.folder <- "../arb-scripts/Marathonas_test/v4_mara/v4_pident_impact/"

# ---- Just go in order ----

import.summary.table <- function(file.name, pident){
  sum.tab <- read.csv(file = file.name)
  sum.tab <- sum.tab[ ,-(2:5)]
  colnames(sum.tab) <- paste(colnames(sum.tab), pident, sep = ".")
  return(sum.tab)
}

sum.96 <- import.summary.table(file.name = file.96, pident = 96)
sum.97 <- import.summary.table(file.name = file.97, pident = 97)
sum.98 <- import.summary.table(file.name = file.98, pident = 98)
sum.99 <- import.summary.table(file.name = file.99, pident = 99)
sum.100 <- import.summary.table(file.name = file.100, pident = 100)

p.compare <- cbind(sum.96[ ,-1], sum.97[ ,-1], sum.98[ ,-1], sum.99[ ,-1], sum.100[ ,-1])
row.names(p.compare) <- sum.96[ ,1]
p.compare <- as.matrix(p.compare)

lin <- cbind(sum.96[ ,2,drop=F], sum.97[ ,2,drop=F], sum.98[ ,2,drop=F], sum.99[ ,2,drop=F], sum.100[ ,2,drop=F])
row.names(lin) <- sum.100[ ,1]
lin <- as.matrix(lin)

cla <- cbind(sum.96[ ,3,drop=F], sum.97[ ,3,drop=F], sum.98[ ,3,drop=F], sum.99[ ,3,drop=F], sum.100[ ,3,drop=F])
row.names(cla) <- sum.100[ ,1]
cla <- as.matrix(cla)

tri <- cbind(sum.96[ ,4,drop=F], sum.97[ ,4,drop=F], sum.98[ ,4,drop=F], sum.99[ ,4,drop=F], sum.100[ ,4,drop=F])
row.names(tri) <- sum.100[ ,1]
tri <- as.matrix(tri)

make.stacked.bar <- function(lev){ # laxy calls from global env
  label.cex <- 1.1
  axis.cex <- 1.2
  title.cex <- 1.4
  
  col.correct <- adjustcolor(col = "darkgreen")
  col.under <- adjustcolor(col = rainbow(n = 20, v = .8)[4])
  col.wrong <- adjustcolor(col = "darkred")
  
  line.loc <- cumsum(lev[ ,ncol(lev)])
  line.loc <- line.loc - (.5 * lev[ ,ncol(lev)])
  label.loc <- line.loc
  for (n in 2:(nrow(lev))){
    if (n < nrow(lev)){
      sep <- .5 * (label.loc[n] - label.loc[n - 1]) + .5 * (label.loc[n + 1] - label.loc[n])
    }else{
      sep <- (label.loc[n] - label.loc[n - 1])
    }
    
    if (sep < 25){
      label.loc[n] <- label.loc[n] + (25 - sep)
    }
  }
  
  max.val <- sum(lev[ ,1])
  
  par(mar = c(2,4,2.4,12))
  bar.loc <- barplot(lev, col = c(col.correct, col.correct, col.correct, col.under, col.under, col.wrong, col.wrong, col.wrong), axisnames = F, axes = F)
  
  mtext(text = colnames(lev), side = 1, line = .5, at = bar.loc, cex = label.cex, las = 1)
  
  axis(side = 2, at = c(0, max.val), labels = c("",max.val), cex = label.cex)
  axis(side = 2, cex = label.cex)
  mtext(text = "Number of Sequences", side = 2, line = 2.3, cex = axis.cex)
  
  text(x = rep.int(x = 6.3, times = nrow(lev)), y = label.loc, labels = row.names(lev), xpd = T, adj = 0, cex = label.cex)
  for (n in 1:nrow(lev)){
    lines(x = c(6.03,6.27), y = c(line.loc[n], label.loc[n]), xpd = T)
  }
  
  mtext(text = "TaxAss Accuracy", side = 3, line = 1.1, cex = title.cex)
  
} 

# ----

pdf(file = paste(output.folder, "1_tribe_detailed.pdf", sep = "/"), width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
# layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
make.stacked.bar(lev = tri)
dev.off()

pdf(file = paste(output.folder, "1_clade_detailed.pdf", sep = "/"), width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
# layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
make.stacked.bar(lev = cla)
dev.off()

pdf(file = paste(output.folder, "1_lineage_detailed.pdf", sep = "/"), width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
# layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
make.stacked.bar(lev = lin)
dev.off()

arb <- data.frame("Greengenes" = sum(lin[c(3,6), 1]), "FreshTrain" = sum(lin[c(1,2,4,5,7,8), 1]))
arb <- t(arb)
colnames(arb) <- "arb"
arb
gg.split <- cbind("96" = sum(lin[c(3,5), 1]), "97" = sum(lin[c(3,5), 2]), "98" = sum(lin[c(3,5), 3]), "99" = sum(lin[c(3,5), 4]), "100" = sum(lin[c(3,5), 5]))
ft.split <- 285 - gg.split                                                                                                  

sum.split <- cbind(arb, rbind(gg.split,ft.split))
sum.split <- as.matrix(sum.split)

sum.errors <- cbind("96" = sum(tri[c(6,7,8), 1]), "97" = sum(tri[c(6,7,8), 2]),"98" = sum(tri[c(6,7,8), 3]) ,"99" = sum(tri[c(6,7,8), 4]),"100" = sum(tri[c(6,7,8), 5]) )

# ----

make.database.split.bar <- function(){
  label.cex <- 1.1
  axis.cex <- 1.2
  title.cex <- 1.4
  
  col.fw <- adjustcolor(col = "darkgreen", alpha.f = 1)
  col.gg <- adjustcolor(col = "darkblue", alpha.f = 1)
  col.error <- adjustcolor(col = "darkred", alpha.f = 1)
  
  max.val <- sum(sum.split[ ,1])
  
  par(mfrow = c(1,2), mar = c(2,1,2.4,1), oma = c(0,3,0,0))
  
  bar.loc <- barplot(sum.split, col = c(col.gg, col.fw), axisnames = F, axes = F)
  
  mtext(text = colnames(sum.split), side = 1, line = .5, at = bar.loc, cex = label.cex, las = 1)
  
  axis(side = 2, at = c(0, max.val), labels = c("",max.val), cex = label.cex)
  axis(side = 2, cex = label.cex)
  mtext(text = "Number of Sequences", side = 2, line = 2.3, cex = axis.cex)
  
  
  barplot(sum.errors, col = "red")
  
  
  
  mtext(text = "TaxAss Accuracy", side = 3, line = -1.3, cex = title.cex, outer = T)
  
  
  
  
  
  
  
  
  
}





