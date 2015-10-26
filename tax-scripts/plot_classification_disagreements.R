# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# Do this for the "take4" results (this used a new database version from Trina)
# compare pidents based on # OTUs and # reads impacted.

#####
# Receive arguments from terminal command line
#####

# userprefs <- commandArgs(trailingOnly = TRUE)
# fw.plus.gg.tax.file.path <- userprefs[1]
# gg.only.tax.file.path <- userprefs[2]
# results.folder.path <- userprefs[3]
# taxonomy.bootstrap.cutoff <- userprefs[4]
# fw.seq.ids.file.path <- userprefs[5]

# Do some sort of loop, the user enters all the folder paths, and then 
# for (p in 1:length(supplied arguments)) {create variables for each path name}

example.user.args <- c("../../take4/compare_percID-93_to_gg-only/conflicts_summary.csv", 93, 
                       "../../take4/compare_percID-94_to_gg-only/conflicts_summary.csv", 94,
                       "../../take4/compare_percID-95_to_gg-only/conflicts_summary.csv", 95,
                       "../../take4/compare_percID-96_to_gg-only/conflicts_summary.csv", 96,
                       "../../take4/compare_percID-97_to_gg-only/conflicts_summary.csv", 97,
                       "../../take4/compare_percID-98_to_gg-only/conflicts_summary.csv", 98,
                       "../../take4/compare_percID-99_to_gg-only/conflicts_summary.csv", 99,
                       "../../take4/compare_percID-100_to_gg-only/conflicts_summary.csv", 100)



#####
# Define functions to import the data
#####

import.conflict.nums.data <- function(FilePath){
  nums <- read.csv(FilePath)
  nums <- nums[,2:3]
  return(nums)
}

import.all.conflict.summaries.into.list <- function(UserArgs){
  example.user.args <- UserArgs
  
  # Import them into a list format
  mismatches.list <- list(NULL)
  counter <- 1
  for (p in seq(from = 1, to = length(example.user.args), by = 2)){
    mismatches.list[[counter]] <- import.conflict.nums.data(FilePath = example.user.args[p])
    names(mismatches.list)[counter] <- example.user.args[p+1]
    counter <- counter + 1
  }
  
  # Rearrange them into a matrix format
  mismatches.matrix <- matrix(0, nrow = nrow(mismatches.list[[1]]), ncol = length(names(mismatches.list)))
  row.names(mismatches.matrix) <- mismatches.list[[1]]$TaxaLevel
  colnames(mismatches.matrix) <- names(mismatches.list)
  for (c in 1:ncol(mismatches.matrix)){
    mismatches.matrix[,c] <- mismatches.list[[c]][,2]
  }
  
  return(mismatches.matrix)
}


#####
# Define functions to analyze the data
#####


#####
# Define functions to plot the data
#####

plot.num.forced.otus <- function(ConflictsSummaryTables, y.axis.limit=0){
  # remove the last row of number FW sequences, that's not needed for this plot.
  mismatches <- ConflictsSummaryTables[1:(nrow(ConflictsSummaryTables)-1),]
  pidents <- colnames(mismatches)
  pidents <- as.numeric(pidents)
  if (y.axis.limit == 0){
    ymax <- max(mismatches)
  }else{
    ymax <- y.axis.limit
  }
  
  # Set up an empty plot
  plot(x = 0, type = "n", ylim = c(0,ymax), xlim = c(min(pidents),max(pidents)),
       main = "How many classification disagreements are there between FW and GG?
         (at different pidents and taxonomic resolutions,\nhow well do gg and fw classifications agree?",
       ylab = "Classification Disagreements- Total number of OTUs", xlab = "\"full length\" pident cutoff (like ANI)")
  
  # Fill Plot with beautiful data
  color <- rainbow(nrow(mismatches))
  for (r in 1:nrow(mismatches)){
    lines(x = pidents, y = mismatches[r,], col = color[r], lwd = 4)
    points(x = pidents, y = mismatches[r,], col = color[r], pch = 19, cex =1.3)
  }
  legend("center",legend = row.names(mismatches), text.col = color, cex=1.3)
}

plot.num.classified.outs <- function(ConflictsSummaryTables){
  num.fw <- ConflictsSummaryTables[nrow(ConflictsSummaryTables),]
  pidents <- colnames(ConflictsSummaryTables)
  pidents <- as.numeric(pidents)
  
  plot(x = pidents, y = num.fw)
}
  
#####
# Use Functions
#####

all.summaries <- import.all.conflict.summaries.into.list(UserArgs = example.user.args)

plot.num.forced.otus(ConflictsSummaryTables = all.summaries)
plot.num.forced.otus(ConflictsSummaryTables = all.summaries, y.axis.limit = 200)
  

  
  






#####
# Make a plot comparing the classification conflicts between pidents
#####

# function to repeat plot easily
plot.forcing <- function(TaxLevel, color){
  t <- TaxLevel
  num.forced <- pid.results[,c(1,3+t)]
  points(num.forced,type = "l", col = color, lwd=5)
}

# Report stats from the old classification pipeline:
old.results <- c(0,34,783,636,1326)
names(old.results) <- colnames(pid.results)[4:8]

# set up plot dimensions, titles
plot(x=0, type="n", xlim= c(75,100), ylim=c(0,500), axes=F, ann=F)
axis(side = 1, lwd=2, lwd.ticks=2, pos=-10, labels =F) #just the axis+ticks
axis(side = 1, padj=0, pos=0, labels =T, tick = F, cex.axis=1.5) #labels separately to scoot them closer
axis(side = 2, lwd=2, lwd.ticks = 2, pos=75, labels=F)
axis(side = 2, tick=F, pos=75.5, labels=T, cex.axis=1.5)

# Set up colors
colors <- rainbow(n = 20)[c(1,3,10,15)]

# Add lines  (for t= 1-kingdom, 2-phylum, 3-class, 4-order, 5-lineage)
for (t in 2:5){
  plot.forcing(TaxLevel = t, color = colors[t-1])
}
lines(x=cbind(c(96,96),c(-20,475)), lwd=3, col= rainbow(20)[17])
lines(x=cbind(c(69,92),c(old.results[2],old.results[2])), lwd=3, col= rainbow(20)[17])

# Do titles with mtext to make them more compact
mtext(expression(bold("Percent ID cutoff (%)")), side = 1, xpd=T, cex=1.5, line= 2, at= 87.5)
mtext(expression(bold("Classification Disagreements (Total #)")), side = 2, xpd=T, cex=1.5, line= 1.8, at= 250)
mtext(expression(bold("Choosing an Appropriate Cutoff")), side=3, xpd=T, cex=2, line= 2, at= 69, adj = 0)
mtext(expression(bold("Phylum")), side=3, , xpd=T, cex=1.5, line=0, at=75, adj=0, col=colors[1])
mtext(expression(bold("Class")), side=3, , xpd=T, cex=1.5, line=-.5, at=86, adj=0, col=colors[2])
mtext(expression(bold("Order")), side=3, , xpd=T, cex=1.5, line=.5, at=89, adj=0, col=colors[3])
mtext(expression(bold("Lineage")), side=3, , xpd=T, cex=1.5, line=0, at=94, adj=0, col=colors[4])
mtext(expression(bold("Chosen")), side=3, , xpd=T, cex=1.3, line=-1, at=96.3, adj=0, padj=1, col=rainbow(20)[17])
mtext(expression(bold("Cutoff:")), side=3, , xpd=T, cex=1.3, line=-2, at=96.3, adj=0, padj=1, col=rainbow(20)[17])
mtext(expression(bold("96%")), side=3, , xpd=T, cex=1.5, line=-3, at=96.6, adj=0, padj=1, col=rainbow(20)[17])
mtext(expression(bold("34 Phyla")), side=1, , xpd=T, cex=1.3, line=-6, at=76, adj=0, padj=1, col=rainbow(20)[17])
mtext(expression(bold("Misclassified")), side=1, , xpd=T, cex=1.3, line=-5, at=76, adj=0, padj=1, col=rainbow(20)[17])
mtext(expression(bold("Previously")), side=1, , xpd=T, cex=1.3, line=-4, at=76, adj=0, padj=1, col=rainbow(20)[17])

#####
# Make a plot showing the proportion of sequences classified in FW database
#####

# Calculate the number classified in FW as a a percent- 
# note there were 21397 OTUs after deblurring + mothur QC
perc.above <- cbind(pid.results$pident, pid.results$above.pident/21397*100)
colnames(perc.above) <- c("pident", "perc.FW") 

# Add in the old pipeline for comparison- 2413 sent to FW out of the same 21397 total
perc.fw.old <- 2413/21397*100

# set up plot dimensions, titles
plot(0, type="n", xlim=c(70,100), ylim=c(0,60), axes=F, ann=F)
axis(side = 1, lwd=2, lwd.ticks=2, pos=-2, labels =F) 
axis(side = 1, padj=0, pos=-1.5, labels =T, tick = F, cex.axis=1.5) 
axis(side = 2, lwd=2, lwd.ticks = 2, pos=68.5, labels=F)
axis(side = 2, tick=F, pos=69, labels=T, cex.axis=1.5)

# Plot the data
points(perc.above, col= rainbow(20)[8], pch=19, cex=1.2)
lines(perc.above, col= rainbow(20)[8], lwd=5)
lines(x=cbind(c(96,96),c(-2,30)), lwd=3, col= rainbow(20)[17])
lines(x=cbind(c(68.5,98),c(11.277,11.277)), lwd=3, col= rainbow(20)[17],xpd=T)

# Do titles with mtext to make them more compact
mtext(expression(bold("Percent ID cutoff (%)")), side = 1, xpd=T, cex=1.5, line= 2, at= 85, lwd=4)
mtext(expression(bold("Sequences Classified with FW (%)")), side = 2, xpd=T, cex=1.5, line= 2.4, at= 30, lwd=4)
mtext(expression(bold("Contribution of FW Taxonomy")), side=3, xpd=T, cex=2, line= 2, at= 65, adj = 0, lwd=5)
mtext(expression(bold("to Total Classifications")), side=3, xpd=T, cex=2, line= 0, at= 68, adj = 0, lwd=4)
mtext(expression(bold("At 96% blast cutoff")), side=1, , xpd=T, cex=1.3, line=-13, at=88, adj=0, padj=1, col=rainbow(20)[17], lwd=4)
mtext(expression(bold("4% classified by FW")), side=1, , xpd=T, cex=1.3, line=-11.9, at=88, adj=0, padj=1, col=rainbow(20)[17], lwd=4)
mtext(expression(bold("In old pipeline")), side=1, , xpd=T, cex=1.3, line=-7.3, at=70, adj=0, padj=1, col=rainbow(20)[17], lwd=4)
mtext(expression(bold("11% classified by FW")), side=1, , xpd=T, cex=1.3, line=-6.2, at=70, adj=0, padj=1, col=rainbow(20)[17], lwd=4)

#####
# Use this to choose colors
plot(1:20, col = rainbow(20)[1:20], cex=5, pch=19)
points(1:20,pch=as.character(c(1:9,0)))
