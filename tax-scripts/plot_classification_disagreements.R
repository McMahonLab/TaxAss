# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# Do this for the "take4" results (this used a new database version from Trina)
# compare pidents based on # OTUs and # reads impacted.

#####
# Import the Data
#####
pid.results <- read.table("data/8-28-15_pident_cutoffs_results", sep="\t", header=T)


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
