# RRR 8-15-16 ----

# fig 3 demonstrates that the workflow is better than using only an 
# ecosystem-specific database by showing the extent of forcing.

# fig 3a shows how forcing can mess up your fine-level taxa analysis
# fig 3b shows how forcing can skew coarse-level taxa analyses

# How am I defining forcing?
# 1. classify everything with custom database
# 2. look at the freshwater classifications of everything that the workflow put in greengenes
# 3. ignore everything that's unclassified in freshwater classifications
# 4. find all the classifications that are different from the workflow classifications-
#    these are the ones I call "forced" into the wrong classification.

# for fig 3a: baseline is freshwater rank abundance
# 1. make a rank abundance plot of the top lineages as defined by the custom-only taxonomy
#    (this leaves out "unclassified" from the taxa names)
# 2.color the bars red for the OTUs that were not classified that way by the workflow

# for fig 3b: baseline is workflow rank abundance
# 1. make a rank abund of the top phyla in freshwater-only classifications
# 2. make a rank abund of the top phyla in workflow classifications
# 3. compare how different freshwater only is from workflow

# from plot class dis:

# make stacked bar
for (t in PlottingLevels){
  # long names go off the plot
  taxa.names <- sub(pattern = ".*__", replacement = "", x = top.taxa[[t]][ ,t])
  taxa.names <- substr(x = taxa.names, start = 1, stop = 20)
  
  # make the y axis not have crazy big numbers on it, put them in rounded percents
  max.bar <- max(top.taxa[[t]][ ,(t+1)])
  y.axis.ticks <- c(0, max.bar * (1/4), max.bar * (1/2), max.bar * (3/4), max.bar)
  y.axis.labels <- round(x = y.axis.ticks / tot.reads * 100, digits = 0)
  
  # generate files of the plots
  png(filename = paste(plots.folder.path, "/", t, "_Forcing_of_top_taxa-", names(top.taxa)[t], ".png", sep = ""), 
      width = 7, height = 5, units = "in", res = 100)
  par(mar = c(10,5,5,2))
  barplot(height = t(as.matrix(top.taxa[[t]][ , (t + 3):(t + 2)])), col = c("grey", "red"),
          main = paste("Forcing at the ", names(top.taxa)[t], " level\nwhen using only the custom database", sep = ""),
          names.arg = taxa.names, legend.text = c("correct","forced"), args.legend = list(x = "topright", bty = "n"),
          las = 2, xpd = TRUE, axes = FALSE, ylab = "Percent of total reads", cex.lab = 1.2)
  axis(side = 2, at = y.axis.ticks, labels = y.axis.labels)
  unnecessary.message <- dev.off()
}