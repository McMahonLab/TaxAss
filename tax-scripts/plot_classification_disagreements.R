# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# compare pidents based on # OTUs and # reads impacted.

#####
# Receive arguments from terminal command line
#####

# userprefs <- commandArgs(trailingOnly = TRUE)
# fw.plus.gg.tax.file.path <- userprefs[1]
# gg.only.tax.file.path <- userprefs[2]
# fw.seq.ids.file.path <- userprefs[5]
# results.folder.path <- userprefs[3]
# taxonomy.bootstrap.cutoff <- userprefs[4]

# Do some sort of loop, the user enters all the folder paths, and then 
# for (p in 1:length(supplied arguments)) {create variables for each path name}

userprefs <- c("../../take5/otus.abund",
               "../../take5/conflicts_70", 70,
               "../../take5/conflicts_80", 80,
               "../../take5/conflicts_85", 85,
               "../../take5/conflicts_90", 90,
               "../../take5/conflicts_91", 91,
               "../../take5/conflicts_92", 92,
               "../../take5/conflicts_93", 93,
               "../../take5/conflicts_94", 94,
               "../../take5/conflicts_95", 95,
               "../../take5/conflicts_96", 96,
               "../../take5/conflicts_97", 97,
               "../../take5/conflicts_98", 98,
               "../../take5/conflicts_99", 99,
               "../../take5/conflicts_100", 100)

#####
# Define functions to import the data
#####

# This funciton is used inside the import.all.conflict.summariers function
import.conflict.nums.data <- function(FilePath){
  nums <- read.csv(FilePath)
  nums <- nums[,2:3]
  return(nums)
}

import.all.conflict.summaries <- function(UserArgs){
  # the first argument is the otu rel abund table, so ignore that here
  user.args <- UserArgs[-1]
  
  # Import them into a list format
  mismatches.list <- list(NULL)
  counter <- 1
  for (p in seq(from = 1, to = length(user.args), by = 2)){
    mismatches.list[[counter]] <- import.conflict.nums.data(FilePath = user.args[p])
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
  mismatches <- ConflictsSummaryTables[1:(nrow(ConflictsSummaryTables)-2),]
  pidents <- colnames(mismatches)
  pidents <- as.numeric(pidents)
  if (y.axis.limit == 0){
    ymax <- max(mismatches)
  }else{
    ymax <- y.axis.limit
  }
  
  # Set up an empty plot
  plot(x = 0, type = "n", ylim = c(0,ymax), xlim = c(min(pidents),max(pidents)),
       main = "Classification Disagreements Between FW and GG\n(Are we forcing OTUs into our favorite FW groups?)",
       ylab = "Classification Disagreements- Total number of OTUs", xlab = "\"full length\" pident cutoff (similar to ANI)")
  
  # Fill Plot with beautiful data
  color <- rainbow(nrow(mismatches))
  for (r in 1:nrow(mismatches)){
    lines(x = pidents, y = mismatches[r,], col = color[r], lwd = 4)
    points(x = pidents, y = mismatches[r,], col = color[r], pch = 19, cex =1.3)
  }
  legend("center",legend = row.names(mismatches), text.col = color, cex=1.3)
}

plot.num.classified.outs <- function(ConflictsSummaryTables){
  num.fw <- ConflictsSummaryTables[nrow(ConflictsSummaryTables)-1,]
  tot.seq <- ConflictsSummaryTables[nrow(ConflictsSummaryTables),1]
  perc.fw <- num.fw / tot.seq * 100
  pidents <- colnames(ConflictsSummaryTables)
  pidents <- as.numeric(pidents)
  
  # Set up and empty plot
  plot(x = pidents, y = perc.fw, type = "n", main = "Contribution of FW Taxonomy to Total Classifications
       (How many of the OTUs are we calling freshwater?)",
       xlab = "\"Full Length\" percent identity (similar to ANI)", ylab = "Percent total OTUs Classified by FW (%)")
  
  # Fill plot with beautiful data
  lines(x = pidents, y = perc.fw, col = "lightsalmon", lwd = 1.5)
  points(x = pidents, y = perc.fw, col = "lightsalmon", pch = 19, cex = 1.3)
}
  

#####
# Use Functions
#####

all.summaries <- import.all.conflict.summaries.into.list(UserArgs = example.user.args)

plot.num.forced.otus(ConflictsSummaryTables = all.summaries)
plot.num.forced.otus(ConflictsSummaryTables = all.summaries, y.axis.limit = 5)

plot.num.classified.outs(ConflictsSummaryTables = all.summaries)








