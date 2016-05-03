# Script written by Robin Rohwer

# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# compare pidents based on # OTUs and # reads impacted.

# This script is step 14 of the taxonomy assignment workflow, used to decide on pident cutoff. It takes arguments from the command line.
# The number of arguments can be variable, but they must be in order
# The variable part is that more pidents can be added, as long as they continue the pattern folder number folder number
# This script is also used in step 16 optionally to plot the forcing that would have occured if you only used FW

# Terminal command line syntax:
# Rscript plot_classification_disagreements.R otus.abund plots conflicts_94 ids.above.94 94 conflicts_96 ids.above.96 96 conflicts_98 ids.above.98 98 ...

#####
# Receive arguments from terminal command line
#####

userprefs <- commandArgs(trailingOnly = TRUE)

# FOR CHOOSING CUTOFF:
userprefs <- c("../../take17/otus.abund",
               "../../take17/plots",
               "regular",
               "regular",
               "../../take17/conflicts_95",
               "../../take17/ids.above.95",
               95,
               "../../take17/conflicts_96",
               "../../take17/ids.above.96",
               96,
               "../../take17/conflicts_97",
               "../../take17/ids.above.97",
               97,
               "../../take17/conflicts_98",
               "../../take17/ids.above.98",
               98,
               "../../take17/conflicts_99",
               "../../take17/ids.above.99",
               99,
               "../../take17/conflicts_100",
               "../../take17/ids.above.100",
               100)
# FOR PLOTTING FORCING (**add this in)


# in case you want to add the db baseline conflict back to the plots, need to specify this path below
# and un-comment the plotting calls that use it at the end of the script.

# db.conflicts.folder.path <- "file path to conflicts_database"
otu.table.path <- userprefs[1]
plots.folder.path <- userprefs[2]
forcing.folder.path <- userprefs[3]
forced.taxonomy.file <-userprefs[4]
rest.of.arguments <- userprefs[-(1:4)]
if (length(rest.of.arguments) > 0){
  pident.folders <- rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)]
  ids.file.paths <- rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)+1]
  pident.values <- as.numeric(rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)+2])
}

# this is automatically exported into the working directory when this script is run normally
seqID.reads.file.path <- "total.reads.per.seqID.csv"


#####
# Define functions to import and process the data
#####

# import all the conflict summary files from each folder and compile them into a matrix
import.all.conflict.summaries <- function(ConflictFolders, PidentsUsed){
  pident.folders <- ConflictFolders
  pident.values <- PidentsUsed
  
  # Import conflict summary files into a list format
  mismatches.list <- list(NULL)
  counter <- 1
  for (p in 1:length(pident.folders)){
    mismatches.list[[counter]] <- read.csv(file = paste(pident.folders[p], "/conflicts_summary.csv", sep = ""))
    names(mismatches.list)[counter] <- as.character(pident.values[p])
    counter <- counter + 1
  }
  
  # Rearrange them into a matrix format
  mismatches.matrix <- matrix(0, nrow = nrow(mismatches.list[[1]]), ncol = length(names(mismatches.list)))
  row.names(mismatches.matrix) <- mismatches.list[[1]]$TaxaLevel
  colnames(mismatches.matrix) <- names(mismatches.list)
  for (c in 1:ncol(mismatches.matrix)){
    mismatches.matrix[ ,c] <- mismatches.list[[c]][ ,2]
  }
  
  return(mismatches.matrix)
}

# import the database conflicts summary file
import.database.conflicts <- function(DatabaseFolder){
  db.conflicts.folder.path <- DatabaseFolder
  db.conflicts <- read.csv(file = paste(db.conflicts.folder.path, "/conflicts_summary.csv", sep =""), stringsAsFactors = FALSE)
  # the sequences totals are not needed here
  return(db.conflicts[1:5,])
}

# import the optional plotting step's forcing conflicts summary file
import.forcing.conflicts <- function(ForcingFolder){
  forcing.folder.path <- ForcingFolder
  forcing.conflicts <- read.csv(file = paste(forcing.folder.path, "/conflicts_summary.csv", sep = ""), stringsAsFactors = FALSE)
  forcing.conflicts.matrix <- as.matrix(forcing.conflicts[ ,2, drop = FALSE])
  row.names(forcing.conflicts.matrix) <- forcing.conflicts[ ,1]
  # remove total that should have been GG-classified (confusingly labelled "numFWseqs")
  forcing.conflicts.matrix <- forcing.conflicts.matrix[-7,1, drop = FALSE]
  return(forcing.conflicts.matrix)
}

# import seqID.reads variable for the forcing plot (it's generated w/ the regular plotting)
import.seqID.reads <- function(FilePath){
  seqID.reads.file.path <- FilePath
  
  seqID.reads <- read.csv(file = seqID.reads.file.path, colClasses = "character")
  seqID.reads[ ,2] <- as.numeric(seqID.reads[ ,2])
  
  return(seqID.reads)
}

# import the OTU table and then pull out just the total reads for each seqID
import.and.reformat.otu.table <- function(OTUtable){
  otu.table.path <- OTUtable
  
  otus <- read.table(file = otu.table.path, sep = "\t", header = TRUE, colClasses = "character")
  samples <- apply(X = otus[ ,-1], MARGIN = 2, FUN = as.numeric)
  seqID.reads <- data.frame(seqID = otus[ ,1], reads = rowSums(samples), stringsAsFactors = FALSE)
  
  return(seqID.reads)
}

# import the ids.above. files that list seqIDs above the pident cutoff
# This is used inside the add.totals.to.read.summaries() function
import.ids.above <- function(FilePaths, PidentsUsed){
  ids.file.paths <- FilePaths
  pident.values <- PidentsUsed
  
  fw.ids <- list(NULL)
  for (p in 1:length(pident.values)){
    # fw.ids[[p]] <- scan(file = ids.file.paths[p])
    
    fw.ids[[p]] <- read.table(file = ids.file.paths[p], colClasses = "character")
    fw.ids[[p]] <- fw.ids[[p]][ ,1]
    
    names(fw.ids)[p] <- pident.values[p]
  }
  
  return(fw.ids)
}

# import files comparing conflicts at each level, record the seqIDs into a list of lists
# structure: outer list- each pident, inner lists- seqIDs at each taxa level
get.conflict.seqIDs <- function(ConflictsFolders, PidentsUsed){
  pident.folders <- ConflictsFolders
  pident.values <- PidentsUsed
  
  # forcing plots go to clade level
  if (pident.values[1] == "forcing"){
    num.taxa.levels <- 6
    conflict.ids.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL, Clade = NULL)
  }else{
    num.taxa.levels <- 5
    conflict.ids.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL)
  }
  
  all.pidents <- list(NULL)
  # for each pident folder
  for (p in 1:length(pident.values)){
    all.files <- list.files(pident.folders[p])
    conflict.ids <- conflict.ids.start
    
    # for each taxonomy file
    for (t in 1:num.taxa.levels){
      conflict.ids[[t]] <- read.csv(file = paste(pident.folders[p], "/", all.files[t], sep = ""), header = TRUE, colClasses = "character")
    }
    
    all.pidents[[p]] <- conflict.ids
    names(all.pidents)[p] <- pident.values[p]
  }
  return(all.pidents)
}

# create a list of total reads that matches the structure of the list of conflicting seqIDs
find.reads.per.seqID <- function(ReadsTable, ConflictsList, Forcing = FALSE){
  otu.reads <- ReadsTable
  conflict.ids <- ConflictsList
  
  # forcing plots go to clade level
  if (Forcing == TRUE){
    num.taxa.levels <- 6
    taxa.list.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL, Clade = NULL)
  }else{
    num.taxa.levels <- 5
    taxa.list.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL)
  }
  
  # for each upper list level (pident folder)
  pidents.list <- list(NULL)
  for (p in 1:length(conflict.ids)){
  
    # for each inner list level (taxonomy level)
    taxa.list <- taxa.list.start
    for (t in 1:num.taxa.levels){
      # for each seqID
      if (length(conflict.ids[[p]][[t]]) > 0){
        for (s in 1:length(conflict.ids[[p]][[t]])){
          index <- which(otu.reads$seqID == conflict.ids[[p]][[t]][s])
          taxa.list[[t]] <- c(taxa.list[[t]], otu.reads[index,2])
        }
      }else taxa.list[[t]] <- 0
    }
    # assign element p of the outer list the inner list it contains
    pidents.list[[p]] <- taxa.list
    names(pidents.list)[p] <- names(conflict.ids)[p]
  }
  return(pidents.list)
}

# collapse the list of total reads per seqID into a table of total reads: taxa level x pident
generate.summary.table.of.reads <- function(ReadsList, Forcing = FALSE){
  reads.list <- ReadsList
  
  # forcing plots go to clade level
  if (Forcing == TRUE){
    num.taxa.levels <- 6
  }else{
    num.taxa.levels <- 5
  }
  
  # Set up empty matrix to fill
  reads.summary <- matrix(0, nrow = num.taxa.levels, ncol = length(reads.list))
  row.names(reads.summary) <- names(reads.list[[1]])
  colnames(reads.summary) <- names(reads.list)
  
  #For each outer list's list (pident)
  for (p in 1:length(reads.list)){
    
    # For each inner list's vector (taxa level)
    taxa.sum <- NULL
    for (t in 1:num.taxa.levels){
      taxa.sum[t] <- sum(reads.list[[p]][[t]])
    }
    reads.summary[,p] <- taxa.sum
  }
  return(reads.summary)
}

# add to the read.summaries table the overall totals to match the otu.summaries structure
# this involves importing the fw_classified_taxonomies.csv file to get the list of fw seqIDs
add.totals.to.read.summaries <- function(ReadSummaryTable, AbundanceTable, PidentsUsed, CustomSeqIDs){
  seqID.reads <- AbundanceTable
  reads.summary <- ReadSummaryTable
  pident.values <- PidentsUsed
  fw.ids <- CustomSeqIDs
  
  fw.reads <- NULL
  # for each pident tried
  for (p in 1:length(pident.values)){
    
    # for each seqID classified with FW
    fw.reads[p] <- 0
    if (length(fw.ids[[p]]) > 0){
      for (s in 1:length(fw.ids[[p]])){
        index <- which(seqID.reads[,1] == fw.ids[[p]][s])
        fw.reads[p] <- fw.reads[p] + seqID.reads[index,2]
      } 
    }
    names(fw.reads)[p] <- paste("pident", pident.values[p], sep = "")
  }
  
  # Add the total reads classified by FW at each pident to the summary table
  reads.summary <- rbind(reads.summary, fw.reads)
  
  # Add the total reads to the summary table
  tot.reads <- sum(seqID.reads[ ,2])
  reads.summary <- rbind(reads.summary, tot.reads)
  
  return(reads.summary)
}

# import bootstrap p-values to compare
import.bootstrap.pvalues <- function(ConflictFolders, PidentsUsed, FW = TRUE){
  pident.folders <- ConflictFolders
  pident.values <- PidentsUsed
  
  if (FW == TRUE){
    db <- "fw"
  }else{
    db <- "gg"
  }
  
  # first import all the files into one list
  pvalues.list <- NULL
  for (p in 1:length(pident.folders)){
    pvalues.list[[p]] <- read.csv(file = paste(pident.folders[p], "/", db, "_classified_bootstraps.csv", sep = ""))
    names(pvalues.list)[p] <- pident.values[p]
  }
  
  # second reformat that list into a new list amenable to box plots
  pident.values <- names(pvalues.list)
  taxa.levels <- c("Kingdom", "Phylum", "Class", "Order", "Lineage", "Clade", "Tribe")
  bootstraps.taxa <- list(NULL)
  for (t in 1:7){
    bootstraps.pidents <- list(NULL)
    for (p in 1:length(pvalues.list)){
      bootstraps.pidents[[p]] <- pvalues.list[[p]][,t]
      names(bootstraps.pidents)[p] <- pident.values[p]
    }
    bootstraps.taxa[[t]] <- bootstraps.pidents
    names(bootstraps.taxa)[t] <- taxa.levels[t]
  }
  
  return(bootstraps.taxa)
}

# import the custom-only classified taxonomy for the "forcing" plots (function borrowed from find_classification_disagreements.R)
import.custom.only.taxonomy <- function(FilePath){
  forced.taxonomy.file.path <- FilePath
  
  # simple imort b/c file format determined by previous script, shouldn't have any weird things
  fw <- read.table(forced.taxonomy.file.path , sep=";", fill=T, colClasses = "character", header = TRUE)
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(fw[,1])
  fw <- fw[index,]
  
  # Convert into a character matrix for faster processing
  fw <- as.matrix(fw)
  
  # Remove row names that will not match between the data tables
  row.names(fw) <- NULL
  
  return(fw)
}

# Group seqIDs at each taxa level- this is used in plot.most.misleading.forced.taxa
group.seqIDs.into.taxa <- function(TaxonomyTable, ReadsPerSeqID){
  forced.taxonomy <- TaxonomyTable
  seqID.reads <- ReadsPerSeqID

  # first put them im the same order (they should be both character seqIDs already, and taxonomy was ordered inside it's import function)
  index <- order(seqID.reads[ ,1])
  seqID.reads.ord <- seqID.reads[index, ]
  if (all.equal(seqID.reads.ord[ ,1], forced.taxonomy[ ,1]) != TRUE){
    cat("Crap something's messed up with the indexing, seqID.reads and forced.taxonomy need to have the same order of seqIDs")
  }
  
  # add total reads to taxonomy table
  tax.reads <- cbind(forced.taxonomy, seqID.reads.ord[ ,2])
  colnames(tax.reads)[9] <- "reads"
  
  # leave out tribe- forced.seqIDs only goes down to clade
  otus.taxa <- tax.reads  # name change is a relic, also makes easier to troubleshoot
  otus.taxa <- otus.taxa[ ,-8]
  
  # Order taxonomy names so that phyla are alphabetical, classes within a phylum are alphabetical, etc
  index <- order(otus.taxa[ ,2], otus.taxa[ ,3], otus.taxa[ ,4], otus.taxa[ ,5], otus.taxa[ ,6], otus.taxa[ ,7])
  otus.taxa.ord <- otus.taxa[index, ]
  
  # Make unique taxonomy table (so that there aren't duplicate names at each level as occurs with "unclassified")
  for (t in 2:7){
    unique.taxa <- unique(otus.taxa.ord[ ,2:t, drop = FALSE])          
    unique.names <- make.unique(unique.taxa[,(t-1)])       
    for (u in 1:nrow(unique.taxa)){                     
      for (r in 1:nrow(otus.taxa.ord)){
        if (all(otus.taxa.ord[r,2:t] == unique.taxa[u,1:(t-1)])){
          otus.taxa.ord[r,t] <- unique.names[u]
        }
      }
    }
  }
  
  # change to being a dataframe so that read numbers can be numeric
  otus.taxa.ord <- as.data.frame(otus.taxa.ord, stringsAsFactors = F)
  otus.taxa.ord[ ,8] <- as.numeric(otus.taxa.ord[ ,8])
  
  # Create a blank list to fill with data
  grouped.taxa <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL)
  
  # Populate list with the unique taxa names at each level and a numeric column for total reads
  for (t in 1:6){
    grouped.taxa[[t]] <- as.data.frame( cbind( unique(otus.taxa.ord[ ,(2:(t+1))]), 
                                               rep( x = 0, times = nrow(unique(otus.taxa.ord[ ,(2:(t+1)), drop=F])) ) ), 
                                        stringsAsFactors = FALSE)
    grouped.taxa[[t]][ ,(t+1)] <- as.numeric(grouped.taxa[[t]][ ,(t+1)])
    colnames(grouped.taxa[[t]]) <- c( colnames(otus.taxa.ord)[(2:(t+1))], "reads")
    rownames(grouped.taxa[[t]])<- NULL
  }
  
  # Fill in the reads columns in each element of the grouped.taxa list
  for (t in 1:6){                                   
    uniquename <- 1
    for (r in 1:nrow(otus.taxa.ord)){
      if ( all( otus.taxa.ord[r, 2:(t + 1)] == grouped.taxa[[t]][uniquename, 1:t]) ){
        grouped.taxa[[t]][uniquename, (t + 1)] <- grouped.taxa[[t]][uniquename, (t + 1)] + otus.taxa.ord[r,8]
      }else{
        uniquename <- uniquename + 1
        grouped.taxa[[t]][uniquename, (t + 1)] <- grouped.taxa[[t]][uniquename, (t + 1)] + otus.taxa.ord[r,8]
      }
    }
    cat("Now you're on taxa level:",t,'\n')
  }
  return(grouped.taxa)
}

# Find the most abundant taxa at each taxonomy level
find.top.taxa.by.total.reads <- function(TaxonomyList, NumberTopTaxa){
  grouped.taxa <- TaxonomyList
  num.taxa <- NumberTopTaxa
  
  # arrange highest to lowest total reads
  grouped.taxa.ord <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL)
  for (t in 1:6){
    index <- order(grouped.taxa[[t]][ ,(t + 1)], decreasing = TRUE)
    grouped.taxa.ord[[t]] <- grouped.taxa[[t]][index, ]
  }
  
  # remove unclassified taxa b/c those really can't be compared on the same taxa level, and likely wouldn't be included in a "top taxa" analysis anyway.
  not.unclassifieds <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL)
  for (t in 1:6){
    index <- grep(x = grouped.taxa.ord[[t]][ ,t], pattern =  "unclassified.*", value = FALSE )
    # this is necessary because you can not use -0 as an index
    if (length(index) != 0){
      not.unclassifieds[[t]] <- grouped.taxa.ord[[t]][-index, ]
    }else{
      not.unclassifieds[[t]] <- grouped.taxa.ord[[t]]
    }
  }
  
  # look just at the top 20 levels
  grouped.taxa.top <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL)
  for (t in 1:6){
    if (nrow(not.unclassifieds[[t]]) < num.taxa){
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]]
    }else{
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]][1:num.taxa, ]
    }
  }
  return(grouped.taxa.top)
}


#####
# Define functions to plot the data
#####

plot.num.forced <- function(ConflictSummaryTable, ResultsFolder, DBconflicts = as.data.frame(FALSE), ByReads = FALSE, AsPercent = FALSE, y.axis.limit = 0){
  sum.table <- ConflictSummaryTable
  results.file.path <- ResultsFolder
  db.conflicts <- DBconflicts
  
  # remove the last 2 rows of number FW sequences- totals info not needed for this plot.
  mismatches <- sum.table[1:(nrow(sum.table)-2),]
  pidents <- colnames(mismatches)
  pidents <- as.numeric(pidents)
  total.seqs.or.reads <- sum.table[nrow(sum.table),1]
  
  # modify plot based on type specified in function calls
  if (ByReads == FALSE){
    plot.of <- "OTU"
  }else{
    plot.of <- "read"
  }
  
  if (AsPercent == FALSE){
    plot.as <- "Total"
  }else{
    plot.as <- "Percent"
    mismatches <- mismatches / total.seqs.or.reads * 100
  }
  
  if (y.axis.limit == 0){
    ymax <- max(mismatches)
    yplotlabel <- ""
  }else{
    ymax <- y.axis.limit
    yplotlabel <- paste("_y-axis_cutoff_",ymax, sep = "")
  }
  
  if (db.conflicts[1,1] == FALSE){
    db.label <- ""
    db.plot.label <- ""
  }else{
    db.label <- "-DBs_included"
    db.plot.label <- "\nHorizonal lines show the conflicts between the databases"
  }
  
  # Save plot as .png file
  png(filename = paste(results.file.path, "/Classification_Disagreements_of_", plot.of, "_", plot.as, "s", yplotlabel, db.label, ".png", sep = ""), 
      width = 5, height = 5, units = "in", res = 100)
  
  # Set up an empty plot
  plot(x = 0, type = "n", ylim = c(0,ymax), xlim = c(min(pidents),max(pidents)),
       main = paste("Disagreements Between Custom and General Taxonomy Classifications\n(Are we forcing OTUs into our favorite groups?)", db.plot.label, sep = ""),
       cex.main = .8,
       ylab = paste("Classification Disagreements (", plot.as, " ", plot.of, "s)", sep = ""), cex.lab = .8, 
       xlab = "\"full length\" pident cutoff (similar to ANI of read to custom database)")
  
  # Fill Plot with beautiful data
  color <- c("seagreen4","violetred4","slateblue4","turquoise4")
  for (r in 1:nrow(mismatches)){
    lines(x = pidents, y = mismatches[r,], col = color[r], lwd = 2)
    points(x = pidents, y = mismatches[r,], col = color[r], pch = 19, cex =1.3)
  }
  
  # Add database conflicts for baseline, if desired:
  if (db.conflicts[1,1] != FALSE){
    for (r in 1:nrow(db.conflicts)){
      abline(h = db.conflicts[r,2], col = color[r], lty =2, lwd =3)
    }
  }
  
  legend("topright", legend = row.names(mismatches), text.col = color, cex=1)
  
  unnecessary.message <- dev.off()
}

plot.num.classified.outs <- function(ConflictSummaryTable, ResultsFolder, ByReads = FALSE, AsPercent = TRUE){
  sum.table <- ConflictSummaryTable
  results.file.path <- ResultsFolder
  
  # pull out data needed for this plot
  num.fw <- sum.table[nrow(sum.table)-1,]
  tot.seq <- sum.table[nrow(sum.table),1]
  pidents <- colnames(sum.table)
  pidents <- as.numeric(pidents)
  
  if (ByReads == FALSE){
    plot.of <- "OTU"
  }else{
    plot.of <- "read"
  }
  
  if (AsPercent == FALSE){
    plot.as <- "Total"
  }else{
    plot.as <- "Percent"
    num.fw <- num.fw / tot.seq * 100
  }
  
  # Save plot as .png file
  png(filename = paste(results.file.path, "/Custom_db_contributions_by_", plot.of, "_", plot.as, "s.png", sep = ""), 
      width = 5, height = 5, units = "in", res = 100)
  
  # Set up and empty plot
  plot(x = pidents, y = num.fw, type = "n", 
       main = "Contribution of Custom Taxonomy \n to Total Classifications", cex.main = 1,
       xlab = "\"Full Length\" percent identity (similar to ANI)", 
       ylab = paste("OTUs Classified by Custom Taxonomy (", plot.as, " ", plot.of, "s)", sep = ""))
  
  # Fill plot with beautiful data
  lines(x = pidents, y = num.fw, col = "lightsalmon", lwd = 1.5)
  points(x = pidents, y = num.fw, col = "lightsalmon", pch = 19, cex = 1.3)
  
  unnecessary.message <- dev.off()
}
  
plot.bootstrap.percents <- function(FWpValues, GGpValues, ResultsFolder){
  fw.pvalues <- FWpValues
  gg.pvalues <- GGpValues
  results.file.path <- ResultsFolder
  
  # set up a new file for the plot
  png(filename = paste(results.file.path, "/Taxonomy_Assignment_Confidences.png", sep = ""), 
      width = 8, height = 10, units = "in", res = 100)
  
  # set up stacked boxplots, left column FW, right column GG
  par(mfcol = c(7,2), omi = c(.6,.9,.6,.1), mai = c(.1,.01,.1,.01))
  
  # plot the FW side
  for (t in 1:7){
    boxplot(fw.pvalues[[t]], range = 0, whisklty = "solid", ylim = c(0,100),axes = F)
    axis(side = 2, at = c(0,100), cex.axis = 1)
    mtext(text = names(fw.pvalues)[t], side = 2, line = 2.2, cex = 1.5)
    # rect(xleft = .5, ybottom = -10, xright = length(fw.pvalues[[1]]) + .5, ytop = 60, col = adjustcolor(col = "grey", alpha.f = .2), border = NA)
  }
  axis(side = 1, at = 1:length(fw.pvalues[[1]]), labels = rep(x = "", times = length(fw.pvalues[[1]])),
       cex.axis = .7, outer = T)
  mtext(side = 1, at = 1:length(fw.pvalues[[1]]), text = names(fw.pvalues[[1]]), cex = .7, line = 1.5)
  
  # plot the GG side
  for (t in 1:7){
    boxplot(gg.pvalues[[t]], range = 0, whisklty = "solid", ylim = c(0,100),axes = F)
    # rect(xleft = .5, ybottom = -10, xright = length(fw.pvalues[[1]]) + .5, ytop = 60, col = adjustcolor(col = "grey", alpha.f = .2), border = NA)
  }
  axis(side = 1, at = 1:length(gg.pvalues[[1]]), labels = rep(x = "", times = length(gg.pvalues[[1]])),
       cex.axis = .7, outer = T)
  mtext(side = 1, at = 1:length(gg.pvalues[[1]]), text = names(gg.pvalues[[1]]), cex = .7, line = 1.5)
  
  # add titles
  mtext("Effect of Cutoff on Assignment Repeatability", side = 3, outer = T, line = 2.5, cex = 1.8)
  mtext(side = 3, at = c(.2,.8), text = c("Custom Classified", "General Classified"), outer = T, line = 0, cex = 1.5)
  mtext("BLAST Full Length pident cutoff (%)", side = 1, outer = T, line = 2.8, cex = 1.2)
  mtext("RDP Classifier Repeatability (that number after your taxonomy name)", side = 2, outer = T, line = 4.9, cex = 1.2)
  
  # finish plotting and create the file
  unnecessary.message <- dev.off()
}

plot.percent.forced <- function(ForcingTable, ResultsFolder, ByReads = FALSE){
  num.forced <- ForcingTable
  plots.folder.path <- ResultsFolder
  
  if (ByReads == FALSE){
    plot.type <- "OTUs"
  }else if (ByReads == TRUE){
    plot.type <- "Reads"
  }
  
  tot <- num.forced[7]
  perc.forced <- num.forced[-7,1,drop = FALSE] / tot * 100
  
  png(filename = paste(plots.folder.path, "/Forcing_by_", plot.type, ".png", sep = ""), 
      width = 7, height = 5, units = "in", res = 100)
  
  par(mar = c(5,6,4,1))
  
  barplot(height = perc.forced, beside = TRUE, names.arg = row.names(perc.forced), ylim = c(0,max(perc.forced)+1),
          space = .2, 
          main = paste("Why the custom database alone should not be used\n(By Percent Total ", plot.type, ")", sep = ""),
          ylab = paste("Percent of total", plot.type, "forced \ninto an incorrect classification"),
          xlab = "Incorrect Classifications at each Taxonomic Level", col = "olivedrab4")
  unnecessary.message <- dev.off()
  
  
}

plot.most.misleading.forced.otus <- function(ReadsPerForcedSeqIDs, ForcedSeqIDs, ReadsPerSeqID, OutputFolder, PlottingLevels = 1:6, NumOTUs = 50){
  forced.seqIDs <- ForcedSeqIDs
  forced.seqID.reads <- ReadsPerForcedSeqIDs
  seqID.reads <- ReadsPerSeqID
  plots.folder.path <- OutputFolder
  
  # reformat the multilevel lists to 1 level, since there's only 1 upper level anyway. 
  forced.reads <- forced.seqID.reads[[1]]
  forced.seqIDs <- forced.seqIDs[[1]]
  
  # find the max OTUs by total abundance
  indeces <- order(seqID.reads[ ,2], decreasing = TRUE)
  max.indeces <- indeces[1:NumOTUs]
  max.seqIDs.reads <- seqID.reads[max.indeces, ]
  max.seqID.reads.perc <- max.seqIDs.reads 
  max.seqID.reads.perc[ ,2] <- max.seqID.reads.perc[ ,2] / sum(seqID.reads[ ,2]) * 100
  
  # Now which ones of those were forced?
  seqID.recorder <- list(kingdom = NULL, phylum = NULL, class = NULL, order = NULL, lineage = NULL, clade = NULL)
  for (t in 1:5){
    for (s in 1:nrow(max.seqID.reads.perc)){
      check.match <- max.seqID.reads.perc[s,1] == forced.seqIDs[[t]]
      index <- which(check.match == 1)
      if (length(index) > 0){
        seqID.recorder[[t]] <- c(seqID.recorder[[t]], s)
      }
    }
  }
  
  # This doesn't show up well because it's by OTU, not by taxonomic assignment which can have several OTUs.
  # plot it and export the plots
  for (t in PlottingLevels){
    color.vector <- rep(x = "grey", times = length(max.indeces))
    color.vector[seqID.recorder[[t]]] <- "red"
    
    png(filename = paste(plots.folder.path, "/Forcing_of_top_OTUs-", names(seqID.recorder)[t], "_level.png", sep = ""), 
        width = 7, height = 5, units = "in", res = 100)
    
    barplot(max.seqID.reads.perc[ ,2], col = color.vector, 
            main = "Do any top OTUs (by total abundance overall) end up\nwith erroneous classifications due to forcing?",
            xlab = paste("Top",length(max.indeces),"OTUs\nRed Bars indicate \"forcing\" at the", names(seqID.recorder)[t], "level"))
    
    unnecessary.message <- dev.off()
  }
}

plot.most.misleading.forced.taxa <- function(TopTaxaList, ForcedTaxonomy, ForcedReadsList, ForcedSeqIDsList, ResultsFolder, PlottingLevels = 1:6, TotalReads){
  top.taxa <- TopTaxaList
  forced.taxonomy <- ForcedTaxonomy
  forced.seqID.reads <- ForcedReadsList
  forced.seqIDs <- ForcedSeqIDsList
  plots.folder.path <- ResultsFolder
  tot.reads <- TotalReads
  
  # reformatting- don't repeat this step twice
  forced.seqID.reads <- forced.seqID.reads[[1]]
  forced.seqIDs <- forced.seqIDs[[1]]
  for (t in 1:6){
    top.taxa[[t]] <- cbind(top.taxa[[t]], forced = 0, correct = 0)
    names(top.taxa[[t]])[(t+1)] <- "total"
  }
  
  # find the contributing forced seqIDs and tally their reads
  taxa.reads.forced <- 0
  for (t in 1:6){ 
    
    for (r in 1:nrow(top.taxa[[t]])){
      index.seqIDs <- which(forced.taxonomy[ ,(t + 1)] == top.taxa[[t]][r,t])
      taxas.seqIDs <- forced.taxonomy[index.seqIDs,1]
      
      for (s in 1:length(taxas.seqIDs)){
        is.forced <- which(forced.seqIDs[[t]] == taxas.seqIDs[s])
        if (length(is.forced) != 0){
          taxa.reads.forced <- taxa.reads.forced + forced.seqID.reads[[t]][is.forced]
        }
      }
      
      top.taxa[[t]][r,(t + 2)] <- taxa.reads.forced
      taxa.reads.forced <- 0
    }
  }
 
  # add in columns of "correct" reads
  for (t in 1:6){
    top.taxa[[t]][ ,(t + 3)] <- top.taxa[[t]][ ,(t + 1)] - top.taxa[[t]][ ,(t + 2)]
  }
  
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
}


#####
# Use Functions
#####

#####
# first check if this is the optional "forcing plot"
#####
if (forcing.folder.path != "regular"){
  otus.forced <- import.forcing.conflicts(ForcingFolder = forcing.folder.path)
  
  seqID.reads <- import.seqID.reads(FilePath = seqID.reads.file.path) # this was exported previously by this script
  
  forced.seqIDs <- get.conflict.seqIDs(ConflictsFolders = forcing.folder.path, PidentsUsed = "forcing")
  
  forced.seqID.reads <- find.reads.per.seqID(ReadsTable = seqID.reads, ConflictsList = forced.seqIDs, Forcing = TRUE)
  
  read.summaries <- generate.summary.table.of.reads(ReadsList = forced.seqID.reads, Forcing = TRUE)
  
  tot.reads <- sum(seqID.reads$reads)
  
  reads.forced <- rbind(read.summaries, tot.reads)
  
  forced.taxonomy <- import.custom.only.taxonomy(FilePath = forced.taxonomy.file)
  
  grouped.taxa <- group.seqIDs.into.taxa(TaxonomyTable = forced.taxonomy, ReadsPerSeqID = seqID.reads)
  
  top.taxa <- find.top.taxa.by.total.reads(TaxonomyList = grouped.taxa, NumberTopTaxa = 20)
  
  plot.percent.forced(ForcingTable = otus.forced, ResultsFolder = plots.folder.path, ByReads = FALSE)
  
  plot.percent.forced(ForcingTable = reads.forced, ResultsFolder = plots.folder.path, ByReads = TRUE)
  
  plot.most.misleading.forced.otus(ReadsPerForcedSeqIDs = forced.seqID.reads, ForcedSeqIDs = forced.seqIDs, 
                                   ReadsPerSeqID = seqID.reads, OutputFolder = plots.folder.path, PlottingLevels = 1:6)
  
  plot.most.misleading.forced.taxa(TopTaxaList = top.taxa, ForcedTaxonomy = forced.taxonomy, 
                                   ForcedReadsList = forced.seqID.reads, ForcedSeqIDsList = forced.seqIDs, 
                                   ResultsFolder = plots.folder.path, PlottingLevels = 4:6, TotalReads = tot.reads)
#####
# If not then do the normal comparison for choosing pident cutoff
#####
}else{
  
  #####
  # examine custom taxonomy disagreements and contribution by number OTUs
  #####
  
  otu.summaries <- import.all.conflict.summaries(ConflictFolders = pident.folders, PidentsUsed = pident.values)
  
  # I removed this from the terminal command also because it is misleading and confusing and doesn't add to the plot. left a commented out input line for use within RStudio at the top.
  # db.summary <- import.database.conflicts(DatabaseFolder = db.conflicts.folder.path)
  
  # ignore lineage- it's way higher than the others doesn't fit on graphs
  otu.summaries <- otu.summaries[-5, ]
  # db.summary <- db.summary[-5, ]
  
  plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path)
  plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, y.axis.limit = 10)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE, y.axis.limit = 1)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, DBconflicts = db.summary, y.axis.limit = max(db.summary[,2]))
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, DBconflicts = db.summary)
  
  # plot.num.classified.outs(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = FALSE)
  plot.num.classified.outs(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE)
  
  #####
  # examine custom taxonomy disagreements and contribution by number reads
  #####
  
  seqID.reads <- import.and.reformat.otu.table(OTUtable = otu.table.path)
  
  conflict.seqIDs <- get.conflict.seqIDs(ConflictsFolders = pident.folders, PidentsUsed = pident.values)
  
  conflict.seqID.reads <- find.reads.per.seqID(ReadsTable = seqID.reads, ConflictsList = conflict.seqIDs)
  
  read.summaries <- generate.summary.table.of.reads(ReadsList = conflict.seqID.reads)
  
  custom.seqIDs <- import.ids.above(FilePaths = ids.file.paths, PidentsUsed = pident.values)
  
  read.summaries <- add.totals.to.read.summaries(ReadSummaryTable = read.summaries, AbundanceTable = seqID.reads, 
                                                 PidentsUsed = pident.values, CustomSeqIDs = custom.seqIDs)
  #leave out lineage on the plots b/c it's too much higher
  read.summaries <- read.summaries[-5, ]
  
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE)
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, y.axis.limit = 10000)
  plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE)
  plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE, y.axis.limit = .5)
  
  # plot.num.classified.outs(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = FALSE)
  plot.num.classified.outs(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE)
  
  # export the reads per seqID for use in the plot_classification_improvement.R script
  write.table(x = seqID.reads, file = "total.reads.per.seqID.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  
  #####
  # examine pident cutoff relationship to bootstrap p-values -this plot is not useful.
  #####
  # fw.pvalues <- import.bootstrap.pvalues(ConflictFolders = pident.folders, PidentsUsed = pident.values, FW = TRUE)
  # gg.pvalues <- import.bootstrap.pvalues(ConflictFolders = pident.folders, PidentsUsed = pident.values, FW = FALSE)
  
  # plot.bootstrap.percents(FWpValues = fw.pvalues, GGpValues = gg.pvalues, ResultsFolder = plots.folder.path)
}



