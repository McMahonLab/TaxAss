# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# compare pidents based on # OTUs and # reads impacted.

#####
# Receive arguments from terminal command line
#####

# userprefs <- commandArgs(trailingOnly = TRUE)

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
# Define functions to import and process the data
#####

# import all the conflict summary files from each folder and compile them into a matrix
import.all.conflict.summaries <- function(UserArgs){
  userprefs <- UserArgs
  
  # the first command line argument is the otu rel abund table, so ignore that here
  user.args <- userprefs[-1]
  
  # Import them into a list format
  mismatches.list <- list(NULL)
  counter <- 1
  for (p in seq(from = 1, to = length(user.args), by = 2)){
    mismatches.list[[counter]] <- read.csv(file = paste(user.args[p], "/conflicts_summary.csv", sep = ""))
    names(mismatches.list)[counter] <- user.args[p+1]
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

# import the OTU table and then pull out just the total reads for each seqID
import.and.reformat.otu.table <- function(UserArgs){
  userprefs <- UserArgs
  
  otu.table.file.path <- userprefs[1]
  
  otus <- read.table(file = otu.table.file.path, header = TRUE, stringsAsFactors = FALSE)
  seqID.reads <- data.frame(seqID = otus[ ,1], reads = rowSums(otus[ ,-1]))
  
  return(seqID.reads)
}

# import files comparing conflicts at each level, record the seqIDs into a list of lists
# structure: outer list- each pident, inner lists- seqIDs at each taxa level
get.conflict.seqIDs <- function(UserArgs){
  userprefs <- UserArgs
  user.args <- userprefs[-1]
  pident.folders <- user.args[seq(from = 1, to = length(user.args), by = 2)]
  pident.values <- user.args[seq(from = 1, to = length(user.args), by = 2)+1]
  
  all.pidents <- list(NULL)
  # for each pident folder
  for (p in 1:length(pident.values)){
    all.files <- list.files(pident.folders[p])
    conflict.ids <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL)
    
    # for each taxonomy file
    for (t in 1:5){
      conflict.ids[[t]] <- read.csv(file = paste(pident.folders[p], "/", all.files[t], sep = ""), header = TRUE, stringsAsFactors = FALSE)
      # first column is the seqID vector
      conflict.ids[[t]] <- conflict.ids[[t]][,1]
    }
    
    all.pidents[[p]] <- conflict.ids
    names(all.pidents)[p] <- pident.values[p]
  }
  return(all.pidents)
}

# create a list of total reads that matches the structure of the list of conflicting seqIDs
find.reads.per.seqID <- function(ReadsTable, ConflictsList){
  otu.reads <- ReadsTable
  conflict.ids <- ConflictsList
  
  # for each upper list level (pident folder)
  pidents.list <- list(NULL)
  for (p in 1:length(conflict.ids)){
  
    # for each inner list level (taxonomy level)
    taxa.list <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL)
    for (t in 1:5){
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
generate.summary.table.of.reads <- function(ReadsList){
  reads.list <- ReadsList
  
  # Set up empty matrix to fill
  reads.summary <- matrix(0, nrow = 5, ncol = length(reads.list))
  row.names(reads.summary) <- names(reads.list[[1]])
  colnames(reads.summary) <- names(reads.list)
  
  #For each outer list's list (pident)
  for (p in 1:length(reads.list)){
    
    # For each inner list's vector (taxa level)
    taxa.sum <- NULL
    for (t in 1:5){
      taxa.sum[t] <- sum(reads.list[[p]][[t]])
    }
    reads.summary[,p] <- taxa.sum
  }
  return(reads.summary)
}

# add to the read.summaries table the overall totals to match the otu.summaries structure
# this involves importing the fw_classified_taxonomies.csv file to get the list of fw seqIDs
add.totals.to.read.summaries <- function(ReadSummaryTable, AbundanceTable, UserArgs){
  seqID.reads <- AbundanceTable
  reads.summary <- ReadSummaryTable
  userprefs <- UserArgs
  user.args <- userprefs[-1]
  pident.folders <- user.args[seq(from = 1, to = length(user.args), by = 2)]
  pident.values <- user.args[seq(from = 1, to = length(user.args), by = 2)+1]
  
  fw.reads <- NULL
  # for each pident tried
  for (p in 1:length(pident.folders)){
    all.files <- list.files(pident.folders[p])
    fw.classified <- read.csv(file = paste(pident.folders[p], "/", all.files[8], sep=""))
    fw.seqIDs <- fw.classified[ ,1]
    
    # for each seqID classified with FW
    fw.reads[p] <- 0
    if (length(fw.seqIDs) > 0){
      for (s in 1:length(fw.seqIDs)){
        index <- which(seqID.reads[,1] == fw.seqIDs[s])
        fw.reads[p] <- fw.reads[p] + seqID.reads[index,2]
      } 
    } # no else b/c it's already zero at start 
    names(fw.reads)[p] <- paste("pident", pident.values[p], sep = "")
  }
  
  # Add the total reads classified by FW at each pident to the summary table
  reads.summary <- rbind(reads.summary, fw.reads)
  
  # Add the total reads to the summary table
  tot.reads <- sum(seqID.reads[ ,2])
  reads.summary <- rbind(reads.summary, tot.reads)
  
  return(reads.summary)
}

#####
# Define functions to plot the data
#####

plot.num.forced.otus <- function(ConflictSummaryTable, ByReads = FALSE, y.axis.limit = 0){
  sum.table <- ConflictSummaryTable
  
  # remove the last row of number FW sequences, that's not needed for this plot.
  mismatches <- sum.table[1:(nrow(sum.table)-2),]
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

otu.summaries <- import.all.conflict.summaries(UserArgs = userprefs)

plot.num.forced.otus(ConflictSummaryTable = otu.summaries)
plot.num.forced.otus(ConflictSummaryTable = otu.summaries, y.axis.limit = 10)

plot.num.classified.outs(ConflictsSummaryTables = otu.summaries)

seqID.reads <- import.and.reformat.otu.table(UserArgs = userprefs)

conflict.seqIDs <- get.conflict.seqIDs(UserArgs = userprefs)

conflict.seqID.reads <- find.reads.per.seqID(ReadsTable = seqID.reads, ConflictsList = conflict.seqIDs)

read.summaries <- generate.summary.table.of.reads(ReadsList = conflict.seqID.reads)

read.summaries <- add.totals.to.read.summaries(ReadSummaryTable = read.summaries, AbundanceTable = seqID.reads, UserArgs = userprefs)
