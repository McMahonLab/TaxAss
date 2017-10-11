# Script written by Robin Rohwer

# Make a plot comparing the classification conflicts between pidents
# Make a plot showing the proportion of sequences classified in FW database
# compare pidents based on # OTUs and # reads impacted.

# This script is step 14 of the taxonomy assignment workflow, used to decide on pident cutoff. It takes arguments from the command line.
# The number of arguments can be variable, but they must be in order
# The variable part is that more pidents can be added, as long as they continue the pattern folder number folder number
# This script is also used in step 16 optionally to plot the forcing that would have occured if you only used FW

# Terminal command line syntax:
# Rscript plot_classification_disagreements.R otus.abund plots regular NA NA conflicts_94 ids.above.94 94 conflicts_96 ids.above.96 96 conflicts_98 ids.above.98 98 ...
# Rscript plot_classification_disagreements.R NA plots conflicts_forcing otus.custom.85.taxonomy otus.98.85.70.taxonomy
# Rscript plot_classification_disagreements.R otus.abund plots conflicts_forcing otus.custom.85.taxonomy otus.98.85.70.taxonomy
# note: the forcing option with otus.abund specified is for if you skipped step 14 so you need to make the seqid.reads file

# ---------------------------------------------------------------------------------------------------------------------
# Receive arguments from terminal command line
# ---------------------------------------------------------------------------------------------------------------------

userprefs <- commandArgs(trailingOnly = TRUE)

# FOR PLOTTING FORCING  **don't forget to change the seqID.reads file path below!!
# cat("fuck you forgot to comment out the file paths in plot_classification_disagreements!\n")
# userprefs <- c(NA, # if seqid.reads exists (i.e. you ran step 14) this is NA, otherwise it's otus.abund file path)
#                "../../ME_plot_test/plots",
#                "../../ME_plot_test/conflicts_forcing",
#                "../../ME_plot_test/otus.custom.85.taxonomy",
#                "../../ME_plot_test/otus.98.85.70.taxonomy")

# # FOR CHOOSING CUTOFF:
# cat("fuck you forgot to comment out the file paths in plot_classification_disagreements!")
# userprefs <- c("../../poster_mend-check/otus.abund",
#                "../../poster_mend-check/plots",
#                "regular",
#                NA,
#                NA,
#                "../../poster_mend-check/conflicts_90",
#                "../../poster_mend-check/ids.above.90",
#                90,
#                "../../poster_mend-check/conflicts_91",
#                "../../poster_mend-check/ids.above.91",
#                91,
#                "../../poster_mend-check/conflicts_92",
#                "../../poster_mend-check/ids.above.92",
#                92,
#                "../../poster_mend-check/conflicts_93",
#                "../../poster_mend-check/ids.above.93",
#                93,
#                "../../poster_mend-check/conflicts_94",
#                "../../poster_mend-check/ids.above.94",
#                94,
#                "../../poster_mend-check/conflicts_95",
#                "../../poster_mend-check/ids.above.95",
#                95,
#                "../../poster_mend-check/conflicts_96",
#                "../../poster_mend-check/ids.above.96",
#                96,
#                "../../poster_mend-check/conflicts_97",
#                "../../poster_mend-check/ids.above.97",
#                97,
#                "../../poster_mend-check/conflicts_98",
#                "../../poster_mend-check/ids.above.98",
#                98,
#                "../../poster_mend-check/conflicts_99",
#                "../../poster_mend-check/ids.above.99",
#                99,
#                "../../poster_mend-check/conflicts_100",
#                "../../poster_mend-check/ids.above.100",
#                100)
# # JUST MAKE SEQID.READS FILE, SKIPPING STEP 14 BUT DOING 15.5.A
# cat("fuck you forgot to comment out the file paths in plot_classification_disagreements!")
# userprefs <- c("../../take18playwith/otus.abund", 
#                "MakeSeqIDReadsOnly")

# in case you want to add the db baseline conflict back to the plots, need to specify this path below
# and un-comment the plotting calls that use it at the end of the script.

# db.conflicts.folder.path <- "file path to conflicts_database"
otu.table.path <- userprefs[1]
plots.folder.path <- userprefs[2]
forcing.folder.path <- userprefs[3]
forced.taxonomy.file <- userprefs[4]
final.taxonomy.file <- userprefs[5]
rest.of.arguments <- userprefs[-(1:5)]
if (length(rest.of.arguments) > 0){
  pident.folders <- rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)]
  ids.file.paths <- rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)+1]
  pident.values <- as.numeric(rest.of.arguments[seq(from = 1, to = length(rest.of.arguments), by = 3)+2])
}

# this is automatically exported into the working directory when this script is run normally
seqID.reads.file.path <- "total.reads.per.seqID.csv"
present.working.directory <- "."

# cat("fuck you forgot to comment out the seqid.reads file path in plot_classification_disagreements!\n")
# seqID.reads.file.path <- "../../ME_plot_test/total.reads.per.seqID.csv"
# present.working.directory <- "../../ME_plot_test/"


# ---------------------------------------------------------------------------------------------------------------------
# Define functions to import and process the data
# ---------------------------------------------------------------------------------------------------------------------

make.plot.directory <- function(FolderPath){
  if (!dir.exists(FolderPath)){
    dir.create(FolderPath)
    cat("made folder: ", FolderPath, "\n")
  }
}

import.all.conflict.summaries <- function(ConflictFolders, PidentsUsed){
  # import all the conflict summary files from each folder and compile them into a matrix
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

import.database.conflicts <- function(DatabaseFolder){
  # import the database conflicts summary file
  db.conflicts.folder.path <- DatabaseFolder
  db.conflicts <- read.csv(file = paste(db.conflicts.folder.path, "/conflicts_summary.csv", sep =""), stringsAsFactors = FALSE)
  # the sequences totals are not needed here
  return(db.conflicts[1:5, ])
}

import.forcing.conflicts <- function(ForcingFolder){
  # import the optional plotting step's forcing conflicts summary file
  forcing.folder.path <- ForcingFolder
  forcing.conflicts <- read.csv(file = paste(forcing.folder.path, "/conflicts_summary.csv", sep = ""), stringsAsFactors = FALSE)
  forcing.conflicts.matrix <- as.matrix(forcing.conflicts[ ,2, drop = FALSE])
  row.names(forcing.conflicts.matrix) <- forcing.conflicts[ ,1]
  # remove total that should have been GG-classified (confusingly labelled "numFWseqs")
  forcing.conflicts.matrix <- forcing.conflicts.matrix[-8,1, drop = FALSE]
  return(forcing.conflicts.matrix)
}

import.seqID.reads <- function(FilePath){
  # import seqID.reads variable for the forcing plot (it's generated w/ the regular plotting)
  seqID.reads.file.path <- FilePath
  
  seqID.reads <- read.csv(file = seqID.reads.file.path, colClasses = "character")
  seqID.reads[ ,2] <- as.numeric(seqID.reads[ ,2])
  
  #alphabetical order matches taxonomy tables
  index <- order(seqID.reads[ ,1])
  seqID.reads <- seqID.reads[index, ]
  
  return(seqID.reads)
}

import.and.reformat.otu.table <- function(OTUtable){
  # import the OTU table and then pull out just the total reads for each seqID
  otu.table.path <- OTUtable
  
  otus <- read.table(file = otu.table.path, sep = "\t", header = TRUE, colClasses = "character")
  seqID <- otus[ ,1]
  otus <- as.matrix(otus[ ,-1])
  otus <- apply(X = otus, MARGIN = 2, FUN = as.numeric)
  reads <- rowSums(otus)
  # otu table is normalized by sample reads- re-normalize by total reads
  reads <- reads / sum(reads) * 100
  
  seqID.reads <- data.frame(seqID, reads, stringsAsFactors = FALSE)
  return(seqID.reads)
}

check.for.seqID.reads <- function(PWDpath){
  # check if you have seqid.reads already or need to make it
  x <- list.files(path = PWDpath)
  check <- any(x == "total.reads.per.seqID.csv")
  return(check)
}

import.ids.above <- function(FilePaths, PidentsUsed){
  # import the ids.above. files that list seqIDs above the pident cutoff
  # This is used inside the add.totals.to.read.summaries() function
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

get.conflict.seqIDs <- function(ConflictsFolders, PidentsUsed, TaxaLevels){
  # import files comparing conflicts at each level, record the seqIDs into a list of lists
  # structure: outer list- each pident, inner lists- seqIDs at each taxa level
  pident.folders <- ConflictsFolders
  pident.values <- PidentsUsed
  num.taxa.levels <- TaxaLevels
  
  # forcing plots go to tribe (7) level, conflicts plots go to lineage (5) level
  conflict.ids.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL, Clade = NULL, Tribe = NULL)
  conflict.ids.start <- conflict.ids.start[1:num.taxa.levels]
  
  all.pidents <- list(NULL)
  # for each pident folder
  for (p in 1:length(pident.values)){
    all.files <- list.files(pident.folders[p])
    conflict.ids <- conflict.ids.start
    
    # for each taxonomy file
    for (t in 1:num.taxa.levels){
      conflict.ids[[t]] <- read.csv(file = paste(pident.folders[p], "/", all.files[t], sep = ""), header = TRUE, colClasses = "character")
      # only want the first column, the seqIDs
      conflict.ids[[t]] <- conflict.ids[[t]][ ,1]
    }
    
    all.pidents[[p]] <- conflict.ids
    names(all.pidents)[p] <- pident.values[p]
  }
  return(all.pidents)
}

find.reads.per.seqID <- function(ReadsTable, ConflictsList, Forcing = FALSE){
  # create a list of total reads that matches the structure of the list of conflicting seqIDs
  otu.reads <- ReadsTable
  conflict.ids <- ConflictsList
  
  taxa.list.start <- list(Kingdom = NULL, Phylum = NULL, Class = NULL, Order = NULL, Lineage = NULL, Clade = NULL, Tribe = NULL)
  num.taxa.levels <- length(conflict.ids[[1]])
  taxa.list.start <- taxa.list.start[1:num.taxa.levels]
  
  # for each upper list level (pident folder)
  pidents.list <- list(NULL)
  for (p in 1:length(conflict.ids)){

    # for each inner list level (taxonomy level)
    taxa.list <- taxa.list.start
    for (t in 1:num.taxa.levels){
      # find the total reads corresponding to each of the conflict seqIDs
      if (length(conflict.ids[[p]][[t]]) > 0){
        conf.ids <- conflict.ids[[p]][[t]]
        reads.ids <- otu.reads$seqID
        combined.ids <- c(conf.ids, reads.ids)
        index <- duplicated(combined.ids)
        index <- index[-c(1:length(conf.ids))]
        taxa.list[[t]] <- otu.reads[index,2]
      }else{
        taxa.list[[t]] <- 0
      }
    }
    # assign element p of the outer list the inner list it contains
    pidents.list[[p]] <- taxa.list
    names(pidents.list)[p] <- names(conflict.ids)[p]
  }
  return(pidents.list)
}

generate.summary.table.of.reads <- function(ReadsList, Forcing = FALSE){
  # collapse the list of total reads per seqID into a table of total reads: taxa level x pident
  reads.list <- ReadsList
 
  # Set up empty matrix to fill
  num.taxa.levels <- length(reads.list[[1]])
  num.pidents <- length(reads.list)
  reads.summary <- matrix(0, nrow = num.taxa.levels, ncol = num.pidents)
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

add.totals.to.read.summaries <- function(ReadSummaryTable, AbundanceTable, PidentsUsed, CustomSeqIDs){
  # add to the read.summaries table the overall totals to match the otu.summaries structure
  # this involves importing the fw_classified_taxonomies.csv file to get the list of fw seqIDs
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
      ids.subset <- fw.ids[[p]]
      reads.ids <- seqID.reads[ ,1]
      index <- duplicated(c(ids.subset, reads.ids))
      index <- index[-c(1:length(ids.subset))]
      fw.reads[p] <- sum(seqID.reads[index,2])
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

import.bootstrap.pvalues <- function(ConflictFolders, PidentsUsed, FW = TRUE){
  # import bootstrap p-values to compare
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

import.taxonomy.file <- function(FilePath, Final = FALSE){
  # import the custom-only classified taxonomy for the "forcing" plots (function borrowed from find_classification_disagreements.R)
  taxonomy.file.path <- FilePath
  if (Final == FALSE){
    delimitor <- ";"
  }else if (Final == TRUE){
    delimitor <- ","
  }
  
  # simple imort b/c file format determined by previous script, shouldn't have any weird things
  tax <- read.table(taxonomy.file.path , sep = delimitor, fill=T, colClasses = "character", header = TRUE)
  
  # Convert into a character matrix for faster processing
  tax <- as.matrix(tax)
  
  # remove percent confidences from final file so that identical names will match
  if (Final == TRUE){
    tax <- apply(X = tax, MARGIN = 2, FUN = remove.parentheses)
  }
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(tax[ ,1])
  tax <- tax[index, ]
  
  # Remove row names that will not match between the data tables
  row.names(tax) <- NULL
  
  return(tax)
}

make.empty.list.structure <- function(ListNames){
  # the list contains zero under each name, not NULL.  NULL is truly empty, but can't have a name
  # this is used in group.taxa()
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

make.unclassifieds.unique <- function(Taxonomy){
  # Taxonomy is a matrix containing only the names, no abundances or OTU numbers
  # this function maintains row order, so those columns can be re-combined after this
  # this is used in group.taxa()
  taxa.levels <- 2:ncol(Taxonomy)
  
  # first add the level-above name to each "unclassified"
  for (t in taxa.levels){
    index <- which(Taxonomy[ ,t] == "unclassified")
    Taxonomy[index,t] <- paste(Taxonomy[index,t], Taxonomy[index,t - 1], sep = ".")
  }
  
  # then remove extra unclassifieds to get only unclassified.last_known_name
  remove.extra.u <- function(x){
    x <- gsub(pattern = ".unclassified", replacement = "", x = x)
    return(x)
  }
  Taxonomy <- apply(X = Taxonomy, MARGIN = 2, FUN = remove.extra.u)
  
  return(Taxonomy)
}

group.taxa <- function(Taxonomy, Abunds, ListNames = c("kingdom","phylum","class","order","lineage","clade","tribe"), UniqueUnclass){
  # Taxonomy is a matrix (or vector) containing only the names, no abundances or OTU numbers
  # Abunds is a matrix (or vector) containing the relative abundances in each sample
  # ListNames is a vector describing grouping level (phylum, class, order, etc- length should match # Tax columns)
  # This function groups the abundances by summing them for each Taxonomy Level
  # It produces a large list, where each element contains an abundance table grouped to that level
  
  # Optionally differentiate the "unclassified" names by their last known name:
  if (UniqueUnclass == TRUE){
    Taxonomy <- make.unclassifieds.unique(Taxonomy = Taxonomy)
  }
  
  # First create the lists of variable categories the function aggregate() requires:
  aggregation.list <- make.empty.list.structure(ListNames = ListNames)
  for(t in 1:length(aggregation.list)){
    aggregation.list[[t]] <- Taxonomy[ ,t]
  }
  
  # Next create the grouped-by-taxa-level list:
  taxa.list <- make.empty.list.structure(ListNames = ListNames)
  for (t in 1:length(taxa.list)){
    taxa.list[[t]] <- aggregate(x = Abunds, by = aggregation.list[1:t], FUN = "sum")
  }
  
  # If don't want "unclassifieds to be differentiated, group them together:
  # note: all unclassifieds at current level t get grouped. so their upper taxa levels are combined too,
  # resulting in a category of "CombinedTaxa" since many diverse things are grouped together potentially.
  if (UniqueUnclass == FALSE){
    for (t in 1:length(taxa.list)){
      index <- which(taxa.list[[t]][ ,t] == "unclassified")
      if (length(index) > 0){
        unclass.abund <- sum(taxa.list[[t]][index, (t + 1)])
        unclass.name <- c(rep.int(x = "CombinedTaxa", times = t - 1), "unclassified")
        unclass.row <- c(unclass.name, unclass.abund)
        taxa.list[[t]] <- taxa.list[[t]][-index, ]
        taxa.list[[t]] <- rbind(taxa.list[[t]], unclass.row)
        taxa.list[[t]][ ,(t + 1)] <- as.numeric(taxa.list[[t]][ ,(t + 1)])
      }
    }
  }
  
  # Change order to by abundance instead of alphabetical:
  for (t in 1:length(taxa.list)){
    index <- order(taxa.list[[t]][ ,(t + 1)], decreasing = TRUE)
    taxa.list[[t]] <- taxa.list[[t]][index, ]
    colnames(taxa.list[[t]])[(t + 1)] <- "abund"
    row.names(taxa.list[[t]]) <- NULL
  }
  
  return(taxa.list)
}

find.top.taxa.by.total.reads <- function(TaxonomyList, NumberTopTaxa = "all", RemoveUnclass){
  # Find the most abundant taxa at each taxonomy level - and decide if you want to include unclassifieds in your rank abund calcs
  grouped.taxa <- TaxonomyList
  
  # trim or don't trim the total number of results
  num.taxa <- NULL
  if (NumberTopTaxa == "all"){
    for (t in 1:length(grouped.taxa)){
      num.taxa[t] <- nrow(grouped.taxa[[t]])
    }
  }else{
    num.taxa <- rep.int(x = NumberTopTaxa, times = length(grouped.taxa))
  }
  
  # arrange highest to lowest total reads
  grouped.taxa.ord <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  for (t in 1:7){
    index <- order(grouped.taxa[[t]][ ,(t + 1)], decreasing = TRUE)
    grouped.taxa.ord[[t]] <- grouped.taxa[[t]][index, ]
  }
  
  if (RemoveUnclass == TRUE){
    # remove unclassified taxa b/c those really can't be compared on the same taxa level, and likely wouldn't be included in a "top taxa" analysis anyway.
    not.unclassifieds <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
    for (t in 1:7){
      index <- grep(x = grouped.taxa.ord[[t]][ ,t], pattern =  "unclassified.*", value = FALSE )
      # this is necessary because you can not use -0 as an index
      if (length(index) != 0){
        not.unclassifieds[[t]] <- grouped.taxa.ord[[t]][-index, ]
      }else{
        not.unclassifieds[[t]] <- grouped.taxa.ord[[t]]
      }
    }
  }else{
    # name is misleading now but this is easier than chaning all the names:
    not.unclassifieds <- grouped.taxa.ord
  }
  
  # look just at the top levels
  grouped.taxa.top <- list("kingdom"=NULL,"phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  for (t in 1:7){
    if (nrow(not.unclassifieds[[t]]) < num.taxa[t]){
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]]
    }else{
      grouped.taxa.top[[t]] <- not.unclassifieds[[t]][1:num.taxa[t], ]
    }
  }
  return(grouped.taxa.top)
}

remove.parentheses <- function(x){
  # Remove parentheses and % confidence- script from find_classification_disagreements.R
  # use with apply()
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

find.forcing.diffs <- function(TopFinalList, AllForcedList){
  # find the differences btwn the workflow and the custom-only classifications
  top.final <- TopFinalList
  all.forced <- AllForcedList
  
  # add a fw.reads column and a difference column to the top final taxa list
  for (t in 1:length(top.final)){
    fw.reads <- NULL
    difference <- NULL
    for (n in 1:nrow(top.final[[t]])){
      index <- which(all.forced[[t]][ ,t] == top.final[[t]][n,t])
      if (length(index) > 0){
        temp.reads <- all.forced[[t]][index,(t+1)]
      }else{
        temp.reads <- 0
      }
      fw.reads <- c(fw.reads, temp.reads)
      difference <- c(difference, temp.reads - top.final[[t]][n,(t+1)])
    }
    top.final[[t]] <- cbind(top.final[[t]], fw.reads, difference)
  }
  
  # add columns for values used to make the stacked grey, red, and blue barplot
  for (t in 1:length(top.final)){
    x <- top.final[[t]]
    grey.bars <- x$abund
    red.bars <- x$difference
    blue.bars <- abs(x$difference)
    
    for (r in 1:nrow(x)){
      if (x$difference[r] > 0){
        blue.bars[r] <- 0
      }else if (x$difference[r] < 0){
        red.bars[r] <- 0
        grey.bars[r] <- grey.bars[r] - blue.bars[r]
      }
    }
    
    top.final[[t]] <- cbind(top.final[[t]], grey.bars, red.bars, blue.bars)
  }
  return(top.final)
}

filter.out.low.abund <- function(TaxaList, CutoffVector){
  # this removes taxa that are low abundance based on their max of red, grey, and blue heights. (so narrows in on interesting bars while maintaining origional ranks)
  for (t in 1:7){
    max.heights <- NULL
    for (r in 1:nrow(TaxaList[[t]])){
      max.heights[r] <- max(TaxaList[[t]][r,(t + 4):(t + 6)])
    }
    index <- which(max.heights < CutoffVector[t])
    if (length(index) > 0){
      TaxaList[[t]] <- TaxaList[[t]][-index, ]
    }
  }
  return(TaxaList)
}

import.pvalues.above.bootstrap.cutoff <- function(FolderPaths, Pidents){
  pident.folders <- FolderPaths
  pident.values <- Pidents
  
  tax.nums.list <- list(NULL)
  for (p in 1:length(pident.values)){
    file.name <- paste(pident.folders[p],"bootstrap_pvalues.csv", sep = "/")
    tax.nums.list[[p]] <- read.csv(file = file.name, header = T, colClasses = "character")
    names(tax.nums.list)[p] <- paste("pident", pident.values[p], sep = "_")
  }
  return(tax.nums.list)
}

convert.to.yes.no.classified <- function(PvaluesList){
  tax.nums <- PvaluesList
  for(p in 1:length(tax.nums)){
    t.n <- apply(X = tax.nums[[p]][ ,-1], MARGIN = 2, FUN = as.numeric)
    t.n <- t.n != 0
    t.n <- apply(X = t.n, MARGIN = 2, FUN = as.numeric)
    tax.nums[[p]][ ,-1] <- t.n
  }
  return(tax.nums)
}

combine.tax.nums.with.seqid.reads <- function(PvaluesList, TotReadsTable){
  seqID.reads <- TotReadsTable
  tax.TF <- PvaluesList
  
  # put seqIDs in same order in all tables
  index.reads <- order(seqID.reads[ ,1])
  seqID.reads <- seqID.reads[index.reads, ]
  
  index.pvals <- NULL
  check.order <- NULL
  for (p in 1:length(tax.TF)){
    index.pvals <- order(tax.TF[[p]][ ,1])
    tax.TF[[p]] <- tax.TF[[p]][index.pvals, ]
    check.order[p] <- all.equal(tax.TF[[p]][ ,1], seqID.reads[ ,1])
  }
  if (all(check.order != TRUE)){
    cat("Fuck there's something messed up w/ the tax.nums indexing in function combine.tax.nums.with.seqid.reads()")
  }
  
  # add read counts to the tax.nums tables
  for (p in 1:length(tax.TF)){
    tax.TF[[p]] <- cbind(seqID.reads, tax.TF[[p]][ ,-1])
  }
  return(tax.TF)
}

calc.reads.per.name <- function(TaxReadsList){
  tax.reads <- TaxReadsList
  for (p in 1:length(tax.reads)){
    tax <- tax.reads[[p]][ ,-c(1:2)]
    tax <- as.matrix(tax)
    reads <- tax.reads[[p]][ ,2]
    t.r <- tax * reads
    tax.reads[[p]] <- data.frame(tax.reads[[p]][ ,1], t.r, stringsAsFactors = FALSE)
  }
  return(tax.reads)
}

create.summary.matrix <- function(TaxList){
  tax.sum <- matrix(0, nrow = length(TaxList), ncol = ncol(TaxList[[1]]) - 1)
  row.names(tax.sum) <- names(TaxList)
  colnames(tax.sum) <- colnames(TaxList[[1]][ ,-1])
  for (p in 1:length(TaxList)){
    tax <- TaxList[[p]][ ,-1]
    tax.sum[p, ] <- colSums(tax)
  }
  return(tax.sum)
}

convert.to.percent <- function(TaxSum, Total){
  TaxSum <- TaxSum / Total * 100
  return(TaxSum)
}

find.alpha.diversity <- function(Taxonomy){
  alpha.div <-NULL
  for (t in 1:ncol(Taxonomy)){
    alpha.div[t] <- length(unique(Taxonomy[ ,t]))
  }
  names(alpha.div) <- colnames(Taxonomy)
  return(alpha.div)
}

find.percent.classified <- function(GroupedTaxaList){
  # note: unclassifieds can't have been made unique
  perc.class <- NULL
  for (t in 1:length(GroupedTaxaList)){
    index <- which(GroupedTaxaList[[t]][ ,t] == "unclassified")
    if (length(index) > 0){
      perc.class[t] <- sum(GroupedTaxaList[[t]][-index, t + 1])
    }else{
      perc.class[t] <- sum(GroupedTaxaList[[t]][ ,t + 1])
    }
  }
  names(perc.class) <- names(GroupedTaxaList)
  return(perc.class)
}

export.fw.stats <- function(AlphaDiv, PercClass, FolderPath){
  file.name <- paste(FolderPath, "Custom-only_classification_stats.csv", sep = "/")
  fw.stats <- data.frame(Percent.Reads.Classified = PercClass, Alpha.Diversity = AlphaDiv)
  write.csv(x = fw.stats, file = file.name, quote = FALSE)
  cat("Made Datafile: ", file.name, "\n")
  return(fw.stats)
}

# ---------------------------------------------------------------------------------------------------------------------
# Define functions to plot the data
# ---------------------------------------------------------------------------------------------------------------------

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
  plot.name <- paste(results.file.path, "/Classification_Disagreements_of_", plot.of, "_", plot.as, "s", yplotlabel, db.label, ".png", sep = "")
  png(filename = plot.name, width = 5, height = 5, units = "in", res = 100)
  
  # Set up an empty plot
  plot.title <- paste("Disagreements Between Custom and General Taxonomy Classifications\n(Are we forcing OTUs into our favorite groups?)", db.plot.label, sep = "")
  y.label <- paste("Classification Disagreements (", plot.as, " ", plot.of, "s)", sep = "")
  x.label <- "\"full length\" pident cutoff (similar to ANI of read to custom database)"
  plot(x = 0, type = "n", ylim = c(0, ymax), xlim = c(min(pidents), max(pidents)),
       main = plot.title, cex.main = .8, ylab = y.label, cex.lab = .8, xlab = x.label)
  
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
  cat("made plot: ", plot.name, "\n")
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
  file.name <- paste(results.file.path, "/Custom_db_contributions_by_", plot.of, "_", plot.as, "s.png", sep = "")
  png(filename = file.name, width = 5, height = 5, units = "in", res = 100)
  
  # Set up and empty plot
  plot.title <- "Contribution of Custom Taxonomy \n to Total Classifications"
  x.label <- "\"Full Length\" percent identity (similar to ANI)"
  y.label <- paste("OTUs Classified by Custom Taxonomy (", plot.as, " ", plot.of, "s)", sep = "")
  plot(x = pidents, y = num.fw, type = "n", main = plot.title, cex.main = 1, xlab = x.label, ylab = y.label)
  
  # Fill plot with beautiful data
  lines(x = pidents, y = num.fw, col = "lightsalmon", lwd = 1.5)
  points(x = pidents, y = num.fw, col = "lightsalmon", pch = 19, cex = 1.3)
  
  unnecessary.message <- dev.off()
  cat("made plot: ", file.name, "\n")
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
  
  tot <- num.forced[8]
  perc.forced <- num.forced[-8,1,drop = FALSE] / tot * 100
  
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

plot.most.misleading.forced.otus <- function(ReadsPerForcedSeqIDs, ForcedSeqIDs, ReadsPerSeqID, OutputFolder, PlottingLevels = 1:7, NumOTUs = 50){
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
  seqID.recorder <- list(kingdom = NULL, phylum = NULL, class = NULL, order = NULL, lineage = NULL, clade = NULL, tribe = NULL)
  for (t in 1:7){
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

plot.most.misleading.forced.taxa <- function(TopTaxaList, ForcedTaxonomy, ForcedReadsList, ForcedSeqIDsList, ResultsFolder, PlottingLevels = 1:7, TotalReads){
  top.taxa <- TopTaxaList
  forced.taxonomy <- ForcedTaxonomy
  forced.seqID.reads <- ForcedReadsList
  forced.seqIDs <- ForcedSeqIDsList
  plots.folder.path <- ResultsFolder
  tot.reads <- TotalReads
  
  # reformatting- don't repeat this step twice
  forced.seqID.reads <- forced.seqID.reads[[1]]
  forced.seqIDs <- forced.seqIDs[[1]]
  for (t in 1:7){
    top.taxa[[t]] <- cbind(top.taxa[[t]], forced = 0, correct = 0)
    names(top.taxa[[t]])[(t+1)] <- "total"
  }
  
  # find the contributing forced seqIDs and tally their reads
  taxa.reads.forced <- 0
  for (t in 1:7){ 
    
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
  for (t in 1:7){
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
    plot.file.name <- paste(plots.folder.path, "/", t, "_Forcing_of_top_taxa-", names(top.taxa)[t], ".png", sep = "")
    png(filename = plot.file.name, width = 7, height = 5, units = "in", res = 100)
    par(mar = c(10,5,5,2))
    barplot(height = t(as.matrix(top.taxa[[t]][ , (t + 3):(t + 2)])), col = c("grey", "red"),
            main = paste("Forcing at the ", names(top.taxa)[t], " level\nwhen using only the custom database", sep = ""),
            names.arg = taxa.names, legend.text = c("correct","forced"), args.legend = list(x = "topright", bty = "n"),
            las = 2, xpd = TRUE, axes = FALSE, ylab = "Percent of total reads", cex.lab = 1.2)
    axis(side = 2, at = y.axis.ticks, labels = y.axis.labels)
    unnecessary.message <- dev.off()
    
    # save data used to make plot as .csv files:
    data.file.name <- paste(plots.folder.path, "/", t, "_Forcing_of_top_taxa-", names(top.taxa)[t], ".csv", sep = "")
    write.csv(x = top.taxa[[t]], file = data.file.name, quote = FALSE, row.names = FALSE)
    
    cat("\ncreated plot file:", plot.file.name, "\nand corresponding data file:", data.file.name, "\n")
  }
}

plot.forcing.diffs <- function(TopTaxaList, NumBars, FolderPath, PlottingLevels = 1:length(TopTaxaList)){
  top.taxa <- TopTaxaList
  
  # pull out and format just the data for the plot
  stacked.data <- list(NULL) 
  for (t in PlottingLevels){
    num.bars <- NumBars
    if (nrow(top.taxa[[t]]) < num.bars){
      num.bars <- nrow(top.taxa[[t]])
    }
    stacked.data[[t]] <- top.taxa[[t]][1:num.bars, c(t + 4, t + 5, t + 6)]
    row.names(stacked.data[[t]]) <- top.taxa[[t]][1:num.bars ,t]
    stacked.data[[t]] <- as.matrix(stacked.data[[t]])
    stacked.data[[t]] <- t(stacked.data[[t]])
    names(stacked.data)[t] <- names(top.taxa)[t]
  }
  
  # export plots and data!
  for (t in PlottingLevels){
    # descriptive file names
    plot.name <- paste(FolderPath, "/", t, "-Forcing_Diffs-", names(stacked.data)[t], ".png", sep = "")
    csv.name <- paste(FolderPath, "/", t, "-Forcing_Diffs-", names(stacked.data)[t], ".csv", sep = "")
    
    # really long names go off the plot
    taxa.names <- sub(pattern = ".*__", replacement = "", x = colnames(stacked.data[[t]]))
    # this ends up looking weird b/c the label are not monospaced font so some of the "short" names get cut off too.
    # for (n in 1:length(taxa.names)){
    #   if (nchar(taxa.names[n], type = "chars") > 20){
    #     taxa.names[n] <- substr(x = taxa.names[n], start = 1, stop = 19)
    #     taxa.names[n] <- paste(taxa.names[n], "-", sep = "")
    #   }
    # }
    colnames(stacked.data[[t]]) <- taxa.names

    # # adjust y-axis to go to the max value
    # max.bar <- max(stacked.data[[t]][1, ] + stacked.data[[t]][2, ]) # grey bar + red bar
    # y.axis.ticks <- c(0, max.bar * (1/4), max.bar * (1/2), max.bar * (3/4), max.bar)
    # y.axis.labels <- round(x = y.axis.ticks / tot.reads * 100, digits = 0)

    # make the plots!
    png(filename = plot.name, width = 7, height = 5, units = "in", res = 100)
    par(mar = c(10,5,5,2))
    barplot(height = stacked.data[[t]], beside = FALSE, col = c("grey","red","blue"), main = names(stacked.data)[t], las = 2, ylab = "Relative Abundance (% reads)", border = NA)
    legend(x = "topright", legend = c("Added Inaccuracy", "Lost Diversity"), fill = c("red", "blue"), border = FALSE, bty = "n", inset = .05)
    unnecessary.message <- dev.off()
    cat("made plot: ", plot.name, "\n")
    
    # export the data! ... do the full data not just stacked data:
    write.csv(x = top.taxa[[t]], file = csv.name, quote = FALSE, row.names = FALSE)
    cat("made datafile: ", csv.name, "\n")
  }
}

export.summary.table <- function(SummaryTable, FolderPath, ByOTU = TRUE, Percents = FALSE){
  folder.path <- FolderPath
  sum.table <- SummaryTable
  
  if (ByOTU == TRUE){
    data.type <- "OTUs"
  }else if (ByOTU == FALSE){
    data.type <- "reads"
  }
  
  make.percent <- function(x){
    perc <- x / tot.read * 100
    return(perc)
  }
  
  if (Percents == TRUE){
    tot.read <- sum.table[nrow(sum.table), 1]
    sum.table <- apply(X = sum.table, MARGIN = 2, FUN = make.percent)
    Percents <- "_percent"
  }else{
    Percents <- ""
  }
  
  file.name <- paste("conflict_summary_by", Percents, "_", data.type, ".csv", sep = "")
  write.csv(x = sum.table, file = paste(FolderPath, "/", file.name, sep = ""), quote = FALSE)
  cat("made datafile: ", file.name, "\n")
}

export.total.forcing.stats <- function(OtuSum, ReadSum, FolderPath){
  # export as a table the total forcing stats (this made the green bar plots that are commented out now) (forcing option)
  # the OTU summary is not in percents, the reads summary is
  otu.tot.row <- nrow(OtuSum)
  OtuSum <- OtuSum / OtuSum[otu.tot.row,1] * 100
  
  both.sum <- cbind(OtuSum, ReadSum)
  colnames(both.sum) <- c("perc.otus.forced", "perc.reads.forced")
  row.names(both.sum)[otu.tot.row] <- "total.reads.or.otus"
  
  file.name <- paste(FolderPath, "/total_percent_forcing_with_only_custom_db.csv", sep = "")
  write.csv(x = both.sum, file = file.name, quote = FALSE, row.names = TRUE)
  cat("Made datafile: ", file.name, "\n")
}

export.grouped.list <- function(Grouped, PlotsPath, FolderName){
  # exported the grouped lists so that each level is another file in a folder (forcing option).
  folder.path <- paste(PlotsPath, "/", FolderName, "/", sep = "")
  dir.create(path = folder.path, showWarnings = FALSE) # the warning is if the folder already exists, but it works regardless.
  for (t in 1:length(Grouped)){
    file.name = paste(folder.path, t, "_grouped_by_", names(Grouped)[t], ".csv", sep = "")
    write.csv(x = Grouped[[t]], file = file.name, quote = FALSE, row.names = FALSE)
  }
  cat("Made datafiles of total reads per taxon in folder: ", folder.path, "\n")
}

plot.total.classified <- function(SummaryMatrix, PidentValues, FolderPath, DataType){
  pidents <- PidentValues
  sum.named <- SummaryMatrix
  
  line.col <- rainbow(n = ncol(sum.named), v = .4)
  x.lim <- c(min(pidents), max(pidents))
  y.lim <- c(40,100)
  y.label <- paste("Percent Classified (", DataType, ")", sep = "")
  x.label <- "Percent Identity Cutoff"
  plot.title <- expression(bold("Choose a percent identity cutoff that maximizes classification"))
  plot.filename <- paste(FolderPath, "/Percent_", DataType, "_Classified_by_Pident.png", sep = "")
  taxa.levels <- sub(pattern = ".fw", replacement = "", x = colnames(sum.named))
  
  png(filename = plot.filename, width = 10, height = 5, units = "in", res = 100)
  
  par(mfrow = c(1,ncol(sum.named)), omi = c(.4,.3,.3,.1), mai = c(.2,.3,.3,0))
  for (t in 1:ncol(sum.named)){
    ass <- sum.named[ ,t]
    # basic plot
    plot(x = pidents, y = ass, col = line.col[t], type = "l", ann = F, lwd = 3, axes = F)
    mtext(text = taxa.levels[t], side = 3, line = .5, outer = F, cex = 1.2, col = line.col[t])
    
    # vertical max line
    index <- which(ass == max(ass))
    max.names <- pidents[index]
    abline(v = max.names, col = adjustcolor(col = line.col[t], alpha.f = .3), lwd = 3)
    
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
    span <- max(ass) - min(ass)
    y.ax <- c(min(ass), min(ass) + (span * 1/3), min(ass) + (span * 2/3), max(ass))
    y.ax.lab <- round(x = y.ax, digits = 0)
    empty.y.labels <- rep("", times = length(y.ax))
    axis(side = 2, at = y.ax, labels = empty.y.labels)
    mtext(text = y.ax.lab, side = 2, line = .7, at = y.ax)
  } 
  mtext(text = plot.title, side = 3, line = .5, outer = T, cex = 1.2)
  mtext(text = x.label, side = 1, line = 1.5, outer = T, cex = 1.2)
  mtext(text = y.label, side = 2, line = .5, outer = T, cex = 1.2)
  
  unnecessary.message <- dev.off()
  cat("made plot: ", plot.filename, "\n")
}

export.total.classified.summary <- function(SummaryMatrix, FolderPath, DataType){
  data.filename <- paste(FolderPath, "/Percent_", DataType, "_Classified_by_Pident.csv", sep = "")
  colnames(SummaryMatrix) <- sub(pattern = ".fw", replacement = "", x = colnames(SummaryMatrix))
  write.csv(x = SummaryMatrix, file = data.filename, quote = F)
  cat("Made Datafile: ", data.filename, "\n")
}

# ---------------------------------------------------------------------------------------------------------------------
# Use Functions
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# first check if all you're doing is making the seqID.reads file
# ---------------------------------------------------------------------------------------------------------------------
if (userprefs[2] == "MakeSeqIDReadsOnly"){
# ---------------------------------------------------------------------------------------------------------------------
  seqID.reads <- import.and.reformat.otu.table(OTUtable = otu.table.path) 
  write.table(x = seqID.reads, file = "total.reads.per.seqID.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("made the datafile: ", "total.reads.per.seqID.csv")

# ---------------------------------------------------------------------------------------------------------------------
# second check if this is the optional "forcing plot"
# ---------------------------------------------------------------------------------------------------------------------
}else if (forcing.folder.path != "regular"){
# ---------------------------------------------------------------------------------------------------------------------

  plots.folder.path <- paste(plots.folder.path, "step_15_5b_Improvement_over_custom-only", sep = "/")
  make.plot.directory(FolderPath = plots.folder.path)
  
  otus.forced <- import.forcing.conflicts(ForcingFolder = forcing.folder.path)
  
  if (check.for.seqID.reads(PWDpath = present.working.directory)){
    seqID.reads <- import.seqID.reads(FilePath = seqID.reads.file.path) # this was exported previously by this script in step 14
  }else{
    seqID.reads <- import.and.reformat.otu.table(OTUtable = otu.table.path) # wasn't exported if you chose pident and skipped step 14
  }
  
  forced.seqIDs <- get.conflict.seqIDs(ConflictsFolders = forcing.folder.path, PidentsUsed = "forcing", TaxaLevels = 7)
  
  forced.seqID.reads <- find.reads.per.seqID(ReadsTable = seqID.reads, ConflictsList = forced.seqIDs, Forcing = TRUE)
  
  read.summaries <- generate.summary.table.of.reads(ReadsList = forced.seqID.reads, Forcing = TRUE)
  
  tot.reads <- sum(seqID.reads$reads)
  
  reads.forced <- rbind(read.summaries, tot.reads)
  
  forced.taxonomy <- import.taxonomy.file(FilePath = forced.taxonomy.file)
  
  grouped.forced.taxa <- group.taxa(Taxonomy = forced.taxonomy[ ,-1], Abunds = seqID.reads[ ,2], UniqueUnclass = FALSE)
  
  # top.taxa <- find.top.taxa.by.total.reads(TaxonomyList = grouped.forced.taxa, NumberTopTaxa = 20, RemoveUnclass = FALSE) # only used for custom-only rank abund order (plot.most.misleading.forced.taxa)                   .                              
  
  final.taxonomy <- import.taxonomy.file(FilePath = final.taxonomy.file, Final = TRUE)

  grouped.final.taxa <- group.taxa(Taxonomy = final.taxonomy[ ,-1], Abunds = seqID.reads[ ,2], UniqueUnclass = FALSE)

  top.final.taxa <- find.top.taxa.by.total.reads(TaxonomyList = grouped.final.taxa, NumberTopTaxa = "all", RemoveUnclass = FALSE) # all here to export all data, numbars (workflow rank abund) is determined in plot function.

  top.final.taxa <- find.forcing.diffs(TopFinalList = top.final.taxa, AllForcedList = grouped.forced.taxa)
  
  top.final.taxa <- filter.out.low.abund(TaxaList = top.final.taxa, CutoffVector = c(0, .5, .5, .5, .5, .5, .5))
  
  # plot.percent.forced(ForcingTable = otus.forced, ResultsFolder = plots.folder.path, ByReads = FALSE)
  
  # plot.percent.forced(ForcingTable = reads.forced, ResultsFolder = plots.folder.path, ByReads = TRUE)
  
  # export.total.forcing.stats(OtuSum = otus.forced, ReadSum = reads.forced, FolderPath = plots.folder.path)
  
  # plot.most.misleading.forced.otus(ReadsPerForcedSeqIDs = forced.seqID.reads, ForcedSeqIDs = forced.seqIDs, 
  #                                  ReadsPerSeqID = seqID.reads, OutputFolder = plots.folder.path, PlottingLevels = 1:7)
  
  # plot.most.misleading.forced.taxa(TopTaxaList = top.taxa, ForcedTaxonomy = forced.taxonomy, 
  #                                  ForcedReadsList = forced.seqID.reads, ForcedSeqIDsList = forced.seqIDs, 
  #                                  ResultsFolder = plots.folder.path, PlottingLevels = 1:7, TotalReads = tot.reads)

  plot.forcing.diffs(TopTaxaList = top.final.taxa, NumBars = 800, FolderPath = plots.folder.path)
  
  # export.grouped.list(Grouped = grouped.forced.taxa, PlotsPath = plots.folder.path, FolderName = "ForcedTaxonomyGroups")
  # export.grouped.list(Grouped = grouped.final.taxa, PlotsPath = plots.folder.path, FolderName = "FinalTaxonomyGroups")
  
  # Export stats on alpha diversity and percent classified using only the FW database:
  fw.unique.unclass <- make.unclassifieds.unique(Taxonomy = forced.taxonomy[ ,-1])
  fw.alpha <- find.alpha.diversity(Taxonomy = fw.unique.unclass)
  
  fw.perc.class <- find.percent.classified(GroupedTaxaList = grouped.forced.taxa)
  
  fw.stats <- export.fw.stats(AlphaDiv = fw.alpha, PercClass = fw.perc.class, FolderPath = plots.folder.path)
  
# ---------------------------------------------------------------------------------------------------------------------
# If not then do the normal comparison for choosing pident cutoff
# ---------------------------------------------------------------------------------------------------------------------
}else{
# ---------------------------------------------------------------------------------------------------------------------
  
# examine custom taxonomy disagreements and contribution by number OTUs -----------------------------------------------
  
  plots.folder.path <- paste(plots.folder.path, "step_14_Choose_pident_cutoff", sep = "/")
  make.plot.directory(FolderPath = plots.folder.path)
  
  otu.summaries <- import.all.conflict.summaries(ConflictFolders = pident.folders, PidentsUsed = pident.values)
  export.summary.table(SummaryTable = otu.summaries, FolderPath = plots.folder.path, ByOTU = TRUE, Percents = FALSE)
  export.summary.table(SummaryTable = otu.summaries, FolderPath = plots.folder.path, ByOTU = TRUE, Percents = TRUE)
  
  # I removed this from the terminal command also because it is misleading and confusing and doesn't add to the plot. left a commented out input line for use within RStudio at the top.
  # db.summary <- import.database.conflicts(DatabaseFolder = db.conflicts.folder.path)
  
  # ignore lineage- it's way higher than the others doesn't fit on graphs
  # ignore order- so many unclassifieds that it looks lower than class and is confusing
  otu.summaries <- otu.summaries[c(-5,-4), ]
  # db.summary <- db.summary[-5, ]
  
  # These plots people generally found confusing, and the message was clearer by looking at the datatables exported above anyway
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, y.axis.limit = 10) #zooms on any phylum forcing
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE, y.axis.limit = 1)
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, DBconflicts = db.summary, y.axis.limit = max(db.summary[,2]))
  # plot.num.forced(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, DBconflicts = db.summary)
  
  # plot.num.classified.outs(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = FALSE)
  plot.num.classified.outs(ConflictSummaryTable = otu.summaries, ResultsFolder = plots.folder.path, AsPercent = TRUE)
  
# examine custom taxonomy disagreements and contribution by number reads ----------------------------------------------
  
  seqID.reads <- import.and.reformat.otu.table(OTUtable = otu.table.path)
  
  conflict.seqIDs <- get.conflict.seqIDs(ConflictsFolders = pident.folders, PidentsUsed = pident.values, TaxaLevels = 5)
  
  conflict.seqID.reads <- find.reads.per.seqID(ReadsTable = seqID.reads, ConflictsList = conflict.seqIDs)
  
  read.summaries <- generate.summary.table.of.reads(ReadsList = conflict.seqID.reads)
  
  custom.seqIDs <- import.ids.above(FilePaths = ids.file.paths, PidentsUsed = pident.values)
  
  read.summaries <- add.totals.to.read.summaries(ReadSummaryTable = read.summaries, AbundanceTable = seqID.reads, 
                                                 PidentsUsed = pident.values, CustomSeqIDs = custom.seqIDs)
  
  # seqID.reads is normalized to total reads, so these are both out of 100 % . An absolute read number is meaningless. 
  # export.summary.table(SummaryTable = read.summaries, FolderPath = plots.folder.path, ByOTU = FALSE)
  export.summary.table(SummaryTable = read.summaries, FolderPath = plots.folder.path, ByOTU = FALSE, Percents = TRUE)
  
  #leave out lineage on the plots b/c it's too much higher
  #leave out order because there's too many unclassifieds, makes it hard to interpret
  read.summaries <- read.summaries[c(-5,-4), ]
  
  # These plots people generally found confusing, and the message was clearer by looking at the datatable exported above anyway
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE)
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, y.axis.limit = 10000)
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE)
  # plot.num.forced(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE, y.axis.limit = .5)
  
  # plot.num.classified.outs(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = FALSE)
  plot.num.classified.outs(ConflictSummaryTable = read.summaries, ResultsFolder = plots.folder.path, ByReads = TRUE, AsPercent = TRUE)
  
  # export the reads per seqID for use in the plot_classification_improvement.R script
  write.table(x = seqID.reads, file = "total.reads.per.seqID.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("made the datafile: ", seqID.reads.file.path, "\n")
  
# examine total percent of dataset that is classified by read and OTU -------------------------------------------------
  
  tax.nums <- import.pvalues.above.bootstrap.cutoff(FolderPaths = pident.folders, Pidents = pident.values)
  tax.TF <- convert.to.yes.no.classified(PvaluesList = tax.nums)
  
  tax.reads <- combine.tax.nums.with.seqid.reads(PvaluesList = tax.TF, TotReadsTable = seqID.reads)
  tax.reads <- calc.reads.per.name(TaxReadsList = tax.reads)
  
  otus.named <- create.summary.matrix(TaxList = tax.TF)
  reads.named <- create.summary.matrix(TaxList = tax.reads)
  
  otus.named <- convert.to.percent(TaxSum = otus.named, Total = nrow(seqID.reads))
  reads.named <- convert.to.percent(TaxSum = reads.named, Total = sum(seqID.reads[ ,2]))
  
  export.total.classified.summary(SummaryMatrix = otus.named, FolderPath = plots.folder.path, DataType = "OTUs")
  export.total.classified.summary(SummaryMatrix = reads.named, FolderPath = plots.folder.path, DataType = "Reads")
  
  # plot.total.classified(SummaryMatrix = reads.named, PidentValues = pident.values, FolderPath = plots.folder.path, DataType = "Reads")
  # k and p mostly all classified, t shows forcing
  plot.total.classified(SummaryMatrix = reads.named[ ,-c(1:2,7)], PidentValues = pident.values, FolderPath = plots.folder.path, DataType = "Reads")
  # trim down the total pidents (for when I did 90:100)
  # plot.total.classified(SummaryMatrix = reads.named[-c(1:5),-c(1:2,7)], PidentValues = pident.values[-c(1:5)], FolderPath = plots.folder.path, DataType = "Reads")

  # The plots by OTU are not very helpful b/c of the strong effects of forcing
  # plot.total.classified(SummaryMatrix = otus.named, PidentValues = pident.values, FolderPath = plots.folder.path, DataType = "OTUs")
  # only see diversity lost with lowering pident at the phylum level- that's probably the cyanos
  plot.total.classified(SummaryMatrix = otus.named[ ,2:4], PidentValues = pident.values, FolderPath = plots.folder.path, DataType = "OTUs")
  # otherwise you actually see forcing at the OTU-level. as more OTUs are put into FW, more OTUs are classified. 
  # these gains must be from stuff forced into a classification in FW, but unclassified in GG.
  
  
# examine pident cutoff relationship to bootstrap p-values -this plot is not useful. ----------------------------------
  
  # fw.pvalues <- import.bootstrap.pvalues(ConflictFolders = pident.folders, PidentsUsed = pident.values, FW = TRUE)
  # gg.pvalues <- import.bootstrap.pvalues(ConflictFolders = pident.folders, PidentsUsed = pident.values, FW = FALSE)
  # plot.bootstrap.percents(FWpValues = fw.pvalues, GGpValues = gg.pvalues, ResultsFolder = plots.folder.path)
}

# ---------------------------------------------------------------------------------------------------------------------


