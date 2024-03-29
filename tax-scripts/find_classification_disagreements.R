# 8-27-15 RRR
# -------------------------------------------------------------

# The goal of this script is to identify the best BLAST cutoff to use
# in my new taxonomy assignment pipeline.
# I will do this by comparing the taxonomy assignments of the seqIDs
# that met different BLAST cutoff percent IDs and were classified by FW
# with the classification of all sequences with GG.
# Then I'll choose a BLAST cutoff level that includes the most sequences
# for classification by FW but still avoids conflicting classifications
# with GG at phylum and class levels.
# More details in 8-26-15 Analysis notes, 10-19-15 Analysis Notes, 10-21-15 Analysis Notes
# This script generates files/data that the next scripts use for plots/stats

# variable names in this script follow this format: fw.percents.fw.only 
  # the first fw refers to which workflow, the one combining fw + gg ("fw") or gg alone "gg"
  # the fw.only means that only the seqIDs that the workflow assigned with fw are included in the table

# This script is also used to identify and fix conflicts between the two databases
# in the create new FreshTrain version pipeline. It finds 
# coarse-level name mismatches that need to be changed in FreshTrain
# a reference list of what freshwater names are called in silva 
# and a list of any fine-level silva names that exist when FreshTrain is unnamed
# (because these might want to be removed from FreshTrain) 

# -------------------------------------------------------------
# Receive arguments from terminal command line
# -------------------------------------------------------------

# Example syntax for pident comparison, final taxonomy generation, final taxonomy w/out comparisons to general-only, database comparison, and forcing analysis, respectively:

# Rscript find_classification_disagreements.R otus.98.taxonomy otus.general.taxonomy ids.above.98 conflicts_98 98 85 70 
# Rscript find_classification_disagreements.R otus.98.taxonomy otus.general.taxonomy ids.above.98 conflicts_98 98 85 70 final
# Rscript find_classification_disagreements.R otus.98.taxonomy quickie ids.above.98 conflicts_98 98 85 70 final
# Rscript find_classification_disagreements.R custom.custom.taxonomy custom.general.taxonomy NA conflicts_database NA NA 70 database
# Rscript find_classification_disagreements.R otus.custom.taxonomy otus.98.85.70.taxonomy ids.above.98 conflicts_forcing NA 85 70 forcing

userprefs <- commandArgs(trailingOnly = TRUE)

# -------------------------------------------------------------
# input for troubleshooting
# -------------------------------------------------------------
# CONFLICT FINDING ONLY:
# cat("fuck you forgot to comment out the file paths in find_classification_disagreements.R!")
# userprefs <- c("~/Desktop/limony/IRD_2020-07-14/limony_Ie.03_R1e2t247_R2e2t195_dmx_dada.98.taxonomy",
#                "~/Desktop/limony/IRD_2020-07-14/limony_Ie.03_R1e2t247_R2e2t195_dmx_dada.silva_nr_v138_taxass.taxonomy",
#                "~/Desktop/limony/IRD_2020-07-14/ids.above.98",
#                "~/Desktop/limony/IRD_2020-07-14/conflicts_98",
#                98,
#                80,
#                80)
# FINAL TABLE GENERATION: note you do need the otus.general.taxonomy file b/c it's used to prep a file for plot_classification_improvement.R in step 15
# cat("fuck you forgot to comment out the file paths in find_classification_disagreements.R!")
# userprefs <- c("~/Desktop/taxass_MEv34_silva138cyanos/MEv34.98.taxonomy",
#                "~/Desktop/taxass_MEv34_silva138cyanos/MEv34.silva_nr_v138_taxass_cyano_edits.taxonomy", # make this one "quickie" if skipping general-only classification
#                "~/Desktop/taxass_MEv34_silva138cyanos/ids.above.98",
#                "~/Desktop/taxass_MEv34_silva138cyanos/conflicts_98",
#                98,
#                80,
#                80,
#                "final")
# DATABASE COMPARISON: part of optional step 11.5
# cat("fuck you forgot to comment out the file paths in find_classification_disagreements.R!")
# userprefs <- c("../../2020-06-02_update_freshtrain/FT_semicol_noprefix_mothur_unnamed_semicol.tax",
#                "../../2020-06-02_update_freshtrain/FT.silva_semicol.taxonomy",
#                "NA",
#                "../../2020-06-02_update_freshtrain/conflicts_database_6-12-20",
#                NA,
#                NA,
#                70,
#                "database")
# FORCING ANALYSIS: part of optional step 15.5
# cat("fuck you forgot to comment out the file paths in find_classification_disagreements.R!")
# userprefs <- c("~/Desktop/taxass_MEv34_silva138cyanos/MEv34.FreshTrain15Jun2020silva138.taxonomy",
#                "~/Desktop/taxass_MEv34_silva138cyanos/MEv34.98.80.80.taxonomy",
#                "~/Desktop/taxass_MEv34_silva138cyanos/ids.above.98",
#                "~/Desktop/taxass_MEv34_silva138cyanos/conflicts_forcing",
#                NA,
#                80,
#                80,
#                "forcing")
# -------------------------------------------------------------

fw.plus.gg.tax.file.path <- userprefs[1]
gg.only.tax.file.path <- userprefs[2]
fw.seq.ids.file.path <- userprefs[3]
results.folder.path <- userprefs[4]
blast.pident.cutoff <- userprefs[5]
taxonomy.pvalue.cutoff.fw <- userprefs[6]
taxonomy.pvalue.cutoff.gg <- userprefs[7]
final.or.database <- userprefs[8]
if (length(userprefs) < 8){final.or.database <- "non-empty string"}


# -------------------------------------------------------------
# Pre-determined Output File Names
# -------------------------------------------------------------
# all of the pre-determined names are used by other scripts and shouldn't be changed. 

# when "final" flag is used, these are exported into the working directory:
otus <- sub(pattern = "^.*/", replacement = "", x = fw.plus.gg.tax.file.path)
otus <- sub(pattern = "\\.\\d+\\.taxonomy$", replacement = "", x = otus)
file.name.final.taxonomy <- paste(otus, ".", blast.pident.cutoff, ".", taxonomy.pvalue.cutoff.fw, ".", taxonomy.pvalue.cutoff.gg, ".taxonomy", sep = "")
file.name.workflow.pvalues <- "final.taxonomy.pvalues"
file.name.general.pvalues <- "final.general.pvalues"
file.name.workflow.names <- "final.taxonomy.names"
file.name.general.names <- "final.general.names"

# when forcing flag is used, these are exported into the working directory:
otus.custom <- sub(pattern = "^.*/", replacement = "", x = fw.plus.gg.tax.file.path)
otus.custom <- sub(pattern = "\\.taxonomy$", replacement = "", x = otus.custom)
file.name.custom.only.taxonomy <- paste(otus.custom, ".", taxonomy.pvalue.cutoff.fw, ".taxonomy", sep = "")

# when script is used for conflict finding, these are exported into the conflicts folder:
file.name.bootstrap.pvalues <- paste(results.folder.path,"bootstrap_pvalues.csv", sep = "/")
file.name.summary.stats <- paste(results.folder.path, "conflicts_summary.csv", sep = "/")
# the conflict file names at each taxa level are created inside of find.conflicting.names()


# -------------------------------------------------------------
# Define Functions for Import and Formatting
# -------------------------------------------------------------

print.poem <- function(){
  # Entertain user with a poem while they wait:
  cat("\nAnd the Days Are Not Full Enough\nby Ezra Pound\n\nAnd the days are not full enough\nAnd the nights are not full enough\nAnd life slips by like a field mouse\n\tNot shaking the grass.\n\n")
}

import.FW.names <- function(FilePath){
  # Import the otu taxonomies assigned with both FW and GG
  fw.plus.gg.tax.file.path <- FilePath
  # Avoid errors from variable row lengths by checking length of all rows (fill=T only checks 1st 5 rows)
  numcol <- max(count.fields(fw.plus.gg.tax.file.path, sep=";"))
  fw <- read.table(fw.plus.gg.tax.file.path, sep=";", fill=T, colClasses = "character", col.names=1:numcol)
  return(fw)
}

import.GG.names <- function(FilePath, final.names = FALSE){
  # Import the otu taxonomies assigned with only GG, or with the workflow's GG+FW combo
  
  if (final.names == FALSE){
    gg.only.tax.file.path <- FilePath
    # Avoid errors from variable row lengths by checking length of all rows (fill=T only checks 1st 5 rows)
    numcol <- max(count.fields(gg.only.tax.file.path, sep=";"))
    gg <- read.table(gg.only.tax.file.path, sep=";", fill=T, colClasses = "character", col.names=1:numcol)

  }else if (final.names == TRUE){   # b/c the final taxonomy file generated is .csv for user conveniance
    workflow.tax.file.path <- FilePath
    numcol <- max(count.fields(workflow.tax.file.path, sep=","))
    gg <- read.table(workflow.tax.file.path, sep=",", fill=T, colClasses = "character", col.names=1:numcol, header = TRUE)
    
  }
  return(gg)
}

reformat.fw <- function(FWtable){
  # Reformat the blast-FW-GG workflow-assigned taxonomy table
  
  fw <- as.matrix(FWtable)
  
  # Remove strain and empty 10th column
  fw <- fw[,-c(9,10)]
  colnames(fw) <- c("seqID.fw","kingdom.fw","phylum.fw","class.fw","order.fw","lineage.fw","clade.fw","tribe.fw")
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(fw[,1])
  fw <- fw[index,]
  
  # Remove row names that will not match between the data tables
  row.names(fw) <- NULL
  
  return(fw)
}

reformat.gg <- function(GGtable){
  # Reformat the green genes taxonomy table
  
  gg <- as.matrix(GGtable)
  
  # Remove empty 9th column
  gg <- gg[,1:8]
  colnames(gg) <- c("seqID.gg","kingdom.gg","phylum.gg","class.gg","order.gg","family.gg","genus.gg","species.gg")
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(gg[,1])
  gg <- gg[index,]
  
  # Remove row names that will not match between the data tables
  row.names(gg) <- NULL
  
  return(gg)
}

check.files.match <- function(FWtable, GGtable){
  # Check that the order of the names is the same in each file, print warning if not.
  
  gg <- GGtable
  fw <- FWtable
  
  # only care about values, but all equal doesn't return true if the names are different (ex. "seqID.fw" != "seqID.gg")
  row.names(fw) <- NULL
  row.names(gg) <- NULL
  colnames(fw) <- NULL
  colnames(gg) <- NULL
  
  order.check <- all.equal(fw[,1], gg[,1])
  if (order.check != TRUE){
    cat("\n\nWARNING!! The Indexing of your files is messed up!!\nYou have to fix this or all your results will be wrong!\n(somehow the seqID order is not the same in both, even though the reformatting functions in this script are supposed to do that. Possibly the seqID names don't match btwn your abund and fasta starting files.)\n\n")
  }
}

import.FW.seq.IDs <- function(FilePath){
  # Import the freshwater sequence IDs as determined by the BLAST cutoff in workflow step 4
  fw.seq.ids.file.path <- FilePath
  fw.seqs <- read.table(file = fw.seq.ids.file.path, colClasses = "character")
  colnames(fw.seqs) <- "seqID.fw"
  return(fw.seqs)
}

remove.parentheses <- function(x){
  # Remove parentheses of % confidence so that names in each table match exactly
  # call this using apply
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

make.degenerates.unique <- function(tax, voldemort){ # modified from reformat_taxonomy_nomenclature.R
  # voldemort is whatever text you're making unique (b/c it cannot be named... voldemort...)
  
  rownames(tax) <- tax[ ,1]
  my.colname <- colnames(tax)[1]
  tax <- tax[ ,-1]
  
  for (t in 1:ncol(tax)){
    index <- which(tax[ ,t] == voldemort)
    tax[index,t] <- paste(tax[index,t], tax[index,t - 1], sep = ".")
    tax[index,t] <- remove.parentheses(x = tax[index,t])
  }
  
  remove.extra.u <- function(x){
    dot.voldemort <- paste0(".", voldemort)
    x <- gsub(pattern = dot.voldemort, replacement = "", x = x)
    return(x)
  }
  tax <- apply(X = tax, MARGIN = 2, FUN = remove.extra.u)
  
  tax <- cbind(rownames(tax), tax)
  colnames(tax)[1] <- my.colname
  row.names(tax) <- NULL
  
  return(tax)
}


# -------------------------------------------------------------
# Define Functions for Data Analysis
# -------------------------------------------------------------

find.fw.seqid.indeces <- function(FullTable, FWids){
  gg.and.fw <- FullTable
  fw.ids <- FWids
  # duplicate marks the 2nd occurance as TRUE
  # combine the ids (fw ids first, all ids second) into one vector with duplicate ids
  # index of FW in all ids is the index of the duplicates - the number of fw ids
  all.ids <- c(fw.ids[ ,1], gg.and.fw[ ,1])
  index <- which(duplicated(all.ids) == TRUE)
  index <- index - nrow(fw.ids)
  return(index)
}

uniform.unclass.names <- function(TaxonomyTable, database = FALSE, final = FALSE){
  # this is used in the function do.bootstrap.cutoff()
  # all unnamed refs are called "unclassified" with no bootstrap value for classification comparisons
  # so for this comparison purpose, unnamed in database is same as unassigned in classification
  # database names were already made uniform in reformat_taxonomy_nomenclature.R
  # and have format unnamed.UpperName
  
  tax <- TaxonomyTable
  
  if(final == TRUE){
    return(tax) # leave "unnamed.Parent" in the final version
  }
  
  index <- grep(pattern = "unnamed", x = tax[ ,-1], value = FALSE, invert = FALSE, ignore.case = TRUE)
  
  if(database == TRUE){
    tax[ ,-1][index] <- "unnamed" # change unnamed.Parent to just be unnamed
  }else{
    tax[ ,-1][index] <- "unclassified" # change all unnamed to be indistinguishable from unclassified
  }
  
  return(tax)
}

pull.out.percent <- function(text){
  # Given a single taxonomy name, pull out the bootstrap percent.
  # text = k__Bacteria(100) returns (as character) 100, But text = "unclassified" returns (as character) "unclassified"
  text <- sub(pattern = '.*\\(', replacement = '\\(', x = text)
  text <- sub(pattern = '\\(', replacement = '', x = text)
  text <- sub(pattern = '\\)', replacement = '', x = text)
  return(text)
}

do.bootstrap.cutoff <- function(TaxonomyTable, BootstrapCutoff, ...){
  # Apply classification bootstrap value cutoff 
  # Everything below the supplied cutoff value is changed to "unclassified"
  # (by bootstrap cutoff I mean the stat that's the % of times it got classified into that taxonomy cluster)
  
  tax <- as.matrix(TaxonomyTable)
  cutoff <- as.numeric(BootstrapCutoff)
  
  # make all unclassified things uniformly called "unclassified"
  tax <- uniform.unclass.names(TaxonomyTable = tax)
  
  # create a matrix of bootstrap numbers and then a T/F matrix
  tax.nums <- apply(tax[ ,2:ncol(tax)], 2, pull.out.percent)
  index <- grep(pattern = "unclassified", x = tax.nums, ignore.case = TRUE, value = F)
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  tax.TF <- tax.nums >= cutoff
  
  # make sure you never get a "classified" under an "unclassified" (this prevents downstream bugs)
  # tax.TF <- t(apply(X = tax.TF, MARGIN = 1, FUN = cummin))
  # took this out for silva, because frequently happens like this in the silva database. hopefully unique unclass will solve bugs.
  
  # make all names in the taxonomy table unclassified if they're below the bootstrap cutoff
  index <- which(tax.TF == 0)
  tax.names <- tax[ ,-1]
  tax.names[index] <- "unclassified"
  tax[ ,-1] <- tax.names
    
  return(tax)
}

find.conflicting.names <- function(FWtable, GGtable, FWtable_percents, GGtable_percents, TaxaLevel, tracker, FolderPath, forcing = FALSE, Database = FALSE){
  # Find seqs misclassified at a given phylogenetic level, t
  fw <- FWtable
  gg <- GGtable
  fw.percents <- FWtable_percents
  gg.percents <- GGtable_percents
  t <- TaxaLevel
  num.mismatches <- tracker
  results.folder.path <- FolderPath
  
  taxa.names <- c("kingdom","phylum","class","order","lineage","clade","tribe")
  
  # compare names in column t+1, because first columns are seqID and t=1 is kingdom, t=7 is tribe
  # ignore names that say unclassified, because that doesn't count as a disagreement
  # except in forcing comparison when fw giving any erroneous name counts as forcing
  # except in database comparison when a missing fw name when silva has a name is a problem
  if (forcing == TRUE){
    index <- which(gg[,t+1] != fw[,t+1] & fw[,t+1] != "unclassified")
  }else if(Database == TRUE){
    index <- which(gg[,t+1] != fw[,t+1] & gg[,t+1] != "unnamed" & gg[,t+1] != "unclassified")
  }else{  
    index <- which(gg[,t+1] != fw[,t+1] & gg[,t+1] != "unclassified" & fw[,t+1] != "unclassified")
  }
  
  # cat("there are ", length(index), " conflicting names at ", taxa.names[t], " level\n")
  num.mismatches[t] <- length(index)
  
  if (Database == TRUE){
    
    # At coarse levels, identify unique mismatched upper-names of lineages (to change in FreshTrain release)
    if (t <= 5){
      check.files.match(GGtable = gg[index,1:(t+1), drop=F], FWtable = fw[index,1:6, drop=F])
      unique.conflicts <- cbind(  gg[index,1:(t+1), drop=F],           fw[index,1:6, drop=F])
    
    # At fine levels, get a list of FreshTrain vs. Silva names (to use as a reference of what non-FreshTrain names mean)
    }else{ 
      check.files.match(GGtable = gg[index,1:(t+1), drop=F], FWtable = fw[index,1:(t+1), drop=F])
      unique.conflicts <- cbind(  gg[index,1:(t+1), drop=F],           fw[index,1:(t+1), drop=F])
    }
    
    unique.conflicts <- unique.conflicts[ ,-c(1,t+2)]
    unique.conflicts <- unique(unique.conflicts)
    unique.conflict.file <- paste(results.folder.path, "/", "unique_conflicts_", t, "_", taxa.names[t],".csv", sep="")
    write.csv(unique.conflicts, file = unique.conflict.file, row.names = FALSE)
    cat("Made file: ", unique.conflict.file, "\n")
    
    if (t == 6){
      # at lineage all fw has names, at species no silva has names, so only check for instances of silva have name when fw is unnamed at clade level
      index <- which(unique.conflicts[ ,t+t] == "unnamed")
      better.in.gg <- unique.conflicts[index, ]
      better.in.gg.file <- paste(results.folder.path, "/", "unnamed_", taxa.names[t], "s_that_have_general_names.csv", sep="")
      write.csv(better.in.gg, file = better.in.gg.file, row.names = FALSE)
      cat("Made file: ", better.in.gg.file, "\n")
    }
  }
  
  # Compare the conflicting tables in entirety, use the original files with percents still in it
  conflicting <- cbind(gg.percents[index, ,drop=F], fw.percents[index, ,drop=F])
  check.files.match(FWtable = conflicting[ ,9:16,drop=F], GGtable = conflicting[ ,1:8,drop=F])
  conflict.file <- paste(results.folder.path, "/", t, "_", taxa.names[t],"_conflicts.csv", sep="")
  write.csv(conflicting, file = conflict.file, row.names = FALSE)
  cat("Made file: ", conflict.file, "\n")
  
  #Track the number of mismatches at each level
  return(num.mismatches)
}

create.summary.vector <- function(Forcing = FALSE){
  # Set up a summary vector to fill
  if(Forcing == TRUE){
    num.mismatches <- vector(mode = "numeric", length = 7)
    names(num.mismatches) <- c("kingdom","phylum","class","order","lineage","clade","tribe")
  }else{
    num.mismatches <- vector(mode = "numeric", length = 5)
    names(num.mismatches) <- c("kingdom","phylum","class","order","lineage")
  }
  return(num.mismatches)
}

export.summary.stats <- function(SummaryVector, FW_seqs, ALL_seqs, FileName){
  # Format and export the summary vector
  num.mismatches <- SummaryVector
  fw.fw.only <- FW_seqs
  fw.percents <- ALL_seqs
  num.mismatches <- c(num.mismatches,"numFWseqs" = nrow(fw.fw.only))
  num.mismatches <- c(num.mismatches, "numALLseqs" = nrow(fw.percents))
  num.mismatches <- data.frame("TaxaLevel" = names(num.mismatches),"NumConflicts" = num.mismatches, row.names = NULL)
  write.csv(num.mismatches, file = FileName, row.names = FALSE)
}

view.bootstraps <- function(TaxonomyTable){
  tax <- TaxonomyTable
  
  # create a matrix of bootstrap numbers, copy from do.bootstrap.cutoff()
  tax.nums <- apply(tax[ ,-1], 2, pull.out.percent)
  index <- which(tax.nums == "unclassified")
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  tax.nums <- cbind(tax[ ,1],tax.nums)
  colnames(tax.nums)[1] <- colnames(tax)[1]
  return(tax.nums)
}


# -------------------------------------------------------------
# Use Functions
# -------------------------------------------------------------

# -------------------------------------------------------------
# Generate a final taxonomy file:
if (final.or.database == "final" | final.or.database == "Final" | final.or.database == "FINAL"){
# -------------------------------------------------------------  
  print.poem()
  
  fw.percents <- import.FW.names(FilePath = fw.plus.gg.tax.file.path)
  fw.percents <- reformat.fw(FWtable = fw.percents)
  
  fw.seq.ids <- import.FW.seq.IDs(FilePath = fw.seq.ids.file.path)
  fw.indeces <- find.fw.seqid.indeces(FullTable = fw.percents, FWids = fw.seq.ids)
  
  fw.percents.fw.only <- fw.percents[fw.indeces, ]
  fw.percents.gg.only <- fw.percents[-fw.indeces, ]
  
  final.taxonomy.fw.only <- do.bootstrap.cutoff(TaxonomyTable = fw.percents.fw.only, BootstrapCutoff = taxonomy.pvalue.cutoff.fw, final = TRUE)
  final.taxonomy.gg.only <- do.bootstrap.cutoff(TaxonomyTable = fw.percents.gg.only, BootstrapCutoff = taxonomy.pvalue.cutoff.gg, final = TRUE)
  
  final.taxonomy <- rbind(final.taxonomy.fw.only, final.taxonomy.gg.only)
  colnames(final.taxonomy) <- c("seqID","kingdom","phylum","class","order","lineage","clade","tribe")
  final.simple <- final.taxonomy # keep unclassified simple for generating tax.names below
  final.taxonomy <- make.degenerates.unique(tax = final.taxonomy, voldemort = "unclassified")
  
  write.table(x = final.taxonomy, file = file.name.final.taxonomy, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Made file: ", file.name.final.taxonomy, "\n")
  
  # the following will be used by the plot_classification_improvement.R script
  if (gg.only.tax.file.path != "quickie"){
    final.simple <- do.bootstrap.cutoff(TaxonomyTable = final.simple, BootstrapCutoff = 1, final = FALSE) # make unnamed say unclassified, unclassified.Parent say unclassified
    
    tax.nums <- view.bootstraps(TaxonomyTable = final.simple)
    write.table(x = tax.nums, file = file.name.workflow.pvalues, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Made file: ", file.name.workflow.pvalues, "\n")
    
    tax.names <- apply(final.simple, 2, remove.parentheses)
    write.table(x = tax.names, file = file.name.workflow.names, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Made file: ", file.name.workflow.names, "\n")
    
    gg.percents <- import.GG.names(FilePath = gg.only.tax.file.path)
    gg.percents <- reformat.gg(GGtable = gg.percents)
    
    gg.taxonomy <- do.bootstrap.cutoff(TaxonomyTable = gg.percents, BootstrapCutoff = taxonomy.pvalue.cutoff.gg)
    colnames(gg.taxonomy) <- c("seqID","kingdom","phylum","class","order","lineage","clade","tribe")
    
    gg.nums <- view.bootstraps(TaxonomyTable = gg.taxonomy)
    write.table(x = gg.nums, file = file.name.general.pvalues, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Made file: ", file.name.general.pvalues, "\n")
    
    gg.names <- apply(X = gg.taxonomy, MARGIN = 2, FUN = remove.parentheses)
    write.table(x = gg.names, file = file.name.general.names, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Made file: ", file.name.general.names, "\n")
  }
  
  
# -------------------------------------------------------------
# Compare databases by looking at how GG classifies the FW representative sequences
}else if (final.or.database == "database"){
# -------------------------------------------------------------
  cat("doing database comparison\n")
  
  fw.percents <- import.FW.names(FilePath = fw.plus.gg.tax.file.path)
  gg.percents <- import.GG.names(FilePath = gg.only.tax.file.path)
  
  gg.percents <- reformat.gg(GGtable = gg.percents)
  fw.percents <- reformat.fw(FWtable = fw.percents)
  
  check.files.match(FWtable = fw.percents, GGtable = gg.percents)
  
  fw <- fw.percents # database only has names
  fw <- uniform.unclass.names(TaxonomyTable = fw, database = TRUE)
  gg.percents <- uniform.unclass.names(TaxonomyTable = gg.percents, database = TRUE)
  gg <- do.bootstrap.cutoff(TaxonomyTable = gg.percents, BootstrapCutoff = taxonomy.pvalue.cutoff.gg) # don't think I need database = TRUE here b/c already used it in prev line... double check
  gg <- apply(gg, 2, remove.parentheses)
  
  check.files.match(FWtable = fw, GGtable = gg)
  
  # Files written in find.conflicting.names() loop: the "TaxaLevel_conflicts.csv" that puts taxonomy tables side by side
  # File written afer loop: the "conflicts_summary.csv" that lists how many conflicts were at each level, and how many seqs were classified by FW
  num.mismatches <- create.summary.vector()
  for (t in 1:7){
    num.mismatches <- find.conflicting.names(FWtable = fw, GGtable = gg,
                                             FWtable_percents = fw.percents, 
                                             GGtable_percents = gg.percents, 
                                             TaxaLevel = t, tracker = num.mismatches,
                                             FolderPath = results.folder.path, Database = TRUE)
  }
  export.summary.stats(SummaryVector = num.mismatches, FW_seqs = fw, ALL_seqs = fw, FileName = file.name.summary.stats)
  cat("Made file: ", file.name.summary.stats, "\n")



# -------------------------------------------------------------
# Look at how bad it'd be if you used only the custom database instead of the workflow
}else if (final.or.database == "forcing"){
# -------------------------------------------------------------  
  cat("examining how the custom database would have classified the dissimilar sequences that didn't belong in it\n")
  
  fw.percents <- import.FW.names(FilePath = fw.plus.gg.tax.file.path) # fw.plus.gg.tax.file.path is actually FW-only here
  gg.percents <- import.GG.names(FilePath = gg.only.tax.file.path, final.names = TRUE) # gg.only.tax.file.path is final workflow taxonomy (fw+gg)
  
  gg.percents <- reformat.gg(GGtable = gg.percents)
  fw.percents <- reformat.fw(FWtable = fw.percents)
  
  check.files.match(FWtable = fw.percents, GGtable = gg.percents)
  
  fw.seq.ids <- import.FW.seq.IDs(FilePath = fw.seq.ids.file.path)
  fw.indeces <- find.fw.seqid.indeces(FullTable = fw.percents, FWids = fw.seq.ids)
  
  fw.percents.gg.only <- fw.percents[-fw.indeces, ]
  gg.percents.gg.only <- gg.percents[-fw.indeces, ]
  
  check.files.match(FWtable = fw.percents.gg.only, GGtable = gg.percents.gg.only)
  
  gg.gg.only <- do.bootstrap.cutoff(TaxonomyTable = gg.percents.gg.only, BootstrapCutoff = taxonomy.pvalue.cutoff.gg)
   
  # Need a custom-only table with bootstrap applied for custom-only forcing plot in plot_classification_disagreements.R with forcing argument
  fw <- do.bootstrap.cutoff(TaxonomyTable = fw.percents, BootstrapCutoff = taxonomy.pvalue.cutoff.fw)
  
  fw.gg.only <- fw[-fw.indeces, ]
  
  check.files.match(FWtable = fw.gg.only, GGtable = gg.gg.only)
  
  # export clean fw-db-only table for forcing plots later
  fw <- apply(fw, 2, remove.parentheses)
  write.table(x = fw, file = file.name.custom.only.taxonomy, sep = ";")
  cat("Made datafile: ", file.name.custom.only.taxonomy, "\n")
  
  fw.gg.only <- apply(fw.gg.only, 2, remove.parentheses)
  gg.gg.only <- apply(gg.gg.only, 2, remove.parentheses)
  
  check.files.match(FWtable = fw.gg.only, GGtable = gg.gg.only)
  
  # Generate the files comparing classifications made by fw to those of gg
  # Generate a summary file listing the total number of classification disagreements at each level
  # Files written in find.conflicting.names() loop: the "TaxaLevel_conflicts.csv" that puts taxonomy tables side by side
  # File written afer loop: the "conflicts_summary.csv" that lists how many conflicts were at each level, and how many seqs were classified by FW
  num.mismatches <- create.summary.vector(Forcing = TRUE)
  for (t in 1:7){
    num.mismatches <- find.conflicting.names(FWtable = fw.gg.only, GGtable = gg.gg.only,
                                             FWtable_percents = fw.percents.gg.only,
                                             GGtable_percents = gg.percents.gg.only, 
                                             TaxaLevel = t, tracker = num.mismatches,
                                             FolderPath = results.folder.path, forcing = TRUE)
  }
  export.summary.stats(SummaryVector = num.mismatches, FW_seqs = fw.gg.only, ALL_seqs = fw.percents, FileName = file.name.summary.stats)
  cat("Made file: ", file.name.summary.stats, "\n")
  
  
  
# ------------------------------------------------------------- 
# Only compare the classifications made by the fw database to the gg classifications, not full tax tables
}else{
# -------------------------------------------------------------  
  cat("comparing seqIDs that TaxAss classified with the custom database to how the general database would have classified them.\n")
  
  fw.percents <- import.FW.names(FilePath = fw.plus.gg.tax.file.path)
  gg.percents <- import.GG.names(FilePath = gg.only.tax.file.path)
  
  gg.percents <- reformat.gg(GGtable = gg.percents)
  fw.percents <- reformat.fw(FWtable = fw.percents)
  
  check.files.match(FWtable = fw.percents, GGtable = gg.percents)
  
  fw.seq.ids <- import.FW.seq.IDs(FilePath = fw.seq.ids.file.path)
  fw.indeces <- find.fw.seqid.indeces(FullTable = fw.percents, FWids = fw.seq.ids)
  
  fw.percents.fw.only <- fw.percents[fw.indeces, ]
  gg.percents.fw.only <- gg.percents[fw.indeces, ]
  
  check.files.match(FWtable = fw.percents.fw.only, GGtable = gg.percents.fw.only)
  
  # save the bootstrap-applied fw-table fw.cutoff.fw.only for creating tax.nums below
  fw.cutoff.fw.only <- do.bootstrap.cutoff(TaxonomyTable = fw.percents.fw.only, BootstrapCutoff = taxonomy.pvalue.cutoff.fw)
  gg.fw.only <- do.bootstrap.cutoff(TaxonomyTable = gg.percents.fw.only, BootstrapCutoff = taxonomy.pvalue.cutoff.gg)
  
  fw.fw.only <- apply(fw.cutoff.fw.only, 2, remove.parentheses)
  gg.fw.only <- apply(gg.fw.only, 2, remove.parentheses)
  
  check.files.match(FWtable = fw.fw.only, GGtable = gg.fw.only)
  
  # Generate the files comparing classifications made by fw to those of gg
  # Generate a summary file listing the total number of classification disagreements at each level
  # Files written in find.conflicting.names() loop: the "TaxaLevel_conflicts.csv" that puts taxonomy tables side by side
  # File written afer loop: the "conflicts_summary.csv" that lists how many conflicts were at each level, and how many seqs were classified by FW
  num.mismatches <- create.summary.vector()
  for (t in 1:5){
    num.mismatches <- find.conflicting.names(FWtable = fw.fw.only, GGtable = gg.fw.only,
                                             FWtable_percents = fw.percents.fw.only,
                                             GGtable_percents = gg.percents.fw.only, 
                                             TaxaLevel = t, tracker = num.mismatches,
                                             FolderPath = results.folder.path)
  }
  export.summary.stats(SummaryVector = num.mismatches, FW_seqs = fw.fw.only, ALL_seqs = fw.percents, FileName = file.name.summary.stats)
  cat("Made file: ", file.name.summary.stats, "\n")
  
  
# ----
  # the following will be used by the plot_classification_disagreements step to help choose a good pident: 
# ----
  fw.percents.gg.only <- fw.percents[-fw.indeces, ]
  fw.percents.gg.only <- do.bootstrap.cutoff(TaxonomyTable = fw.percents.gg.only, BootstrapCutoff = taxonomy.pvalue.cutoff.gg)
  
  workflow.taxonomy <- rbind(fw.cutoff.fw.only, fw.percents.gg.only)
  
  tax.nums <- view.bootstraps(TaxonomyTable = workflow.taxonomy)
  write.csv(x = tax.nums, file = file.name.bootstrap.pvalues, quote = FALSE, row.names = FALSE)
  cat("Made file: ", file.name.bootstrap.pvalues, "\n")
}









# ---- Old stuff ----

# # Generate a file of the fw-assigned taxonomies, and a matching table of just their bootstrap values
#   # File written: the "fw_classified_bootstraps.csv" that lists the bootstrap of all sequences classified by freshwater
#   # File written: the "fw_classified_taxonomies.csv" that lists the taxonomy of all sequences classified by freshwater
# fw.bootstraps <- view.bootstraps(TaxonomyTable = fw.percents.fw.only)
# write.csv(fw.bootstraps, file = paste(results.folder.path, "/", "fw_classified_bootstraps.csv", sep=""), row.names = FALSE)
# write.csv(fw.percents.fw.only, file = paste(results.folder.path, "/", "fw_classified_taxonomies.csv", sep=""), row.names = FALSE)
# 
# # Generate a file of the gg taxonomies for the fw-assigned sequences, and a matching table of gg bootstraps
#   # File written: the "gg_classified_bootstraps.csv" that lists the bootstrap of all sequences classified by freshwater
#   # File written: the "gg_classified_taxonomies.csv" that lists the taxonomy of all sequences classified by freshwater
# gg.bootstraps <- view.bootstraps(TaxonomyTable = gg.percents.fw.only)
# write.csv(gg.bootstraps, file = paste(results.folder.path, "/", "gg_classified_bootstraps.csv", sep=""), row.names = FALSE)
# write.csv(gg.percents.fw.only, file = paste(results.folder.path, "/", "gg_classified_taxonomies.csv", sep=""), row.names = FALSE)
# ----
