# ---- Convert mothur output into a Workflow-compatible OTU table (Step 0) ----

# Script written by Ben Peterson, edited by Robin Rohwer

# Q: What are all these f****** mothur formats???  

# A: a count_table is made early on, and is used if you do not go on to clustering steps.
# it is tab delimited, with column names. 
# rows are seqIDs, columns are samples, except col 1 = the seqID and col 2 = the total (summed) for each seqID
# the values are raw sequencing counts.

# a shared table ***help Ben!  I stopped before making one...**
# it is tab delimited
# Remove the first three columns
# First column is a label of the clustering cutoff
# Second column is sample name. 
# Never fear, these will be added in as row names! 
# Third column is the total number of OTUs

# Syntax in Command Line:

# Rscript StupidLongMothurName.count_table count_table otus.abund

# ---- Accept Arguments from Terminal Command Line ----

userprefs <- commandArgs(trailingOnly = TRUE)
path.to.mothur.file <- userprefs[1] 
mothur.extension <- userprefs[2]
path.to.output.file <- userprefs[3]

# path.to.mothur.file <- "~/Desktop/TaxonomyTrainingSets/BLASTing/StartFiles19/danube-10.count_table"
# mothur.extension <- "count_table" # because I think..... several ones will work??? clearly .shared worked for ben!
# path.to.output.file <- "~/Desktop/TaxonomyTrainingSets/BLASTing/StartFiles19/danube-10.otus"

# ---- Define Functions ----

import.and.format.count_table <- function(FilePath){
  count.table <- read.table(file = FilePath, header = TRUE, sep = "\t", colClasses = "character")
  # col 2 is "total", col 1 is seqID ("Representative_Sequence")
  seq.ids <- count.table[ ,1]
  count.table <- count.table[ ,-c(1:2)]
  # matrices are faster
  count.table <- apply(X = count.table, MARGIN = 2, FUN = as.numeric)
  row.names(count.table) <- seq.ids
  # for some reason I can only get apply to work over rows so transpose for later
  count.table <- t(count.table)
  return(count.table)
}

convert.to.rel.abund <- function(OTUs){
  sample.totals <- rowSums(OTUs)
  take.perc <- function(x){
    x <- x / sample.totals * 100
    return(x)
  }
  rel.otus <- apply(X = OTUs, MARGIN = 2, FUN = take.perc)
  return(rel.otus)
}

reformat.for.workflow <- function(OTUs){
  # rows = OTUs, cols = samples, col 1 = seqIDs
  otu.table <- t(OTUs)
  otu.table <- cbind(row.names(otu.table), otu.table)
  colnames(otu.table)[1] <- "seqID"
  return(otu.table)
}

export.workflow.table <- function(OTUs, FilePath){
  write.table(x = OTUs, file = FilePath, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("Reformatted OTU table saved as", FilePath)
}

import.and.format.shared <- function(FilePath){
  # This reads in the shared file that was generated in mothur
  OTU.table <- read.table("dataEdited/mothur/filterTest.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE)
  
  # Save sample names in separate vector
  sampleID <- OTU.table[, 2]
  
  # Remove the first three columns
  # First column is a label of the clustering cutoff
  # Second column is sample name. 
  # Never fear, these will be added in as row names! 
  # Third column is the total number of OTUs
  OTU.table <- OTU.table[, -c(1:3)]
  
  # Add sample names as row names
  row.names(OTU.table) <- sampleID
}

# ---- Use Functions ----

if (mothur.extension == "count_table"){
  otu.matrix <- import.and.format.count_table(FilePath = path.to.mothur.file)
}else if (mothur.extension == "shared"){
  otu.matrix <- import.and.format.shared(FilePath = path.to.mothur.file)
}

otu.matrix <- convert.to.rel.abund(OTUs = otu.matrix)

otu.matrix <- reformat.for.workflow(OTUs = otu.matrix)

export.workflow.table(OTUs = otu.matrix, FilePath = path.to.output.file)



