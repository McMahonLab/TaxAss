# 8-27-15 RRR

# The goal of this script is to identify the best BLAST cutoff to use 
# in my new taxonomy assignment pipeline.
# I will do this by comparing the taxonomy assignments of the seqIDs
# that met different BLAST cutoff percent IDs and were classified by FW
# with the classification of all sequences with GG.
# Then I'll choose a BLAST cutoff level that includes the most sequences
# for classification by FW but still avoids conflicting classifications
# with GG at phylum and class levels.
# More details in 8-26-15 Analysis notes.

#####
# Receive arguments from terminal command line
#####

# userprefs <- commandArgs(trailingOnly = TRUE)
# fw.plus.gg.tax.file.path <- userprefs[1]
# gg.only.tax.file.path <- userprefs[2]
# results.folder.path <- userprefs[3]
# taxonomy.bootstrap.cutoff <- userprefs[4]

fw.plus.gg.tax.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take4/otus.94.taxonomy"
gg.only.tax.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take4/otus.gg.taxonomy"
results.folder.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take4/compare_percID-94_to_gg-only/"
taxonomy.bootstrap.cutoff <- 60
fw.seq.ids.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take4/ids.above.94"

#####
# Define Functions for Import and Formatting
#####

# Import the freshwater sequence IDs as determined by the BLAST cutoff in workflow step 4
import.FW.seq.IDs <- function(){
  fw.seqs <- scan(file = fw.seq.ids.file.path)
  fw.seqs <- as.character(fw.seqs)
  return(fw.seqs)
}

# Import the otu taxonomies assigned with both FW and GG
import.FW.names <- function(){
  # Avoid errors from variable row lengths by checking length of all rows (fill=T only checks 1st 5 rows)
  numcol <- max(count.fields(fw.plus.gg.tax.file.path, sep=";"))
  fw <- read.table(fw.plus.gg.tax.file.path, sep=";", fill=T, stringsAsFactors = F, col.names=1:numcol)
  return(fw)
}

# Import the otu taxonomies assigned with only GG
import.GG.names <- function(){
  # Avoid errors from variable row lengths by checking length of all rows (fill=T only checks 1st 5 rows)
  numcol <- max(count.fields(gg.only.tax.file.path, sep=";"))
  gg <- read.table(gg.only.tax.file.path, sep=";", fill=T, stringsAsFactors = F, col.names=1:numcol)
  return(gg)
}

# Reformat the blast-FW-GG workflow-assigned taxonomy table
reformat.fw <- function(FWtable){
  fw <- FWtable
  
  # Remove strain and empty 10th column
  fw <- fw[,-c(9,10)]
  
  # Rename columns
  colnames(fw) <- c("seqID.fw","kingdom.fw","phylum.fw","class.fw","order.fw","linege.fw","clade.fw","tribe.fw")
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(fw[,1])
  fw <- fw[index,]
  
  # Convert into a character matrix (from dataframe w/ seqID's integer) for faster processing
  fw <- as.matrix(fw)
  
  # Remove row names that will not match between the data tables
  row.names(fw) <- NULL
  
  return(fw)
  
}

# Reformat the green genes taxonomy table
reformat.gg <- function(GGtable){
  gg <- GGtable
  
  # Remove empty 9th column
  gg <- gg[,1:8]
  
  # Rename columns
  colnames(gg) <- c("seqID.gg","kingdom.gg","phylum.gg","class.gg","order.gg","family.gg","genus.gg","species.gg")
  
  # Reorder sequence IDs so can match them to the other file
  index <- order(gg[,1])
  gg <- gg[index,]
  
  # Convert into a character matrix (from dataframe w/ seqID's integer) instead of a dataframe for faster processing
  gg <- as.matrix(gg)
  
  # Remove row names that will not match between the data tables
  row.names(gg) <- NULL
  
  return(gg)
}

# Check that the order of the names is the same in each file:
check.files.match <- function(FWtable, GGtable){
  gg <- GGtable
  fw <- FWtable
  
  order.check <- all.equal(fw[,1], gg[,1])
  if (order.check == FALSE){
    cat("\n\nWARNING!! The Indexing of your files is messed up!!\n\n")
  }
}

# Remove parentheses of % confidence so that names in each table match exactly
remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

# Entertain user with a poem while they wait:
print.poem <- function(){
  cat("\nAnd the Days Are Not Full Enough\nby Ezra Pound\n\nAnd the days are not full enough\nAnd the nights are not full enough\nAnd life slips by like a field mouse\n\tNot shaking the grass.\n\n")
}

#####
# Define Functions for Data Analysis
#####

# Find index of sequences that were classified by the freshwater database, not green genes
# note this returns indeces that match the gg and fw tables, so a vector of indeces, not a tax table
find.fw.indeces <- function(TaxonomyTable, SeqIDs){
  tax <- TaxonomyTable
  ids <- SeqIDs
  
  index <- NULL
  for (e in 1:length(ids)){
    index <- c(index, which(tax[,1] == ids[e]))
  }
  return(index)
}

# this makes all the unclassified/unknown/any other word for it names be uniformly called "unclassified"
# finds them b/c those names do not have bootstrap percents in parentheses, i.e. the (70)
uniform.unclass.names <- function(TaxonomyTable){
  tax <- TaxonomyTable
  
  # Warn user the names you are changing
  odd.entries <- unique(grep(pattern = '.*\\(', x <- tax[,2:8], value = TRUE, invert=T))
  if (length(odd.entries) > 0){
  cat("\nWarning: These names in your taxonomy table are missing a bootstrap taxonomy assignment value:\n\n",
      odd.entries,
      "\n\nThese names will be renamed as \"unclassified\". If that seems incorrect",
      "then you have to figure out why the parentheses are missing from them.", 
      "\nHave ALL your names here? Check that in the step 7 mothur command probs=T\n")
  }
  
  # Change all those names to unclassified (sometimes, for example, they might be "unknown")
  index <- grep(pattern = '.*\\(', x <- tax[,2:8], value = FALSE, invert=T)
  tax[,2:8][index] <- "unclassified"
  
  return(tax)
}

# Given a single taxonomy name, pull out the bootstrap percent.
# For example, text = k__Bacteria(100) would return (as character) 100
# But also, text = "unclassified" would return (as character) "unclassified"
pull.out.percent <- function(text){
  text <- sub(pattern = '.*\\(', replacement = '\\(', x = text)
  text <- sub(pattern = '\\(', replacement = '', x = text)
  text <- sub(pattern = '\\)', replacement = '', x = text)
  return(text)
}

# When the value supplied is FALSE, the text is changed to "unclassified"
# This is used in do.bootstrap.cutoff()
make.unclassified <- function(text,val){
  if (val == F){
    text <- "unclassified"
  }
  return(text)
}

# Apply classification bootstrap value cutoff 
# Everything below the supplied cutoff value is changed to "unclassified"
# (this stat is the % of times it got classified into that taxonomy cluster)
do.bootstrap.cutoff <- function(TaxonomyTable, BootstrapCutoff){
  tax <- as.matrix(TaxonomyTable)
  cutoff <- BootstrapCutoff
  
  # parentheses and are named something that is not "unclassified" or "unknown"
  tax.nums <- apply(tax[,2:ncol(tax)],2,pull.out.percent)
  index <- which(tax.nums == "unclassified")
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  tax.TF <- tax.nums >= cutoff
  
  # could you do this faster with sweep() ?  I think so...
  # ehh, but sweep requires stats calculated from itself, and this is a separate matrix
  for (r in 1:nrow(tax)){
    for (c in 2:ncol(tax)){
      tax[r,c] <- make.unclassified(text = tax[r,c], val = tax.TF[r,c-1])
    }
  }
  
  return(tax)
}

# Find seqs misclassified at a given phylogenetic level, t
find.conflicting.names <- function(FWtable, GGtable, FWtable_percents, GGtable_percents, TaxaLevel){
  fw <- FWtable
  gg <- GGtable
  fw.percents <- FWtable_percents
  gg.percents <- GGtable_percents
  t <- TaxaLevel
  
  taxa.names <- c("kingdom","phylum","class","order","lineage","clade","tribe")
  
  # compare names in column t+1, because first columns are seqID and t=1 is kingdom, t=7 is tribe
  # ignore names that say unclassified
  index <- which(gg[,t+1] != fw[,t+1] & gg[,t+1] != "unclassified" & fw[,t+1] != "unclassified")
  cat("there are ", length(index), " conflicting names at ", taxa.names[t], " level\n")
  
  # Compare the conflicting tables in entirety, use the original files with percents still in it
  conflicting <- cbind(gg.percents[index,,drop=F], fw.percents[index,,drop=F])
  
  # Check that the files still line up correctly
  check.files.match(FWtable = conflicting[,9:16,drop=F], GGtable = conflicting[,1:8,drop=F])
  
  # Export a file with the conflicting rows side by side.
  write.csv(conflicting, file = paste(results.folder.path, "/", taxa.names[t],"_conflicts.csv", sep=""))
}


#####
# Use Functions
#####

print.poem()

fw.percents <- import.FW.names()
gg.percents <- import.GG.names()

gg.percents <- reformat.gg(GGtable = gg.percents)
fw.percents <- reformat.fw(FWtable = fw.percents)

check.files.match(FWtable = fw.percents, GGtable = gg.percents)

# Only compare the classifications made by the fw database to the gg classifications, not full tax tables
fw.seq.ids <- import.FW.seq.IDs()
fw.indeces <- find.fw.indeces(TaxonomyTable = fw.percents, SeqIDs = fw.seq.ids)
fw.percents.fw.only <- fw.percents[fw.indeces,]
gg.percents.fw.only <- gg.percents[fw.indeces,]

check.files.match(FWtable = fw.percents.fw.only, GGtable = gg.percents.fw.only)

fw.fw.only <- do.bootstrap.cutoff(TaxonomyTable = fw.percents.fw.only, BootstrapCutoff = taxonomy.bootstrap.cutoff)
cat("\nFinished bootstrap value cutoff on workflow's taxonomy file.\n")
gg.fw.only <- do.bootstrap.cutoff(TaxonomyTable = gg.percents.fw.only, BootstrapCutoff = taxonomy.bootstrap.cutoff)
cat("Finished bootstrap cutoff on comparison taxonomy file.\n\n")

check.files.match(FWtable = fw.fw.only, GGtable = gg.fw.only)

fw.fw.only <- apply(fw.fw.only, 2, remove.parentheses)
gg.fw.only <- apply(gg.fw.only, 2, remove.parentheses)

check.files.match(FWtable = fw.fw.only, GGtable = gg.fw.only)

for (t in 1:5){
  find.conflicting.names(FWtable = fw.fw.only, GGtable = gg.fw.only, FWtable_percents = fw.percents.fw.only, 
                         GGtable_percents = gg.percents.fw.only, TaxaLevel = t)
}


