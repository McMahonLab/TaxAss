# 8-27-15 RRR

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

#####
# Receive arguments from terminal command line
#####

userprefs <- commandArgs(trailingOnly = TRUE)
fw.plus.gg.tax.file.path <- userprefs[1]
gg.only.tax.file.path <- userprefs[2]
fw.seq.ids.file.path <- userprefs[3]
results.folder.path <- userprefs[4]
blast.pident.cutoff <- userprefs[5]
taxonomy.bootstrap.cutoff <- userprefs[6]
check.or.final <- userprefs[7]
if (length(userprefs) < 7){check.or.final <- "check"}


# fw.plus.gg.tax.file.path <- "../../take7/otus.98.taxonomy"
# gg.only.tax.file.path <- "../../take7/otus.general.taxonomy"
# fw.seq.ids.file.path <- "../../take7/ids.above.98"
# results.folder.path <- "../../take7/conflicts_98/"
# blast.pident.cutoff <- 98
# taxonomy.bootstrap.cutoff <- 70
# check.or.final <- "final"

#####
# Define Functions for Import and Formatting
#####

# Entertain user with a poem while they wait:
print.poem <- function(){
  cat("\nAnd the Days Are Not Full Enough\nby Ezra Pound\n\nAnd the days are not full enough\nAnd the nights are not full enough\nAnd life slips by like a field mouse\n\tNot shaking the grass.\n\n")
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
  
  # convert seqIDs to characters in case they are numeric b/c as.matrix on numbers adds spaces but as.character doesn't
  fw[,1] <- as.character(fw[,1])
  
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
  
  # convert seqIDs to characters in case they are numeric b/c as.matrix() on numbers adds spaces but as.character() doesn't
  gg[,1] <- as.character(gg[,1])
  
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

# Import the freshwater sequence IDs as determined by the BLAST cutoff in workflow step 4
import.FW.seq.IDs <- function(){
  fw.seqs <- scan(file = fw.seq.ids.file.path)
  fw.seqs <- as.character(fw.seqs)
  return(fw.seqs)
}

# Remove parentheses of % confidence so that names in each table match exactly
remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}


#####
# Define Functions for Data Analysis
#####

# Find index of sequences that were classified by the freshwater database. returns a vector of indeces.
find.fw.indeces <- function(TaxonomyTable, SeqIDs){
  tax <- TaxonomyTable
  ids <- SeqIDs
  
  # when seqIDs are numbers written as characters, sometimes one file could have whitespace placeholders for
  # for the shorter characters.  But can't just change to numeric, in case some files have non-number seqIDs
  
  index <- NULL
  for (e in 1:length(ids)){
    index <- c(index, which(tax[,1] == ids[e]))
  }
  return(index)
}

# makes all the unclassified/unknown/any other word for it names be uniformly called "unclassified"
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
  # text = k__Bacteria(100) returns (as character) 100, But text = "unclassified" returns (as character) "unclassified"
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
  # (by bootstrap cutoff I mean the stat that's the % of times it got classified into that taxonomy cluster)
do.bootstrap.cutoff <- function(TaxonomyTable, BootstrapCutoff){
  tax <- as.matrix(TaxonomyTable)
  cutoff <- BootstrapCutoff
  
  # create a matrix of bootstrap numbers and then a T/F matrix
  tax.nums <- apply(tax[,2:ncol(tax)],2,pull.out.percent)
  index <- which(tax.nums == "unclassified")
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  tax.TF <- tax.nums >= cutoff
  
  # make all names in the taxonomy table unclassified if they're below the bootstrap cutoff
  for (r in 1:nrow(tax)){
    for (c in 2:ncol(tax)){
      tax[r,c] <- make.unclassified(text = tax[r,c], val = tax.TF[r,c-1])
    }
  }
  
  return(tax)
}

# Find seqs misclassified at a given phylogenetic level, t
find.conflicting.names <- function(FWtable, GGtable, GGtable_percents, TaxaLevel, tracker){
  fw <- FWtable
  gg <- GGtable
  gg.percents <- GGtable_percents
  t <- TaxaLevel
  num.mismatches <- tracker
  
  taxa.names <- c("kingdom","phylum","class","order","lineage","clade","tribe")
  
  # compare names in column t+1, because first columns are seqID and t=1 is kingdom, t=7 is tribe
  # ignore names that say unclassified
  index <- which(gg[,t+1] != fw[,t+1] & gg[,t+1] != "unclassified" & fw[,t+1] != "unclassified")
  cat("there are ", length(index), " conflicting names at ", taxa.names[t], " level\n")
  num.mismatches[t] <- length(index)
  
  # Compare the conflicting tables in entirety, use the original files with percents still in it
  conflicting <- cbind(gg.percents[index,,drop=F], fw[index,,drop=F])
  
  # Check that the files still line up correctly
  check.files.match(FWtable = conflicting[,9:16,drop=F], GGtable = conflicting[,1:8,drop=F])
  
  # Export a file with the conflicting rows side by side.
  write.csv(conflicting, file = paste(results.folder.path, "/", t, "_", taxa.names[t],"_conflicts.csv", sep=""), row.names = FALSE)
  
  #Track the number of mismatches at each level
  return(num.mismatches)
}

# Set up a summary vector to fill
create.summary.vector <- function(){
  num.mismatches <- vector(mode = "numeric", length = 5)
  names(num.mismatches) <- c("kingdom","phylum","class","order","lineage")
  return(num.mismatches)
}

# Format and export the summary vector
export.summary.stats <- function(SummaryVector){
  num.mismatches <- SummaryVector
  num.mismatches <- c(num.mismatches,"numFWseqs" = nrow(fw.fw.only))
  num.mismatches <- c(num.mismatches, "numALLseqs" = nrow(fw.percents))
  num.mismatches <- data.frame("TaxaLevel" = names(num.mismatches),"NumConflicts" = num.mismatches, row.names = NULL)
  write.csv(num.mismatches, file = paste(results.folder.path, "/", "conflicts_summary.csv", sep=""), row.names = FALSE)
}

# Check how high the freshwater bootstrap values end up given your cutoff.
view.bootstraps <- function(TaxonomyTable){
  tax <- TaxonomyTable
  
  # create a matrix of bootstrap numbers, copy from do.bootstrap.cutoff()
  tax.nums <- apply(tax[,2:ncol(tax)],2,pull.out.percent)
  index <- which(tax.nums == "unclassified")
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  
  return(tax.nums)
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

# Generate a final taxonomy file:
if (check.or.final == "final" | check.or.final == "Final" | check.or.final == "FINAL"){
  final.taxonomy <- do.bootstrap.cutoff(TaxonomyTable = fw.percents, BootstrapCutoff = taxonomy.bootstrap.cutoff)
  colnames(final.taxonomy) <- c("seqID","kingdom","phylum","class","order","lineage","clade","tribe")
  write.table(x = final.taxonomy, file = paste("otus.", blast.pident.cutoff, ".", taxonomy.bootstrap.cutoff, ".taxonomy", sep = ""), 
              sep = ";", row.names = FALSE)
}

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

# Generate the files comparing classifications made by fw to those of gg
# Generate a summary file listing the total number of classification disagreements at each level
  # Files written in find.conflicting.names() loop: the "TaxaLevel_conflicts.csv" that puts taxonomy tables side by side
  # File written afer loop: the "conflicts_summary.csv" that lists how many conflicts were at each level, and how many seqs were classified by FW
num.mismatches <- create.summary.vector()
for (t in 1:5){
  num.mismatches <- find.conflicting.names(FWtable = fw.fw.only, GGtable = gg.fw.only,
                                           GGtable_percents = gg.percents.fw.only, 
                                           TaxaLevel = t, tracker = num.mismatches)
}
export.summary.stats(SummaryVector = num.mismatches)

# Generate a file of the fw-assigned taxonomies, and a matching table of just their bootstrap values
  # File written: the "fw_classified_bootstraps.csv" that lists the bootstrap of all sequences classified by freshwater
  # File written: the "fw_classified_taxonomies.csv" that lists the taxonomy of all sequences classified by freshwater
fw.bootstraps <- view.bootstraps(TaxonomyTable = fw.percents.fw.only)
write.csv(fw.bootstraps, file = paste(results.folder.path, "/", "fw_classified_bootstraps.csv", sep=""), row.names = FALSE)
write.csv(fw.percents.fw.only, file = paste(results.folder.path, "/", "fw_classified_taxonomies.csv", sep=""), row.names = FALSE)

# Generate a file of the gg taxonomies for the fw-assigned sequences, and a matching table of gg bootstraps
  # File written: the "gg_classified_bootstraps.csv" that lists the bootstrap of all sequences classified by freshwater
  # File written: the "gg_classified_taxonomies.csv" that lists the taxonomy of all sequences classified by freshwater
gg.bootstraps <- view.bootstraps(TaxonomyTable = gg.percents.fw.only)
write.csv(gg.bootstraps, file = paste(results.folder.path, "/", "gg_classified_bootstraps.csv", sep=""), row.names = FALSE)
write.csv(gg.percents.fw.only, file = paste(results.folder.path, "/", "gg_classified_taxonomies.csv", sep=""), row.names = FALSE)
