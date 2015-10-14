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

userprefs <- commandArgs(trailingOnly = TRUE)
fw.plus.gg.tax.file.path <- userprefs[1]
gg.only.tax.file.path <- userprefs[2]
results.folder.path <- userprefs[3]
taxonomy.bootstrap.cutoff <- userprefs[4]


#####
# Define Functions
#####

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
  cat("\n\nDo the sequence ID numbers in each table match?", as.character(order.check), "\n\n")
}

# Remove parentheses of % confidence so that names in each table match exactly
remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

# Apply classification bootstrap value cutoff 
# (this stat is the % of times it got classified into that taxonomy cluster)
do.bootstrap.cutoff <- function(TaxonomyTable, BootstrapCutoff){
  tax <- as.matrix(TaxonomyTable)
  cutoff <- BootstrapCutoff
  
  # define two internal functions:
  pull.out.percent <- function(text){
    text <- sub(pattern = '.*\\(', replacement = '\\(', x = text)
    text <- sub(pattern = '\\(', replacement = '', x = text)
    text <- sub(pattern = '\\)', replacement = '', x = text)
    return(text)
  }
  make.unclassified <- function(text,val){
    if (val == F){
      text <- "unclassified"
    }
    return(text)
  }
  
  # Do calculations
  tax.nums <- apply(tax[,2:ncol(tax)],2,pull.out.percent)
  index <- which(tax.nums == "unclassified" | tax.nums == "unknown")
  tax.nums[index] <- 0
  tax.nums <- apply(X = tax.nums, MARGIN = 2, FUN = as.numeric)
  tax.TF <- tax.nums >= cutoff
  
  # could you do this faster with sweep() ?  I think so...
  for (r in 1:nrow(tax)){
    for (c in 2:ncol(tax)){
      tax[r,c] <- make.unclassified(text = tax[r,c], val = tax.TF[r,c-1])
    }
  }
  
  return(tax)
}


# Find clones misclassified at a given phylogenetic level, t
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
  cat("there are ", length(index), " conflicting names at ", taxa.names[t], " level")
  
  # Compare the conflicting tables in entirety, use the original files with percents still in it
  conflicting <- cbind(gg.percents[index,], fw.percents[index,])
  
  # Check that the files still line up correctly
  check.files.match(FWtable = conflicting[,9:16], GGtable = conflicting[,1:8])
  
  # Export a file with the conflicting rows side by side.
  write.csv(conflicting, file = paste(results.folder.path, "/", taxa.names[t],"_conflicts.csv", sep=""))
}


#####
# Use Functions
#####

fw.percents <- import.FW.names()
gg.percents <- import.GG.names()

gg.percents <- reformat.gg(GGtable = gg.percents)
fw.percents <- reformat.fw(FWtable = fw.percents)

check.files.match(FWtable = fw.percents, GGtable = gg.percents)

fw <- do.bootstrap.cutoff(TaxonomyTable = fw.percents, BootstrapCutoff = taxonomy.bootstrap.cutoff)
cat("\nOHHHHHHH, we're half-way the-ere\n")
gg <- do.bootstrap.cutoff(TaxonomyTable = gg.percents, BootstrapCutoff = taxonomy.bootstrap.cutoff)
cat("OOOOOO-OOOH, livin' on a prayer\n")

check.files.match(FWtable = fw, GGtable = gg)

fw <- apply(fw, 2, remove.parentheses)
gg <- apply(gg, 2, remove.parentheses)

check.files.match(FWtable = fw, GGtable = gg)

for (t in 1:5){
  find.conflicting.names(FWtable = fw, GGtable = gg, FWtable_percents = fw.percents, GGtable_percents = gg.percents, TaxaLevel = t)
}


