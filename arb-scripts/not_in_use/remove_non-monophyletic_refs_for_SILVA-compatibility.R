# RRR 2018-4-30
# This script removes references from the FreshTrain that we manually looked at first.
# These references end up conflicted in SILVA because the FreshTrain gives them the same lineage
# name while silva has 2 different order names. To fix the problem, we remove half the lineage 
# from the freshtrain so that it only has one order name. This is a simple approach based on
# that fact that there are very few references for the lineages in question so they are not
# well enough flushed out to justify keeping them and changing them in the FreshTrain.
# I used the Database_Improvement workflow and classified FreshTrain against silva. Then I
# manually looked at the unique conflicts in excel and then in the taxonomy files to identify
# these sequences to remove.

# ---- define paths ----

userprefs <- commandArgs(trailingOnly = TRUE)

start.fasta.file <- userprefs[1]
start.taxonomy.file <- userprefs[2]
end.fasta.file <- userprefs[3]
end.taxonomy.file <- userprefs[4]
silva.version <- userprefs[5]

cat("\n\ncrap re-comment file paths\n\n")
start.fasta.file <- "../FreshTrain-files/FreshTrain25Jan2018Greengenes13_5/FreshTrain25Jan2018Greengenes13_5.fasta"
start.taxonomy.file <- "../FreshTrain-files/FreshTrain25Jan2018Greengenes13_5/FreshTrain25Jan2018Greengenes13_5.taxonomy"
end.fasta.file <- "~/Desktop/test.fasta"
end.taxonomy.file <- "~/Desktop/test.taxonomy"
silva.version <- "v132"


# ---- functions ----

import.mothur.formatted.tax <- function(filepath){
  # seqID delimited from taxonomy by tab, taxa levels delimited by semicolon
  numcol <- max(count.fields(filepath, sep=";"))
  silva <- read.table(file = filepath, sep=";", fill=T, colClasses = "character", col.names=1:numcol)
  silva <- as.matrix(silva)
  kingdom <- gsub(pattern = ".*\t", replacement = "", x = silva[ ,1])
  seqid <- gsub(pattern = "\t.*$", replacement = "", x = silva[ ,1])
  silva <- silva[ ,-1]
  silva <- cbind(seqid, kingdom, silva)
  if(ncol(silva) == 8){
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported with a species column, the unique species names are: ", unique(silva[ ,8]),"\n")
  }else if(ncol(silva) == 7){
    silva <- cbind(silva, "unnamed")
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported without a species column, added one filled with \"unnamed\"\n")
  }else if(ncol(silva) == 9){
    extracol <- unique(silva[ ,9])
    silva <- silva[ ,1:8]
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported with an extra 9th column containing:", extracol, ". This column was removed.\n")
  }else{
    cat("something is wrong with the database import, it has ", ncol(silva), "columns in it.\n")
  }
  return(silva)
}

revert.to.mothur.format <- function(silva){
  seqid.kingdom <- paste(silva[ ,1], silva[ ,2], sep = "\t")
  silva <- silva[ ,-(1:2)]
  # also add a blank column to get semicolons at the end of the line also
  silva <- cbind(seqid.kingdom, silva, "")
  return(silva)
}

# ---- go ----

# this is so much more work than it's worth. doing it manually. 
