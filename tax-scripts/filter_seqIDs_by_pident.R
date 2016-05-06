# RRR 8-25-15 wrote script
# RRR 1/12/16 split script into two scripts- this and calc_full_length_pident.R

# This script generates a \n separated file of OTU (query) sequence IDs.
# The user specifies if they want seqIDs >= cutoff (TRUE), or seqIDs < cutoff (FALSE)
# The output of this script can be fed into the fetch_fastas_with_seqIDs.py python script (step 7).

# The script is called from the command line with the following 4 arguments in order:
# 1st		the path to the reformatted blast file. This is the otus.taxonomy.blast.table.modified file
#       that is like the blast output table with the calculated full length pident included
# 2nd		the path the to file you are creating. This is the list of sequence ID's matching your criteria.
# 3rd		the cutoff percent identity to use. Note the cutoff is applied to the calculated "full length" percent identity
# 4th		want matches? TRUE or FALSE. TRUE: return seqID's >= cutoff, FALSE: return seqID's < cutoff
# Example syntax in terminal:
# Rscript scriptname.R blastfilename outputfilename cutoff# TRUE/FALSE

# The format of the blast file is qseqid, pident, length, qlen, q.align, true.pids, hit.num.best.ids
# The format of the output file is \n separated qseqids


# ####
# Receive arguments from terminal command line
# ####

userprefs <- commandArgs(trailingOnly = TRUE)
blast.file.path <- userprefs[1]
output.file.path <- userprefs[2]
users.cutoff <- as.numeric(userprefs[3])
user.wants.matches <- as.logical(userprefs[4])

# blast.file.path <- "../../take16/otus.custom.blast.table.modified"
# output.file.path <- "../../take16/ids.above.98"
# users.cutoff <- as.numeric(98)
# user.wants.matches <- as.logical(TRUE)

# ####
# Define Functions
# ####

# Import the modified blast table
import.BLAST.data <- function(File){
  blast.file.path <- File
  blast <- read.table(file = blast.file.path, sep = "\t", header = F, stringsAsFactors = F, colClasses = "character")
  colnames(blast) <- c("qseqid","pident","length","qlen","q.align","true.pids","hit.num.best.ids")
  return(blast)
}

# select qseqids that are above or below the true.percent.id cutoff
filter.seqIDs <- function(BlastTable, cutoff.perc.id, want.matches, File){
  blast <- BlastTable
  cutoff <- cutoff.perc.id
  hits <- want.matches
  output.file.path <- File
  
  # had to import everything as character because very long numbers for seqIDs got changed to scientific notation on import
  
  
  
  if (hits == TRUE){
    index <- which(as.numeric(blast$true.pids) >= cutoff)
  }else{
    index <- which(as.numeric(blast$true.pids) < cutoff)
  }
  
  seq.ids <- blast[index,1]
  
  write(x = seq.ids, file = output.file.path, sep="\n")
  
  if (hits == TRUE){
    cat("\nAdded",length(seq.ids),"sequences matching the custom database with >= ", 
        cutoff, "% sequence identity to",output.file.path,
        "\nAssign taxonomy to these sequences using the Small Custom database\n\n")
  }else{
    cat("\nAdded",length(seq.ids),"sequences matching the custom database with < ", 
    cutoff, "% sequence identity to",output.file.path,
    "\nAssign taxonomy to these sequences using the Large General database\n\n")
  }
  
}

# ####
# Use Functions
# ####

blast <- import.BLAST.data(File = blast.file.path)

filter.seqIDs(BlastTable = blast, cutoff.perc.id = users.cutoff, want.matches = user.wants.matches, File = output.file.path)
