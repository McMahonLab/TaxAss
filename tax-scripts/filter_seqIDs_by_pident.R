# RRR 8-25-15

# This script generates a \n separated file of query sequence IDs that match a criteria.

# The script is called from the command line with the following 5 arguments in order:
# 1st		the path to the r script.  The script is called 8-24-15_Filter_BLAST_for_tax_assignment.R
# 2nd		the path to the formatted blast file. This is the otus.taxonomy.blast.table file from step 3
# 3rd		the path the to file you are creating. This is the list of sequence ID's matching your criteria.
# 4th   the path to the other file you're creating.  This is used to assess the appropriateness of BLAST settings.
# 5th		the cutoff percent identity to use. Note this is the calculated "true" full length percent identity
# 6th		want matches? TRUE or FALSE. TRUE: return seqID's >= cutoff, FALSE: return seqID's < cutoff
# Example syntax in terminal:
# Rscript scriptname.R blastfilename outputfilename statsfilename cutoff# TRUE/FALSE

# The format of the blast file is -outfmt "6 qseqid pident length qlen qstart qend" 
# The format of the output file is \n separated qseqids
# The format of the stats file is very similar to the blast.table file, but the headings are 
# qseqid, pident, length, qlen, q.align, true.pids, hit.num.best.ids

# This script caclulates the full-length query's percent ID (BLAST may leave gaps at the ends).
# Then it finds all the query IDs >= OR < the specified (full-length) percent ID cutoff, depending which you ask for.
# The output of this script can be fed into the fetch_fastas_with_seqIDs.py python script.
# The stats file output can be fed into the plot_blast_hit_stats.R R script.

#####
# Receive arguments from terminal command line
#####

userprefs <- commandArgs(trailingOnly = TRUE)
blast.file.path <- userprefs[1]
hit.stats.path <- userprefs[2]
output.file.path <- userprefs[3]
users.cutoff <- as.numeric(userprefs[4])
user.wants.matches <- as.logical(userprefs[5])

# blast.file.path <- "../../take6/otus.custom.blast.table"
# output.file.path <- "../../take6/ids.above.98"
# hit.stats.path <- "../../take6/hit.stats.98"
# users.cutoff <- 98
# user.wants.matches <- TRUE

#####
# Define Functions
#####

# Import the blast output
import.BLAST.data <- function(File){
  blast.file.path <- File
  blast <- read.table(file = blast.file.path, sep = "\t", stringsAsFactors = F)
  colnames(blast) <- c("qseqid","pident","length","qlen","qstart","qend")
  return(blast)
}

# Format the blast data for downstream analyses
format.BLAST.data <- function(BlastTable){
  blast <- BlastTable
  
  # format numbers as numberic instead of integer
  blast[,2:6] <- apply(blast[,2:6], 2, as.numeric)
  
  # add a query alignment length column 
  q.align <- blast$qend - blast$qstart + 1   # add 1 b/c start is 1 instead of 0.  for 4 aligned pairs, start=1, end=4, length=4-1+1=4
  blast <- cbind(blast,q.align)              
  
  # remove the now unnecessary qend/qstart columns
  blast <- blast[,-(5:6)]
  
  # are most of the alignments full length?
#   cat("\nalignment length \tmean:", mean(blast$length), "\t\tmin:", min(blast$length), "\t\tmax:",max(blast$length),
#       "\nthe query length \tmean:", mean(blast$qlen), "\t\tmin:", min(blast$qlen), "\t\tmax:",max(blast$qlen),"\n")
  
  return(blast)
}

# Calculate pident for full length queries (account for any gaps at ends of HSPs) for the whole table
calc.full.pIDs <- function(BlastTable){
  blast <- BlastTable
  
  # define function to apply to each row
  pid.fun <- function(NumericBlastMatrixRow){
    b <- NumericBlastMatrixRow
    pid <- (b[1] * b[2]) / (b[2] - b[4] + b[3]) # (pident * length) / (length - q.align + qlen)
    return(pid)
  }
  
  # create a numeric blast matrix to use the apply funciton
  b <- blast[,-1]
  b <- as.matrix(b)
  
  # find the "true" pid for the full length sequence of each blast hit
  true.pids <- apply(b, 1, pid.fun)
  
  # add this to the blast table
  blast <- cbind(blast, true.pids)
  
  return(blast)
}

# choose the hit with the best "full length" pident- report which # hit it is.
choose.best.hit <- function(BlastTable, OutputFile){
  blast <- BlastTable
  hit.stats.path <- OutputFile
  
  # get a list of all the unique id names, and set up blank vectors to record which are best
  unique.qids <- unique(blast$qseqid)
  index.best.ids <- NULL
  hit.num.best.ids <- NULL
  
  # for each unique query id name, report the index and hit number of the best "full length" hit
  for (q in 1:length(unique.qids)){
    index.same.ids <- which(blast$qseqid == unique.qids[q])
    
    # start with the first hit as the default best one
    index.best <- index.same.ids[1]
    hit.num <- 1
    
    # if there are multiple hits, replace if another is better
    if (length(index.same.ids) > 1){
      for (i in 2:length(index.same.ids)){
        if (blast$true.pids[index.same.ids][i] > blast$true.pids[index.best]){
          index.best <- index.same.ids[i]
          hit.num <- i
        }
      }
    }
    
    # record the blast table index of the best hit
    index.best.ids <- c(index.best.ids, index.best)
    # record which number hit that was.
    hit.num.best.ids <- c(hit.num.best.ids, hit.num)
  }
  
  # change blast table to only contain the best hits
  blast <- blast[index.best.ids,]
  # add a column with the hit number
  blast <- cbind(blast, hit.num.best.ids)
  
  # save this table to feed into a plotting/analysis script
  write.csv(x = blast, file = hit.stats.path, row.names = F)
  
  return(blast)
}

# select qseqids that are above or below the true.percent.id cutoff
filter.seqIDs <- function(BlastTable, cutoff.perc.id, want.matches, File){
  blast <- BlastTable
  cutoff <- cutoff.perc.id
  hits <- want.matches
  output.file.path <- File
  
  if (hits == TRUE){
    index <- which(blast$true.pids >= cutoff)
  }else{
    index <- which(blast$true.pids < cutoff)
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

# print poem to enrich the user experience
print.poem <- function(){
  # Tomas Tranströmer won the nobel prizer in literature in 2011
  # I like this poem because of it's central metaphor of a snow covered lake being like a blank page
  cat("\nFrom March 1979\n",
      "   by Tomas Tranströmer\n",
      "   translated to english by _______\n",
      "\nTired of all who come with words, words but no language\n",
      "I went to the snow-covered island.\n",
      "The wild does not have words.\n",
      "The unwritten pages spread out on all sides!\n",
      "I come upon the tracks of roe deer in the snow.\n",
      "Language but no words.\n\n", sep = "")
}

#####
# Use Functions
#####

print.poem()

blast <- import.BLAST.data(File = blast.file.path)

blast <- format.BLAST.data(BlastTable = blast)

blast <- calc.full.pIDs(BlastTable = blast)

blast <- choose.best.hit(BlastTable = blast, OutputFile = hit.stats.path)

filter.seqIDs(BlastTable = blast, cutoff.perc.id = users.cutoff, want.matches = user.wants.matches, File = output.file.path)
