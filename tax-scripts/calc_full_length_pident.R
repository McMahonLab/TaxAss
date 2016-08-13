# RRR 8-25-15 wrote
# RRR 1/12/16 split into two separate scripts

# This script caclulates the full-length query's percent ID (BLAST may leave gaps at the ends).
# This script generates a new blast output table that has this "full length" pident in it.
# This script calculates which blast hit number was best after correction to "full length"
# The generated table feeds into two scripts: filter_seqIDs_by_pident.R and plot_blast_hit_stats.R 

# The script is called from the command line with the following 2 arguments in order:
# 1st		the path to the formatted blast file. This is the otus.taxonomy.blast.table file from step 3
# 2nd		the path the to file you are creating. This is the otus.taxonomy.blast.table.modified file
#       that contains the blast table with full length pident calculated and which hit was best specified.
#       This is used later both to filter the seqIDs based on their full length pident, in filter_seqIDs_by_pident.R
#       and to examine the appropriateness of the blast setting in plot_blast_hit_stats.R
# Example syntax in terminal:
# Rscript scriptname.R blastfilename outputfilename

# The format of the blast table is -outfmt "6 qseqid pident length qlen qstart qend" 
# The format of the output file is very similar to the blast.table file: tab delimited and no column names, 
# but the columns are qseqid, pident, length, qlen, q.align, true.pids, hit.num.best.ids


#####
# Receive arguments from terminal command line
#####

userprefs <- commandArgs(trailingOnly = TRUE)
blast.file.path <- userprefs[1]
hit.stats.path <- userprefs[2]

# blast.file.path <- "../../take19/otus.custom.blast.table"
# hit.stats.path <- "../../take16/otus.custom.blast.table.modified"

#####
# Define Functions
#####

# Import the blast output
import.BLAST.data <- function(File){
  blast.file.path <- File
  blast <- read.table(file = blast.file.path, sep = "\t", stringsAsFactors = F, colClasses = "character")
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
  write.table(x = blast, file = hit.stats.path, row.names = F, col.names = F, sep = "\t", quote =  FALSE)
  
  return(blast)
}

# print poem to enrich the user experience
print.poem <- function(){
  # Tomas Tranströmer won the nobel prizer in literature in 2011
  # I like this poem because of it's central metaphor of a snow covered lake being like a blank page
  cat("\nFrom March 1979\n",
      "   by Tomas Tranströmer\n",
      "   translated to english by Robin Fulton\n",
      "\nWeary of all who come with words, words but no language\n",
      "I make my way to the snow-covered island.\n",
      "The untamed has no words.\n",
      "The unwritten pages spread out on every side!\n",
      "I come upon the tracks of deer in the snow.\n",
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


