# RRR 2018-7-20

# dada2 has an internal implementation of the Wang classifier/ RDP classifier / mothur default / TaxAss choice
# BUT, just keep using mothur with TaxAss because formatting is painful.
# this script takes the dada2 output file (seqtab_nochim) and cretes otus.fasta and otus.abund
# it also creates an file otus.count which lists the total reads in each sample.

# command line syntax:

# Rscript reformat_dada2_seqtabs.R seqtab_nochim.rds otus.fasta otus.abund otus.count

# ---- Accept Arguments from Terminal Command Line ----

userprefs <- commandArgs(trailingOnly = TRUE)
path.to.seqtab <- userprefs[1] 
fasta.output <- userprefs[2] 
abund.output <- userprefs[3]
count.output <- userprefs[4]

# cat("\nfuck you forgot to comment out the file paths in reformat_dada2_seqtabs.R!\n")
# path.to.seqtab <- "/Users/athena/Desktop/dada2-meV34/dada2/seqtab_nochim.rds"
# fasta.output <- "/Users/athena/Desktop/dada2-meV34/taxass/otus.fasta"
# abund.output <- "/Users/athena/Desktop/dada2-meV34/taxass/otus.abund"
# count.output <- "/Users/athena/Desktop/dada2-meV34/taxass/otus.count"

# ---- define functions ----

import.dada2.file <- function(path){
  is.rds <- grepl(pattern = "*.rds", x = path)
  if (is.rds){
    seqtab_nochim <- readRDS(file = path)
    return(seqtab_nochim)
  }else{
    cat("input must be rds file ending in \".rds\"\nSave it in this format at the end of dada2 pipeline using\nsaveRDS(object = seqtab_nochim, file = \"filename\")")
  }  
}

make.fasta.file <- function(dadatable, fasta.path){
  fasta.seqs <- colnames(dadatable)
  otu.names <- paste("otu", 1:length(fasta.seqs), sep = "_")
  fasta.names <- paste(">", otu.names, sep = "")
  fasta.file <- paste(fasta.names, fasta.seqs, sep = "\n")
  write.table(x = fasta.file, file = fasta.path, quote = F, row.names = F, sep = "\n", col.names = F)
  cat("Made file: ", fasta.path, "\n")
  return(otu.names)
}

find.zero.samples <- function(tot.reads, otu.table){
  read.stats <- boxplot(x = tot.reads, plot = F)
  if(length(read.stats$out) < 1){
    cat("No samples have an outlier number of total reads. Your samples have ", 
        round(mean(tot.reads), 0),"+/-", round(sd(tot.reads), 0), "(+/-", round(sd(tot.reads)/mean(tot.reads) * 100), "%) total reads.\n")
  }else{
    cat("Your samples have ",
        round(mean(tot.reads), 0),"+/-", round(sd(tot.reads), 0), "(+/-", round(sd(tot.reads)/mean(tot.reads) * 100), "%) total reads.\n")
    cat("\nThese samples have outlier read counts (You might consider removing them later):\n", paste(names(read.stats$out), " : ", read.stats$out, sep = "", "\n"))
  }
  index <- which(tot.reads == 0)
  if (length(index) > 0){
    cat("These samples with ZERO reads are being removed now to avoid errors in TaxAss:\n",
        names(tot.reads)[index], "\n")
  }
  return(index)
}

convert.to.rel.abund <- function(OTUs){
  sample.totals <- rowSums(OTUs)
  # vectors are applied to matrices by stepping down rows in a column
  norm.otus <- OTUs / sample.totals * 100
  return(norm.otus)
}

reformat.for.taxass <- function(OTUs){
  # rows = OTUs, cols = samples, col 1 = seqIDs
  otu.table <- t(OTUs)
  otu.table <- cbind(row.names(otu.table), otu.table)
  colnames(otu.table)[1] <- "seqID"
  return(otu.table)
}

make.abund.file <- function(OTUs, FilePath){
  write.table(x = OTUs, file = FilePath, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("Made file: ", FilePath, "\n")
}

make.count.file <- function(Count, FilePath){
  Count <- cbind(names(Count), Count)
  colnames(Count) <- c("Sample.Name", "Total.Reads")
  write.table(x = Count, file = FilePath, quote = F, sep = "\t", row.names = F)
  cat("Made file: ", FilePath, "\n")
}

# ---- use functions ----

seqtab_nochim <- import.dada2.file(path = path.to.seqtab)

seqIDs <- make.fasta.file(dadatable = seqtab_nochim, fasta.path = fasta.output)

colnames(seqtab_nochim) <- seqIDs

read.count <- rowSums(seqtab_nochim)

index <- find.zero.samples(tot.reads = read.count, otu.table = seqtab_nochim)
if(length(index) > 0){
  read.count <- read.count[-index]
  seqtab_nochim <- seqtab_nochim[-index, ]
}

seqtab_nochim <- convert.to.rel.abund(OTUs = seqtab_nochim)

seqtab_nochim <- reformat.for.taxass(OTUs = seqtab_nochim)

make.abund.file(OTUs = seqtab_nochim, FilePath = abund.output)

make.count.file(Count = read.count, FilePath = count.output)




