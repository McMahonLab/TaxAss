# RRR 2018-7-20

# dada2 has an internal implementation of the Wang classifier/ RDP classifier / mothur default / TaxAss choice
# BUT, just keep using mothur with TaxAss because formatting is painful.
# this script takes the dada2 output file (seqtab_nochim) and cretes otus.fasta and otus.abund

# command line syntax:

# Rscript reformat.dada2.seqtabs.R seqtab_nochim.rds otus.fasta otus.abund

# ---- Accept Arguments from Terminal Command Line ----

userprefs <- commandArgs(trailingOnly = TRUE)
path.to.seqtab <- userprefs[1] 
fasta.output <- userprefs[2] 
abund.output <- userprefs[3]

# cat("fuck you forgot to comment out the file paths in reformat_dada2_seqtabs.R!")
# path.to.seqtab <- "~/Desktop/dada2-meV45/seqtab_nochim.rds"
# fasta.output <- "~/Desktop/dada2-meV45/taxass/otus.fasta"
# abund.output <- "~/Desktop/dada2-meV45/taxass/otus.abund"

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

# ---- use functions ----

seqtab_nochim <- import.dada2.file(path = path.to.seqtab)

seqIDs <- make.fasta.file(dadatable = seqtab_nochim, fasta.path = fasta.output)

colnames(seqtab_nochim) <- seqIDs

seqtab_nochim <- convert.to.rel.abund(OTUs = seqtab_nochim)

seqtab_nochim <- reformat.for.taxass(OTUs = seqtab_nochim)

make.abund.file(OTUs = seqtab_nochim, FilePath = abund.output)
