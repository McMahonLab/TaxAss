# RRR 6/3/20
# make a smarter way of unaligning that doesn't touch the sequence IDs
# Remove .
# Remove -
# IF SPECIFIED Replace U with T (only used in creating FreshTrain releases)
# IF SPECIFIED remove white space and everything after it from the seqID line

# ---- input ----

userprefs <- commandArgs(trailingOnly = TRUE)
aligned.fasta <- userprefs[1]
created.fasta <- userprefs[2]
format.ids <- userprefs[3]
replace.U <- userprefs[4]

# # Manual Troubleshooting  
# cat("\n\nForgot to comment out file paths!!!\n")
# aligned.fasta <- userprefs[1]
# created.fasta <- userprefs[2]
# format.ids <- userprefs[3]
# replace.U <- userprefs[4]


if (is.na(userprefs[3])){
  format.ids <- FALSE
}else{
  format.ids <- TRUE
}
if (is.na(userprefs[4])){
  replace.U <- FALSE
}else{
  replace.U <- TRUE
}

# ---- do things ----

fasta <- scan(file = aligned.fasta, what = "character", sep = "\n", quote = NULL)

ids <- seq.int(from = 1, to = length(fasta) - 1, by = 2)

seqs <- seq.int(from = 2, to = length(fasta), by = 2)

fasta[seqs] <- gsub(pattern = "-", replacement = "", x = fasta[seqs])
cat("Removed alignment dashes\n")

fasta[seqs] <- gsub(pattern = "\\.", replacement = "", x = fasta[seqs])
cat("Removed alignment dots\n")

if (replace.U){
  fasta[seqs] <- gsub(pattern = "U", replacement = "T", x = fasta[seqs], ignore.case = TRUE)
  cat("Replaced U and u with T.\n")
}

if (format.ids){
  fasta[ids] <- gsub(pattern = "[[:blank:]].*$", replacement = "", x = fasta[ids])
  cat("Deleted any white space and anything after it from the seqIDs.\n")
}


# ---- export ----

write(x = fasta, file = created.fasta, ncolumns = 1, append = F)
cat("Made file: ", created.fasta, "\n")

