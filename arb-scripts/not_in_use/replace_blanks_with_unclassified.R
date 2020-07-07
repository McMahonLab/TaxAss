# RRR 2/11/16
# reformat the taxonomy training set so that there are no empty spots.
# replace empty levels with the word "unnamed"
# empty spots at genus level yield incorrect p-values in RDP classifier

# syntax from command line:

# Rscript replace_blanks_with_unclassified.R inputfile outputfile

# ---- Define Input ----

# # cat("\nfuck you forgot to comment out the file paths!!!\n")
# file.with.blanks.semicolon.delim <- "~/Desktop/TaxonomyTrainingSets/BLASTing/SILVA_downloads/12-7-16_format_mothur_silva/semicolons.tax"
# reformatted.file.semicolon.delim <- "~/Desktop/TaxonomyTrainingSets/BLASTing/SILVA_downloads/12-7-16_format_mothur_silva/allranks.tax"

userprefs <- commandArgs(trailingOnly = TRUE)
file.with.blanks.semicolon.delim <- userprefs[1] 
reformatted.file.semicolon.delim <- userprefs[2]

# ---- Go ----

# import taxonomy file:
numcol <- max(count.fields(file.with.blanks.semicolon.delim, sep=";", quote=""))
taxonomy.unformatted <- read.table(file = file.with.blanks.semicolon.delim, header = FALSE, sep = ";", quote = "", 
                                   col.names = 1:numcol, fill = TRUE, stringsAsFactors = FALSE)

# format it:
taxonomy <- taxonomy.unformatted[ ,1:8]
taxonomy <- as.matrix(taxonomy)
index <- which(taxonomy[ ,-1] == "")
taxonomy[ ,-1][index] <- "unnamed"
index <- which(is.na(taxonomy[ ,-1]))
taxonomy[ ,-1][index] <- "unnamed"

# export new version- this is semicolon delim.

write.table(x = taxonomy, file = reformatted.file.semicolon.delim, quote = FALSE, sep = ";", row.names = FALSE, col.names = FALSE)
