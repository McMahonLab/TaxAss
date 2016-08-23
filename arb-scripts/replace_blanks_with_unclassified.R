# RRR 2/11/16
# reformat the taxonomy training set so that there are no empty spots.
# replace empty levels with the word "unclassified"
# empty spots at genus level yield incorrect p-values in RDP classifier

file.with.blanks.semicolon.delim <- "~/Desktop/TaxonomyTrainingSets/BLASTing/Primers/v4_all_FreshTrain/FW_al"
reformatted.file.semicolon.delim <- "~/Desktop/TaxonomyTrainingSets/BLASTing/Primers/v4_all_FreshTrain/cleaned4"


# import taxonomy file:
numcol <- max(count.fields(file.with.blanks.semicolon.delim, sep=";", quote=""))
taxonomy.unformatted <- read.table(file = file.with.blanks.semicolon.delim, header = FALSE, sep = ";", quote = "", 
                                   col.names = 1:numcol, fill = TRUE, stringsAsFactors = FALSE)

# format it:
taxonomy <- taxonomy.unformatted[ ,1:8]
taxonomy <- as.matrix(taxonomy)
index <- which(taxonomy[ ,-1] == "")
taxonomy[ ,-1][index] <- "unclassified"

# export new version- this is semicolon delim.

write.table(x = taxonomy, file = reformatted.file.semicolon.delim, quote = FALSE, sep = ";", row.names = FALSE, col.names = FALSE)
