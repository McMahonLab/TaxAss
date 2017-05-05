# Generate the assorted stats I list in the paper text:

# Percent Mendota dataset that is cyanobacteria:

file.path.tax <- "../../ME_GG/data/otus.98.85.70.taxonomy"
file.path.abund <- "../../ME_GG/data/otus.abund"

tax <- read.csv(file = file.path.tax, colClasses = "character")
abund <- read.table(file = file.path.abund, header = T, sep = "\t", colClasses = "character")
abund[ ,-1] <- apply(X = as.matrix(abund[ ,-1]), MARGIN = 2, FUN = as.numeric)
index.tax <- order(tax[ ,1])
index.abund <- order(abund[ ,1])
tax <- tax[index.tax, ]
abund <- abund[index.abund, ]
all.equal(abund[ ,1], tax[ ,1])
tot.abund <- rowSums(abund[ ,-1])
phyla <- tax[ ,3]
names(phyla) <- tax[ ,1]
remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}
phyla <- remove.parentheses(x = phyla)
all.equal(names(phyla), names(tot.abund))
phyla.tots <- data.frame(phylum = phyla, abundance = tot.abund)
phyla.tots <- aggregate(x = phyla.tots[ ,2], by = list(phyla.tots[ ,1]), FUN = sum)
index.cyano <- which(phyla.tots[ ,1] == "p__Cyanobacteria")
perc.cyanos <- phyla.tots[index.cyano,2] / sum(phyla.tots[ ,2]) * 100
perc.cyanos
