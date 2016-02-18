# RRR 2/12/16

# plot Actinobacteria to show importance of lower levels.

# import taxonomy file

otu.table.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/otus.abund"
taxonomy.table.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take13/otus.98.85.70.taxonomy"

otus <- read.table(file = otu.table.path, sep = "\t", header = TRUE)
tax <- read.table(file = taxonomy.table.path, sep = ",", header = TRUE)

remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

tax <- apply(tax, 2, remove.parentheses)

otus <- otus[ ,-97]
seqIDs <- otus[ ,1]
otus <- otus[ ,-1]
otus <- apply(otus, 2, as.numeric)
row.names(otus) <- seqIDs

# put dates in order!
samples <- colnames(otus)
samples <- substr(samples, start=3, stop=9)
samples <- as.Date(samples, format="%d%b%y")
colnames(otus) <- as.character(format(samples, "%m-%d-%y") ) 
index <- order(samples)
otus <- otus[ ,index]

# pull out all actinobacteria from taxonomy, find their seqIDs, find them in the OTU table, group them together, plot reads/time
index <- which(tax[ ,3] == "p__Actinobacteria")
p.actino <- tax[index, ]
index <- which(p.actino[ ,4] == "c__Actinobacteria")
c.actino <- p.actino[index, ]
index <- which(c.actino[ ,5] == "o__Actinomycetales")
o.actino <- c.actino[index, ]
index <- which(o.actino[ ,6] == "acI")
l.actino <- o.actino[index, ]
index <- which(l.actino[ ,7] == "acI-A")
A.actino <- l.actino[index, ]
index <- which(l.actino[ ,7] == "acI-B")
B.actino <- l.actino[index, ]
index <- which(l.actino[ ,7] == "acI-C")
C.actino <- l.actino[index, ]


