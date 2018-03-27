# RRR 3-27-18 ----

# SILVA database has extremely inconsistent names for unclassified entries.
# This is a problem.

# Why?
# 1. Names need to be unique for downstream analysis.
# Example: 
  # phylum-A;nonunique-class;order-A;...
  # phylum-B;nonunique-class;order-B;...
# Analysis scripts would group order-A and order-B as belonging to the same class, 
# when really those are different classes because they belong to different phylums.
# 2. TaxAss needs names to be unique to tally-up unclassified entries.
# Example: this is done for Figure 2 in the manuscript, in step 15.5.a.

# How does this script fix it?

# 1. This script takes all things meaning "we don't know the name" and changes them to a few versions of "we don't know". 
# Example: unknown
  # phylum-A;Unknown_Class;order-A;... becomes
  # phylum-A;unknown;order-A;...
# Example: unnamed
  # phylum-A;phylum-A_cl;order-A;... becomes
  # phylum-A;unnamed;order-A;...
# Example: uncultured
  # phylum-A;uncultured;order-A;... becomes
  # phylum-A;uncultured;order-A;...
# Example: blank
  # order-A;family-A;genus-A;; becomes
  # order-A;family-A;genus-A;blank;
# Example: unclassified
  # there are no unclassifieds in silva right now, but anything with unclassified in it will be changed to simply "unclassified"

# 2. This script makes all the unclassified entries unique based on their upper-level taxonomy.
# Example:
  # phylum-A;uncultured;order-A;... becomes
  # phylum-A;uncultured.phylum-A;order-A;...
# Example:
  # phylum-A;uncultured;uncultured;... becomes
# phylum-A;uncultured.phylum-A;uncultured.phylum-A;...

# Why does that fix it?
# 1. Each daughter name will end up with unique parent names.
# 2. The nomenclature is consistent so TaxAss can still tally up the unclassifieds.
# Just call anything unclassified that begins with (blank)|(uncultured)|(unnamed)|(unknown)|(unclassified)  

# What doesn't it fix?
# 1. I don't know if "uncultured" and "unknown" should actually be distinguished or not. 
# Now they are being distinguished, but possibly that's incorrect.
# 2. There are still daughter names below unknown names. 
# This somewhat breaks logic, but with the parent names made unique hopefully it will not break scripts.


# ---- file paths/commandline input ----

mothur.formatted.silva <- "../../StartFiles-Databases/Silva.nr_v132/silva_132_noperiod.taxonomy"
mothur.formatted.silva <- "../../StartFiles-Databases/withGG/general.taxonomy"

new.silva.file <- "~/Desktop/cleansilva.taxonomy"

# ---- functions ----

import.mothur.formatted.tax <- function(filepath){
  # seqID delimited from taxonomy by tab, taxa levels delimited by semicolon
  numcol <- max(count.fields(filepath, sep=";"))
  silva <- read.table(file = filepath, sep=";", fill=T, colClasses = "character", col.names=1:numcol)
  silva <- as.matrix(silva)
  kingdom <- gsub(pattern = ".*\t", replacement = "", x = silva[ ,1])
  seqid <- gsub(pattern = "\t.*$", replacement = "", x = silva[ ,1])
  silva <- silva[ ,-1]
  silva <- cbind(seqid, kingdom, silva)
  if(ncol(silva) == 8){
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("silva imported with a species column, the unique species names are: ", unique(silva[ ,8]),"\n")
  }else if(ncol(silva) == 7){
    silva <- cbind(silva, "unclassified")
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("silva imported without a species column, added one filled with \"unclassified\"\n")
  }else{
    cat("something is wrong with the silva import, it has ", ncol(silva), "columns in it.\n")
  }
  return(silva)
}

fill.blanks <- function(tax){
  # tax is only the taxonomy columns, seqID not included
  index <- which(tax == "")
  tax[index] <- "blank"
  cat("changing", length(index), "blanks to say \"blank\"\n")
  return(tax)
}

uniform.uncultured <- function(tax){
  # anything with uncultured in the name 
  index <- grep(pattern = 'uncultured', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"uncultured\"")
  tax[index] <- "uncultured"
  return(tax)
}

uniform.unknown <- function(tax){
  # anything with unknown in the name, esp like Unknown_Family
  index <- grep(pattern = 'unknown', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unknown\"")
  tax[index] <- "unknown"
  return(tax)
}

uniform.unnamed <- function(tax){
  # anything in front, _fa _ge with nothing after at end
  index <- grep(pattern = '.*_((ph)|(cl)|(or)|(fa)|(ge))$', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  if (length(index) > 0){
    cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\nIt might seem like some of these are real names, but if you search for them in silva you'll see they are just listing the name from an upper level. So don't worry. They'll get that upper level name back, just with unnamed as the prefix instead.")
  }
  tax[index] <- "unnamed"
  
  # p__ c__ at front with anything after
  index <- grep(pattern = '.{1}__$', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  if (length(index) > 0){
    cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  }
  tax[index] <- "unnamed"
  
  return(tax)
}

uniform.unclassified <- function(tax){
  # anything with unclassified in the name
  index <- grep(pattern = 'unclassified', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unclassified\"")
  tax[index] <- "unclassified"
  return(tax)
}

make.degenerates.unique <- function(tax, voldemort){
  # tax is the taxonomy table without seqid column
  # voldemort is whatever text you're making unique (b/c it cannot be named... voldemort...)
  for (t in 1:ncol(tax)){
    index <- which(tax[ ,t] == voldemort)
    tax[index,t] <- paste(tax[index,t], tax[index,t - 1], sep = ".")
  }
  
  remove.extra.u <- function(x){
    dot.voldemort <- paste(".", voldemort, sep = "")
    x <- gsub(pattern = dot.voldemort, replacement = "", x = x)
    return(x)
  }
  tax <- apply(X = tax, MARGIN = 2, FUN = remove.extra.u)
  
  return(tax)
}

revert.to.mothur.format <- function(silva){
  seqid.kingdom <- paste(silva[ ,1], silva[ ,2], sep = "\t")
  silva <- silva[ ,-(1:2)]
  silva <- cbind(seqid.kingdom, silva)
  return(silva)
}

# ---- go ----

stupid.silva <- import.mothur.formatted.tax(filepath = mothur.formatted.silva)
stupid <- stupid.silva
stupid.silva <- stupid

stupid.silva[ ,-1] <- fill.blanks(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.uncultured(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unknown(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unnamed(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unclassified(tax = stupid.silva[ ,-1])

x <- stupid.silva
stupid.silva <- x

stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "blank")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "uncultured")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unknown")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unnamed")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unclassified")

uniform.silva <- revert.to.mothur.format(silva = stupid.silva)

write.table(x = uniform.silva, file = new.silva.file, quote = F, col.names = FALSE, row.names = FALSE)

