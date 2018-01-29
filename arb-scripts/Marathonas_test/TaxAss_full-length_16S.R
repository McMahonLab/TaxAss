# RRR 1/17/18
# To satisfy reviewer #2, and also I think Pat Schloss bioRxiv review:
# Concern was: in blue/red improvement plot, "truth" is TaxAss to support using TaxAss over GG alone.
# The issue is that I'm testing TaxAss, so it shouldn't be used as the reference even though
# this example is just meant to demonstrate obvious errors in the GG-only method
# Approach for a better "truth":
# Ryan aligned some new full-length 16S sequences to the FreshTrain
# Take their assignment from using full length and arb as the "truth"
# I'll chop them into a variable region and then assign them using TaxAss
# This should validate that TaxAss works, and address concerns about overclassifying missing references

# 1. use mothur to trim full length seqs: trim_full_length_to_primer_region.sh

# 2. use TaxAss to assign taxonomy: results in otus.98.80.taxonomy

# 3. Compare to arb classifications: this custom script because differentiating unclassifieds more than find_classification_disagreements.R

# ---- File paths ----

file.arb.tax <- "../arb-scripts/Marathonas_test/mara.taxonomy"

file.v4.tax <- "../arb-scripts/Marathonas_test/mara_v4.98.80.80.taxonomy"
file.v4.ids <- "../arb-scripts/Marathonas_test/add-tax-scripts-and-databases/ids.above.98"


# ---- Functions taken from find_classification_disagreements.R ----

remove.parentheses <- function(x){
  # Remove parentheses of % confidence so that names in each table match exactly
  # call this using apply
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

find.fw.seqid.indeces <- function(FullTable, FWids){
  gg.and.fw <- FullTable
  fw.ids <- FWids
  # duplicate marks the 2nd occurance as TRUE
  # combine the ids (fw ids first, all ids second) into one vector with duplicate ids
  # index of FW in all ids is the index of the duplicates - the number of fw ids
  all.ids <- c(fw.ids[ ,1], gg.and.fw[ ,1])
  index <- which(duplicated(all.ids) == TRUE)
  index <- index - nrow(fw.ids)
  return(index)
}

# ---- Functions taken from plot_classification_improvement.R ----

make.unclassifieds.unique <- function(Taxonomy){
  # Taxonomy is a matrix containing only the names, no abundances or OTU numbers
  # this function maintains row order, so those columns can be re-combined after this
  # this is used in group.taxa()
  taxa.levels <- 2:ncol(Taxonomy)
  
  # first add the level-above name to each "unclassified"
  for (t in taxa.levels){
    index <- which(Taxonomy[ ,t] == "unclassified")
    Taxonomy[index,t] <- paste(Taxonomy[index,t], Taxonomy[index,t - 1], sep = ".")
  }
  
  # then remove extra unclassifieds to get only unclassified.last_known_name
  remove.extra.u <- function(x){
    x <- gsub(pattern = ".unclassified", replacement = "", x = x)
    return(x)
  }
  Taxonomy <- apply(X = Taxonomy, MARGIN = 2, FUN = remove.extra.u)
  
  return(Taxonomy)
}

# ---- Functions made/modified for this analysis specifically ----




find.correct.classifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  correct.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  for (t in 1:length(correct.class)){
    # subset to only classified names
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    if(length(index) > 0){
      fw.arb <- fw.arb[-index, ]
    }
    
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    if(length(index) > 0){
      fw.tag <- fw.tag[-index, ]
    }
    
    # subset to only shared ID's:
    arb <- fw.arb[ ,1]
    tag <- fw.tag[ ,1]
    shared <- intersect(x = arb, y = tag)
    
    fw.arb <- as.data.frame(fw.arb, stringsAsFactors = F)
    fw.arb <- merge(x = fw.arb, y = shared, by = 1, sort = T)
    
    fw.tag <- as.data.frame(fw.tag, stringsAsFactors = F)
    fw.tag <- merge(x = fw.tag, y = shared, by = 1, sort = T)
    
    # subset to only matching names:
    all.equal(fw.tag[ ,1],fw.arb[ ,1])
    index <- which(fw.tag[ ,t + 1] == fw.arb[ ,t + 1])
    correct.class[[t]] <- fw.tag[index, ]
  }
  return(correct.class)
}

correct.unclass.list <- find.correct.classifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)

find.correct.unclassifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  correct.unclass <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  for (t in 1:length(correct.unclass)){
    # subset to only unclassified names
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    fw.arb <- fw.arb[index, ]
    
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    fw.tag <- fw.tag[index, ]
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      correct.unclass[[t]] <- fw.tag
      next
    }
    
    # subset to only shared ID's:
    arb <- fw.arb[ ,1]
    tag <- fw.tag[ ,1]
    shared <- intersect(x = arb, y = tag)
    
    fw.arb <- as.data.frame(fw.arb, stringsAsFactors = F)
    fw.arb <- merge(x = fw.arb, y = shared, by = 1, sort = T)
    
    fw.tag <- as.data.frame(fw.tag, stringsAsFactors = F)
    fw.tag <- merge(x = fw.tag, y = shared, by = 1, sort = T)
    
    # subset to only matching upper-level names:
    for (u in 1:(t - 1)){
      all.equal(fw.tag[ ,1],fw.arb[ ,1])
      index <- which(fw.tag[ ,(t + 1)- u] == fw.arb[ ,(t + 1) - u])
      fw.tag <- fw.tag[index, ]
      fw.arb <- fw.arb[index, ]
    }
    correct.unclass[[t]] <- fw.tag[index, ]
  }
  return(correct.unclass)
}

find.underclassifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  correct.unclass <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  for (t in 1:length(correct.unclass)){
    # subset to only unclassified names
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    fw.arb <- fw.arb[index, ]
    
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    fw.tag <- fw.tag[index, ]
    
    # subset to only shared ID's:
    arb <- fw.arb[ ,1]
    tag <- fw.tag[ ,1]
    shared <- intersect(x = arb, y = tag)
    
    fw.arb <- as.data.frame(fw.arb, stringsAsFactors = F)
    fw.arb <- merge(x = fw.arb, y = shared, by = 1, sort = T)
    
    fw.tag <- as.data.frame(fw.tag, stringsAsFactors = F)
    fw.tag <- merge(x = fw.tag, y = shared, by = 1, sort = T)
    
    # subset to only matching upper-level names:
    for (u in 1:(t - 1)){
      all.equal(fw.tag[ ,1],fw.arb[ ,1])
      index <- which(fw.tag[ ,(t + 1)- u] == fw.arb[ ,(t + 1) - u])
      fw.tag <- fw.tag[index, ]
      fw.arb <- fw.arb[index, ]
    }
    correct.unclass[[t]] <- fw.tag[index, ]
  }
  return(correct.unclass)
}

find.incorrect.classifications 

find.incorrect.unclassifications 

find.incorrect.classifications <- function(taxass, arb, FWtable_percents, GGtable_percents, TaxaLevel, tracker, FolderPath, forcing = FALSE){
  # Find seqs misclassified at a given phylogenetic level, t
  fw <- FWtable
  gg <- GGtable
  fw.percents <- FWtable_percents
  gg.percents <- GGtable_percents
  t <- TaxaLevel
  num.mismatches <- tracker
  results.folder.path <- FolderPath
  
  taxa.names <- c("kingdom","phylum","class","order","lineage","clade","tribe")
  
  # compare names in column t+1, because first columns are seqID and t=1 is kingdom, t=7 is tribe
  # ignore names that say unclassified, except in forcing comparison when fw giving any erroneous name counts as forcing
  if (forcing == TRUE){
    index <- which(gg[,t+1] != fw[,t+1] & fw[,t+1] != "unclassified")
  }else{
    index <- which(gg[,t+1] != fw[,t+1] & gg[,t+1] != "unclassified" & fw[,t+1] != "unclassified")
  }
  
  cat("there are ", length(index), " conflicting names at ", taxa.names[t], " level\n")
  num.mismatches[t] <- length(index)
  
  # Compare the conflicting tables in entirety, use the original files with percents still in it
  conflicting <- cbind(gg.percents[index,,drop=F], fw.percents[index,,drop=F])
  
  # Check that the files still line up correctly
  check.files.match(FWtable = conflicting[,9:16,drop=F], GGtable = conflicting[,1:8,drop=F])
  
  # Export a file with the conflicting rows side by side.
  write.csv(conflicting, file = paste(results.folder.path, "/", t, "_", taxa.names[t],"_conflicts.csv", sep=""), row.names = FALSE)
  
  #Track the number of mismatches at each level
  return(num.mismatches)
}

# # example of logic:
# arb <- letters[1:5]
# tag <- letters[2:6]
# index.arb <- duplicated(c(tag,arb))
# index.arb <- index.arb[-(1:length(tag))]
# index.tag <- duplicated(c(arb,tag))
# index.tag <- index.tag[-(1:length(arb))]
# arb <- arb[index.arb]
# tag <- tag[index.tag]


# ---- Import & Format Data ----

arb.tax <- read.csv(file.arb.tax, header = T, colClasses = "character")
arb.tax <- as.matrix(arb.tax)

index <- which(arb.tax[ ,6] == "unclassified")
fw.arb.tax <- arb.tax[-index, ]

v4.tax <- read.csv(file = file.v4.tax, header = T, colClasses = "character")
v4.tax <- as.matrix(v4.tax)
v4.tax <- apply(X = v4.tax, MARGIN = 2, FUN = remove.parentheses)

v4.ids <- read.csv(file.v4.ids, colClasses = "character", header = F)

index <- find.fw.seqid.indeces(FullTable = v4.tax, FWids = v4.ids)
fw.v4.tax <- v4.tax[index, ]

# ---- Analyze Data ----

correct.class.list <- find.correct.classifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)



