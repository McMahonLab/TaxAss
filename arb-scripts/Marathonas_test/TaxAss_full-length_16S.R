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

# ---- Functions made for this analysis specifically ----

find.correct.classifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  # correct classification defined as: taxa names agree and exist
  
  correct.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  for (t in 1:length(correct.class)){
    # subset to only classified names
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    if(length(index) > 0){
      fw.arb <- fw.arb[-index, ,drop = F]
    }
    
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    if(length(index) > 0){
      fw.tag <- fw.tag[-index, ,drop = F]
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
    correct.class[[t]] <- fw.tag[index, ,drop = F]
  }
  return(correct.class)
}

find.correct.unclassifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  # correct unclassification defined as: taxa names agree at upper levels and don't exist at lower levels
  
  correct.unclass <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  # save to re-start each loop step-down
  fwa <- fw.arb
  fwt <- fw.tag
  
  for (t in 1:length(correct.unclass)){
    fw.arb <- fwa
    fw.tag <- fwt
    
    # subset to only unclassified names
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    fw.arb <- fw.arb[index, ,drop = F]
    
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    fw.tag <- fw.tag[index, ,drop = F]
    
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
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      correct.unclass[[t]] <- fw.tag
      next
    }
    
    # subset to unclassifieds with same level unclassification as arb
    fw.tag[ ,-1] <- make.unclassifieds.unique(Taxonomy = fw.tag[ ,-1])
    fw.arb[ ,-1] <- make.unclassifieds.unique(Taxonomy = fw.arb[ ,-1])
    index <- which(fw.tag[ ,t + 1] == fw.arb[ ,t + 1])
    correct.unclass[[t]] <- fw.tag[index, ,drop = F]
  }
  return(correct.unclass)
}

find.underclassifications <- function(fw.arb, fw.tag){
  # input data is already subsetted to FreshTrain-classified seqs
  
  # underclassification defined as: taxa name exists in arb and unclassified in tag, but all existing upper names agree
  # only includes FreshTrain seqs because don't know "correct" class of others
  
  under.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  # save to re-start each loop step-down
  fwa <- fw.arb
  fwt <- fw.tag
  
  for (t in 1:length(under.class)){
    fw.arb <- fwa
    fw.tag <- fwt
    
    # subset to only unclassified names in tag data
    tag <- fw.tag[ ,t + 1]
    index <- which(tag == "unclassified")
    fw.tag <- fw.tag[index, ,drop = F]
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      under.class[[t]] <- fw.tag
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
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      under.class[[t]] <- fw.tag
      next
    }
    
    # subset to only classified names in arb data
    cat(all.equal(fw.tag[ ,1],fw.arb[ ,1]))
    index <- which(fw.arb[ ,t + 1] != "unclassified")
    fw.arb <- fw.arb[index, ]
    fw.tag <- fw.tag[index, ,drop = F]
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      under.class[[t]] <- fw.tag
      next
    }
    
    # subset to only matching upper-level names:
    # step up taxa levels, keeping only matching seqIDs as go
    for (u in 1:(t - 1)){
      cat(all.equal(fw.tag[ ,1],fw.arb[ ,1]))
      index <- which(fw.tag[ ,(t + 1)- u] == fw.arb[ ,(t + 1) - u] | fw.tag[ ,(t + 1)- u] == "unclassified")
      fw.tag <- fw.tag[index, ,drop = F]
      fw.arb <- fw.arb[index, ]
    }
    under.class[[t]] <- fw.tag
  }
  return(under.class)
}

find.over.classifications <- function(fw.arb, fw.tag){
  # input data is subsetted to FreshTrain-classified seqs
  
  # over-classification defined as: tag has name while arb is unclass, all existing upper names match
  # only includes FreshTrain arb b/c no way to know if upper-levels match for the non-FreshTrain arb classifications
  # this is like a special type of misclassification
  
  over.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  # save to re-start each loop step-down
  fwa <- fw.arb
  fwt <- fw.tag
  
  for (t in 1:length(over.class)){
    fw.arb <- fwa
    fw.tag <- fwt
    
    # subset to only unclassified names in arb data
    arb <- fw.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    fw.arb <- fw.arb[index, ,drop = F]
    
    if (nrow(fw.arb) < 1){ # skip ahead if there are none
      over.class[[t]] <- fw.arb
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
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      over.class[[t]] <- fw.tag
      next
    }
    
    # subset to only classified names in tag data
    cat(all.equal(fw.tag[ ,1],fw.arb[ ,1]))
    index <- which(fw.tag[ ,t + 1] != "unclassified")
    fw.arb <- fw.arb[index, ]
    fw.tag <- fw.tag[index, ,drop = F]
    
    if (nrow(fw.tag) < 1){ # skip ahead if there are none
      over.class[[t]] <- fw.tag
      next
    }
    
    # subset to only matching upper-level names:
    # step up taxa levels, keeping only matching seqIDs as go
    for (u in 1:(t - 1)){
      cat(all.equal(fw.tag[ ,1],fw.arb[ ,1]))
      index <- which(fw.tag[ ,(t + 1)- u] == fw.arb[ ,(t + 1) - u] | fw.arb[ ,(t + 1)- u] == "unclassified")
      fw.tag <- fw.tag[index, ,drop = F]
      fw.arb <- fw.arb[index, ]
    }
    over.class[[t]] <- fw.tag
  }
  return(over.class)
}

find.correct.GG.classifications <- function(all.arb, fw.arb, all.tag, fw.tag){
  # subset to GG-classified only
  gg.arb.ids <- setdiff(x = all.arb[ ,1], y = fw.arb[ ,1])
  gg.tag.ids <- setdiff(x = all.tag[ ,1], fw.tag[ ,1])
  
  all.arb <- as.data.frame(all.arb, stringsAsFactors = F)
  gg.arb <- merge(x = all.arb, y = gg.arb.ids, by = 1, sort = T)
  
  all.tag <- as.data.frame(all.tag, stringsAsFactors = F)
  gg.tag <- merge(x = all.tag, y = gg.tag.ids, by = 1, sort = T)
  
  # subset to shared ID's
  shared <- intersect(gg.arb.ids, gg.tag.ids)
  gg.arb <- merge(x = gg.arb, y = shared, by = 1, sort = T)
  gg.tag <- merge(x = gg.tag, y = shared, by = 1, sort = T)
  
  # no results list because "true" GG name unknown so all that GG handled both times are "correct"
  return(gg.tag)
}

# ---- In progress ----





find.mis.classifications <- function(all.arb, all.tag){
  # input data is all classified seqs
  
  # mis-classification defined as: 
  # - mismatched names, 
  # - unclass with mismatched upper names
  # - GG name when should be FreshTrain
  # - FreshTrain name when should be GG
  # excludes: unclass with matching upper names (over & under classifications)
  
  mis.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  # save to re-start each loop step-down
  aa <- all.arb
  at <- all.tag
  
  for (t in 1:length(mis.class)){
    all.arb <- aa
    all.tag <- at
    
    # subset to only unclassified names in arb data
    arb <- all.arb[ ,t + 1]
    index <- which(arb == "unclassified")
    all.arb <- all.arb[index, ,drop = F]
    
    if (nrow(all.arb) < 1){ # skip ahead if there are none
      mis.class[[t]] <- all.arb
      next
    }
    
    # subset to only shared ID's:
    arb <- all.arb[ ,1]
    tag <- all.tag[ ,1]
    shared <- intersect(x = arb, y = tag)
    
    all.arb <- as.data.frame(all.arb, stringsAsFactors = F)
    all.arb <- merge(x = all.arb, y = shared, by = 1, sort = T)
    
    all.tag <- as.data.frame(all.tag, stringsAsFactors = F)
    all.tag <- merge(x = all.tag, y = shared, by = 1, sort = T)
    
    if (nrow(all.tag) < 1){ # skip ahead if there are none
      mis.class[[t]] <- all.tag
      next
    }
    
    # subset to only classified names in tag data
    cat(all.equal(all.tag[ ,1],all.arb[ ,1]))
    index <- which(all.tag[ ,t + 1] != "unclassified")
    all.arb <- all.arb[index, ]
    all.tag <- all.tag[index, ,drop = F]
    
    if (nrow(all.tag) < 1){ # skip ahead if there are none
      mis.class[[t]] <- all.tag
      next
    }
    
    # subset to only matching upper-level names:
    # step up taxa levels, keeping only matching seqIDs as go
    for (u in 1:(t - 1)){
      cat(all.equal(all.tag[ ,1],all.arb[ ,1]))
      index <- which(all.tag[ ,(t + 1)- u] == all.arb[ ,(t + 1) - u] | all.arb[ ,(t + 1)- u] == "unclassified")
      all.tag <- all.tag[index, ,drop = F]
      all.arb <- all.arb[index, ]
    }
    mis.class[[t]] <- all.tag
  }
  return(mis.class)
  
  
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

correct.unclass.list <- find.correct.unclassifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)

under.class.list <- find.underclassifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)

over.class.list <- find.over.classifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)

gg.class.table <- find.correct.GG.classifications(all.arb = arb.tax, all.tag = v4.tax, fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
