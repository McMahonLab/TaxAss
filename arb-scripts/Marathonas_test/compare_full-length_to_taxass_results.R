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

userprefs <- commandArgs(trailingOnly = TRUE)

# for import
file.arb.tax <- userprefs[1]

file.v4.tax <- userprefs[2]
file.v4.ids <- userprefs[3]

# for export
folder.v4 <- userprefs[4]

# # troubleshoot: ----
# cat("CRAP COMMENT OUT FILE PATHS!!!")
# # for import
# file.arb.tax <- "../arb-scripts/Marathonas_test/mara.taxonomy"
# 
# file.v4.tax <- "../arb-scripts/Marathonas_test/v4_mara/otus.100.80.80.taxonomy"
# file.v4.ids <- "../arb-scripts/Marathonas_test/v4_mara/ids.above.100"
# 
# # for export
# folder.v4 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_100/"
# # end troubleshoot ----

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
      fw.arb <- fw.arb[index, ,drop = F]
    }
    under.class[[t]] <- cbind(fw.arb, fw.tag)
    colnames(under.class[[t]]) <- paste(colnames(under.class[[t]]), c(rep.int(x = "arb", times = 8), rep.int(x = "tag", times = 8)), sep = ".")
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
      fw.arb <- fw.arb[index, ,drop = F]
    }
    over.class[[t]] <- cbind(fw.arb, fw.tag)
    colnames(over.class[[t]]) <- paste(colnames(over.class[[t]]), c(rep.int(x = "arb", times = 8), rep.int(x = "tag", times = 8)), sep = ".")
  }
  return(over.class)
}

find.correct.GG.classifications <- function(all.arb, fw.arb, all.tag, fw.tag){
  # wasn't in FreshTrain in arb alignment, was classified by gg in TaxAss
  
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
  
  # no results list because "true" GG name unknown so all that GG handled both times are "correct" at every level
  return(gg.tag)
}

find.incorrect.GG.classifications <- function(fw.arb, all.tag, fw.tag){
  # was in the FreshTrain in arb alignment, was classified by GG in TaxAss
  
  # subset to GG-classified tags
  gg.tag.ids <- setdiff(x = all.tag[ ,1], fw.tag[ ,1])
  all.tag <- as.data.frame(all.tag, stringsAsFactors = F)
  gg.tag <- merge(x = all.tag, y = gg.tag.ids, by = 1, sort = T)
  
  # subset to shared IDs with freshwater arb IDs
  shared <- intersect(gg.tag[ ,1], fw.arb[ ,1])
  fw.arb <- as.data.frame(fw.arb, stringsAsFactors = F)
  gg.arb <- merge(x = fw.arb, y = shared, by = 1, sort = T)
  gg.tag <- merge(x = gg.tag, y = shared, by = 1, sort = T)
  
  # combine b/c interesting to compare
  all.equal(gg.arb[ ,1], gg.tag[ ,1])
  incorrectly.gg <- cbind(gg.arb, gg.tag)
  colnames(incorrectly.gg) <- paste(colnames(incorrectly.gg), c(rep.int(x = "arb", times = 8), rep.int(x = "tag", times = 8)), sep = ".")
  return(incorrectly.gg)
}

find.incorrect.FT.classifications <- function(all.arb, fw.arb, fw.tag){
  # was not in FreshTrain in arb alignment, was classified by FreshTrain in TaxAss
  
  # subset to non-FreshTrain arb
  gg.arb.ids <- setdiff(x = all.arb[ ,1], y = fw.arb[ ,1])
  all.arb <- as.data.frame(all.arb, stringsAsFactors = F)
  gg.arb <- merge(x = all.arb, y = gg.arb.ids, by = 1, sort = T)
  
  # subset to shared IDs with FreshTrain tags
  shared <- intersect(x = gg.arb[ ,1], fw.tag[ ,1])
  fw.tag <- as.data.frame(fw.tag, stringsAsFactors = F)
  gg.arb <- merge(x = gg.arb, y = shared, by = 1, sort = T)
  fw.tag <- merge(x = fw.tag, y = shared, by = 1, sort = T)
  
  # combine b/c interesting to compare
  all.equal(gg.arb[ ,1], fw.tag[ ,1])
  incorrectly.ft <- cbind(gg.arb, fw.tag)
  colnames(incorrectly.ft) <- paste(colnames(incorrectly.ft), c(rep.int(x = "arb", times = 8), rep.int(x = "tag", times = 8)), sep = ".")
  return(incorrectly.ft)
}

find.mis.classifications <- function(c.class, c.unclass, under, over, c.gen, i.gen, i.fresh, all.arb, all.tag){
  # input data is all the other summaries because this is basically just all the rest
  
  # mis-classification defined as: mismatched names and unclass with mismatched upper names withing FreshTrain orgs
  # excludes: unclass with matching upper names (over & under classifications)
  # excludes: seqs classified in GG in one and FT in other & vice versa (incorrect gg and ft classifications)
  
  mis.class <- list("kingdom" = NULL,"phylum" = NULL,"class" = NULL,"order" = NULL,"lineage" = NULL,"clade" = NULL,"tribe" = NULL)
  
  all.tag <- as.data.frame(all.tag, stringsAsFactors = F)
  all.arb <- as.data.frame(all.arb, stringsAsFactors = F)
  
  for (t in 1:length(mis.class)){
    cc <- c.class[[t]][ ,1]
    cu <- c.unclass[[t]][ ,1]
    un <- under[[t]][ ,1]
    ov <- over[[t]][ ,1]
    cg <- c.gen[ ,1]
    ig <- i.gen[ ,1]
    ft <- i.fresh[ ,1]
    
    # Check that no duplicate IDs
    taken.ids <- c(cc,cu,un,ov,cg,ig,ft)
    cat(length(unique(x = taken.ids)) == length(cc) + length(cu) + length(un) + length(ov) + length(cg) + length(ig) + length(ft))
    
    # subset to the not-taken seqIDs
    mis.ids <- setdiff(x = all.tag[ ,1], y = taken.ids)
    mis.tag <- merge(x = all.tag, y = mis.ids, by = 1, sort = T)
    mis.arb <- merge(x = all.arb, y = mis.ids, by = 1, sort = T)
    
    if (nrow(all.arb) < 1){ # skip ahead if there are none
      mis.class[[t]] <- all.arb
      next
    }
    
    mis.class[[t]] <- cbind(mis.arb, mis.tag)
    colnames(mis.class[[t]]) <- paste(colnames(mis.class[[t]]), c(rep.int(x = "arb", times = 8), rep.int(x = "tag", times = 8)), sep = ".")
  }
  return(mis.class)
}

# ---- Functions to visualize results ----

list.to.csv <- function(the.folder, the.list){
  dir.create(path = the.folder)
  for (n in 1:length(the.list)){
    file.name <- paste(the.folder, "/", n, "-", names(the.list)[n], ".csv", sep = "")
    write.csv(x = the.list[[n]], file = file.name)
    cat("made file:", file.name, "\n")
  }
}

create.summary.table <- function(){ # lazy calls to global environment...
  results.v4 <- data.frame(matrix(data = "bla", nrow = 9, ncol = 7))
  for( t in 1:7){
    temp <- rbind(names(correct.class.list)[t], # taxa level
                  nrow(over.class.list[[t]]),      # arb is FW unclass, tag is FW name, upper-names match
                  nrow(mis.class.list[[t]]),       # arb FW name != tag FW name
                  nrow(incorrect.ft.table),        # arb is GG, tag is FW
                  nrow(incorrect.gg.table),        # arb is FW, tag is GG
                  nrow(under.class.list[[t]]),     # arb is FW name, tag is FW unclass
                  nrow(correct.gg.table),          # arb is GG, tag is GG
                  nrow(correct.unclass.list[[t]]), # arb is FW unclass, tag is FW unclass
                  nrow(correct.class.list[[t]])    # arb FW name == tag FW name
                  )
    results.v4[ ,t] <- temp
  }
  colnames(results.v4) <- results.v4[1, ]
  results.v4 <- results.v4[-1, ]
  row.names(results.v4) <- c("overclassifications",
                             "misclassifications", 
                             "incorrectly in FreshTrain",
                             "incorrectly in greengenes",
                             "underclassifications",
                             "correctly in greengenes",
                             "correct unclassifications",
                             "correct classifications")
  results.v4 <- results.v4[8:1, ] # flip table so correct on bottom in bar plots
  temp <- row.names(results.v4)
  results.v4 <- as.matrix(results.v4)
  results.v4 <- apply(X = results.v4, MARGIN = 2, FUN = as.numeric)
  row.names(results.v4) <- temp
  return(results.v4)
}

make.stacked.bar <- function(){ # laxy calls from global env
  label.cex <- 1.2
  axis.cex <- 1.3
  title.cex <- 1.5
  
  col.correct <- adjustcolor(col = "darkgreen")
  col.under <- adjustcolor(col = rainbow(n = 20, v = .8)[4])
  col.wrong <- adjustcolor(col = "darkred")
  
  line.loc <- cumsum(results.v4[ ,7])
  line.loc <- line.loc - (.5 * results.v4[ ,7])
  label.loc <- line.loc
  for (n in 2:(nrow(results.v4))){
    if (n < nrow(results.v4)){
      sep <- .5 * (label.loc[n] - label.loc[n - 1]) + .5 * (label.loc[n + 1] - label.loc[n])
    }else{
      sep <- (label.loc[n] - label.loc[n - 1])
    }
    
    if (sep < 15){
      label.loc[n] <- label.loc[n] + (15 - sep)
    }
  }
  
  max.val <- sum(results.v4[ ,1])
  
  par(mar = c(2,4,3,12))
  bar.loc <- barplot(results.v4[ ,-(1:4)], col = c(col.correct, col.correct, col.correct, col.under, col.under, col.wrong, col.wrong, col.wrong), axisnames = F, axes = F)
  
  mtext(text = colnames(results.v4)[-(1:4)], side = 1, line = .5, at = bar.loc, cex = label.cex)
  
  axis(side = 2, at = c(0, max.val), labels = c("",max.val), cex = label.cex)
  axis(side = 2, cex = label.cex)
  mtext(text = "Number of Sequences", side = 2, line = 2.3, cex = axis.cex)
  
  text(x = rep.int(x = 3.8, times = nrow(results.v4)), y = label.loc, labels = row.names(results.v4), xpd = T, adj = 0, cex = label.cex)
  for (n in 1:nrow(results.v4)){
    lines(x = c(3.62,3.78), y = c(line.loc[n], label.loc[n]), xpd = T)
  }
  
  mtext(text = "TaxAss Accuracy", side = 3, line = 1.3, cex = title.cex)
  
} 




# ---- Import & Format Data ----

arb.tax <- read.csv(file.arb.tax, header = T, colClasses = "character")
arb.tax <- as.matrix(arb.tax)

index <- which(arb.tax[ ,6] == "unclassified") # which lineage == unclassified
fw.arb.tax <- arb.tax[-index, ]

v4.tax <- read.csv(file = file.v4.tax, header = T, colClasses = "character")
v4.tax <- as.matrix(v4.tax)
v4.tax <- apply(X = v4.tax, MARGIN = 2, FUN = remove.parentheses)

v4.ids <- read.csv(file.v4.ids, colClasses = "character", header = F)

index <- find.fw.seqid.indeces(FullTable = v4.tax, FWids = v4.ids)
fw.v4.tax <- v4.tax[index, ]

# ---- Compare Classifications ----

# identify classification types ----
correct.class.list <- find.correct.classifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
correct.unclass.list <- find.correct.unclassifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
under.class.list <- find.underclassifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
over.class.list <- find.over.classifications(fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
correct.gg.table <- find.correct.GG.classifications(all.arb = arb.tax, all.tag = v4.tax, fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
incorrect.gg.table <- find.incorrect.GG.classifications(fw.arb = fw.arb.tax, all.tag = v4.tax, fw.tag = fw.v4.tax)
incorrect.ft.table <- find.incorrect.FT.classifications(all.arb = arb.tax, fw.arb = fw.arb.tax, fw.tag = fw.v4.tax)
mis.class.list <- find.mis.classifications(c.class = correct.class.list, c.unclass = correct.unclass.list, under = under.class.list, over = over.class.list, c.gen = correct.gg.table, i.gen = incorrect.gg.table, i.fresh = incorrect.ft.table, all.arb = arb.tax, all.tag = v4.tax)

# save tables ----
the.folder <- paste(folder.v4, "correct_class", sep = "/")
list.to.csv(the.folder = the.folder, the.list = correct.class.list)

the.folder <- paste(folder.v4, "correct_unclass", sep = "/")
list.to.csv(the.folder = the.folder, the.list = correct.unclass.list)

the.folder <- paste(folder.v4, "under_class", sep = "/")
list.to.csv(the.folder = the.folder, the.list = under.class.list)

the.folder <- paste(folder.v4, "over_class", sep = "/")
list.to.csv(the.folder = the.folder, the.list = over.class.list)

the.file <- paste(folder.v4, "correct_gg.csv", sep = "/")
write.csv(x = correct.gg.table, file = the.file)

the.file <- paste(folder.v4, "incorrect_gg.csv", sep = "/")
write.csv(x = incorrect.gg.table, file = the.file)

the.file <- paste(folder.v4, "incorrect_ft.csv", sep = "/")
write.csv(x = incorrect.ft.table, file = the.file)

the.folder <- paste(folder.v4, "mis_class", sep = "/")
list.to.csv(the.folder = the.folder, the.list = mis.class.list)


# explore results ----
results.v4 <- create.summary.table()
write.csv(x = results.v4, file = paste(folder.v4, "summary_table.csv", sep = "/"), quote = F)

pdf(file = paste(folder.v4, "stacked_bar.pdf", sep = "/"), width = 6.875, height = 3, family = "Helvetica", title = "Marathonas Validation", colormodel = "srgb")
layout(mat = matrix(data = c(1,1,1,2,2), nrow = 1, ncol = 5))
make.stacked.bar()
dev.off()



# x <- incorrect.gg.table[ ,c(1,6:8,14:16)]
# # how many that went to GG are unclassified at clade/tribe
# y <- x[ ,-1] == "unclassified"
# y <- colSums(y) / nrow(y) * 100
# barplot(y, main = "Percent Unclassified")
# # which ones are classified in each?
# index <- which(x[ ,3] != "unclassified" & x[ ,6] != "unclassified") # clade
# x[index, ]
# index <- which(x[ ,2] != "unclassified" & x[ ,5] != "unclassified") # lineage
# x[index, ]













