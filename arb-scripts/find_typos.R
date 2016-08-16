# RRR 2-8-16
# finding typos in FW taxonomy by looking for duplicate names with different phylogenies above them

#####
# Define File Paths
#####

comma.delim.trainingset <- "~/Desktop/custom.taxonomy-semicolon-check"
results.folder <- "~/Desktop/tax_typos/fix_check"


#####
# Define Functions
#####

# import FW training set (semicolon delimited)
import.taxonomy <- function(FilePath){
  numcol <- max(count.fields(FilePath, sep=";"))
  taxonomy <- read.table(FilePath, sep=";", fill=T, stringsAsFactors = F, col.names=1:numcol)
  taxonomy<- taxonomy[1:8]
  colnames(taxonomy) <- c("seqID","kingdom","phylum","class","order","lineage","clade","tribe")
  taxonomy <- as.matrix(taxonomy)
  return(taxonomy)
}

# fix the blank spots that cause p-value errors in RDP-classifier
make.empty.names.unclassified <- function(TaxonomyTable){
  taxonomy <- TaxonomyTable
  index <- which(taxonomy == "")
  taxonomy[index] <- "unclassified"
  return(taxonomy)
}

# Make a unique taxonomy table (so that duplicate names are indexed)
make.names.unique <- function(TaxonomyTable){
  taxonomy <- TaxonomyTable
  index <- order(taxonomy[ ,2], taxonomy[ ,3], taxonomy[ ,4], taxonomy[ ,5], taxonomy[ ,6], taxonomy[ ,7], taxonomy[ ,8])
  taxonomy.ord <- taxonomy[index, ]
  
  for (t in 2:8){
    unique.taxa <- unique(taxonomy.ord[ ,2:t, drop = FALSE])          
    unique.names <- make.unique(unique.taxa[,(t-1)])       
    for (u in 1:nrow(unique.taxa)){                     
      for (r in 1:nrow(taxonomy.ord)){
        if (all(taxonomy.ord[r,2:t] == unique.taxa[u,1:(t-1)])){
          taxonomy.ord[r,t] <- unique.names[u]
        }
      }
    }
  }
  return(taxonomy.ord)
}

# Create a list with unique names at each taxonomic level 
list.names.by.tax.level <- function(TaxonomyTable){
  taxonomy.ord <- TaxonomyTable
  grouped.taxa <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  
  for (t in 1:7){
    grouped.taxa[[t]] <- unique(taxonomy.ord[ ,(2:(t+1)), drop = FALSE])
    colnames(grouped.taxa[[t]]) <- colnames(taxonomy.ord)[(2:(t+1))]
    rownames(grouped.taxa[[t]])<- NULL
  }
  return(grouped.taxa)
}

# pull out all of the names with periods in them from indexing, then look at only ones that weren't unclassified
find.typos.from.repeat.names <- function(TaxonomyList){
  grouped.taxa <- TaxonomyList
  typos <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  
  for (t in 1:7){
    typos[[t]] <- grep(x = grouped.taxa[[t]][ ,t], pattern = ".*\\..*", value = TRUE)
    index.unclassifieds <- grep(x = typos[[t]], pattern =  "unclassified.*", value = FALSE )
    typos[[t]] <- typos[[t]][-index.unclassifieds]
  }
  return(typos)
}

# find the corresponding seqIDs to check in ARB
find.taxonomy.of.typos <- function(TyposList){
  typos <- TyposList
  typo.taxonomies <- list("kingdom"=NULL, "phylum"=NULL, "class"=NULL, "order"=NULL, "lineage"=NULL, "clade"=NULL, "tribe"=NULL)
  
  for (t in 1:7){
    typo.names <- sub(x = typos[[t]], pattern = "\\..*", replacement = "")
    for (n in 1:length(typo.names)){
      seqID.indeces <- which(taxonomy[ ,(t+1)] == typo.names[n])
      typo.taxonomies[[t]] <- rbind(typo.taxonomies[[t]], taxonomy[seqID.indeces, ])
    }
  }
  return(typo.taxonomies)
}

# write these names to excel sheets for Trina to fix in ARB :) thanks trina!
export.typos.into.excel <- function(TyposList, FilePath){
  typo.taxonomies <- TyposList
  results.folder <- FilePath
  for (t in 1:7){
    write.csv(x = typo.taxonomies[[t]], file = paste(results.folder, "/", t, "_dups_at_", names(typos)[t], "_level.csv", sep = ""), row.names = FALSE)
  }
  cat("Files written in folder ", results.folder)
}

# Find the specific seqID names of the typos to make Trina's life easier
pull.out.incorrect.seqIDs <- function(TaxonomyList){
  typo.taxonomies <- TaxonomyList
  
  typo.seqIDs <- NULL
  seqID.counter <- 1
  # loop counters are t, d, n, s, o in that order- don't re-use these variables! :)
  for (t in 1:7){
    # ignore the taxa levels that have no duplicate names
    if (nrow(typo.taxonomies[[t]]) > 0){
      
      # look at the taxa for each duplicate name separately
      dup.names <- unique(typo.taxonomies[[t]][ ,(t + 1)])
      for (d in 1:length(dup.names)){
        dup.indeces <- which(typo.taxonomies[[t]][ ,(t + 1)] == dup.names[d])
        dup.taxonomies <- typo.taxonomies[[t]][dup.indeces, ]
        typo.names <- unique(dup.taxonomies[ ,t])
        
        # ignore the duplicates that arrise from typos more than one taxa level up
        if (length(typo.names) > 1){
          
          # assume the name that occurs most is correct
          num.names <- 0
          for (n in 1:length(typo.names)){
            num.names[n] <- length(which(dup.taxonomies[ ,t] == typo.names[n]))
          }
          name.index <- which(num.names == max(num.names))
          the.typo.is <- typo.names[-(name.index)]
          the.typo.isnt <- typo.names[name.index]
          
          # get the seqIDs for all of the names the typo is (i)
          seqID.index <- NULL
          for (i in 1:length(the.typo.is)){
            seqID.index <- which(dup.taxonomies[ ,t] == the.typo.is[i])
            the.seqID.is <- dup.taxonomies[seqID.index,1]
            
            # save the seqIDs of interest, allowing for there to be multiple seqIDs (s) for each different name the typo is (i) 
            for (s in 1:length(the.seqID.is)){
              typo.seqIDs[seqID.counter] <- paste("the seqID ", the.seqID.is[s], " has ", names(typo.taxonomies)[t-1], " ", the.typo.is[i],  
                                                  " while most ", dup.names[d], " are ", the.typo.isnt, sep = "")
              cat(typo.seqIDs[seqID.counter],"\n")
              seqID.counter <- seqID.counter + 1
            }
          }
        }
      }
    }
  }
  return(typo.seqIDs)
}

export.seqIDs.to.fix <- function(FilePath, SeqIDsWithTypos){
  results.folder <- FilePath
  typo.seqIDs <- SeqIDsWithTypos
  
  write.csv(x = typo.seqIDs, file = paste(results.folder, "/SeqIDs_to_Fix.csv", sep = ""), row.names = FALSE)
}

#####
# Use Functions
#####

taxonomy <- import.taxonomy(FilePath = comma.delim.trainingset)

taxonomy <- make.empty.names.unclassified(TaxonomyTable = taxonomy)

taxonomy.ord <- make.names.unique(TaxonomyTable = taxonomy)

grouped.taxa <- list.names.by.tax.level(TaxonomyTable = taxonomy.ord)

typos <- find.typos.from.repeat.names(TaxonomyList = grouped.taxa)

typo.taxonomies <- find.taxonomy.of.typos(TyposList = typos)

export.typos.into.excel(TyposList = typo.taxonomies, FilePath = results.folder)

typo.seqIDs <- pull.out.incorrect.seqIDs(TaxonomyList = typo.taxonomies)

export.seqIDs.to.fix(FilePath = results.folder, SeqIDsWithTypos = typo.seqIDs)




