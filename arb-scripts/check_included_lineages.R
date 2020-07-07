# RRR 6/5/20
# Before we used a python script to only include specific FreshTrain lineages,
# But now Trina has figured out how to mark which ones are part of the FreshTrain
# and export only those references from Arb. 
# So I am replacing that filtering script with a checking script. 
# This lists all the lineages included in the export, to double check that they are both
# all the ones you want and 
# only the ones you want!
# Having poor references means they will have high pident in blast and non-FreshTrain seqs
# will get recruited into the FreshTrain-classification, ending up with worse classification
# than they would have gotten in silva (either unclassified or forced into incorrect name)

# ---- input ----

userprefs <- commandArgs(trailingOnly = TRUE)
input.file <- userprefs[1] 
created.file.folder <- userprefs[2]


# # Manual Troubleshooting
# cat("\n\nForgot to comment out file paths!!!\n\n")
# input.file <- "../../2020-06-02_update_freshtrain/FT_semicol_noprefix_mothur_unnamed.tax"
# created.file.folder <- "../../2020-06-02_update_freshtrain/included_lineages"

# ---- functions ----

# this function copied from reformat_taxonomy_nomenclature.R
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
    cat("database imported with a species column, the unique species names are: ", unique(silva[ ,8]),"\n")
  }else if(ncol(silva) == 7){
    silva <- cbind(silva, "unnamed")
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported without a species column, added one filled with \"unnamed\"\n")
  }else if(ncol(silva) == 9){
    extracol <- unique(silva[ ,9])
    silva <- silva[ ,1:8]
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported with an extra 9th column containing:", extracol, ". This column was removed.\n")
  }else{
    cat("something is wrong with the database import, it has ", ncol(silva), "columns in it.\n")
  }
  return(silva)
}

# ---- do stuff ----

ft <- import.mothur.formatted.tax(filepath = input.file)

index <- grep(pattern = 'unnamed', x = ft[ ,6], value = FALSE, invert = FALSE, ignore.case = TRUE)

if (length(index) > 0){
  cat("There are ", length(index), "FreshTrain references with no lineage-level classification. These definitely need to be removed!\n")
  
  nft <- ft[index, ]
  i <- order(nft[ ,2], nft[ ,3], nft[ ,4], nft[ ,5], nft[ ,6], nft[ ,7], nft[ ,8])
  nft <- nft[i, ]
  
  ft <- ft[-index, ]
  
}else{
  cat("There are no unclassified lineages in the FreshTrain export- good.\n")
}


ft.lineages <- unique(ft[ ,2:6])
i <- order(ft.lineages[ ,1], ft.lineages[ ,2], ft.lineages[ ,3], ft.lineages[ ,4], ft.lineages[ ,5])
ft.lineages <- ft.lineages[i, ]

ft.clades <- unique(ft[ ,2:7])
i <- order(ft.clades[ ,1], ft.clades[ ,2], ft.clades[ ,3], ft.clades[ ,4], ft.clades[ ,5], ft.clades[ ,6])
ft.clades <- ft.clades[i, ]

ft.tribes <- unique(ft[ ,2:8])
i <- order(ft.tribes[ ,1], ft.tribes[ ,2], ft.tribes[ ,3], ft.tribes[ ,4], ft.tribes[ ,5], ft.tribes[ ,6], ft.tribes[ ,7])
ft.tribes <- ft.tribes[i, ]


# ---- export ----

if (length(index) > 0){
  created.file <- file.path(created.file.folder,"STOP-remove_these_from_FreshTrain_export.csv")
  cat("Making file: ", created.file, "\n")
  write.csv(x = nft, file = created.file, quote = F, row.names = F)
}

created.file <- file.path(created.file.folder,"Check_this_list_of_unique_lineages_included_in_FreshTrain_export.csv")
cat("Making file: ", created.file, "\n")
write.csv(x = ft.lineages, file = created.file, quote = F, row.names = F)

created.file <- file.path(created.file.folder,"Check_this_list_of_unique_clades_included_in_FreshTrain_export.csv")
cat("Making file: ", created.file, "\n")
write.csv(x = ft.clades, file = created.file, quote = F, row.names = F)

created.file <- file.path(created.file.folder,"Check_this_list_of_unique_tribes_included_in_FreshTrain_export.csv")
cat("Making file: ", created.file, "\n")
write.csv(x = ft.tribes, file = created.file, quote = F, row.names = F)

