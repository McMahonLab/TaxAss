# RRR 6/6/20
# Make an R script that changes corse-level database conflicts by directly 
# reading in the output of the find_classifications_disagreements.R -database
# script. I don't feel like manually changing the python script each time. 
# Esp b/c now the phylum actinobacteria name changed and that means there's 
# a lot more lineages who need a phylum name adjustment. 

# ---- input ----

userprefs <- commandArgs(trailingOnly = TRUE)
folder.path <- userprefs[1]
old.FreshTrain <- userprefs[2]
created.noconflict.FreshTrain <- userprefs[3]



# # Troubleshoot with local paths
# cat("\n\n\nYOU FORGOT TO COMMENT OUT FILE PATHS!!\n\n\n")
# folder.path <- "../../2020-06-02_update_freshtrain/run_taxass/conflicts_database/"
# old.FreshTrain <- "../../2020-06-02_update_freshtrain/run_taxass/FT.taxonomy"
# created.noconflict.FreshTrain <- "../../2020-06-02_update_freshtrain/run_taxass/FT_noconflicts.taxonomy"

# ---- functions ----

change.conflicts <- function(t, conflicts, ft){
  if(nrow(conflicts) > 0){
    lineages <- conflicts[ ,c(t,ncol(conflicts))]
    for (l in 1:nrow(lineages)){
      index <- which(ft[ ,6] == lineages[l,2])
      cat("Changing ", lineages[l,2], " from ", ft[index[1],t+1], " to ", lineages[l,1], "\n")
      ft[index,t+1] <- lineages[l,1]
    }
  }else{
    cat("No conflicts\n")
  }
  return(ft)
}

# ---- go! ----

kingdom.conflicts <- read.csv(file = file.path(folder.path,"unique_conflicts_1_kingdom.csv"))
phylum.conflicts <- read.csv(file = file.path(folder.path,"unique_conflicts_2_phylum.csv"))
class.conflicts <- read.csv(file = file.path(folder.path,"unique_conflicts_3_class.csv"))
order.conflicts <- read.csv(file = file.path(folder.path,"unique_conflicts_4_order.csv"))

FreshTrain <- read.delim(file = old.FreshTrain, header = F, sep = ";", colClasses = "character")

cat("\nKingdom conflicts\n")
FreshTrain <- change.conflicts(t = 1, conflicts = kingdom.conflicts, ft = FreshTrain)

cat("\nPhylum conflicts\n")
FreshTrain <- change.conflicts(t = 2, conflicts = phylum.conflicts, ft = FreshTrain)

cat("\nClass conflicts\n")
FreshTrain <- change.conflicts(t = 3, conflicts = class.conflicts, ft = FreshTrain)

cat("\nOrder conflicts\n")
FreshTrain <- change.conflicts(t = 4, conflicts = order.conflicts, ft = FreshTrain)

# ---- export ----

cat("Making file: ", created.noconflict.FreshTrain, "\n")
write.table(x = FreshTrain, file = created.noconflict.FreshTrain, append = F, quote = F, sep = ";", row.names = F)
