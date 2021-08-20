# RRR
# make a key to match old names to new names
# save a version as rds for making the changes
# save a version with changes noted, to put in excel and make a reference


library("readxl")

edited <- read_excel(path = "silva138_cyanos_manually_edited.xlsx")
edited <- edited[ ,c(1,10)]

mono <- read.table(file = "new_names_monophyletic_semicol.taxonomy", sep = ";")
mono <- mono[ ,-9]
colnames(mono) <- c("Key.Ref.Num","Edited.Kingdom","Edited.Phylum","Edited.Class","Edited.Order","Edited.Family","Edited.Genus","Edited.Species")

silva <- read.table(file = "silva138_semicol.taxonomy", header = F, sep = ";", colClasses = "character")

old <- silva[ ,-c(1,9)]
index <- grep(pattern = "Cyanobacteria", x = old[ ,2], value = F)
old <- old[index, ]
old <- unique(old)
colnames(old) <- c("Silva138.Kingdom","Silva138.Phylum","Silva138.Class","Silva138.Order","Silva138.Family","Silva138.Genus","Silva138.Species")

old <- merge(y = edited, x = old, by.y = "Old.Species", by.x = "Silva138.Species")
old <- old[ ,c(8,2:7,1)]
old[1:5, ]

key <- merge(x = old, y = mono, by = "Key.Ref.Num")
key[1:5, ]

index <- order(key[ ,15], key[ ,14], key[ ,13], key[, 12], key[ ,11], key[ ,10], key[ ,9])
key <- key[index, ]

# Note changes for more manual edits:

get.change.index <- function(sil, ed){
  index <- sil != ed
  return(index)
}

which(get.change.index(sil = key$Silva138.Kingdom, ed = key$Edited.Kingdom))
which(get.change.index(sil = key$Silva138.Phylum, ed = key$Edited.Phylum))
which(get.change.index(sil = key$Silva138.Class, ed = key$Edited.Class))
which(get.change.index(sil = key$Silva138.Order, ed = key$Edited.Order))
which(get.change.index(sil = key$Silva138.Family, ed = key$Edited.Family))
which(get.change.index(sil = key$Silva138.Genus, ed = key$Edited.Genus))
which(get.change.index(sil = key$Silva138.Species, ed = key$Edited.Species))

# Only genus and species names have changes
genus.changed <- get.change.index(sil = key$Silva138.Genus, ed = key$Edited.Genus)
species.changed <- get.change.index(sil = key$Silva138.Species, ed = key$Edited.Species)

colnames(key)

key <- data.frame("Key.Ref.Num" = key[ ,1], "genus.changed" = genus.changed, "species.changed" = species.changed, key[ ,2:15])


write.csv(x = key, file = "Cyanobacteria_silva138_Name_Edits.csv", quote = F, row.names = F)
