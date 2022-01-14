# RRR

library(readxl)

key <- read_excel(path = "Cyanobacteria_silva138_Manual_Name_Edits.xlsx")

silva <- read.delim(file = "silva138_semicol.taxonomy", sep = ";", header = F)

for (r in 1:nrow(key)){
  species <- key$Silva138.Species[r]
  index <- which(silva[ ,8] == species)
  silva[index,2:8] <- key[r,12:18]
}

silva[ ,9] <- ""
colnames(silva) <- NULL

write.table(x = silva, file = "silva138_EditedCyanos_semicol.taxonomy", quote = F, sep = ";", row.names = F)
