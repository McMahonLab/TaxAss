# RRR

library(readxl)

key <- read_excel(path = "silva138_cyanos_manually_edited.xlsx")

key <- key[ ,1:9]

key[ ,9] <- ""

colnames(key) <- NULL

key

write.table(x = key, file = "new_names_semicoldelim.taxonomy", quote = F, row.names = F, sep = ";")
