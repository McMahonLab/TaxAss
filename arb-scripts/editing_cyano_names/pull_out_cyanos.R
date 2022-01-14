# RRR editing silva cyanos 8/19/21

silva <- read.table(file = "data/2021-08-19_edit-cyano-taxonomy/silva138_semicol.taxonomy", header = F, sep = ";", colClasses = "character")

cyano.index <- grep(pattern = "Cyanobacteria", x = silva[ ,3], value = F)

cyanos <- silva[cyano.index, ]

length(unique(cyanos[ ,7])) # 241 is do-able manually

cyanos <- cyanos[, -c(1,9)] # multiple references for some names, 4842 refs but only 241 names to edit

cyanos <- unique(cyanos)

colnames(cyanos) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

index <- order(cyanos[ ,7], cyanos[ ,6], cyanos[ ,5], cyanos[, 4], cyanos[ ,3], cyanos[ ,2], cyanos[ ,1])

cyanos <- cyanos[index, ]

write.csv(x = cyanos, file = "data/2021-08-19_edit-cyano-taxonomy/silva138_cyanos.csv", quote = F, row.names = F)
