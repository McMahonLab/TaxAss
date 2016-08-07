# crap there are 850 ACK-M1 in my dataset!  how are they not acI?
# taxonomy, OTU, and blast tables from Ananke

# ----- import -----
otu.table.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take18/data/otus.abund"
taxa.table.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take18/data/otus.98.80.70.taxonomy"
blast.file.path <- "~/Desktop/TaxonomyTrainingSets/BLASTing/take18 copy/otus.custom.blast.table.modified"

otus <- read.table(otu.table.path, header = T, colClasses = "character")

blast <- read.table(blast.file.path, colClasses = "character")
colnames(blast) <- c("qseqid","pident","length","qlen","q.align","true.pids","hit.num")

taxa <- read.table(taxa.table.path, header = T, colClasses = "character", sep = ",")
remove.parentheses <- function(x){
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}
taxa <- apply(taxa, 2, remove.parentheses)

# ----- which ones are f__ACK-M1 -----
# there are 850 of them
index <- which(taxa[ ,6] == "f__ACK-M1")
length(index)
ack.seqIDs <- taxa[index,1]

# their blast results are mostly low, but some are very close to 98% cutoff
index <- NULL
for (s in 1:length(ack.seqIDs)){
  index <- c(index, which(blast[ ,1] == ack.seqIDs[s]))
}
blast.ack <- blast[index, ]
all(blast.ack[ ,1] == ack.seqIDs)
blast.ack[ ,-1] <- apply(blast.ack[ ,-1], 2, as.numeric)
summary(blast.ack$pident)
summary(blast.ack$true.pids)
# most are first hit
summary(blast.ack$hit.num)
# 40 % have 0-1 overhang, 45 % are 2-5 overhangs, 15 % are longer
# overhang = query length - query bp's in alignment = qlen - q.align 
overhangs <- (blast.ack$qlen - blast.ack$q.align)
summary(overhangs)
hist(x = overhangs)
sum(overhangs < 2) / length(overhangs)
sum(overhangs > 1 & overhangs < 6) / length(overhangs)
sum(overhangs > 5) / length(overhangs)
# of the overhangs longer than 1, how many could have been included if there was a single match in the overhang
index <- which(overhangs > 1)
blast.ack.hangs <- blast.ack[index, ]
summary(blast.ack.hangs$true.pids)
# a new pid function gives a freebie match in the overhang:
pid.fun <- function(NumericBlastMatrixRow){
  b <- NumericBlastMatrixRow
  pid <- (b[1] * b[2] + 1) / (b[2] - b[4] + b[3]) # (pident * length + 1) / (length - q.align + qlen)
  return(pid)
}
b <- blast.ack.hangs[ ,-1]
b <- as.matrix(b)
freebie.pids <- apply(b, 1, pid.fun)
summary(freebie.pids)
# the max is still 97.93, so they wouldn't have made the cutoff even if it wasn't conservative!

# you can't really see a difference when an overhang match is included
plot(density(freebie.pids), col = "red")
par(new = T)
lines(density(blast.ack.hangs$true.pids), col = "black")
x <- NULL
for (s in 1:length(freebie.pids)){
  x <- c(x, freebie.pids[s] - blast.ack.hangs$true.pids[s])
}
summary(x)
# so the increase there is miniscule- btwn .009 and .013 .  think this is an OK estimate then, certainly better than just the HSP

# generally how are things changing with correction?
plot(density(blast.ack$pident), col = "black")
lines(density(blast.ack$true.pids), col = "red")

# ----- this is a great demo of it's importance! show with everything! -----
blast[ ,-1] <- apply(blast[ ,-1], 2, as.numeric)
summary(blast$pident)
summary(blast$true.pids)

# the corrected tail makes it hard to see
plot(density(c(blast$pident, blast$true.pids)), type = "n", main = "Corrected vs. Uncorrected Blast Results")
lines(density(blast$pident), col = "black")
lines(density(blast$true.pids), col = "red")

plot(density(blast$pident), col = "black")
lines(density(blast$true.pids), col = "red")

# focus only on pident choices that are reasonable
plot(density(c(blast$pident, blast$true.pids)), type = "n", main = "Corrected vs. Uncorrected Blast Results", xlim = c(90,max(density(blast$pident)$x)), ylim = c(0,max(density(blast$pident)$y)))
lines(density(blast$pident), col = "black")
lines(density(blast$true.pids), col = "red")

plot(blast$true.pids, cex = .2, type = "n", main = "black: pident, red: corrected")
points(blast$pident, cex = .05, pch = 19, col = "black")
points(blast$true.pids, cex = .05, pch = 19, col = "red")

par(mfrow = c(1,2))
plot(blast$pident, cex = .05, pch = 19, col = "black", main = "pident")
plot(blast$true.pids, cex = .05, pch = 19, col = "red", main = "corrected")

par(mfrow = c(1,1))

# what about cyanos- are they  great example?
index <- which(taxa[ ,3] == "p__Cyanobacteria")
cyano.seqIDs <- taxa[index,1]
index <- NULL
for (s in 1:length(cyano.seqIDs)){
  index <- c(index, which(blast[ ,1] == cyano.seqIDs[s]))
}
blast.cyano <- blast[index, ]
# some were probably left out of BLAST, that's why they don't match
all(blast.cyano[ ,1] == cyano.seqIDs)
blast.cyano[ ,-1] <- apply(blast.cyano[ ,-1], 2, as.numeric)

# wow how is there a bunp at 100! wonder if those are the ones that were called acIV...
par(mfrow = c(1,2))
plot(density(blast.cyano$pident))
plot(density(blast.cyano$true.pids), col = "red")

par(mfrow = c(1,1))
plot(density(c(blast.cyano$pident, blast.cyano$true.pids)), type = "n", 
     main = "Cyanos: black = pident, red = corrected",
     xlim = c(min(density(blast.cyano$true.pids)$x), max(density(blast.cyano$pident)$x)),
     ylim = c(0, max(density(blast.cyano$pident)$y)))
lines(density(blast.cyano$pident))
lines(density(blast.cyano$true.pids), col = "red")

cyano.hist <- hist(c(blast.cyano$pident,blast.cyano$true.pids), plot = F, breaks = 100)
hist(blast.cyano$pident, breaks = cyano.hist$breaks, ylim = c(0, max(cyano.hist$counts)), border = "black", col = adjustcolor("black", alpha.f = .2), main = "Cyanos: black = pident, red = corrected")
par(new = T)
hist(blast.cyano$true.pids, breaks = cyano.hist$breaks, ylim = c(0, max(cyano.hist$counts)), border = "red", col = adjustcolor("red", alpha.f = .2), main = "")

plot(blast.cyano$true.pids, cex = .2, type = "n", main = "cyanos. black: pident, red: corrected")
points(blast.cyano$pident, cex = .05, pch = 19, col = "black")
points(blast.cyano$true.pids, cex = .05, pch = 19, col = "red")



# ----- how does f__ACK-M1 behave over time? -----

# each row is an OTU, each column is a date
index <- NULL
for (s in 1:length(ack.seqIDs)){
  index <- c(index, which(otus[ ,1] == ack.seqIDs[s]))
}
otus.ack <- otus[index, ]
ack.seqIDs <- otus.ack[ ,1]
otus.ack <- otus.ack[ ,-1]
otus.ack <- apply(otus.ack, 2, as.numeric)
row.names(otus.ack) <- ack.seqIDs
seqID.tots <- rowSums(otus.ack)
for (s in 1:length(seqID.tots)){
  for (d in 1:ncol(otus.ack)){
    otus.ack[s,d] <- otus.ack[s,d] / seqID.tots[s]
  }
}

plot(x = 0, type = "n", xlim = c(0,95), ylim = c(0, max(otus.ack)))
for (s in 1:nrow(otus.ack)){
  lines(otus.ack[s, ], col = adjustcolor("black", alpha.f = .2))
}

# ----- how does that compare to acI? -----

# there are 5023 of them
index <- which(taxa[ ,6] == "acI")
length(index)
acI.seqIDs <- taxa[index,1]

# their blast results seem pretty even btwn 98 and 100
index <- NULL
for (s in 1:length(acI.seqIDs)){
  index <- c(index, which(blast[ ,1] == acI.seqIDs[s]))
}
blast.acI <- blast[index, ]
all(blast.acI[ ,1] == acI.seqIDs)
blast.acI[ ,-1] <- apply(blast.acI[ ,-1], 2, as.numeric)
summary(blast.acI$pident)
summary(blast.acI$true.pids)
# this is actually pretty cool. can see that some of the 100% hits moved to 98 after correction
plot(density(blast.acI$pident), col = "black")
par(new = T)
lines(density(blast.acI$true.pids), col = "red")
# Check that this kernel-density thing isn't showing a smoothing artifact:
plot(blast.acI$pident, cex = .2)
# most are first hit
summary(blast.acI$hit.num)
# the max overhang is 2, prob the max mismatches also. most have no overhang
overhangs <- (blast.acI$qlen - blast.acI$q.align)
summary(overhangs)
hist(x = overhangs)
sum(overhangs == 0) / length(overhangs)
sum(overhangs == 1) / length(overhangs)
sum(overhangs == 2) / length(overhangs)



