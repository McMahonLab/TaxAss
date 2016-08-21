# acI clades across ecosystems:

mendota.unclust.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_mend_unclust/plots/FinalTaxonomyGroups/"

bogs.epi.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_epi/plots/FinalTaxonomyGroups/"

bogs.hypo.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_bogs_hypo/plots/FinalTaxonomyGroups/"

michigan.hiseq.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_michigan_hiseq/plots/FinalTaxonomyGroups/"

danube.10.final <- "~/Desktop/TaxonomyTrainingSets/BLASTing/poster_danube_10/plots/FinalTaxonomyGroups/"


# ---- define functions ----

import.grouped.folder <- function(FolderPath){
  grouped <- list(kingdom=NULL, phylum=NULL, class=NULL, order=NULL, lineage=NULL, clade=NULL, tribe=NULL)
  files <- list.files(FolderPath)
  for (t in 1:7){
    grouped[[t]] <- read.csv(file = paste(FolderPath, files[t], sep = ""), colClasses = "character")
    grouped[[t]][ ,t+1] <- as.numeric(grouped[[t]][ ,t+1])
  }
  return(grouped)
}

pull.out.taxon <- function(BigList, Tax, Level){
  big.table <- BigList[[Level]]
  index <- which(big.table[ ,Level] == Tax)
  taxon.tot <- big.table[index, Level + 1]
  return(taxon.tot)
}


mend <- import.grouped.folder(FolderPath = mendota.unclust.final)
mich <- import.grouped.folder(FolderPath = michigan.hiseq.final)
epi <- import.grouped.folder(FolderPath = bogs.epi.final)
hypo <- import.grouped.folder(FolderPath = bogs.hypo.final)
dan <- import.grouped.folder(FolderPath = danube.10.final)


a.mend <- pull.out.taxon(BigList = mend, Tax = "acI-A", Level = 6)
a.epi <- pull.out.taxon(BigList = epi, Tax = "acI-A", Level = 6)
a.hypo <- pull.out.taxon(BigList = hypo, Tax = "acI-A", Level = 6)
a.mich <- pull.out.taxon(BigList = mich, Tax = "acI-A", Level = 6)
a.dan <- pull.out.taxon(BigList = dan, Tax = "acI-A", Level = 6)

a <- c(a.mend, a.mich, a.epi, a.hypo, a.dan)
names(a) <- c("Eutrophic", "Oligotrophic", "Bog Epi", "Bog Hypo","River")

b.mend <- pull.out.taxon(BigList = mend, Tax = "acI-B", Level = 6)
b.epi <- pull.out.taxon(BigList = epi, Tax = "acI-B", Level = 6)
b.hypo <- pull.out.taxon(BigList = hypo, Tax = "acI-B", Level = 6)
b.mich <- pull.out.taxon(BigList = mich, Tax = "acI-B", Level = 6)
b.dan <- pull.out.taxon(BigList = dan, Tax = "acI-B", Level = 6)

b <- c(b.mend, b.mich, b.epi, b.hypo, b.dan)
names(b) <- c("Eutrophic", "Oligotrophic", "Bog Epi", "Bog Hypo","River")

c.mend <- pull.out.taxon(BigList = mend, Tax = "acI-C", Level = 6)
c.epi <- pull.out.taxon(BigList = epi, Tax = "acI-C", Level = 6)
c.hypo <- pull.out.taxon(BigList = hypo, Tax = "acI-C", Level = 6)
c.mich <- pull.out.taxon(BigList = mich, Tax = "acI-C", Level = 6)
c.dan <- pull.out.taxon(BigList = dan, Tax = "acI-C", Level = 6)

c <- c(c.mend, c.mich, c.epi, c.hypo, c.dan)
names(c) <- c("Eutrophic", "Oligotrophic", "Bog Epi", "Bog Hypo","River")


barplot(a)
barplot(b)
barplot(c)
max(a)
max(b)
max(c)


a.title <- expression(bold("acI-A"))
b.title <- expression(bold("acI-B"))
c.title <- expression(bold("acI-C"))

a.axis <- c(0,3,6,9,12)
a.ax.lab <- expression(bold("Rel. Abundance (% Reads)"))
file.name.a <- "~/Dropbox/Trina/8-20-16_ISME16_figures/next_steps_acI-A.png"
file.name.b <- "~/Dropbox/Trina/8-20-16_ISME16_figures/next_steps_acI-B.png"
file.name.c <- "~/Dropbox/Trina/8-20-16_ISME16_figures/next_steps_acI-C.png"
file.name.lab <- "~/Dropbox/Trina/8-20-16_ISME16_figures/next_steps_ylab.png"


# a

png(filename = file.name.a, width = 4.65, height = 5.54, units = "in", res = 300)
par(mar=c(5,2,2.5,0))
loc.labels <- barplot(height = a, names.arg = c("","","","",""), axes = F, col = "purple", ylim = c(0,13.5))
text(x = loc.labels, y = -.4, labels = names(a), adj = 1, cex = 1.5, srt = 40, xpd = T)
axis(side = 2, at = a.axis, labels = F, line = 0, lwd = 3, xpd=T)
mtext(text = a.axis, side = 2, line = .5, at = a.axis, cex = 1.5)
mtext(text = a.title, side = 3, line = 1, cex = 2)
dev.off()

png(filename = file.name.b, width = 4.65, height = 5.54, units = "in", res = 300)
par(mar=c(5,2,2.5,0))
loc.labels <- barplot(height = b, names.arg = c("","","","",""), axes = F, col = "red", ylim = c(0,13.5))
text(x = loc.labels, y = -.4, labels = names(a), adj = 1, cex = 1.5, srt = 40, xpd = T)
axis(side = 2, at = a.axis, labels = F, line = 0, lwd = 3, xpd=T)
mtext(text = a.axis, side = 2, line = .5, at = a.axis, cex = 1.5)
mtext(text = b.title, side = 3, line = 1, cex = 2)
dev.off()

c.axis <- c(0,.5,1,1.5)

png(filename = file.name.c, width = 4.65, height = 5.54, units = "in", res = 300)
par(mar=c(5,2,2.5,0))
loc.labels <- barplot(height = c, names.arg = c("","","","",""), axes = F, col = "blue")
text(x = loc.labels, y = -.06, labels = names(a), adj = 1, cex = 1.5, srt = 40, xpd = T)
axis(side = 2, at = c.axis, labels = F, line = 0, lwd = 3, xpd=T)
mtext(text = c.axis, side = 2, line = .5, at = c.axis, cex = 1.5)
mtext(text = c.title, side = 3, line = 1, cex = 2)
dev.off()



# get a label:

png(filename = file.name.lab, width = 4.65, height = 5.54, units = "in", res = 300)
plot.new()
mtext(text = a.ax.lab, side = 2, line = 2, cex = 1.5)
dev.off()

