# Generate the assorted stats I list in the paper text:


# (Abstract) ----
# Range of improvement in tribe/clade new classifications.
# run beginning of fig_2 script

stacked.eco.list
newly.named <- data.frame(matrix(data = 0, nrow = 5, ncol = 3))
colnames(newly.named) <- colnames(stacked.eco.list[[1]])[5:7]
for (e in 1:5){
  newly.named[e,] <- stacked.eco.list[[e]][3,5:7]
  row.names(newly.named)[e] <- names(stacked.eco.list)[e]
}

newly.named

min(newly.named$lineage)
max(newly.named$lineage)
mean(newly.named$lineage)

min(newly.named$clade)
max(newly.named$clade)
mean(newly.named$clade)

min(newly.named$tribe)
max(newly.named$tribe)
mean(newly.named$tribe)

# Marathonas error rates
# run fig_marathonas script to import data.
v3
v4
v5

incorrect <- cbind(v4[6:8,5:7], v5[6:8,5:7], v3[6:8,5:7])
incorrect <- rbind(incorrect, "tot" = colSums(incorrect))
incorrect # goes v4, v4-5, v3-4 

incorrect.perc <- incorrect / sum(v4[,1]) * 100
incorrect.perc # goes v4, v4-5, v3-4  (285 tot clones)

correct.perc <- 100 - incorrect.perc[4, ]
correct.perc

# (Results) ----

# Mara validation (fig 2) ----
# run fig_marathonas script to import data.

v4

correct.tribe.perc <- sum(v4[1:4,7]) / sum(v4[ ,7]) * 100 # correct, correct unclass, correct in Silva, underclass
correct.tribe.perc # all green / tot

perc.incorrect.incorrectly.in.silva <- v4[5,7] / sum(v4[5:8,7]) * 100
perc.incorrect.incorrectly.in.silva # yellow / red + yellow

really.bad.tribe <- c(v4[6:8,7], sum(v4[6:8,7])) / sum(v4[ ,7]) * 100 
really.bad.tribe # red and sum of red / tot


# Fine Res Class Increased (fig 3) ----
# run beginning of red-yellow fig script

# to list the increases in classification percent:
beside.eco.list$Mendota

# to list the increases in renamed reads:
stacked.eco.list$Mendota

# percent of GG family-names re-named by TaxAss:
beside.eco.list$Mendota[1,5] # named by GG
stacked.eco.list$Mendota[2,5] # renamed by TaxAss
stacked.eco.list$Mendota[2,5] / beside.eco.list$Mendota[1,5] * 100

# consistent across ecosystems 
# to list the genus-level improvements by ecosystem:
stacked.tax.list$`Genus/Clade`
stacked.tax.list$`Genus/Clade`[2, ] + stacked.tax.list$`Genus/Clade`[3, ] # sum "improved" reads
stacked.tax.list$`Species/Tribe`[2, ] + stacked.tax.list$`Species/Tribe`[3, ]
stacked.tax.list$`Family/Lineage`[2, ] + stacked.tax.list$`Family/Lineage`[3, ]

# Richness maintained (fig 3) ----

# go to supp_fig_1_pident_recalc_needed.R, run the "quick look" section ** skip to the "==> cyanos" comment to save time
# Percent Mendota dataset that is cyanobacteria:
sum(cyanos.all.reads[ ,2])

# go to supp_fig_1_pident_recalc_needed.R, run the "quick look" section
# percent cyanos classified as something else:
perc.cyano.classifications

# amount of focing at fine-level (fig 3b)
# source beginning of fig 3 script (select for Mendota)
# bacI looks worst:
forcing.data$lineage[ ,1:10]
forcing.data$lineage[2,5] # 2 = red bar, 5 = bacI
forcing.data$lineage[1,5] # 1 = grey bar
forcing.data$lineage[2,5] / forcing.data$lineage[1,5] * 100



