# RRR
# Table S1
# Coarse-level forcing check
# generated directly by find_classification_disagreements.R
# found after clean-up in analysis/conflicts_98/conflicts_summary.csv
# note: in draft had all the different pident choices so could see increase in forcing at lower levels.
# but the increase not noticeable until like 92% or so. So better to have just one simple example as used pident.

file.path <- "~/Desktop/TaxAss-BatchFiles-go/Mendota/TaxAss-Mendota/analysis/conflicts_98/conflicts_summary.csv"

table.path <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Supplemental_Table_1.csv"

coarse.conflicts <- read.csv(file = file.path, stringsAsFactors = F)

coarse.conflicts <- cbind(coarse.conflicts, PercConflicts = coarse.conflicts[ ,2] / coarse.conflicts[7,2] * 100)

coarse.conflicts <- coarse.conflicts[-(6:7), ]

colnames(coarse.conflicts) <- c("Taxa Level", "Total Conflicts", "Percent Conflicts")

coarse.conflicts[ ,1] <- c("Kingdom","Phylum","Class","Order","Family/Lineage")

coarse.conflicts$`Percent Conflicts` <- signif(x = coarse.conflicts$`Percent Conflicts`, digits = 2)

coarse.conflicts

write.csv(x = coarse.conflicts, file = table.path, row.names = F, quote = F)

# in excel, make the boxes and column names pretty. 
# change font to helvetica, size 12, bold headers and labels, italicize and unbold superscripts
# goal is width 6.875 inches. 1 page = 8.5 x 11, so set side margins to 0.82 and fit to page layout. 8.5 - 6.875 = 1.625, 1.625 / 2 = 0.8125
# page size doesn't matter too much because this is supplemental. Also leaving as vertical instead of horizontal since supplemental.
# add superscripts a b c and d

# print to pdf from excel, open in ppt (word resizes), click "reset size" under format picture, crop to size, select cropped image, save as image, as pdf