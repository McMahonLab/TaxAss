# RRR
# Table 3 ****NOW TABLE 2****
# BLAST hit chosen. This is generated directly in plot_blast_hit_stats.R

# Decided to remove the pident 0 from it though for the paper's table
# Decided including it just made it a little confusing.
# Also removing the lower pidents, 90 - 94

taxass.file <- "~/Desktop/TaxAss-BatchFiles-go/Mendota/TaxAss-Mendota/analysis/plots/step_4_check_BLAST_settings/BLAST_hits_used_for_pidents_0-100.csv"

paper.file <- "~/Dropbox/PhD/Write It/draft 6/new_figs/Table_2.csv"

paper.table <- read.csv(file = taxass.file)

paper.table

paper.table[ ,1] <- c("Hit 1", "Hit 2", "Hit 3", "Hit 4", "Hit 5")

paper.table <- paper.table[ ,-(2:7)]

colnames(paper.table) <- c("BLAST result", paste(95:100, "% ID", sep = ""))

options(scipen = 999)
paper.table[ ,2:7] <- signif(x = paper.table[ ,2:7], digits = 2) 

paper.table

write.csv(x = paper.table, file = paper.file, quote = F, row.names = F)
cat("Made file: ", output.file.path)

# in excel, make the boxes and column names pretty. 
# add column headers: "Percent Results with Highest Percent Identity"
# add superscripts a b and c to explain blast result, % ID, and title.
# change font to helvetica, size 11, bold headers and labels, italicize and unbold superscripts
# fix significant digits to be by significant digits- don't understand wtf R thinks that means...
# goal is width 6.875 inches. 1 page = 8.5 x 11, so set side margins to 0.82 and fit to page layout. 8.5 - 6.875 = 1.625, 1.625 / 2 = 0.8125

# not sure how to format, but can save as .ps from excel. 
# for now did this: 
# print to pdf from excel, open in ppt (word resizes), crop to size (unclick lock aspect ratio), select cropped image, save as image, as pdf
