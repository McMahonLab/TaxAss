# RRR
# Table 3 ****NOW TABLE 2****
# BLAST hit chosen. This is generated directly in plot_blast_hit_stats.R

# Decided to remove the pident 0 from it though for the paper's table
# Decided including it just made it a little confusing.
# Also removing the lower pidents, 90 - 94

taxass.file <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/step_4_check_BLAST_settings/BLAST_hits_used_for_pidents_0-100.csv"

# paper.file <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/table_blast-hit.csv"

paper.table <- read.csv(file = taxass.file)

paper.table

paper.table[ ,1] <- c("Hit 1", "Hit 2", "Hit 3", "Hit 4", "Hit 5")

paper.table <- paper.table[ ,-(2:7)]

colnames(paper.table) <- c("BLAST result", paste(95:100, "% ID", sep = ""))

options(scipen = 999)
paper.table[ ,2:7] <- round(x = paper.table[ ,2:7], digits = 3) 

paper.table

write.csv(x = paper.table, file = paper.file, quote = F, row.names = F)
cat("Made file: ", paper.file)

# in excel, make the boxes and column names pretty. 
# add column headers: "Percent Identity Cutoff Applied
# merge Blast result header into 2 vertical cells.
# add superscript a after Percent Identity Cutoff Applied
# add superscript b after the top left value
# change font to helvetica, size 11, bold headers and labels, italicize and unbold superscripts, center text
# show only max digits needed per column to display smallest value
# goal is width 6.875 inches. 1 page = 8.5 x 11, so set side margins to 0.82 and fit to page layout. 8.5 - 6.875 = 1.625, 1.625 / 2 = 0.8125
# save as xlxs not csv or formatting is lost

# print to pdf from excel
# open in ppt,unclick lock aspect ratio box, click on "reset size" to get original size back,
# crop to size, select cropped image, save as image, as pdf
