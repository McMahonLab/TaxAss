# RRR
# Table 1 compares the database performances at each taxa level for both metrics:
# by % reads classified i.e. percent classified i.e. red/orange plot
# by % otus classified i.e. alpha diversity i.e. blue/red plot

# ---- Define File Paths ----

run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Mendota/TaxAss-Mendota/plots/"

# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/table_best-of-both-databases.csv"

# these filenames shouldn't change b/c scripts output as named:
taxass.and.gg.alpha <- "step_15_5a_Improvement_over_general-only/alpha_diversity_TaxAss_vs_General.csv"
taxass.and.gg.perc <- "step_15_5a_Improvement_over_general-only/WorkflowImprovement-BesideData-Reads.csv"
fw.alpha.and.perc <- "step_15_5b_Improvement_over_custom-only/Custom-only_classification_stats.csv"

# ---- Create Table 1 ----

a <- read.csv(file = paste(run.folder, taxass.and.gg.alpha, sep = "/"))
p <- read.csv(file = paste(run.folder, taxass.and.gg.perc, sep = "/"))
f <- read.csv(file = paste(run.folder, fw.alpha.and.perc, sep = "/"))

table.1 <- data.frame(TaxaLevel = c("Kingdom", "Phylum", "Class", "Order", "Family/Lineage", "Genus/Clade", "Species/Tribe"),
                      Silva.perc = unlist(p[1,-1]), TaxAss.perc = unlist(p[2,-1]), FreshTrain.perc = f[ ,2],
                      Silva.alpha = a[ ,3], TaxAss.alpha = a[ ,2], FreshTrain.alpha = f[ ,3])

table.1 <- table.1[-1, ]

table.1[ ,2:4] <- round(x = as.matrix(table.1[ ,2:4]), digits = 0)

table.1

write.csv(x = table.1, file = output.file.path, quote = F, row.names = F)

cat("Made file: ", output.file.path)
# in excel, make the boxes and column names pretty. 
# Change to simply "Greengenes TaxAss FreshTrain" as repeated headers
# add 2 column headers: "Percent Classified" over 1st 3 columns, "Taxonomic Richness" over second 3 columns
# change font to helvetica, size 11, bold headers and labels, italicize and unbold superscripts
# goal is width 6.875 inches. 1 page = 8.5 x 11, so set side margins to 0.82 and fit to page layout. 8.5 - 6.875 = 1.625, 1.625 / 2 = 0.8125
# column widths are 10.33 and taxa label column is 14

# not sure how to format, but can save as .ps from excel. 
# for now did this: 
# print to pdf from excel, open in ppt (word resizes), crop to size, select cropped image, save as image, as pdf

# ---- repeat with this input for supplemental other datasets views ----
run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Michigan/TaxAss-Michigan/plots/"
# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/Supplemental_table_best-of-both-Michigan.csv"

run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/Danube/TaxAss-Danube/plots/"
# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/Supplemental_table_best-of-both-Danube.csv"

run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/TroutBogEpi/TaxAss-TroutBogEpi/plots/"
# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/Supplemental_table_best-of-both-BogEpi.csv"

run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/TroutBogHypo/TaxAss-TroutBogHypo/plots/"
# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/Supplemental_table_best-of-both-BogHypo.csv"

run.folder <- "~/Desktop/2018-05-10_taxass_server_results_for_resubmission/MouseGut/TaxAss-MouseGut/plots/"
# output.file.path <- "~/Dropbox/PhD/Write It/draft 7/re-submission_figures/Supplemental_table_best-of-both-MouseGut.csv"


# ----