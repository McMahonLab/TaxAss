# 2018-3-8 RRR

# Look at how pident choice changes database assignment errors

# Re-Run step 15 on the v4_mara or whichever folder using different pidents.

# Run through main R script to get summary tables or error types exported

# Those tables are the input into this script.

# Do number of seqs incorrectly-assigned to GG/FT change?

# ---- File Paths ----

file.96 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_96/summary_table.csv"
file.97 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_97/summary_table.csv"
file.98 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_98/summary_table.csv"
file.99 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_99/summary_table.csv"
file.100 <- "../arb-scripts/Marathonas_test/v4_mara/v4_results_100/summary_table.csv"

# ---- Functions ----

import.summary.table <- function(file.name, pident){
  sum.tab <- read.csv(file = file.name)
  sum.tab <- sum.tab[ ,-(2:5)]
  colnames(sum.tab) <- paste(colnames(sum.tab), pident, sep = ".")
  return(sum.tab)
}

sum.96 <- import.summary.table(file.name = file.96, pident = 96)
sum.97 <- import.summary.table(file.name = file.97, pident = 97)
sum.98 <- import.summary.table(file.name = file.98, pident = 98)
sum.99 <- import.summary.table(file.name = file.99, pident = 99)
sum.100 <- import.summary.table(file.name = file.100, pident = 100)

p.compare <- cbind(sum.96[ ,-1], sum.97[ ,-1], sum.98[ ,-1], sum.99[ ,-1], sum.100[ ,-1])
row.names(p.compare) <- sum.96[ ,1]
p.compare <- as.matrix(p.compare)

barplot(p.compare)
