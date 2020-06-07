# RRR 6/7/2020
# After fixing conflicts, re-run the reformat_taxonomy_nomenclature.R 
# to check that nothing was reformatted. If it was, dig into why manually

userprefs <- commandArgs(trailingOnly = TRUE)
x <- userprefs[1]
y <- userprefs[2]

# cat("\n\nComment out file paths!!\n\n")
# x <- "../../2020-06-02_update_freshtrain/FT_semicol_noprefix_mothur_unnamed_semicol_noconflicts_mothur.tax"
# y <- "../../2020-06-02_update_freshtrain/FT_semicol_noprefix_mothur_unnamed_semicol_noconflicts_mothur_unnamed.tax"


x <- read.delim(file = x, sep = ";")
y <- read.delim(file = y, sep = ";")

my.check <- all.equal(x,y)

if (my.check == TRUE){
  cat("Files are all equal.\n")
}else{
  cat("Something is different between the files. On your own to figure out why... \n")
  cat(my.check)
}