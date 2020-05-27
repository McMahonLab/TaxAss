# RRR 3-27-18 ----

# SILVA database has extremely inconsistent names for unclassified entries.
# This is a problem.

# Why?
# 1. Names need to be unique for downstream analysis.
# Example: 
  # phylum-A;nonunique-class;order-A;...
  # phylum-B;nonunique-class;order-B;...
# Analysis scripts would group order-A and order-B as belonging to the same class, 
# when really those are different classes because they belong to different phylums.
# 2. TaxAss needs names to be unique to tally-up unclassified entries.
# Example: this is done for Figure 2 in the manuscript, in step 15.5.a.

# How does this script fix it?

# 1. This script takes all things meaning "we don't know the name" and changes them to a few versions of "we don't know the name". 
# Example: unknown
  # phylum-A;Unknown_Class;order-A;... becomes
  # phylum-A;unknown;order-A;...
# Example: unnamed
  # phylum-A;phylum-A_cl;order-A;... becomes
  # phylum-A;unnamed;order-A;...
# Example: uncultured
  # phylum-A;uncultured;order-A;... becomes
  # phylum-A;uncultured;order-A;...
# Example: uncertain
  # phylum-A;Incertae_Sedis_phylum-A;order-A;... becomes
  # phylum-A;uncertain;order-A;...
# Example: blank
  # order-A;family-A;genus-A;; becomes
  # order-A;family-A;genus-A;unnamed;
# Example: unclassified
  # there are no unclassifieds in silva right now, but anything with unclassified in it will be changed to simply "unclassified"


# 2. This script makes all the unclassified entries unique based on their upper-level taxonomy.
# Example:
  # phylum-A;uncultured;order-A;... becomes
  # phylum-A;uncultured.phylum-A_C;order-A;...
# Example:
  # phylum-A;uncultured;uncultured;... becomes
# phylum-A;uncultured.phylum-A_C;uncultured.phylum-A_O;...

# 3. This script makes all repeated names unique based on their taxon-level.
# Example: Actinobacteria is a phylum name and a class name
  # Actinobacteria;Actinobacteria;order-A;... becomes
  # Actinobacteria_P;Actinobacteria_C;order-A;...

# Why does that fix it?
# 1. Each daughter name will end up with unique parent names.
# 2. The nomenclature is consistent so TaxAss can still tally up the unclassifieds.
# By calling anything unclassified that begins with (uncultured)|(unnamed)|(unknown)|(unclassified)  

# What doesn't it fix?
# 1. I don't know if "uncultured" and "unknown" should actually be distinguished or not. 
# Now they are being distinguished, but possibly that's incorrect and would be simpler to group them all.
# 2. There are still daughter names below unknown names. 
# This somewhat breaks logic, but with the parent names made unique hopefully it will not break scripts.

# Note: this works on greengenes too, but greengenes is more uniform. 
# all it changes is k__ p__ c__ etc that are blank afterwards.
# Also, greengenes reads in with 9 columns, so 9th column must be manually deleted after import.

# Note: when taxon-level identifiers are added, uppercase letters are used for general taxonomy, 
# and lowercase letters are used for FreshTrain taxonomy. That way you know capitol C is Class 
# and lowercase c is clade.  

# ---- file paths/commandline input ----

userprefs <- commandArgs(trailingOnly = TRUE)
mothur.formatted.silva <- userprefs[1]
new.silva.file <- userprefs[2]
if (!is.na(userprefs[3])){
  file.type <- userprefs[3]
}else{
  file.type <- "General"
}
 

# Troubleshoot with local paths:
cat("crap forgot to comment out file paths!!!")
mothur.formatted.silva <- "../../StartFiles-Databases/Silva.nr_v132/silva.nr_v132.tax"
mothur.formatted.silva <- "../../StartFiles-Databases/withGG/general.taxonomy"
mothur.formatted.silva <- "../../StartFiles-Databases/with132/custom.taxonomy"
new.silva.file <- "~/Desktop/cleansilva.taxonomy"
file.type <- "FreshTrain"
file.type <- "General"



if(file.type == "FreshTrain"){
  level.abbreviations <- c("K","P","C","O","l","c","t")
}else{
  level.abbreviations <- c("K","P","C","O","F","G","S")
}


# ---- functions ----

import.mothur.formatted.tax <- function(filepath){
  # seqID delimited from taxonomy by tab, taxa levels delimited by semicolon
  numcol <- max(count.fields(filepath, sep=";"))
  silva <- read.table(file = filepath, sep=";", fill=T, colClasses = "character", col.names=1:numcol)
  silva <- as.matrix(silva)
  kingdom <- gsub(pattern = ".*\t", replacement = "", x = silva[ ,1])
  seqid <- gsub(pattern = "\t.*$", replacement = "", x = silva[ ,1])
  silva <- silva[ ,-1]
  silva <- cbind(seqid, kingdom, silva)
  if(ncol(silva) == 8){
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported with a species column, the unique species names are: ", unique(silva[ ,8]),"\n")
  }else if(ncol(silva) == 7){
    silva <- cbind(silva, "unnamed")
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported without a species column, added one filled with \"unnamed\"\n")
  }else if(ncol(silva) == 9){
    extracol <- unique(silva[ ,9])
    silva <- silva[ ,1:8]
    colnames(silva) <- c("seqid","kingdom","phylum","class","order","family","genus","species")
    cat("database imported with an extra 9th column containing:", extracol, ". This column was removed.\n")
  }else{
    cat("something is wrong with the database import, it has ", ncol(silva), "columns in it.\n")
  }
  return(silva)
}

make.taxon.names.unique.across.levels <- function(tax, abbr){
  all.names <- NULL
  for (t in 1:ncol(tax)){
    level.names <- unique(tax[ ,t])
    all.names <- c(all.names, level.names)
  }
  non.uniq.names <- all.names[duplicated(all.names)]
  non.uniq.names <- unique(non.uniq.names)
  for (t in 1:ncol(tax)){
    for (n in 1:length(non.uniq.names)){
      index <- which(tax[ ,t] == non.uniq.names[n])
      if (length(index) > 0){
        new.name <- paste0(non.uniq.names[n], ".", abbr[t])
        tax[index,t] <- new.name
        cat("Replacing  ",  non.uniq.names[n], "  with  ", new.name, "\n")
      }
    }
  }
  return(tax)
}

fill.blanks.with.unnamed <- function(tax){
  # tax is only the taxonomy columns, seqID not included
  index <- which(tax == "")
  tax[index] <- "unnamed"
  cat("changing", length(index), "blanks to say \"unnamed\"\n")
  return(tax)
}

uniform.uncultured <- function(tax){
  # anything with uncultured in the name 
  index <- grep(pattern = 'uncultured', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"uncultured\"\n")
  tax[index] <- "uncultured"
  return(tax)
}

uniform.unknown <- function(tax){
  # anything with unknown in the name, esp like Unknown_Family
  index <- grep(pattern = 'unknown', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unknown\"\n")
  tax[index] <- "unknown"
  return(tax)
}

uniform.uncertain <- function(tax){
  # anything with Incertae_Sedis in the name, esp like Unknown_Family
  index <- grep(pattern = 'Incertae_Sedis', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"uncertain\"\n")
  tax[index] <- "uncertain"
  return(tax)
}

uniform.unnamed.silva <- function(tax){
  # anything in front, _fa _ge with nothing after at end
  index <- grep(pattern = '.*_((ph)|(cl)|(or)|(fa)|(ge))$', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  if (length(index) > 20){
    cat("changing",length(index), "names like \n", unique(tax[index])[1:20], "\nto say \"unnamed\"\n")
    cat("It might seem like some of these are real names, but if you search for them in silva you'll see they are just listing the name from an upper level. So don't worry. They'll get that upper level name back, just with unnamed as the prefix instead.\n")
  }else{
    cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  }
  tax[index] <- "unnamed"
  return(tax)
}

uniform.unnamed.gg <- function(tax){
  # p__ c__ at front with anything after (this is for greengenes)
  index <- grep(pattern = '.{1}__$', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
  return(tax)
}

uniform.unclassified <- function(tax){
  # anything with unclassified in the name
  index <- grep(pattern = 'unclassified', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unclassified\"\n")
  tax[index] <- "unclassified"
  return(tax)
}

make.degenerates.unique <- function(tax, voldemort){
  # tax is the taxonomy table without seqid column
  # voldemort is whatever text you're making unique (b/c it cannot be named... voldemort...)
  for (t in 1:ncol(tax)){
    index <- which(tax[ ,t] == voldemort)
    tax[index,t] <- paste(tax[index,t], tax[index,t - 1], sep = ".")
  }
  
  remove.extra.u <- function(x){
    dot.voldemort <- paste0(".", voldemort)
    x <- gsub(pattern = dot.voldemort, replacement = "", x = x)
    return(x)
  }
  tax <- apply(X = tax, MARGIN = 2, FUN = remove.extra.u)
  
  return(tax)
}

remove.extra.voldemorts <- function(tax){
  # since there's multiple way things must not be named, make sure each thing just has one not-named label
  # these are the only 2 I saw in silva 132, may need to re-check in future versions for other combos
  
  # # Check for additional ones to fix with new versions:
  # grep(pattern = "uncertain.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncertain.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncertain.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncertain.unclassified", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "uncultured.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncultured.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncultured.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncultured.unclassified", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "unknown.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unknown.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unknown.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unknown.unclassified", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "unnamed.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unnamed.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unnamed.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unnamed.unclassified", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "unclassified", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  
  index <- grep(pattern = "unnamed.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing", length(index), " \"unnamed.unknown\" to lead with just \"unknown\"\n")
  tax <- gsub(pattern = "unnamed.unknown", replacement = "unknown", x = tax, ignore.case = TRUE)
  
  index <- grep(pattern = "unnamed.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing", length(index), "\"unnamed.uncultured\" to lead with just \"uncultured\"\n")
  tax <- gsub(pattern = "unnamed.uncultured", replacement = "uncultured", x = tax, ignore.case = TRUE)
  
  index <- grep(pattern = "unnamed.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing", length(index), " \"unnamed.uncertain\" to lead with just \"uncertain\"\n")
  tax <- gsub(pattern = "unnamed.uncertain", replacement = "uncertain", x = tax, ignore.case = TRUE)
  
  return(tax)
}

revert.to.mothur.format <- function(silva){
  seqid.kingdom <- paste(silva[ ,1], silva[ ,2], sep = "\t")
  silva <- silva[ ,-(1:2)]
  # also add a blank column to get semicolons at the end of the line also
  silva <- cbind(seqid.kingdom, silva, "")
  return(silva)
}

# ---- go ----

stupid.silva <- import.mothur.formatted.tax(filepath = mothur.formatted.silva)

stupid.silva[ ,-1] <- make.taxon.names.unique.across.levels(tax = stupid.silva[ ,-1], abbr = level.abbreviations)

stupid.silva[ ,-1] <- fill.blanks.with.unnamed(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.uncultured(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unknown(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.uncertain(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unnamed.silva(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unnamed.gg(tax = stupid.silva[ ,-1])
stupid.silva[ ,-1] <- uniform.unclassified(tax = stupid.silva[ ,-1])


stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "uncultured")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unknown")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "uncertain")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unnamed")
stupid.silva[ ,-1] <- make.degenerates.unique(tax = stupid.silva[ ,-1], voldemort = "unclassified")
stupid.silva[ ,-1] <- remove.extra.voldemorts(tax = stupid.silva[ ,-1])

uniform.silva <- revert.to.mothur.format(silva = stupid.silva)

write.table(x = uniform.silva, file = new.silva.file, quote = F, col.names = FALSE, row.names = FALSE, sep = ";")
cat("Made file: ", new.silva.file, "\n")
