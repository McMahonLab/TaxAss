# RRR 3-27-18 ----

# SILVA database has extremely inconsistent names for unclassified entries.
# This is a problem because can't group your abundances by taxonomy when the names aren't unique,
# and the whole point of TaxAss and good taxa names is to explore by taxa-groups instead of OTUs!
# (Also, TaxAss needs names to be unique to tally-up unclassified entries, as done in step 15.5.a.)

# How does this script adjust unclassified names?

# 1. This script takes all things meaning "we don't know the name" and changes them to "unnamed". 
     # This looses any info carried by the specific starting word, but I can't figure out if they 
     # mean anything and it saves SOOO much time and confusion to just make everything the same.
     # I chose the "word" unnamed because it has no implications except that there wasn't a name.
     # I was tempted to call everything that couldn't be named Voldemort, but I want deatheaters to cite me too  
  # unknown
  # unnamed
  # uncultured
  # Incertae_Sedis
  # blank
  # p__ etc (what blanks look like in greengenes)
  # unclassified
  # known-name_ge etc (what unknown genus of known-name coarser level looks like in silva)


# 2. This script changes all the created "unnamed" entries to "unnamed.nearest-known-coarser-name"
     # Note that sometimes silva has gaps with no name followed by a fine-resolution name
  # Example:
    # phylum-A    unnamed             unnamed             unnamed             ...   becomes
    # phylum-A    unnamed.phylum-A    unnamed.phylum-A    unnamed.phylum-A    ...
  # Example:
    # phylum-A    unnamed             order-B             unnamed             ...   becomes
    # phylum-A    unnamed.phylum-A    order-B             unnamed.order-B    ...


# How does this script enforce unique names?

# 3. This script finds any names that exist at more than one taxa-level and then 
     # adds to the end a period followed by a the first letter of the taxa-level.
     # The letter is uppercase for general and lowercase letters for FreshTrain. 
     # That way you know capitol C is Class and lowercase c is clade.
  # Example: in silva 132
    # Actinobacteria        Actinobacteria        order-A     ...     becomes
    # Actinobacteria_P      Actinobacteria_C      order-A     ...     
  # Example: in FreshTrain 
    # alfI-A                alfI-A                tribe-A     ...     becomes
    # alfI-A_l              alfI-A_c              tribe-A     ...     

# 4. This script finds any daughter names that have multiple parent-names and then 
    # adds to the end a period followed by the parent name that it belongs to. 
    # This could potentially become unweildy, but hopefully this hardly ever happens!
  # Example:
    # phylum-A      class-A         order-C               family-D          ...   
    # phylum-A      class-B         order-C               family-E          ...
    #                                                                         becomes
    # phylum-A      class-A         order-C.class-A       family-D          ...   
    # phylum-A      class-B         order-C.class-B       family-E          ...


# Why does that fix it?
# 1. Each daughter name will end up with unique parent names.
# 2. The nomenclature is consistent so TaxAss can still tally up the unclassifieds.
# By calling anything unclassified that begins with (uncultured)|(unnamed)|(unknown)|(unclassified)  

# What doesn't it fix?
# 1. I don't know if "uncultured" and "unknown" should actually be distinguished or not. 
# Now they are not being distinguished, but possibly that's incorrect and I'm losing information.
# 2. There are still daughter names below unknown names. 
# This somewhat breaks logic, but with the parent names made unique hopefully it will not break scripts.

# Note: this works on greengenes too, but greengenes is more uniform. 
# all it changes is k__ p__ c__ etc that are blank afterwards.
# Also, greengenes reads in with 9 columns, so 9th column must be manually deleted after import.


# ---- file paths/commandline input ----

userprefs <- commandArgs(trailingOnly = TRUE)
mothur.formatted.tax <- userprefs[1]
new.tax.file <- userprefs[2]
if (!is.na(userprefs[3])){
  file.type <- userprefs[3]
}else{
  file.type <- "General"
}
 

# # Troubleshoot with local paths:
# cat("\n\nyou forgot to comment out file paths!!!\n\n")
# # mothur.formatted.tax <- "../../StartFiles-Databases/Silva.nr_v138/silva.nr_v138.tax"
# # mothur.formatted.tax <- "../../StartFiles-Databases/withGG/general.taxonomy"
# mothur.formatted.tax <- "../../2020-06-02_update_freshtrain/FT_noconflicts_mothur.taxonomy"
# new.tax.file <- "../../2020-06-02_update_freshtrain/FT_noconflicts_unnamed.taxonomy"
# file.type <- "FreshTrain"
# # file.type <- "General"



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
  non.uniq.names <- grep(pattern = "unnamed", x = non.uniq.names, ignore.case = F, value = T, invert = T)
  if(length(non.uniq.names) == 0){
    cat("No duplicated names across taxa-levels.\n")
  }
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

uniform.unnamed <- function(tax){
  # anything with unnamed in the name 
  index <- grep(pattern = 'unnamed', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
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
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
  return(tax)
}

uniform.unknown <- function(tax){
  # anything with unknown in the name, esp like Unknown_Family
  index <- grep(pattern = 'unknown', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
  return(tax)
}

uniform.uncertain <- function(tax){
  # anything with Incertae_Sedis in the name, esp like Unknown_Family
  index <- grep(pattern = 'Incertae_Sedis', x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
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
  cat("changing",length(index), "\n", unique(tax[index]), "\nto say \"unnamed\"\n")
  tax[index] <- "unnamed"
  return(tax)
}

make.degenerates.unique <- function(tax, voldemort){
  # tax is the taxonomy table without seqid column
  # voldemort is whatever text you're making unique (b/c it cannot be named... voldemort...)
  # originally this script maintained all the different words, but now they all just become unnamed
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
  # Not used anymore, got too complicated so just all voldemorts "unnamed"
  
  # since there's multiple way things must not be named, make sure each thing just has one not-named label
  # these are the only 2 I saw in silva 132, may need to re-check in future versions for other combos
  
  # # Check for additional ones to fix with new versions:
  # grep(pattern = "uncertain.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncertain.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncertain.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "uncultured.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncultured.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "uncultured.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "unknown.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unknown.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unknown.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # 
  # grep(pattern = "unnamed.uncertain", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unnamed.unknown", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  # grep(pattern = "unnamed.uncultured", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
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
  
  index <- grep(pattern = "uncertain.unnamed", x = tax, value = FALSE, invert = FALSE, ignore.case = TRUE)
  cat("changing", length(index), " \"uncertain.unnamed\" to lead with just \"uncertain\"\n")
  tax <- gsub(pattern = "uncertain.unnamed", replacement = "uncertain", x = tax, ignore.case = TRUE)
  
  return(tax)
}

make.taxon.names.unique.within.level <- function(tax){
  # # test data
  # tax <- matrix(data = c(letters[1:3], letters[1:3], letters[c(1,26,3)], letters[c(1,2,4)], letters[5:7], letters[c(5,6,8)], letters[c(24,6,16)], letters[9:11], letters[9:11], letters[c(9,10,12)], letters[c(9,13,14)], letters[c(9,13,14)], letters[c(9,13,15)] ), ncol = 3, byrow = T)
  
  for (t in 2:ncol(tax)){
    temp <- tax[ ,1:t]
    temp <- unique(temp)
    index <- duplicated(temp[ ,t])
    dups <- unique(temp[index,t])
    cat("at level ", t, "\n")
    if (length(dups) > 0){
      for (d in 1:length(dups)){
        index <- which(tax[ ,t] == dups[d])
        cat(unique(tax[index,t]), " is not unique and is being replaced by:  ", unique(paste0(tax[index,t], ".", tax[index,t-1])), "\n")
        tax[index,t] <- paste0(tax[index,t], ".", tax[index,t-1])
      }
    }
  }
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

stupid <- import.mothur.formatted.tax(filepath = mothur.formatted.tax)

stupid[ ,-1] <- fill.blanks.with.unnamed(tax = stupid[ ,-1])

stupid[ ,-1] <- make.taxon.names.unique.across.levels(tax = stupid[ ,-1], abbr = level.abbreviations)

stupid[ ,-1] <- uniform.unnamed(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.uncultured(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.unknown(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.uncertain(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.unnamed.silva(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.unnamed.gg(tax = stupid[ ,-1])
stupid[ ,-1] <- uniform.unclassified(tax = stupid[ ,-1])

stupid[ ,-1] <- make.degenerates.unique(tax = stupid[ ,-1], voldemort = "unnamed")

stupid[ ,-1] <- make.taxon.names.unique.within.level(tax = stupid[ ,-1])

less.stupid <- revert.to.mothur.format(silva = stupid)

write.table(x = less.stupid, file = new.tax.file, quote = F, col.names = FALSE, row.names = FALSE, sep = ";")
cat("Made file: ", new.tax.file, "\n")
