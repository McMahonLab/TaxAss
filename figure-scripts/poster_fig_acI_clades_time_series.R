


# ---- file paths ----

grouped.taxa.folder.path.workflow <- "~/Desktop/Exploring_taxa_in_10yr_TAGs/data/"
library("lubridate")


# ---- functions ----

import.grouped.data <- function(FolderPath){
  if (!exists("grouped")){
    grouped.taxa.folder.path <- FolderPath
    file.names <- list.files(path = grouped.taxa.folder.path)
    grouped <- list(NULL)
    for (t in 1:length(file.names)){
      grouped[[t]] <- read.csv(file = paste(grouped.taxa.folder.path, "/", file.names[t], sep = ""), header = FALSE, stringsAsFactors = FALSE)
      colnames(grouped[[t]]) <- grouped[[t]][1, ]
      grouped[[t]] <- grouped[[t]][-1, ]
      row.names(grouped[[t]]) <- NULL
      tot.col <- ncol(grouped[[t]])
      grouped[[t]][ ,(t + 1):tot.col] <- apply(X = grouped[[t]][ ,(t + 1):tot.col], MARGIN = 2, FUN = as.numeric)
    }
    names(grouped) <- c("kingdom","phylum","class","order","lineage","clade","tribe","OTU")
  }else{
    cat("note: did not re-load \"grouped\" because it already exists.")
  }
  return(grouped)
}

find.name.abund <- function(name.taxa.level, name, abund.cutoff, Grouped){
  grouped <- Grouped
  t <- name.taxa.level
  index <- which(grouped[[t]][ ,t] == name)
  abund <- as.vector(grouped[[t]][index, (t + 1):ncol(grouped[[t]])])
  abund <- unlist(abund)
  return(abund)
}

find.unique.daughter.names <- function(parent.t.level, name, Grouped){
  grouped <- Grouped
  t <- parent.t.level
  index <- which(grouped[[t + 1]][ ,t] == name)
  daughter.names <- unique(grouped[[t + 1]][index, t + 1])
  index <- which(substr(x = daughter.names, start = 1, stop = 12) == "unclassified")
  if (length(index) != 0){
    daughter.names <- daughter.names[-index]
  }
  return(daughter.names)
}

collect.names.and.abunds.into.nested.lists <- function(top.name, top.level, low.level, abund.cutoff, Grouped){
  grouped <- Grouped
  # the taxa list will have an element for each taxa level, this element will be the names list
  taxa.list <- list(NULL)
  daughter.names <- NULL
  for (t in top.level:low.level){
    # the names list will have an element for each name, this element will be an abundance vector
    names.list <- list(NULL)
    if (is.null(daughter.names)){
      daughter.names <- top.name
    }
    low.abund.indeces <- NULL
    for (n in 1:length(daughter.names)){
      names.list[[n]] <- find.name.abund(name.taxa.level = t, name = daughter.names[n], Grouped = grouped)
      names(names.list)[n] <- daughter.names[n]
      names.list[[n]] <- names.list[[n]] / tot.reads * 100
      if (max(names.list[[n]]) < abund.cutoff){
        low.abund.indeces <- c(low.abund.indeces, n)
      }
    }
    if (is.null(low.abund.indeces) != TRUE){
      names.list <- names.list[-(low.abund.indeces)]
    }
    
    taxa.list[[t]] <- names.list
    names(taxa.list)[t] <- paste("level", t, sep = "")
    
    #set up the daughter names for the next run through the loop
    name.compiler <- NULL
    for (n in 1:length(daughter.names)){
      name.compiler <- c(name.compiler, find.unique.daughter.names(parent.t.level = t, name = daughter.names[n], Grouped = grouped))
    }
    daughter.names <- name.compiler
  }
  names(taxa.list)[top.level:low.level] <- names(grouped)[top.level:low.level]
  return(taxa.list)
}

make.polygon <- function(x, y, col, alpha){
  poly.x <- c( min(x), x, max(x))
  poly.y <- c(0, y, 0)
  poly.col <- adjustcolor(col = col, alpha.f = alpha)
  polygon(x = poly.x, y = poly.y, col = poly.col, border = NA)
}



# ---- go go go ----

grouped <- import.grouped.data(FolderPath = grouped.taxa.folder.path.workflow)
tot.reads <- sum(grouped[[8]][ ,9])
sample.dates <- parse_date_time(colnames(grouped[[2]][3:ncol(grouped[[2]])]), orders = "mdy", tz = "Etc/GMT-5")

parent.taxon.name <- "acI"
highest.taxa.level <- 5
lowest.taxa.level <- 6
plot.taxa.levels <- c(5,6)
perc.abund.cutoff <- 1                # if your cutoff leaves nothing in a plotted level it throws an error
taxa.list <- collect.names.and.abunds.into.nested.lists(top.name = parent.taxon.name, top.level = highest.taxa.level, 
                                                        low.level = lowest.taxa.level, abund.cutoff = perc.abund.cutoff, 
                                                        Grouped = grouped)
# ---- ugh ----

a <- taxa.list[[6]][[1]]

b <- taxa.list[[6]][[2]]

c <- taxa.list[[6]][[3]]


plot.title <- expression(bold("Clades Under acI Lineage in Lake Mendota"))
y.label <- expression(bold("Rel. Abundance (% Reads)"))
x.ax <- parse_date_time(x = unique(c(year(sample.dates),2012)), orders = "y", tz = "Etc/GMT-5")
x.ax.lab <- unique(c(year(sample.dates),"2012"))
y.ax <- c(0, 10, 20, 30)
legend.txt <- expression(bold("acI-A"), bold("acI-B"), bold("acI-C"))
file.name <- "~/Dropbox/Trina/8-20-16_ISME16_figures/ME_time_series_acI.png"
x.label <- expression(bold("Year"))
legend.loc <- parse_date_time(x = "2013", orders = "y", tz = "Etc/GMT-5")
x.label.loc <- parse_date_time(x = "7-20-05", orders = "mdy", tz = "Etc/GMT-5")

png(filename = file.name, width = 11.27, height = 5.41, units = "in", res = 300)
par(mar = c(4,3,2,8))
plot(y = a, x = sample.dates, type = "n", axes = F, ann = F)
axis(side = 1, at = x.ax, labels = F, line = -.5, lwd = 3, xpd = T)
mtext(side = 1, text = x.ax.lab, line = 0, at = x.ax, cex = 1.5, las =2)
mtext(text = x.label, side = 1, line = 3, cex = 1.5, at = x.label.loc)
axis(side = 2, at = y.ax, labels = F, line = 0, lwd = 3, xpd = T)
mtext(text = y.ax, side = 2, line = .5, cex = 1.5, at = y.ax)
mtext(text = plot.title, side = 3, line = 0, cex = 2)
mtext(text = y.label, side = 2, line = 1.5, cex = 1.5)
mtext(text = legend.txt, side = 1, line = c(-5,-3,-1), at = legend.loc, cex = 1.5, col = c("purple", "red", "blue"))

lines(y = a, x = sample.dates, col = "purple", lwd = 3)
make.polygon(x = sample.dates, y = a, col = "purple", alpha = .5)
lines(y = b, x = sample.dates, col = "red", lwd = 3)
make.polygon(x = sample.dates, y = b, col = "red", alpha = .3)
lines(y = c, x = sample.dates, col = "blue", lwd = 3)
make.polygon(x = sample.dates, y = c, col = "blue", alpha = .5)
dev.off()
