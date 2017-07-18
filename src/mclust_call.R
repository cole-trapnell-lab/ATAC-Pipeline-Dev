#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("no args")
} else {
  prefix <- args[1]
}

suppressMessages(library(mclust, quietly=TRUE))

celldf <- read.table(paste0(prefix, ".cell_read_counts.txt"))
celltotals <- celldf$V2
cellcall <-  Mclust(data.frame(log10(celltotals)),G=2)
cellfloor <- min(celltotals[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05)])

print(cellfloor)
