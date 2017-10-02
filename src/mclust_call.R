#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("no args")
} else {
  prefix <- args[1]
}

suppressMessages(library(mclust, quietly=TRUE))
suppressMessages(library(ggplot2, quietly=TRUE))

path <- gsub("(.*)/.*","\\1",prefix)

celldf <- read.table(paste0(prefix, ".cell_read_counts.txt"))
celltotals <- celldf$V2
cellcall <-  Mclust(data.frame(log10(celltotals)),G=2)
cellfloor <- min(celltotals[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05)])

pdf(paste0(path, "/qc_info/read_count_density.pdf"), height = 3, width = 3)
qplot(log10(celltotals), geom="density") + geom_vline(xintercept = log10(cellfloor)) + labs(x = "Log10 of cell read count") + theme_bw()
dev.off()

print(cellfloor)
