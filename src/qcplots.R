#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("no args")
} else {
  prefix <- args[1]
}

suppressMessages(library(ggplot2, quietly=TRUE))

path <- gsub("(.*)/.*","\\1",prefix)

insertsdf <- read.table(paste0(prefix, ".cell_read_counts.txt"))

pdf(paste0(path, "/qc_info/insert_sizes.pdf"), height = 3, width = 3)
qplot(insertsdf$V1, geom="density") + labs(x = "Insert size") + theme_bw()
dev.off()

