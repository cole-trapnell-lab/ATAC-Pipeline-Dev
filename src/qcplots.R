#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("no args")
} else {
  dir <- args[1]
  prefix <- args[2]
}

suppressMessages(library(ggplot2, quietly=TRUE))

#path <- gsub("(.*)/.*","\\1",prefix)

insertsdf <- read.table(paste0(dir, "/", prefix, ".insertsize.txt"))

pdf(paste0(dir, "/insert_sizes.pdf"), height = 3, width = 3)
qplot(insertsdf$V1, geom="density") + labs(x = "Insert size") + theme_bw()
graphics.off()

