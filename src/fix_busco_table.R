#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

busco_table <- fread(snakemake@input[[1]],
                     fill = TRUE,
                     header = TRUE,
                     sep = "\t",
                     skip = 2,
                     verbose = TRUE)

head(busco_table)

fwrite(busco_table,
       snakemake@output[[1]])

sessionInfo()
