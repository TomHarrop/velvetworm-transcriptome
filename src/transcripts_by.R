#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)

expr_file <- snakemake@input[["expr"]]
quant_sf_file <- snakemake@input[["qf"]][[1]]
gtm_file <- snakemake@input[["gtm"]]

# dev
# expr_file <- "output/040_trinity-abundance/raw/salmon.isoform.TMM.EXPR.matrix"
# quant_sf_file <- "output/040_trinity-abundance/raw/D1/quant.sf"
# gtm_file <- "output/030_trinity/trinity_raw/Trinity.fasta.gene_trans_map"

# use the trinity gene_trans_map with the salmon quant
expr <- fread(expr_file)
quant_sf <- fread(quant_sf_file)
gtm <- fread(gtm_file, header = FALSE, col.names = c("gene", "trans"))

# merge length and isoform
gtm_length <- merge(gtm, quant_sf, by.x = "trans", by.y = "Name")
longest_i <- gtm_length[, .I[which.max(Length)], by = gene][, V1]
longest_trans <- gtm_length[longest_i, unique(trans)]

# merge expr and isoform
expr_sum <- melt(expr, id.vars = "V1")[, .(tmm_sum = sum(value)), by = V1]
gtm_expr <- merge(gtm, expr_sum, by.x = "trans", by.y = "V1")
expr_i <- gtm_expr[, .I[which.max(tmm_sum)], by = gene][, V1]
expr_trans <- gtm_expr[expr_i, unique(trans)]

# save lists
fwrite(list(longest_trans),
       snakemake@output[["ln"]],
       col.names = FALSE)
fwrite(list(expr_trans),
       snakemake@output[["expr"]],
       col.names = FALSE)

sessionInfo()

