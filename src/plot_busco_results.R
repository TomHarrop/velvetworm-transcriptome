#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########


busco_files <- snakemake@input[["busco_files"]]
plot_file <- snakemake@output[["plot"]]

# dev
# plot_file <- "test/busco_plot.pdf"
# busco_files <- list.files("output/099_busco",
#                           recursive = TRUE,
#                           pattern = "fixed_full_table.csv",
#                           full.names = TRUE)

########
# MAIN #
########

# read data
names(busco_files) <- sapply(busco_files, function(x)
    unlist(strsplit(x, "/", fixed = TRUE))[[3]])
# print(busco_files)
busco_list <- lapply(busco_files, fread)
busco_data <- rbindlist(busco_list, idcol = "Transcriptome")

# mung
busco_data[, c("Type",
               "Filtering") :=
               tstrsplit(Transcriptome, ".", fixed = TRUE)]

# summarise
busco_counts <- busco_data[, .(busco_count = length(unique(`# Busco id`))),
                           by = .(Type, Filtering, Status)]
busco_counts[, total_buscos := sum(busco_count),
             by = .(Type, Filtering)]
busco_counts[, busco_percent := busco_count * 100 / total_buscos]

# fill blanks
setkey(busco_counts, Type, Filtering, Status)
pd <- busco_counts[CJ(unique(Type),
                      unique(Filtering),
                      unique(Status))]

# order the plot
type_order <- c("raw" = "Raw",
                     "merged" = "Merged")
filt_order <- c("expr" = "Expression",
                "length" = "Length")
status_order <- c("Complete" = "Single copy",
                  "Duplicated" = "Multi copy",
                  "Fragmented" = "Fragmented",
                  "Missing" = "Missing")

pd[, Type := factor(plyr::revalue(Type, type_order),
                        levels = type_order)]
pd[, Filtering := factor(plyr::revalue(Filtering, filt_order),
                    levels = filt_order)]
pd[, Status := factor(plyr::revalue(Status, status_order),
                    levels = status_order)]

# set up colours
fill_colours = viridis::viridis_pal()(5)

# winning assembly completeness
max_pct <- pd[Status == "Single copy", max(busco_percent, na.rm = TRUE)]

# plot
gp <- ggplot(pd, aes(x = Status,
                     y = busco_percent,
                     fill = Status)) +
    coord_flip() +
    ylab("%") + xlab(NULL) +
    scale_fill_manual(values = rev(fill_colours[2:5]),
                      guide = FALSE) +
    facet_grid(Filtering ~ Type) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = max_pct, colour = fill_colours[1])

# write output
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
