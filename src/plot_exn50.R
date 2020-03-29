#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(ggplot2)

exn50_file <- snakemake@input[["exn50"]]
plot_file <- snakemake@output[["plot"]]

exn50 <- fread(exn50_file)

gt <- parse(text = paste0(
  'italic("Ex")[90] * italic("N")[50] ==',
  exn50[Ex == 90, ExN50]))
st <- parse(text = paste0(
  'italic(N)[50] ~ "of the" ~ ',
  exn50[Ex == 90, num_transcripts],
  ' ~ "most abundant transcripts that account for 90% of total reads"'
))

vd <- viridisLite::viridis(3)
gp <- ggplot(exn50, aes(x = Ex, y = ExN50)) +
  theme_grey(base_size = 12) +
  geom_path(colour = vd[[2]], size = 1) +
  geom_point(colour = vd[[1]], size = 3, shape = 16) +
  ggtitle(gt,
          subtitle = st) +
  ylab(expression(italic("ExN")[50])) +
  xlab(expression(italic("Ex")))

# write output
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
