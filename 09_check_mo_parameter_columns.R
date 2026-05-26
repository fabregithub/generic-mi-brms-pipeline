#!/usr/bin/env Rscript
# Inspect monotonic-effect parameter columns in results/parameter_draws.rds

suppressPackageStartupMessages({
  library(tidyverse)
})

f <- "results/parameter_draws.rds"

if (!file.exists(f)) {
  stop("File not found: ", f)
}

draws <- readRDS(f)
nms <- names(draws)

cat("\nColumns matching ^b_|^bsp_|^simo_:\n")
print(nms[grepl("^(b_|bsp_|simo_)", nms)])

cat("\nColumns containing C6yincome, medu, or mo:\n")
print(nms[grepl("C6yincome|medu|mo", nms)])

cat("\nRecommended parameter_draw_regex for mo() models:\n")
cat('analysis_spec$model$parameter_draw_regex <- "^(b_|bsp_|sd_|simo_)"\n')
