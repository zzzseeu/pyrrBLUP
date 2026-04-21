#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
location_cols <- strsplit(args[[1]], ",", fixed = TRUE)[[1]]
runs_root <- if (length(args) >= 2) args[[2]] else "runs"

for (location_col in location_cols) {
  out_dir <- file.path(runs_root, location_col, "r")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2("Rscript", c("scripts/build_r_baseline.R", location_col, out_dir))
  if (status != 0) {
    quit(status = status)
  }
}
