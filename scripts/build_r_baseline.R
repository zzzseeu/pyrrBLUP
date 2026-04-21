#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rrBLUP)
})

args <- commandArgs(trailingOnly = TRUE)
location_col <- if (length(args) >= 1) args[[1]] else "loc_BeiJ"
out_dir <- if (length(args) >= 2) args[[2]] else "baselines/current"

repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
geno_path <- normalizePath(file.path(repo_root, "..", "..", "project", "genomic_selection_project", "data", "processed", "model_dataset", "agront_ds_202603", "genotype_012.csv"), mustWork = TRUE)
pheno_path <- normalizePath(file.path(repo_root, "..", "..", "project", "genomic_selection_project", "data", "processed", "phenotype", "phenotype_ETN_processed_20260227.csv"), mustWork = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

geno <- read.csv(geno_path, check.names = FALSE)
pheno <- read.csv(pheno_path, check.names = FALSE)

stopifnot("sample" %in% names(geno))
stopifnot(all(c("sample", "value", location_col) %in% names(pheno)))

pheno_loc <- pheno[pheno[[location_col]] == 1, c("sample", "value")]
common <- intersect(pheno_loc$sample, geno$sample)
stopifnot(length(common) > 1)

pheno_aligned <- pheno_loc[match(common, pheno_loc$sample), , drop = FALSE]
geno_aligned <- geno[match(common, geno$sample), , drop = FALSE]
rownames(geno_aligned) <- geno_aligned$sample
M <- as.matrix(geno_aligned[, setdiff(names(geno_aligned), "sample"), drop = FALSE])
A <- A.mat(M)
ms <- mixed.solve(y = pheno_aligned$value, K = A)
kb <- kin.blup(data = pheno_aligned, geno = "sample", pheno = "value", K = A)

write.csv(data.frame(sample = common), file.path(out_dir, "samples.csv"), row.names = FALSE)
write.csv(data.frame(sample = pheno_aligned$sample, value = pheno_aligned$value), file.path(out_dir, "phenotype.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(A), A, check.names = FALSE), file.path(out_dir, "amat.csv"), row.names = FALSE)
write.csv(data.frame(beta = ms$beta), file.path(out_dir, "mixed_solve_beta.csv"), row.names = FALSE)
write.csv(data.frame(sample = names(ms$u), u = as.numeric(ms$u)), file.path(out_dir, "mixed_solve_u.csv"), row.names = FALSE)
write.csv(data.frame(vu = ms$Vu, ve = ms$Ve), file.path(out_dir, "mixed_solve_varcomp.csv"), row.names = FALSE)
write.csv(data.frame(sample = names(kb$pred), pred = as.numeric(kb$pred)), file.path(out_dir, "kin_blup_predictions.csv"), row.names = FALSE)
writeLines(location_col, file.path(out_dir, "location.txt"))
