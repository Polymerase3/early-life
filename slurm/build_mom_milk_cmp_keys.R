#!/usr/bin/env Rscript

source("R/zzz.R")

is_new_mom_milk_meta_cmp <- function(cmp) {
  length(cmp) == 2 &&
    all(grepl("^mom_milk_T\\d+_.+", cmp))
}

mask <- vapply(comparisons2, is_new_mom_milk_meta_cmp, logical(1))
keys <- vapply(comparisons2[mask], paste, collapse = "|", FUN.VALUE = character(1))

out_path <- "scripts/mom_milk_cmp_keys.txt"
writeLines(keys, out_path)
cat(sprintf("Wrote %d comparison keys to %s\n", length(keys), out_path))
