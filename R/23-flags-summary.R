#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

summary_dir_default <- file.path("results", "alex_flags", "summary")
in_path_default <- file.path(summary_dir_default, "delta_table_all_microbiome_recoded.csv")

in_path <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else in_path_default
out_dir <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else summary_dir_default

nominal_threshold <- 0.05
bh_threshold      <- 0.10
timepoint_levels  <- c("m6", "m9", "m12", "combined")

if (!file.exists(in_path)) {
  stop("Missing input table: ", in_path)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")

required_cols <- c("timepoint", "p_perm", "p_adj_bh")
missing_cols  <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "))
}

d <- d %>%
  mutate(
    timepoint  = as.character(.data$timepoint),
    p_perm     = as.numeric(.data$p_perm),
    p_adj_bh   = as.numeric(.data$p_adj_bh)
  ) %>%
  filter(.data$timepoint %in% timepoint_levels)

summary_tbl <- d %>%
  group_by(timepoint) %>%
  summarise(
    n_total   = n(),
    n_nominal = sum(!is.na(.data$p_perm) & .data$p_perm < nominal_threshold),
    n_bh      = sum(!is.na(.data$p_adj_bh) & .data$p_adj_bh <= bh_threshold),
    .groups   = "drop"
  ) %>%
  mutate(timepoint = factor(.data$timepoint, levels = timepoint_levels)) %>%
  arrange(.data$timepoint)

cat("\n=== DELTA results summary ===\n")
cat(sprintf("Input:  %s\n", in_path))
cat(sprintf("Nominal threshold: p_perm < %.2f\n", nominal_threshold))
cat(sprintf("BH threshold:      p_adj_bh <= %.2f\n\n", bh_threshold))
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\n")

out_csv <- file.path(out_dir, "flags_summary_by_timepoint.csv")
readr::write_csv(summary_tbl, out_csv)
cat("Summary table written to:", out_csv, "\n")
