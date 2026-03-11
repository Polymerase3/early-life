#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

default_in_path <- file.path(
  "results", "results",
  "kid_serum_T2_preterm_yes_vs_kid_serum_T2_preterm_no",
  "DELTA_framework", "delta_table.csv"
)
in_path <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_in_path

if (!file.exists(in_path)) {
  stop("Missing input file: ", in_path)
}

out_dir <- file.path("results", "BH_selection_metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_delta_path <- file.path(out_dir, "preterm_T2_delta_table_bh.csv")
out_summary_path <- file.path(out_dir, "preterm_T2_bh_significance_summary.csv")

delta_df <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")
required_cols <- c("T_obs", "m_eff", "p_perm")
missing_cols <- setdiff(required_cols, names(delta_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "))
}

delta_full <- delta_df %>%
  dplyr::mutate(
    T_obs = as.numeric(.data$T_obs),
    m_eff = as.numeric(.data$m_eff),
    p_perm = as.numeric(.data$p_perm)
  ) %>%
  dplyr::group_by(.data$T_obs, .data$m_eff) %>%
  dplyr::mutate(.test_id = dplyr::cur_group_id()) %>%
  dplyr::ungroup()

delta_unique <- delta_full %>%
  dplyr::group_by(.data$.test_id) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_bh = p.adjust(.data$p_perm, method = "BH"))

delta_out <- delta_full %>%
  dplyr::left_join(
    delta_unique %>% dplyr::select(".test_id", "p_adj_bh"),
    by = ".test_id"
  ) %>%
  dplyr::select(-".test_id")

alpha_levels <- c(0.01, 0.05, 0.10)

summary_row <- tibble::tibble(alpha = alpha_levels) %>%
  dplyr::mutate(
    scope = "row_level",
    n_total = sum(!is.na(delta_out$p_adj_bh)),
    n_significant = vapply(
      alpha,
      function(a) sum(delta_out$p_adj_bh < a, na.rm = TRUE),
      integer(1)
    )
  )

summary_unique <- tibble::tibble(alpha = alpha_levels) %>%
  dplyr::mutate(
    scope = "unique_test_level",
    n_total = sum(!is.na(delta_unique$p_adj_bh)),
    n_significant = vapply(
      alpha,
      function(a) sum(delta_unique$p_adj_bh < a, na.rm = TRUE),
      integer(1)
    )
  )

summary_tbl <- dplyr::bind_rows(summary_row, summary_unique) %>%
  dplyr::mutate(prop_significant = .data$n_significant / .data$n_total) %>%
  dplyr::select("scope", "alpha", "n_total", "n_significant", "prop_significant")

readr::write_csv(delta_out, out_delta_path)
readr::write_csv(summary_tbl, out_summary_path)

cat("Input:", in_path, "\n")
cat("Saved BH-adjusted table:", out_delta_path, "\n")
cat("Saved significance summary:", out_summary_path, "\n\n")
print(summary_tbl, n = nrow(summary_tbl))
