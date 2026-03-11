#!/usr/bin/env Rscript

# FDR computation across all comparisons under a results root:
# - Discover comparison IDs in results/results (or CLI arg)
# - Load DELTA_framework/delta_table.csv for each comparison
# - Compute BH q-values (p_perm) and empirical-null FDR (locfdr on T_obs_stand)
# - Save per-comparison table as DELTA_framework/delta_table_fdr.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(locfdr)
})

bh_cutoff <- 0.1
emp_fdr_cutoff <- 0.2

args <- commandArgs(trailingOnly = TRUE)
base_dir_default <- file.path("results", "results")
base_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else base_dir_default

if (!dir.exists(base_dir)) {
  stop("Missing results directory: ", base_dir)
}

all_comparison_ids <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
has_delta_table <- file.exists(file.path(base_dir, all_comparison_ids, "DELTA_framework", "delta_table.csv"))
comparison_ids <- sort(all_comparison_ids[has_delta_table])

if (length(comparison_ids) == 0) {
  stop("No DELTA tables found under: ", base_dir)
}

missing_delta <- sort(all_comparison_ids[!has_delta_table])
if (length(missing_delta) > 0) {
  warning(
    "Skipping ", length(missing_delta),
    " comparison directories without DELTA_framework/delta_table.csv"
  )
}

cat("Results directory:", base_dir, "\n")
cat("Comparisons to process:", length(comparison_ids), "\n")

run_one <- function(comparison) {
  delta_dir  <- file.path(base_dir, comparison, "DELTA_framework")
  delta_path <- file.path(delta_dir, "delta_table.csv")
  out_path   <- file.path(delta_dir, "delta_table_fdr.csv")

  if (!file.exists(delta_path)) {
    warning("Missing DELTA table for ", comparison, ": ", delta_path)
    return(NULL)
  }

  delta_df <- readr::read_csv(delta_path, show_col_types = FALSE)
  required_cols <- c("T_obs", "m_eff", "p_perm", "T_obs_stand")
  missing_cols <- setdiff(required_cols, names(delta_df))
  if (length(missing_cols) > 0) {
    warning(
      "Skipping ", comparison, "; missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }

  # Match 11-normal-meta.R: BH and locfdr are computed on distinct(T_obs, m_eff).
  cmp_df_uniq <- delta_df %>%
    dplyr::distinct(T_obs, m_eff, .keep_all = TRUE) %>%
    dplyr::mutate(p_adj_bh = p.adjust(p_perm, method = "BH"))

  z_vals <- cmp_df_uniq$T_obs_stand
  keep <- is.finite(z_vals)
  cmp_df_uniq$fdr_emp_null <- NA_real_

  if (sum(keep) > 0) {
    ef <- tryCatch(
      locfdr::locfdr(z_vals[keep], plot = 0, df = 7, nulltype = 1),
      error = function(e) {
        warning("locfdr failed for ", comparison, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(ef) && !is.null(ef$fdr) && length(ef$fdr) == sum(keep)) {
      cmp_df_uniq$fdr_emp_null[keep] <- ef$fdr
    } else if (!is.null(ef)) {
      warning("locfdr did not return per-observation fdr for ", comparison)
    }
  } else {
    warning("No finite T_obs_stand values for ", comparison)
  }

  fdr_by_test <- cmp_df_uniq %>%
    dplyr::select(T_obs, m_eff, p_adj_bh, fdr_emp_null)

  cmp_df_out <- delta_df %>%
    dplyr::select(-dplyr::any_of(c("p_adj_bh", "fdr_emp_null"))) %>%
    dplyr::left_join(fdr_by_test, by = c("T_obs", "m_eff"))

  readr::write_csv(cmp_df_out, out_path)

  bh_hits <- sum(cmp_df_uniq$p_adj_bh < bh_cutoff, na.rm = TRUE)
  emp_hits <- sum(cmp_df_uniq$fdr_emp_null < emp_fdr_cutoff, na.rm = TRUE)
  cat(
    sprintf(
      "Saved: %s | source: %s | rows: %d | unique_tests: %d | BH q<%.2f: %d | locfdr<%.2f: %d\n",
      out_path, base_dir, nrow(cmp_df_out), nrow(cmp_df_uniq), bh_cutoff, bh_hits, emp_fdr_cutoff, emp_hits
    )
  )

  invisible(cmp_df_out)
}

results <- lapply(comparison_ids, run_one)
names(results) <- comparison_ids

n_ok <- sum(!vapply(results, is.null, logical(1)))
n_fail <- length(results) - n_ok
cat("\nDone. Processed:", n_ok, " | Failed:", n_fail, "\n")

# Build a joined BH-selection table for metadata comparisons requested by the report workflow.
metadata_pairs <- list(
  c("kid_serum_T2_siblings", "kid_serum_T2_no_siblings"),
  c("kid_serum_T6_siblings", "kid_serum_T6_no_siblings"),
  c("kid_serum_T8_siblings", "kid_serum_T8_no_siblings"),
  c("kid_serum_T2_delmode_VG", "kid_serum_T2_delmode_CS"),
  c("kid_serum_T6_delmode_VG", "kid_serum_T6_delmode_CS"),
  c("kid_serum_T8_delmode_VG", "kid_serum_T8_delmode_CS"),
  c("kid_serum_T2_delplace_home", "kid_serum_T2_delplace_hospital"),
  c("kid_serum_T6_delplace_home", "kid_serum_T6_delplace_hospital"),
  c("kid_serum_T8_delplace_home", "kid_serum_T8_delplace_hospital"),
  c("kid_serum_T2_feeding_BF", "kid_serum_T2_feeding_MF"),
  c("kid_serum_T6_feeding_BF", "kid_serum_T6_feeding_MF"),
  c("kid_serum_T8_feeding_BF", "kid_serum_T8_feeding_MF"),
  c("kid_serum_T2_preterm_yes", "kid_serum_T2_preterm_no"),
  c("kid_serum_T6_preterm_yes", "kid_serum_T6_preterm_no"),
  c("kid_serum_T8_preterm_yes", "kid_serum_T8_preterm_no"),
  c("kid_serum_T2_CDrisk_yes", "kid_serum_T2_CDrisk_no"),
  c("kid_serum_T6_CDrisk_yes", "kid_serum_T6_CDrisk_no"),
  c("kid_serum_T8_CDrisk_yes", "kid_serum_T8_CDrisk_no"),
  c("kid_serum_T2_lockdown_before", "kid_serum_T2_lockdown_after"),
  c("kid_serum_T6_lockdown_before", "kid_serum_T6_lockdown_after"),
  c("kid_serum_T8_lockdown_before", "kid_serum_T8_lockdown_after"),
  c("kid_serum_T2_smoking_yes", "kid_serum_T2_smoking_no"),
  c("kid_serum_T6_smoking_yes", "kid_serum_T6_smoking_no"),
  c("kid_serum_T8_smoking_yes", "kid_serum_T8_smoking_no"),
  c("kid_serum_T2_sex_male", "kid_serum_T2_sex_female"),
  c("kid_serum_T6_sex_male", "kid_serum_T6_sex_female"),
  c("kid_serum_T8_sex_male", "kid_serum_T8_sex_female")
)

build_comparison_id <- function(pair) {
  paste0(pair[[1]], "_vs_", pair[[2]])
}

load_selected_fdr <- function(pair) {
  lhs <- pair[[1]]
  rhs <- pair[[2]]
  cmp_id <- build_comparison_id(pair)
  fdr_path <- file.path(base_dir, cmp_id, "DELTA_framework", "delta_table_fdr.csv")

  if (!file.exists(fdr_path)) {
    warning("Missing delta_table_fdr.csv for selected comparison: ", cmp_id)
    return(NULL)
  }

  cmp_df <- readr::read_csv(fdr_path, show_col_types = FALSE, name_repair = "unique_quiet")
  if (!("p_adj_bh" %in% names(cmp_df))) {
    warning("Missing p_adj_bh in selected comparison: ", cmp_id)
    return(NULL)
  }
  if (!("fdr_emp_null" %in% names(cmp_df))) {
    warning("Missing fdr_emp_null in selected comparison: ", cmp_id)
    return(NULL)
  }
  if (!("p_perm" %in% names(cmp_df))) {
    warning("Missing p_perm in selected comparison: ", cmp_id)
    return(NULL)
  }

  cmp_df %>%
    dplyr::distinct(T_obs, m_eff, .keep_all = TRUE) %>%
    dplyr::mutate(
      comparison_id = cmp_id,
      group_left = lhs,
      group_right = rhs
    )
}

selected_joined <- dplyr::bind_rows(lapply(metadata_pairs, load_selected_fdr))
selected_bh <- selected_joined %>%
  dplyr::filter(!is.na(.data$p_perm), .data$p_perm < 0.05, !is.na(.data$p_adj_bh), .data$p_adj_bh < bh_cutoff)
selected_fdr20 <- selected_joined %>%
  dplyr::filter(!is.na(.data$p_perm), .data$p_perm < 0.05, !is.na(.data$fdr_emp_null), .data$fdr_emp_null < emp_fdr_cutoff)

bh_out_dir <- file.path("results", "BH_selection_metadata")
dir.create(bh_out_dir, recursive = TRUE, showWarnings = FALSE)
bh_out_path <- file.path(bh_out_dir, "BH_selection_metadata.csv")
readr::write_csv(selected_bh, bh_out_path)
fdr20_out_path <- file.path(bh_out_dir, "FDR_selection_metadata_fdr_lt_0.2.csv")
readr::write_csv(selected_fdr20, fdr20_out_path)

cat(
  "Saved BH-selected metadata table:", bh_out_path,
  " | rows:", nrow(selected_bh),
  " | comparisons_requested:", length(metadata_pairs),
  " | filter: p_perm<0.05 & p_adj_bh<", bh_cutoff,
  "\n"
)
cat(
  "Saved empirical-FDR-selected metadata table:", fdr20_out_path,
  " | rows:", nrow(selected_fdr20),
  " | filter: p_perm<0.05 & fdr_emp_null<", emp_fdr_cutoff,
  "\n"
)

# Build one global DELTA table by stacking all per-comparison delta_table_fdr.csv files.
fdr_paths <- file.path(base_dir, comparison_ids, "DELTA_framework", "delta_table_fdr.csv")
has_fdr <- file.exists(fdr_paths)

if (any(!has_fdr)) {
  warning(
    "Skipping ", sum(!has_fdr),
    " comparison directories without DELTA_framework/delta_table_fdr.csv"
  )
}

global_fdr_list <- lapply(which(has_fdr), function(i) {
  cmp_id <- comparison_ids[[i]]
  fdr_path <- fdr_paths[[i]]

  cmp_df <- readr::read_csv(fdr_path, show_col_types = FALSE, name_repair = "unique_quiet")
  cmp_df %>%
    dplyr::mutate(source_folder = file.path(base_dir, cmp_id), .before = 1)
})

global_fdr <- dplyr::bind_rows(global_fdr_list)
global_out_path <- file.path(base_dir, "delta_table.csv")
readr::write_csv(global_fdr, global_out_path)

cat(
  "Saved global stacked DELTA table:", global_out_path,
  " | rows:", nrow(global_fdr),
  " | source_folders:", sum(has_fdr),
  "\n"
)

