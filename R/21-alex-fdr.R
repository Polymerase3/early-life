#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(locfdr)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- file.path("results", "alex_flags")
base_dir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else base_dir_default

summary_dir_default <- file.path(base_dir, "summary")
summary_dir <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else summary_dir_default

out_csv_default <- file.path(summary_dir, "delta_table_all_microbiome_recoded.csv")
out_rds_default <- file.path(summary_dir, "delta_table_all_microbiome_recoded.rds")
out_csv <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else out_csv_default
out_rds <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else out_rds_default

if (!dir.exists(base_dir)) {
  stop("Missing input directory: ", base_dir)
}

delta_paths <- list.files(
  base_dir,
  pattern = "^delta_table\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(delta_paths) == 0L) {
  stop("No delta_table.csv files found under: ", base_dir)
}

extract_group_parts <- function(x) {
  m <- stringr::str_match(x, "^(.*)=([01])$")
  list(
    microbiome = m[, 2],
    group_bin = suppressWarnings(as.integer(m[, 3]))
  )
}

read_one <- function(path) {
  rel <- sub(paste0("^", base_dir, "/?"), "", path)
  rel_parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  timepoint <- if (length(rel_parts) >= 1) rel_parts[[1]] else NA_character_
  comparison_id <- if (length(rel_parts) >= 2) rel_parts[[2]] else NA_character_
  micro_from_path <- sub("_1_vs_0$", "", comparison_id)

  d <- readr::read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet") %>%
    dplyr::select(-dplyr::any_of(c("...1", "X1", "X")))

  if (!all(c("group1", "group2") %in% names(d))) {
    warning("Skipping file without group1/group2: ", path)
    return(NULL)
  }

  g1 <- extract_group_parts(d$group1)
  g2 <- extract_group_parts(d$group2)

  mismatch <- !is.na(g1$microbiome) & !is.na(g2$microbiome) & g1$microbiome != g2$microbiome
  if (any(mismatch)) {
    warning("Mismatched group1/group2 microbiome labels in: ", path)
  }

  d %>%
    dplyr::mutate(
      timepoint = timepoint,
      comparison_id = comparison_id,
      microbiome_tested = dplyr::coalesce(g1$microbiome, g2$microbiome, micro_from_path),
      group1 = g1$group_bin,
      group2 = g2$group_bin
    )
}

delta_all <- purrr::map(delta_paths, read_one) %>%
  purrr::compact() %>%
  dplyr::bind_rows()

if (nrow(delta_all) == 0L) {
  stop("No valid rows were loaded from DELTA tables.")
}

required_cols <- c(
  "microbiome_tested", "timepoint", "rank", "group1", "group2", "design",
  "n_peptides_used", "T_obs", "T_null_mean", "T_null_sd", "T_obs_stand",
  "Z_from_p", "p_perm", "b", "max_delta", "frac_delta_pos", "comparison_id"
)
missing_required <- setdiff(required_cols, names(delta_all))
if (length(missing_required) > 0L) {
  for (nm in missing_required) {
    delta_all[[nm]] <- NA
  }
}

delta_all <- delta_all %>%
  dplyr::mutate(
    group1 = as.integer(.data$group1),
    group2 = as.integer(.data$group2)
  ) %>%
  dplyr::rename(
    D_unstandardized = T_obs,
    D_null_mean = T_null_mean,
    D_null_sd = T_null_sd,
    D_standardized = T_obs_stand
  ) %>%
  dplyr::mutate(
    p_adj_bh = p.adjust(.data$p_perm, method = "BH"),
    .after = "p_perm"
  )

# Efron's empirical null over the whole delta_all table.
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

delta_all$fdr_emp_null <- NA_real_
z_vals <- delta_all$D_standardized
keep <- is.finite(z_vals)
locfdr_plot_path <- file.path(dirname(out_csv), "delta_all_locfdr_plot.png")

if (sum(keep) > 0L) {
  grDevices::png(locfdr_plot_path, width = 1800, height = 1400, res = 200)
  ef <- tryCatch(
    locfdr::locfdr(z_vals[keep], plot = 1, df = 7, nulltype = 1),
    error = function(e) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste("locfdr failed:", conditionMessage(e)))
      warning("locfdr failed on delta_all: ", conditionMessage(e))
      NULL
    }
  )
  grDevices::dev.off()

  if (!is.null(ef) && !is.null(ef$fdr) && length(ef$fdr) == sum(keep)) {
    delta_all$fdr_emp_null[keep] <- ef$fdr
  } else if (!is.null(ef)) {
    warning("locfdr did not return per-observation fdr for delta_all.")
  }
} else {
  warning("No finite D_standardized values in delta_all; skipping locfdr.")
}

delta_all <- delta_all %>%
  dplyr::select(
    "microbiome_tested", "timepoint", "rank",
    "group1", "group2", "design", "n_peptides_used",
    "D_unstandardized", "D_null_mean", "D_null_sd", "D_standardized",
    "Z_from_p", "p_perm", "p_adj_bh", "fdr_emp_null",
    "b", "max_delta", "frac_delta_pos", "comparison_id"
  )

readr::write_csv(delta_all, out_csv)
saveRDS(delta_all, out_rds)

# Collect "interesting" POP plots for nominally significant comparisons.
interesting_dir <- file.path(summary_dir, "interesting_plots")
dir.create(interesting_dir, recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

sig_comparisons <- delta_all %>%
  dplyr::filter(!is.na(.data$p_perm), .data$p_perm < 0.05) %>%
  dplyr::distinct(.data$timepoint, .data$comparison_id, .data$microbiome_tested, .data$p_perm)

copied_n <- 0L
missing_n <- 0L

if (nrow(sig_comparisons) > 0L) {
  for (i in seq_len(nrow(sig_comparisons))) {
    tp <- sig_comparisons$timepoint[[i]]
    cmp <- sig_comparisons$comparison_id[[i]]
    micro <- sig_comparisons$microbiome_tested[[i]]
    pval <- sig_comparisons$p_perm[[i]]

    src_plot <- file.path(
      base_dir, tp, cmp, "POP_framework", "plots", "peptide_id_static.svg"
    )

    dst_name <- paste0(
      sanitize(tp), "__", sanitize(micro), "__", sanitize(cmp),
      "__pperm_", formatC(pval, format = "f", digits = 6),
      "__peptide_id_static.svg"
    )
    dst_plot <- file.path(interesting_dir, dst_name)

    if (file.exists(src_plot)) {
      file.copy(src_plot, dst_plot, overwrite = TRUE)
      copied_n <- copied_n + 1L
    } else {
      missing_n <- missing_n + 1L
    }
  }
}

cat("Input directory:", base_dir, "\n")
cat("Delta files loaded:", length(delta_paths), "\n")
cat("Rows in combined table:", nrow(delta_all), "\n")
cat("Saved CSV:", out_csv, "\n")
cat("Saved RDS:", out_rds, "\n")
cat("Saved locfdr plot:", locfdr_plot_path, "\n")
cat("Interesting plots directory:", interesting_dir, "\n")
cat("Nominally significant comparisons (p_perm < 0.05):", nrow(sig_comparisons), "\n")
cat("Interesting plots copied:", copied_n, " | missing source plots:", missing_n, "\n")
