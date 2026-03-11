#!/usr/bin/env Rscript

# Statistical summary for selected non-milk contrasts
# - Load DELTA_framework/delta_table.csv for selected comparisons under results/results
# - Compute BH q-values (p_perm) for reference
# - Compute empirical-null FDR via locfdr on standardized test stats
# - Build the final hit table from raw significance (p_perm < 0.05)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(locfdr)
  library(png)
  library(grid)
})

base_dir <- file.path("results", "results")
emp_fdr_cutoff <- 0.2
bh_cutoff <- 0.1
plot_dir <- file.path("results", "other_plots", "locfdr_normal")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Requested comparison pairs
comparisons <- list(
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

# Extended comparisons for analyses: add mom_kid_serum_T2_* fallbacks alongside kid_serum_T2_*.
kid_t2_mask <- vapply(comparisons, function(x) any(grepl("^kid_serum_T2", x)), logical(1))
comparisons2 <- comparisons
if (any(kid_t2_mask)) {
  comparisons2_extra <- lapply(
    comparisons[kid_t2_mask],
    function(x) sub("^kid_serum_T2", "mom_kid_serum_T2", x)
  )
  drop_mask <- vapply(
    comparisons2_extra,
    function(x) identical(x, c("mom_serum_T2", "mom_kid_serum_T2")),
    logical(1)
  )
  comparisons2_extra <- comparisons2_extra[!drop_mask]
  comparisons2 <- c(comparisons, comparisons2_extra)
}

pair_to_dir <- function(x) paste(x, collapse = "_vs_")
requested_comparisons <- unique(vapply(comparisons2, pair_to_dir, character(1)))

# Keep only requested comparisons that exist on disk
all_delta_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
all_delta_dirs <- all_delta_dirs[
  file.exists(file.path(all_delta_dirs, "DELTA_framework", "delta_table.csv"))
]
available_comparisons <- basename(all_delta_dirs)

missing_comparisons <- setdiff(requested_comparisons, available_comparisons)
if (length(missing_comparisons) > 0) {
  warning(
    "Requested comparisons not found under ", base_dir, ": ",
    paste(missing_comparisons, collapse = ", ")
  )
}

comparisons_to_run <- requested_comparisons[requested_comparisons %in% available_comparisons]

if (length(comparisons_to_run) == 0) {
  stop("None of the requested DELTA tables were found under ", base_dir)
}

run_one <- function(comparison) {
  delta_path <- file.path(base_dir, comparison, "DELTA_framework", "delta_table.csv")

  if (!file.exists(delta_path)) {
    stop("Missing DELTA_framework/delta_table.csv for ", comparison, " at ", delta_path)
  }

  delta_df <- readr::read_csv(delta_path, show_col_types = FALSE)

  # Drop duplicates on test/m_eff to mimic the quick-look logic
  cmp_df <- delta_df |>
    distinct(T_obs, m_eff, .keep_all = TRUE)

  cat("\nComparison:", comparison, "\n")
  cat("Tests after dedup (distinct T_obs, m_eff): ", nrow(cmp_df), "\n", sep = "")

  # Benjamini-Hochberg on permutation p-values (reference only)
  cmp_df <- cmp_df |>
    mutate(p_adj_bh = p.adjust(p_perm, method = "BH"))

  # Efron's empirical null on standardized test statistic
  z_vals <- cmp_df$T_obs_stand
  keep <- is.finite(z_vals)
  cmp_df$fdr_emp_null <- NA_real_
  locfdr_plot_path <- file.path(
    plot_dir,
    paste0(gsub("[^A-Za-z0-9_-]+", "_", comparison), "_locfdr.png")
  )

  if (sum(keep) > 0) {
    ef <- NULL
    grDevices::png(
      filename = locfdr_plot_path,
      width = 1600,
      height = 1200,
      res = 150
    )
    ef <- tryCatch(
      locfdr::locfdr(z_vals[keep], plot = 1, df = 7, nulltype = 1),
      error = function(e) {
        warning("locfdr failed for ", comparison, ": ", conditionMessage(e))
        NULL
      }
    )
    grDevices::dev.off()

    if (!is.null(ef) && !is.null(ef$fdr) && length(ef$fdr) == sum(keep)) {
      cmp_df$fdr_emp_null[keep] <- ef$fdr
    } else {
      warning("locfdr did not return per-observation fdr for ", comparison)
      if (file.exists(locfdr_plot_path)) file.remove(locfdr_plot_path)
      locfdr_plot_path <- NA_character_
    }
  } else {
    warning("No finite T_obs_stand values for ", comparison)
    if (file.exists(locfdr_plot_path)) file.remove(locfdr_plot_path)
    locfdr_plot_path <- NA_character_
  }

  bh_hits <- sum(cmp_df$p_adj_bh < bh_cutoff, na.rm = TRUE)
  emp_hits <- sum(cmp_df$fdr_emp_null < emp_fdr_cutoff, na.rm = TRUE)

  cat("BH (p_perm, q<", bh_cutoff, ") hits: ", bh_hits, "\n", sep = "")
  cat("Empirical null (locfdr, fdr<", emp_fdr_cutoff, ") hits: ", emp_hits, "\n", sep = "")

  neg_nom_count <- sum(cmp_df$T_obs_stand < 0 & cmp_df$p_perm < 0.05, na.rm = TRUE)
  pos_nom_count <- sum(cmp_df$T_obs_stand > 0 & cmp_df$p_perm < 0.05, na.rm = TRUE)

  cat("Counts by sign:\n")
  cat("  Nominal p_perm<0.05  - : ", neg_nom_count, "  + : ", pos_nom_count, "\n", sep = "")
  cat("  Empirical fdr<", emp_fdr_cutoff, " total : ", emp_hits, "\n", sep = "")

  anno_hits <- cmp_df %>%
    filter(p_perm < 0.05, grepl("^anno_", feature)) %>%
    arrange(fdr_emp_null)

  if (nrow(anno_hits) > 0) {
    cat("\nSignificant anno_* features (p_perm<0.05):\n")
    print(anno_hits %>% select(feature, rank, p_perm, p_adj_bh, T_obs_stand, fdr_emp_null))
  }

  invisible(list(
    all = cmp_df,
    locfdr_plot_path = locfdr_plot_path
  ))
}

results <- lapply(comparisons_to_run, run_one)
names(results) <- comparisons_to_run

summary_all <- dplyr::bind_rows(
  lapply(seq_along(results), function(i) {
    cmp <- comparisons_to_run[[i]]
    all_tbl <- results[[i]]$all

    tibble(
      comparison = cmp,
      n_tests = nrow(all_tbl),
      n_nominal_p_lt_0_05 = sum(all_tbl$p_perm < 0.05, na.rm = TRUE),
      n_bh_q_lt_0_1 = sum(all_tbl$p_adj_bh < bh_cutoff, na.rm = TRUE),
      n_emp_fdr_lt_0_2 = sum(all_tbl$fdr_emp_null < emp_fdr_cutoff, na.rm = TRUE),
      n_emp_neg = sum(all_tbl$fdr_emp_null < emp_fdr_cutoff & all_tbl$T_obs_stand < 0, na.rm = TRUE),
      n_emp_pos = sum(all_tbl$fdr_emp_null < emp_fdr_cutoff & all_tbl$T_obs_stand > 0, na.rm = TRUE)
    )
  })
) %>%
  arrange(desc(n_emp_fdr_lt_0_2), desc(n_bh_q_lt_0_1), comparison)

empirical_hits_all <- dplyr::bind_rows(
  lapply(seq_along(results), function(i) {
    cmp <- comparisons_to_run[[i]]
    tbl <- results[[i]]$all
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    dplyr::mutate(tbl, comparison = cmp, .before = 1)
  })
) %>%
  filter(p_perm < 0.05) %>%
  arrange(desc(abs(T_obs_stand)))

locfdr_plot_info <- tibble(
  comparison = comparisons_to_run,
  plot_path = vapply(
    results,
    function(x) if (!is.null(x$locfdr_plot_path)) x$locfdr_plot_path else NA_character_,
    character(1)
  )
) %>%
  filter(!is.na(plot_path), file.exists(plot_path))

save_locfdr_grid <- function(plot_paths, titles, out_path, ncol = 8) {
  n <- length(plot_paths)
  if (n == 0) return(invisible(NULL))

  ncol <- max(1, min(ncol, n))
  nrow <- ceiling(n / ncol)

  grDevices::png(
    filename = out_path,
    width = 1400 * ncol,
    height = 1050 * nrow,
    res = 150
  )
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = nrow, ncol = ncol)
  grid::pushViewport(grid::viewport(layout = lay))

  for (i in seq_len(n)) {
    row_i <- ceiling(i / ncol)
    col_i <- i - (row_i - 1) * ncol
    img <- png::readPNG(plot_paths[[i]])

    grid::pushViewport(grid::viewport(layout.pos.row = row_i, layout.pos.col = col_i))
    grid::grid.rect(gp = grid::gpar(fill = "white", col = "grey85"))
    grid::grid.text(
      label = titles[[i]],
      x = 0.5, y = 0.98,
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    grid::grid.raster(img, x = 0.5, y = 0.45, width = 0.95, height = 0.86, interpolate = TRUE)
    grid::popViewport()
  }

  grid::popViewport()
  grDevices::dev.off()
  invisible(out_path)
}

locfdr_grid_path <- file.path(plot_dir, "locfdr_all_comparisons_grid.png")
if (nrow(locfdr_plot_info) > 0) {
  save_locfdr_grid(
    plot_paths = locfdr_plot_info$plot_path,
    titles = locfdr_plot_info$comparison,
    out_path = locfdr_grid_path,
    ncol = 8
  )
  cat("\nSaved locfdr grid to: ", locfdr_grid_path, "\n", sep = "")
} else {
  cat("\nNo locfdr plot files were generated.\n")
}

cat("\nRan ", length(comparisons_to_run), " requested comparisons (", length(missing_comparisons), " missing).\n", sep = "")
cat("\n\n=== Summary across selected results comparisons ===\n")
print(summary_all, n = nrow(summary_all))

cat("\n=== Combined significant hits (p_perm<0.05) ===\n")
if (nrow(empirical_hits_all) > 0) {
  print(
    empirical_hits_all %>%
      transmute(
        comparison,
        rank,
        feature,
        n_peptides_used,
        `Z-score` = T_obs_stand,
        p_perm,
        max_delta,
        p_adj_bh,
        fdr_emp_null
      ),
    n = min(500, nrow(empirical_hits_all))
  )
} else {
  cat("None\n")
}

if (interactive()) {
  if (nrow(empirical_hits_all) > 0) {
    View(empirical_hits_all %>%
           filter(p_adj_bh < 0.1) %>%
           #filter(fdr_emp_null < 0.2) %>%
      select(comparison, rank, feature, n_peptides_used, T_obs_stand, p_perm, max_delta, p_adj_bh, fdr_emp_null)
    )
  }
  if (nrow(locfdr_plot_info) > 0 && file.exists(locfdr_grid_path)) {
    utils::browseURL(locfdr_grid_path)
  }
}
