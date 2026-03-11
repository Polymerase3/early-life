#!/usr/bin/env Rscript

# Combine all DELTA_framework/delta_table.csv files under results/results
# into a single DuckDB table with the directory name recorded as "comparison".

library(duckdb)
library(DBI)
library(dplyr)
library(purrr)
library(readr)
library(phiper)
library(Cairo)
library(plotly)
library(htmlwidgets)

base_dir <- file.path("results", "results")
db_path <- file.path("results", "delta_tables.duckdb")

# locate all comparison directories that contain the expected file
delta_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
delta_dirs <- delta_dirs[
  file.exists(file.path(delta_dirs, "DELTA_framework", "delta_table.csv"))
]

if (length(delta_dirs) == 0) {
  stop("No DELTA_framework/delta_table.csv files found under ", base_dir)
}

read_delta <- function(dir_path) {
  comparison <- basename(dir_path)
  csv_path <- file.path(dir_path, "DELTA_framework", "delta_table.csv")
  readr::read_csv(csv_path, show_col_types = FALSE) |>
    mutate(comparison = comparison) |>
    relocate(comparison, .before = 1)
}

delta_tbl <- map_dfr(delta_dirs, read_delta)

# -------------------------------------------------------------------------
# Quick comparison: BH vs empirical-null (Efron) for one contrast
# -------------------------------------------------------------------------

target_cmp <- "kid_serum_T2_vs_kid_serum_T6"

cmp_df <- delta_tbl |>
  filter(comparison == target_cmp) |>
  distinct(T_obs, m_eff, .keep_all = TRUE)

# Benjamini-Hochberg on permutation p-values
cmp_df <- cmp_df |>
  mutate(p_adj_bh = p.adjust(p_perm, method = "BH"))

bh_hits <- sum(cmp_df$p_adj_bh < 0.1, na.rm = TRUE)

# Efron's empirical null on standardized test statistic
# Use locfdr directly on finite z-scores; flag hits at fdr < 0.15
z_vals <- cmp_df$T_obs_stand
keep <- is.finite(z_vals)
cmp_df$fdr_emp_null <- NA_real_

if (sum(keep) > 0) {
  ef <- locfdr::locfdr(z_vals[keep], plot = 1, df = 7, nulltype = 0)
  if (!is.null(ef$fdr) && length(ef$fdr) == sum(keep)) {
    cmp_df$fdr_emp_null[keep] <- ef$fdr
  } else {
    warning("locfdr did not return per-observation fdr; empirical-null FDR set to NA")
  }
} else {
  warning("No finite T_obs_stand values; empirical-null FDR set to NA")
}

emp_hits <- sum(cmp_df$fdr_emp_null < 0.15, na.rm = TRUE)

cat("Comparison:", target_cmp, "\n")
cat("BH (p_perm, q<0.1) hits:", bh_hits, "\n")
cat("Empirical null (locfdr on T_obs_stand, fdr<0.15) hits:", emp_hits, "\n")

# define delta hits: empirical fdr < 0.15
delta_hits <- cmp_df %>%
  filter(fdr_emp_null < 0.15) %>%
  arrange(fdr_emp_null)

# -------------------------------------------------------------------------
# Load POP single-peptide table for the target comparison
# -------------------------------------------------------------------------
single_pep_path <- file.path(base_dir, target_cmp, "POP_framework", "single_peptide.csv")

pep_tbl <- readr::read_csv(single_pep_path, show_col_types = FALSE)

# peptide library to map features -> peptides
peplib_path <- file.path("results", "peptide_library.rds")

peplib <- readRDS(peplib_path)

tax_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
tax_cols <- intersect(tax_ranks, names(peplib))
special_features <- c(
  "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
  "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
  "is_probio", "is_IgA", "is_flagellum", "signalp6_slow",
  "is_topgraph_new", "is_allergens",
  "anno_is_fungi", "anno_is_food", "anno_is_homo_sapiens", "anno_is_lacto_phage"
)

# helper: map feature to peptide_ids (same logic as 03-analysis)
get_binary_and_ids <- function(feature, peplib, tax_cols, peptide_col = "peptide_id", rank_name = NULL) {
  if (feature %in% special_features) {
    vals <- peplib[[feature]]
    present <- !is.na(vals) & as.logical(vals)
  } else if (feature %in% names(peplib)) {
    vals <- peplib[[feature]]
    present <- as.logical(vals)
    present[is.na(present)] <- FALSE
  } else if (!is.null(rank_name) && rank_name %in% names(peplib)) {
    vals <- peplib[[rank_name]]
    present <- !is.na(vals) & vals == feature
  } else {
    if (length(tax_cols) == 0L) {
      stop("No taxonomic columns found in peptide library.")
    }
    matches <- lapply(tax_cols, function(col) {
      vals <- peplib[[col]]
      !is.na(vals) & vals == feature
    })
    present <- Reduce(`|`, matches)
  }

  peptide_ids <- as.character(peplib[[peptide_col]][present])
  peptide_ids <- unique(peptide_ids)
  list(peptide_ids = peptide_ids)
}

downsample_for_static <- function(df, prop = 0.1, seed = 1L) {
  if (is.null(df) || !nrow(df)) return(df)
  n <- nrow(df)
  size <- max(1L, floor(n * prop))
  if (size >= n) return(df)
  set.seed(seed)
  dplyr::slice_sample(df, n = size)
}

add_background_static <- function(p, bg,
                                  size = 0.8,
                                  alpha = 0.12,
                                  color = "#808080") {
  if (is.null(bg) || !nrow(bg)) return(p)
  bg_layer <- ggplot2::geom_point(
    data = bg,
    mapping = ggplot2::aes(x = percent1, y = percent2),
    inherit.aes = FALSE,
    color = color,
    size = size,
    alpha = alpha,
    show.legend = FALSE
  )
  p$layers <- c(list(bg_layer), p$layers)
  p
}

# output dirs for quick scatter plots
hit_out_dir <- file.path(base_dir, target_cmp, "DELTA_framework", "interesting_features_quick")
scatter_dir <- file.path(hit_out_dir, "scatter")
dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

library(Cairo)
for (i in seq_len(nrow(delta_hits))) {
  feat <- delta_hits$feature[i]
  grp1 <- delta_hits$group1[i]
  grp2 <- delta_hits$group2[i]
  safe_name <- gsub("[^A-Za-z0-9_-]+", "_", as.character(feat))
  file_prefix <- file.path(hit_out_dir, safe_name)
  scatter_prefix <- file.path(scatter_dir, safe_name)

  ids <- get_binary_and_ids(feat, peplib, tax_cols, rank_name = delta_hits$rank[i])$peptide_ids
  feature_data <- pep_tbl %>% filter(feature %in% ids)

  if (!nrow(feature_data)) {
    feature_data <- pep_tbl %>%
      dplyr::filter(rank == delta_hits$rank[i], feature == feat)
  }

  if (!nrow(feature_data)) {
    message("No peptide rows for feature: ", feat)
    next
  }

  bg_df <- pep_tbl %>% filter(!(feature %in% ids))
  feature_data_static <- downsample_for_static(feature_data, prop = 1, seed = 1L)
  bg_df_static <- downsample_for_static(bg_df, prop = 1, seed = 1L)

  # scatter static
  CairoSVG(paste0(scatter_prefix, "_scatter_static.svg"),
    dpi = 300,
    height = 30, width = 30, unit = "cm", bg = "white"
  )
  p_scatter <- phiper::scatter_static(
    df = feature_data_static,
    xlab = grp1,
    ylab = grp2,
    point_size = 2,
    point_alpha = 0.85,
    jitter_width_pp = 0.15,
    jitter_height_pp = 0.15,
    font_family = "Montserrat",
    font_size = 12
  ) +
    ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
    ggplot2::theme(
      plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
      text        = ggplot2::element_text(family = "Montserrat")
    ) +
    ggplot2::ggtitle(feat)

  p_scatter <- add_background_static(
    p_scatter,
    bg = bg_df_static,
    size = 1,
    alpha = 0.35
  )
  print(p_scatter)
  dev.off()

  # scatter interactive
  bg_df <- bg_df %>%
    dplyr::mutate(category = "all peptides")
  cat_cols <- c(
    "significant (wBH, per rank)" = "#FF1744",
    "nominal only"                = "#00E676",
    "not significant"             = "#2979FF",
    "all peptides"                = "#7A7A7A"
  )
  p_inter <- phiper::scatter_interactive(
    df = feature_data,
    xlab = grp1,
    ylab = grp2,
    peplib = peplib,
    show_background = TRUE,
    background_df = bg_df,
    background_name = "all peptides",
    background_color = "#808080",
    background_alpha = 0.40,
    background_size = 7,
    background_max_n = Inf,
    category_colors = cat_cols,
    point_size = 16,
    point_alpha = 0.95,
    jitter_width_pp = 0.05,
    jitter_height_pp = 0.05,
    font_size = 12
  )
  p_inter <- plotly::layout(
    p_inter,
    autosize = TRUE,
    margin   = list(l = 70, r = 30, t = 10, b = 70),
    xaxis    = list(range = c(-2, 102), automargin = TRUE),
    yaxis    = list(range = c(-2, 102), automargin = TRUE)
  )
  htmlwidgets::saveWidget(
    p_inter,
    file = paste0(file_prefix, "_scatter_interactive.html"),
    selfcontained = FALSE
  )
}

# -------------------------------------------------------------------------
# Breastmilk comparison counts: raw-sig and BH-sig by Z-score direction
# -------------------------------------------------------------------------

breastmilk_comparisons <- c(
  "kid_serum_T6_vs_mom_milk_T6",
  "mom_serum_T2_vs_mom_milk_T4",
  "mom_milk_T4_vs_mom_milk_T7"
)

bm_counts <- purrr::map_dfr(breastmilk_comparisons, function(cmp) {
  path <- file.path(base_dir, cmp, "DELTA_framework", "delta_table_fdr.csv")
  if (!file.exists(path)) {
    warning("Missing delta_table_fdr.csv for ", cmp)
    return(NULL)
  }
  df <- readr::read_csv(path, show_col_types = FALSE) |>
    distinct(T_obs, m_eff, .keep_all = TRUE)
  tibble::tibble(
    comparison  = cmp,
    total_tests = nrow(df),
    raw_sig_pos = sum(df$p_perm   < 0.05 & df$T_obs_stand > 0, na.rm = TRUE),
    raw_sig_neg = sum(df$p_perm   < 0.05 & df$T_obs_stand < 0, na.rm = TRUE),
    bh_sig_pos  = sum(df$p_adj_bh < 0.1  & df$T_obs_stand > 0, na.rm = TRUE),
    bh_sig_neg  = sum(df$p_adj_bh < 0.1  & df$T_obs_stand < 0, na.rm = TRUE)
  )
})

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("BREASTMILK COMPARISON SIGNIFICANCE COUNTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")
cat(sprintf("%-42s  %12s  %11s  %11s  %10s  %10s\n",
            "comparison", "total_tests", "raw_sig_pos", "raw_sig_neg", "bh_sig_pos", "bh_sig_neg"))
cat(strrep("-", 102), "\n")
for (i in seq_len(nrow(bm_counts))) {
  cat(sprintf("%-42s  %12d  %11d  %11d  %10d  %10d\n",
              bm_counts$comparison[i],
              bm_counts$total_tests[i],
              bm_counts$raw_sig_pos[i],
              bm_counts$raw_sig_neg[i],
              bm_counts$bh_sig_pos[i],
              bm_counts$bh_sig_neg[i]))
}
cat("\n")
