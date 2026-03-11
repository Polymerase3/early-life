#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phiper)
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(Cairo)
  library(htmlwidgets)
  library(grid)
})

base_dir <- file.path("results", "results")
base_dir_milk <- file.path("results", "results_milk")
peplib <- readRDS(file.path("results", "peptide_library.rds"))

selected_pairs <- list(
  # c("mom_serum_T2", "mom_milk_T4"),
  # c("mom_serum_T2", "mom_milk_T7"),
  # c("kid_serum_T6", "mom_milk_T4"),
  # c("kid_serum_T6", "mom_milk_T7"),
  # c("kid_serum_T2_siblings", "kid_serum_T2_no_siblings"),
  # c("kid_serum_T6_siblings", "kid_serum_T6_no_siblings"),
  # c("kid_serum_T8_siblings", "kid_serum_T8_no_siblings"),
  # c("kid_serum_T2_delmode_VG", "kid_serum_T2_delmode_CS"),
  # c("kid_serum_T6_delmode_VG", "kid_serum_T6_delmode_CS"),
  # c("kid_serum_T8_delmode_VG", "kid_serum_T8_delmode_CS"),
  # c("kid_serum_T2_delplace_home", "kid_serum_T2_delplace_hospital"),
  # c("kid_serum_T6_delplace_home", "kid_serum_T6_delplace_hospital"),
  # c("kid_serum_T8_delplace_home", "kid_serum_T8_delplace_hospital"),
  # c("kid_serum_T2_feeding_BF", "kid_serum_T2_feeding_MF"),
  # c("kid_serum_T6_feeding_BF", "kid_serum_T6_feeding_MF"),
  # c("kid_serum_T8_feeding_BF", "kid_serum_T8_feeding_MF"),
  # c("kid_serum_T2_preterm_yes", "kid_serum_T2_preterm_no"),
  # c("kid_serum_T6_preterm_yes", "kid_serum_T6_preterm_no"),
  # c("kid_serum_T8_preterm_yes", "kid_serum_T8_preterm_no"),
  # c("kid_serum_T2_CDrisk_yes", "kid_serum_T2_CDrisk_no"),
  # c("kid_serum_T6_CDrisk_yes", "kid_serum_T6_CDrisk_no"),
  # c("kid_serum_T8_CDrisk_yes", "kid_serum_T8_CDrisk_no"),
  # c("kid_serum_T2_lockdown_before", "kid_serum_T2_lockdown_after"),
  # c("kid_serum_T6_lockdown_before", "kid_serum_T6_lockdown_after"),
  # c("kid_serum_T8_lockdown_before", "kid_serum_T8_lockdown_after"),
  # c("kid_serum_T2_smoking_yes", "kid_serum_T2_smoking_no"),
  # c("kid_serum_T6_smoking_yes", "kid_serum_T6_smoking_no"),
  # c("kid_serum_T8_smoking_yes", "kid_serum_T8_smoking_no"),
  # c("kid_serum_T2_sex_male", "kid_serum_T2_sex_female"),
  # c("kid_serum_T6_sex_male", "kid_serum_T6_sex_female"),
  # c("kid_serum_T8_sex_male", "kid_serum_T8_sex_female")
)

selected_comparisons <- unique(vapply(selected_pairs, paste, collapse = "_vs_", FUN.VALUE = character(1)))

# Include all comparisons from results/results_milk as an additional set.
milk_comparisons <- character(0)
if (dir.exists(base_dir_milk)) {
  milk_dirs <- list.dirs(base_dir_milk, recursive = FALSE, full.names = FALSE)
  milk_comparisons <- sort(unique(milk_dirs[nzchar(milk_dirs)]))
}

tax_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
tax_cols <- intersect(tax_ranks, names(peplib))
special_features <- c(
  "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
  "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
  "is_probio", "is_IgA", "is_flagellum", "signalp6_slow",
  "is_topgraph_new", "is_allergens",
  "anno_is_fungi", "anno_is_food", "anno_is_homo_sapiens", "anno_is_lacto_phage"
)

get_binary_and_ids <- function(feature, peplib, tax_cols, peptide_col = "peptide_id") {
  if (feature %in% special_features) {
    vals <- peplib[[feature]]
    present <- !is.na(vals) & as.logical(vals)
  } else if (feature %in% names(peplib)) {
    vals <- peplib[[feature]]
    present <- as.logical(vals)
    present[is.na(present)] <- FALSE
  } else {
    if (length(tax_cols) == 0L) stop("No taxonomic columns found in peptide library.")
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

downsample_for_static <- function(df, prop = 1, seed = 1L) {
  if (is.null(df) || !nrow(df)) return(df)
  n <- nrow(df)
  size <- max(1L, floor(n * prop))
  if (size >= n) return(df)
  set.seed(seed)
  dplyr::slice_sample(df, n = size)
}

add_background_static <- function(p, bg, size = 0.8, alpha = 0.12, color = "#808080") {
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

plot_feature_all <- function(feature_name, group1, group2, feature_data, out_dir, pep_tbl, peplib) {
  if (is.null(feature_data) || nrow(feature_data) == 0L) {
    message("Skipping ", feature_name, " (no peptide data)")
    return(invisible(NULL))
  }

  safe_name <- gsub("[^A-Za-z0-9_-]+", "_", as.character(feature_name))
  base_plot_dir <- file.path(out_dir, "DELTA_framework", "interesting_features")
  dir.create(base_plot_dir, recursive = TRUE, showWarnings = FALSE)
  scatter_dir <- file.path(base_plot_dir, "scatter")
  dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- file.path(base_plot_dir, safe_name)
  scatter_prefix <- file.path(scatter_dir, safe_name)

  SHOW_BG <- TRUE
  BG_MAX_N_INTERACTIVE <- Inf
  BG_SEED <- 1L

  bg_df <- NULL
  if (isTRUE(SHOW_BG)) {
    bg_df <- pep_tbl
    if ("feature" %in% names(bg_df) && "feature" %in% names(feature_data)) {
      bg_df <- bg_df %>% dplyr::filter(!(feature %in% feature_data$feature))
    }
    if (all(c("percent1", "percent2") %in% names(bg_df))) {
      set.seed(BG_SEED)
      bg_df <- bg_df %>% dplyr::slice_sample(prop = 1) %>% dplyr::distinct(percent1, percent2, .keep_all = TRUE)
    }
    if (nrow(bg_df) > BG_MAX_N_INTERACTIVE) {
      set.seed(BG_SEED)
      bg_df <- bg_df %>% dplyr::slice_sample(n = BG_MAX_N_INTERACTIVE)
    }
  }

  feature_data_static <- downsample_for_static(feature_data, prop = 1, seed = 1L)
  bg_df_static <- downsample_for_static(bg_df, prop = 1, seed = BG_SEED)

  CairoSVG(paste0(scatter_prefix, "_scatter_static.svg"), dpi = 300, height = 30, width = 30, unit = "cm", bg = "white")
  p_scatter <- scatter_static(
    df = feature_data_static,
    xlab = group1,
    ylab = group2,
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
      text = ggplot2::element_text(family = "Montserrat")
    ) +
    ggplot2::ggtitle(feature_name)

  p_scatter <- add_background_static(p_scatter, bg = bg_df_static, size = 1, alpha = 0.35)
  print(p_scatter)
  grDevices::dev.off()

  if (!is.null(bg_df) && nrow(bg_df)) {
    bg_df <- bg_df %>%
      dplyr::select(-dplyr::any_of(c("category"))) %>%
      dplyr::mutate(category = "all peptides")
  }

  cat_cols <- c(
    "significant (wBH, per rank)" = "#FF1744",
    "nominal only" = "#00E676",
    "not significant" = "#2979FF",
    "all peptides" = "#7A7A7A"
  )

  p_inter <- scatter_interactive(
    df = feature_data,
    xlab = group1,
    ylab = group2,
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
    margin = list(l = 70, r = 30, t = 10, b = 70),
    xaxis = list(range = c(-2, 102), automargin = TRUE),
    yaxis = list(range = c(-2, 102), automargin = TRUE)
  )

  htmlwidgets::saveWidget(
    p_inter,
    file = paste0(file_prefix, "_scatter_interactive.html"),
    selfcontained = FALSE
  )

  use_smooth <- nrow(feature_data) >= 7
  smooth_k <- if (use_smooth) 3L else 1L

  CairoSVG(paste0(file_prefix, "_deltaplot_static.svg"), dpi = 300, height = 30, width = 30, unit = "cm", bg = "white")
  p_delta_static <- tryCatch(
    deltaplot(
      prev_tbl = feature_data,
      group_pair_values = c(group1, group2),
      group_labels = c(group1, group2),
      point_jitter_width = 0.01,
      point_jitter_height = 0.01,
      point_alpha = 0.6,
      point_size = 6,
      add_smooth = use_smooth,
      smooth_k = smooth_k,
      arrow_head_length_mm = 4
    ) + ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12)),
    error = function(e) {
      ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::ggtitle(paste0("No delta plot for ", feature_name, "\n(", group1, " vs ", group2, ")"))
    }
  )
  print(p_delta_static)
  grDevices::dev.off()

  CairoSVG(paste0(file_prefix, "_ecdfplot_static.svg"), dpi = 300, height = 30, width = 30, unit = "cm", bg = "white")
  p_ecdf_static <- tryCatch(
    ecdf_plot(
      prev_tbl = feature_data,
      group_pair_values = c(group1, group2),
      group_labels = c(group1, group2),
      line_width_pt = 1,
      line_alpha = 1,
      show_median_lines = TRUE,
      show_ks_test = TRUE
    ) + ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12)),
    error = function(e) {
      ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::ggtitle(paste0("No ECDF plot for ", feature_name, "\n(", group1, " vs ", group2, ")"))
    }
  )
  print(p_ecdf_static)
  grDevices::dev.off()

  invisible(NULL)
}

replot_tasks <- dplyr::bind_rows(
  tibble::tibble(
    comparison = selected_comparisons,
    base_dir = base_dir,
    source = "results"
  ),
  tibble::tibble(
    comparison = milk_comparisons,
    base_dir = base_dir_milk,
    source = "results_milk"
  )
) %>%
  dplyr::distinct(source, comparison, .keep_all = TRUE)

for (task_idx in seq_len(nrow(replot_tasks))) {
  cmp <- replot_tasks$comparison[[task_idx]]
  cmp_base_dir <- replot_tasks$base_dir[[task_idx]]
  cmp_source <- replot_tasks$source[[task_idx]]
  out_dir <- file.path(cmp_base_dir, cmp)
  delta_path <- file.path(out_dir, "DELTA_framework", "delta_table_interesting.csv")
  pep_path <- file.path(out_dir, "POP_framework", "single_peptide.csv")

  if (!file.exists(delta_path) || !file.exists(pep_path)) {
    message("Skipping comparison (missing file) [", cmp_source, "]: ", cmp)
    next
  }

  message("\n=== Replotting interesting features for [", cmp_source, "]: ", cmp, " ===")
  res_filtered <- readr::read_csv(delta_path, show_col_types = FALSE)
  pep_tbl <- readr::read_csv(pep_path, show_col_types = FALSE)

  if (nrow(res_filtered) == 0 || nrow(pep_tbl) == 0) {
    message("Skipping ", cmp, " (empty delta or peptide table)")
    next
  }

  res_with_pep <- res_filtered %>%
    mutate(
      peptide_ids = map(feature, ~ get_binary_and_ids(.x, peplib, tax_cols)$peptide_ids),
      pep_tbl_subset = map(peptide_ids, ~ pep_tbl %>% filter(feature %in% .x))
    )

  for (i in seq_len(nrow(res_with_pep))) {
    feature_name <- res_with_pep$feature[[i]]
    feature_data <- res_with_pep$pep_tbl_subset[[i]]

    if (is.null(feature_data) || nrow(feature_data) == 0) {
      message("Skipping ", feature_name, " in ", cmp, " (no mapped peptide rows)")
      next
    }

    group1 <- feature_data$group1[[1]]
    group2 <- feature_data$group2[[1]]
    message("  Plotting [", i, "/", nrow(res_with_pep), "]: ", feature_name)
    plot_feature_all(feature_name, group1, group2, feature_data, out_dir, pep_tbl, peplib)
  }
}

message("\nDone.")
