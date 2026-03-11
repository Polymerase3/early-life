# load required packages
library(phiper)
library(rlang)
library(ggplot2)
library(Cairo)
library(openxlsx)
library(dplyr)
library(purrr)
library(locfdr)

set.seed(632961)

# parse command-line arguments for optional parameters
N_CORES <- 30
LOG <- TRUE
LOG_FILE <- NULL
MAX_GB <- 40
CMP_INDEX <- NA_integer_
CMP_KEY <- NULL
DATA_RDS_PATH <- file.path("data", "alex_flagellins", "phipseq_with_microbiome_all.rds")
RESULTS_ROOT <- file.path("results", "alex_flags")
args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("=", arg)) {
    parts <- strsplit(arg, "=")[[1]]
    key <- parts[1]
    value <- parts[2]
    if (key == "N_CORES") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) N_CORES <- val_num
    } else if (key == "MAX_GB") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) MAX_GB <- val_num
    } else if (key == "LOG") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        LOG <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        LOG <- FALSE
      }
    } else if (key == "LOG_FILE") {
      LOG_FILE <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key %in% c("CMP_INDEX", "COMPARISON_INDEX")) {
      val_num <- suppressWarnings(as.integer(value))
      if (!is.na(val_num)) CMP_INDEX <- val_num
    } else if (key %in% c("CMP_KEY", "COMPARISON_KEY")) {
      CMP_KEY <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key %in% c("DATA_RDS", "DATA_RDS_PATH")) {
      DATA_RDS_PATH <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key %in% c("OUT_DIR", "RESULTS_ROOT")) {
      RESULTS_ROOT <- sub("^['\\\"]|['\\\"]$", "", value)
    }
  }
}

if (!file.exists(DATA_RDS_PATH)) {
  stop("Input RDS not found: ", DATA_RDS_PATH)
}

analysis_input <- readRDS(DATA_RDS_PATH)
analysis_input <- as.data.frame(analysis_input, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "timepoint", "peptide_id", "exist")
missing_required <- setdiff(required_cols, names(analysis_input))
if (length(missing_required) > 0) {
  stop("Missing required columns in input: ", paste(missing_required, collapse = ", "))
}

if (!("fold_change" %in% names(analysis_input))) {
  analysis_input$fold_change <- analysis_input$exist
}
if (!("subject_id" %in% names(analysis_input))) {
  analysis_input$subject_id <- sub("-.*$", "", analysis_input$sample_id)
}

analysis_input$sample_id <- as.character(analysis_input$sample_id)
analysis_input$subject_id <- as.character(analysis_input$subject_id)
analysis_input$timepoint <- as.character(analysis_input$timepoint)
analysis_input$peptide_id <- as.character(analysis_input$peptide_id)
analysis_input$exist <- as.integer(analysis_input$exist)
analysis_input$fold_change <- as.numeric(analysis_input$fold_change)

exist_col_idx <- match("exist", names(analysis_input))
if (is.na(exist_col_idx) || exist_col_idx >= ncol(analysis_input)) {
  stop("Could not identify microbiome columns after 'exist'.")
}
microbiome_cols <- names(analysis_input)[(exist_col_idx + 1):ncol(analysis_input)]
if (length(microbiome_cols) == 0L) {
  stop("No microbiome columns found after 'exist'.")
}

timepoints <- unique(analysis_input$timepoint)
comparison_grid <- expand.grid(
  timepoint = timepoints,
  microbe_col = microbiome_cols,
  stringsAsFactors = FALSE
)
comparison_grid$cmp_key <- paste(comparison_grid$timepoint, comparison_grid$microbe_col, sep = "|")

if (is.na(CMP_INDEX)) {
  env_idx <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (nzchar(env_idx)) {
    val_num <- suppressWarnings(as.integer(env_idx))
    if (!is.na(val_num)) CMP_INDEX <- val_num
  }
}

if (!is.null(CMP_KEY) && !is.na(CMP_INDEX)) {
  stop("Provide only one of CMP_KEY or CMP_INDEX (or SLURM_ARRAY_TASK_ID).")
}

if (!is.null(CMP_KEY)) {
  match_idx <- match(CMP_KEY, comparison_grid$cmp_key)
  if (is.na(match_idx)) {
    stop("CMP_KEY not found. Example valid key: ", comparison_grid$cmp_key[[1]])
  }
  comparison_grid <- comparison_grid[match_idx, , drop = FALSE]
} else if (!is.na(CMP_INDEX)) {
  if (CMP_INDEX < 1L || CMP_INDEX > nrow(comparison_grid)) {
    stop("CMP_INDEX out of range: ", CMP_INDEX, " (1..", nrow(comparison_grid), ")")
  }
  comparison_grid <- comparison_grid[CMP_INDEX, , drop = FALSE]
}

# use the merged long table directly (no phip_convert / DuckDB CSV import)
ps <- analysis_input %>%
  dplyr::mutate(dplyr::across(where(is.logical), as.integer))

# create base results directory and save peptide library
dir.create(RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)
peplib_candidates <- c(
  file.path("results", "peptide_library.rds"),
  file.path("data", "peptide_library.rds")
)
peplib_path <- peplib_candidates[file.exists(peplib_candidates)][1]
if (is.na(peplib_path)) {
  stop("No peptide library RDS found. Expected one of: ", paste(peplib_candidates, collapse = ", "))
}
peplib_global <- readRDS(peplib_path)
saveRDS(peplib_global, file = file.path(RESULTS_ROOT, "peptide_library.rds"))

safe_component <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

as_phip_data_view <- function(df, peptide_library = NULL) {
  structure(
    list(
      data_long = df,
      comparisons = tibble::tibble(
        comparison = character(),
        group1 = character(),
        group2 = character(),
        variable = character()
      ),
      peptide_library = peptide_library,
      meta = list(full_cross = FALSE, exist_prop = NA_real_)
    ),
    class = "phip_data"
  )
}

# columns always to keep in saved data
base_cols <- c("sample_id", "peptide_id", "group_char", "exist", "fold_change")

# helper function to safely save an RDS file after creating directories
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}

# helper function: select base columns (and any extras), drop NA rows, collect and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))
  avail <- names(data)
  missing <- setdiff(cols_to_select, avail)
  if (length(missing)) {
    message("Skipping missing columns: ", paste(missing, collapse = ", "))
  }
  df <- data %>% dplyr::select(dplyr::any_of(cols_to_select))
  if (isTRUE(drop_na)) {
    df <- df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(intersect(cols_to_select, names(.))), ~ !is.na(.)))
  }
  save_rds_safe(df, out_path)
  invisible(df)
}

# set up parallel processing for DELTA framework analysis
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
options(future.globals.maxSize = MAX_GB * 1024^3, future.scheduling = 1)
original_plan <- future::plan()
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
} else {
  if (N_CORES > 1L) {
    future::plan(future::multicore, workers = N_CORES)
  } else {
    future::plan(future::sequential)
  }
}

# -------------------------------------------------------------------
# internal: per-plot PHIP palette mapping (recycled & named)
# -------------------------------------------------------------------
.phip_palette_map <- function(levels_vec, palette = phip_palette) {
  if (length(levels_vec) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }
  stats::setNames(rep_len(palette, length(levels_vec)), levels_vec)
}

# loop through each comparison and perform the same analysis
for (cmp_idx in seq_len(nrow(comparison_grid))) {
  cmp_timepoint <- comparison_grid$timepoint[[cmp_idx]]
  microbe_col <- comparison_grid$microbe_col[[cmp_idx]]
  var1 <- paste0(microbe_col, "=1")
  var2 <- paste0(microbe_col, "=0")

  message("Running comparison: [", cmp_timepoint, "] ", microbe_col, " (1 vs 0)")

  label_dir <- paste0(safe_component(microbe_col), "_1_vs_0")
  out_dir <- file.path(RESULTS_ROOT, safe_component(cmp_timepoint), label_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # filter to one timepoint and split by microbiome variable (1 vs 0)
  ps_cmp <- ps %>%
    filter(timepoint == cmp_timepoint) %>%
    filter(!is.na(!!sym(microbe_col)) & ((!!sym(microbe_col)) %in% c(0, 1))) %>%
    mutate(
      group_char  = if_else(!!sym(microbe_col) == 1L, var1, var2),
      group_dummy = if_else(!!sym(microbe_col) == 1L, 1L, 0L)
    )

  group_sample_counts <- ps_cmp %>%
    dplyr::select(sample_id, group_char) %>%
    dplyr::distinct() %>%
    dplyr::count(group_char, name = "n_samples")

  if (nrow(group_sample_counts) < 2L) {
    message("Skipping [", cmp_timepoint, "] ", microbe_col, " (need both 0 and 1 groups).")
    next
  }

  # save the filtered dataset for this comparison
  make_and_save(
    data       = ps_cmp,
    out_path   = file.path(out_dir, paste0(label_dir, "_data.rds")),
    extra_vars = unique(c("subject_id", "timepoint", microbe_col))
  )

  # ------------------ enrichment counts ------------------
  CairoSVG(file.path(out_dir, "enrichment_counts.svg"),
    dpi = 300,
    height = 30, width = 40, unit = "cm", bg = "white"
  )
  ps_cmp_phip <- as_phip_data_view(ps_cmp, peptide_library = peplib_global)
  p_enrich <- plot_enrichment_counts(ps_cmp_phip, group_cols = "group_char") +
    theme(text = element_text(family = "DejaVu Sans"))
  print(p_enrich)
  dev.off()

  # ------------------ alpha diversity --------------------
  carry_cols_alpha <- intersect(c("sex", "age", "subject_id", "timepoint"), names(ps_cmp))
  if (!length(carry_cols_alpha)) {
    carry_cols_alpha <- NULL
  }
  alpha_div <- compute_alpha_diversity(ps_cmp, group_cols = "group_char", carry_cols = carry_cols_alpha)
  dir.create(file.path(out_dir, "alpha_diversity"), recursive = TRUE, showWarnings = FALSE)
  write.xlsx(alpha_div, file.path(out_dir, "alpha_diversity", "table.xlsx"))

  CairoSVG(file.path(out_dir, "alpha_diversity", "plot.svg"),
    dpi = 300,
    height = 30, width = 40, unit = "cm", bg = "white"
  )
  p_alpha <- plot_alpha_diversity(alpha_div, metric = "richness", group_col = "group_char")
  print(p_alpha)
  dev.off()

  # ------------------ beta diversity ---------------------
  ps_cmp <- ps_cmp %>%
    mutate(sample_id = as.character(sample_id))

  dist_bc <- phiper:::compute_distance(ps_cmp,
    value_col = "exist",
    method_normalization = "hellinger",
    distance = "kulczynski", n_threads = 10
  )
  dir.create(file.path(out_dir, "beta_diversity"), recursive = TRUE, showWarnings = FALSE)
  dist_mat <- as.matrix(dist_bc)
  openxlsx::write.xlsx(dist_mat, file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"), rowNames = TRUE)

  pcoa_res <- phiper:::compute_pcoa(dist_bc, neg_correction = "none", n_axes = 109)
  saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))

  cap_res <- phiper:::compute_capscale(dist_bc, ps = ps_cmp, formula = ~group_char)
  saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))

  permanova_res <- phiper:::compute_permanova(dist_bc, ps = ps_cmp, group_col = "group_char")
  saveRDS(permanova_res, file.path(out_dir, "beta_diversity", "permanova_results.rds"))

  disp_res <- phiper:::compute_dispersion(dist_bc, ps = ps_cmp, group_col = "group_char")
  saveRDS(disp_res, file.path(out_dir, "beta_diversity", "dispersion_results.rds"))
  print(disp_res)

  n_samples <- nrow(dist_mat)
  tsne_perp <- function(default_p, n) {
    if (is.na(n) || n < 3L) return(NA_integer_)
    p <- min(as.integer(default_p), as.integer(n) - 1L)
    p <- max(2L, p)
    if (p >= n) return(NA_integer_)
    p
  }

  # 2D t-SNE (adaptive perplexity; skip gracefully for tiny n)
  p2 <- tsne_perp(15L, n_samples)
  if (!is.na(p2)) {
    tsne_res2d <- tryCatch(
      phiper:::compute_tsne(
        ps = ps_cmp, dist_obj = dist_bc, dims = 2L,
        perplexity = p2, meta_cols = c("group_char")
      ),
      error = function(e) {
        message("Skipping 2D t-SNE for ", var1, " vs ", var2, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(tsne_res2d)) {
      openxlsx::write.xlsx(tsne_res2d, file = file.path(out_dir, "beta_diversity", "tsne2d_results.xlsx"), rowNames = TRUE)
      CairoSVG(file.path(out_dir, "beta_diversity", "tsne2d_plot.svg"),
        dpi = 300,
        height = 30, width = 30, unit = "cm", bg = "white"
      )
      p_tsne2d <- phiper:::plot_tsne(tsne_res2d, view = "2d", colour = "group_char", palette = c("blue", "green"))
      print(p_tsne2d)
      dev.off()
    }
  } else {
    message("Skipping 2D t-SNE for ", var1, " vs ", var2, ": too few samples (n=", n_samples, ").")
  }

  # 3D t-SNE (adaptive perplexity; requires at least 4 samples)
  p3 <- if (n_samples >= 4L) tsne_perp(20L, n_samples) else NA_integer_
  if (!is.na(p3)) {
    tsne_res3d <- tryCatch(
      phiper:::compute_tsne(
        ps = ps_cmp, dist_obj = dist_bc, dims = 3L,
        perplexity = p3, meta_cols = c("group_char")
      ),
      error = function(e) {
        message("Skipping 3D t-SNE for ", var1, " vs ", var2, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(tsne_res3d)) {
      openxlsx::write.xlsx(tsne_res3d, file = file.path(out_dir, "beta_diversity", "tsne3d_results.xlsx"), rowNames = TRUE)
      p3d <- phiper:::plot_tsne(tsne_res3d, view = "3d", colour = "group_char", palette = c("blue", "green"))
      htmlwidgets::saveWidget(p3d, file = file.path(out_dir, "beta_diversity", "tsne3d_plot.html"), selfcontained = FALSE)
    }
  } else {
    message("Skipping 3D t-SNE for ", var1, " vs ", var2, ": too few samples (n=", n_samples, ").")
  }

  # add group information to PCoA sample coordinates
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
      ps_cmp %>% dplyr::select(sample_id, group_char) %>% dplyr::distinct(),
      by = "sample_id"
    )

  # PCoA plot with group centroids and ellipses
  CairoSVG(file.path(out_dir, "beta_diversity", "pcoa_plot.svg"),
    dpi = 300,
    height = 30, width = 30, unit = "cm", bg = "white"
  )
  p_pcoa <- phiper:::plot_pcoa(pcoa_res,
    axes = c(1, 2), group_col = "group_char",
    ellipse_by = "group", show_centroids = TRUE
  )
  print(p_pcoa)
  dev.off()

  # scree plot for first 15 axes of PCoA
  CairoSVG(file.path(out_dir, "beta_diversity", "scree_plot.svg"),
    dpi = 300,
    height = 30, width = 30, unit = "cm", bg = "white"
  )
  p_scree <- phiper:::plot_scree(pcoa_res, n_axes = 15, type = "line") + theme()
  print(p_scree)
  dev.off()

  # determine which contrast label actually exists in the dispersion object
  available_contrasts <- unique(disp_res$distances$contrast)

  # preferred pairwise label: var1 vs var2 (e.g. "control vs MCI")
  pair_contrast <- paste(var1, "vs", var2)

  if (pair_contrast %in% available_contrasts) {
    contrast_to_use <- pair_contrast
  } else if ("<global>" %in% available_contrasts) {
    # fallback: use global dispersion if pairwise distances are not stored
    contrast_to_use <- "<global>"
  } else {
    # last resort: just take the first available contrast and warn
    contrast_to_use <- available_contrasts[1]
    message(
      "Warning: requested contrast '", pair_contrast,
      "' not found in disp_res$distances$contrast. Using '",
      contrast_to_use, "' instead."
    )
  }

  CairoSVG(file.path(out_dir, "beta_diversity", "dispersion_plot.svg"),
    dpi = 300,
    height = 30, width = 30, unit = "cm", bg = "white"
  )
  p_disp <- phiper:::plot_dispersion(
    disp_res,
    scope        = "group",
    contrast     = contrast_to_use,
    show_violin  = TRUE,
    show_box     = TRUE,
    show_points  = TRUE
  )
  print(p_disp)
  dev.off()

  # ------------------ POP framework ----------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))
  peplib <- peplib_global

  extract_tbl <- function(obj) {
    if (is.data.frame(obj)) {
      return(tibble::as_tibble(obj))
    }
    for (nm in c("data", "table", "tbl", "df", "result", "results")) {
      if (!is.null(obj[[nm]])) {
        return(tibble::as_tibble(obj[[nm]]))
      }
    }
    out <- try(tibble::as_tibble(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
    out <- try(as.data.frame(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(tibble::as_tibble(out))
    }
    stop("Cannot extract a data table from the POP result object.")
  }

  downsample_for_static <- function(df, prop = 1, seed = 1L) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    n <- nrow(df)
    size <- max(1L, floor(n * prop))
    if (size >= n) {
      return(df)
    }
    set.seed(seed)
    dplyr::slice_sample(df, n = size)
  }

  dir.create(file.path(out_dir, "POP_framework"), recursive = TRUE, showWarnings = FALSE)

  paired_prev <- FALSE

  prev_res_pep <- phiper::ph_prevalence_compare(
    x                 = data_frameworks,
    group_cols        = "group_char",
    rank_cols         = "peptide_id",
    compute_ratios_db = TRUE,
    paired            = paired_prev,
    parallel          = TRUE,
    collect           = TRUE
  )
  pep_tbl <- extract_tbl(prev_res_pep)
  write.csv(pep_tbl, file.path(out_dir, "POP_framework", "single_peptide.csv"))

  ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  prev_res_rank <- phiper::ph_prevalence_compare(
    x                 = data_frameworks,
    group_cols        = "group_char",
    rank_cols         = ranks_tax,
    compute_ratios_db = FALSE,
    paired            = paired_prev,
    parallel          = TRUE,
    peptide_library   = peplib,
    collect           = TRUE
  )
  rank_tbl <- extract_tbl(prev_res_rank)
  write.csv(rank_tbl, file.path(out_dir, "POP_framework", "taxa_ranks.csv"))

  ranks_combined <- c(ranks_tax, "peptide_id")
  plots_dir <- file.path(out_dir, "POP_framework", "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  for (rank_name in ranks_combined) {
    rank_chr <- as.character(rank_name)
    out_name <- file.path(plots_dir, rank_chr)
    df_rank <- if (rank_chr == "peptide_id") {
      pep_tbl
    } else {
      rank_tbl %>% filter(rank == rank_chr)
    }
    df_rank_static <- downsample_for_static(df_rank, prop = 1, seed = 1L)
    p_static <- scatter_static(
      df = df_rank_static,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      point_size = 2,
      jitter_width_pp = 0.15,
      jitter_height_pp = 0.15,
      point_alpha = 0.85,
      font_size = 12
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat")
      )
    ggsave(paste0(out_name, "_static.svg"), p_static,
      dpi = 300,
      height = 30, width = 30, unit = "cm", bg = "white",
      device = Cairo::CairoSVG
    )

    p_inter <- scatter_interactive(
      df = df_rank,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      peplib = peplib,
      point_size = 10,
      jitter_width_pp = 0.25,
      jitter_height_pp = 0.25,
      point_alpha = 0.85,
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
      file = paste0(out_name, "_interactive.html"),
      selfcontained = FALSE
    )
  }

  # ------------------ DELTA framework ---------------------
  dir.create(file.path(out_dir, "DELTA_framework"), recursive = TRUE, showWarnings = FALSE)

  peptides_in_test <- data_frameworks %>%
    dplyr::pull(peptide_id) %>%
    as.character() %>%
    unique() %>%
    stats::na.omit() %>%
    as.character()

  peplib_delta <- tibble::tibble(
    peptide_id = peptides_in_test,
    is_flagellin = "1"
  )

  # keep full metadata peplib for plotting, but add a matching flag column
  peplib$is_flagellin <- ifelse(as.character(peplib$peptide_id) %in% peptides_in_test, "1", "0")
  peplib[] <- lapply(peplib, as.character)
  peplib_delta[] <- lapply(peplib_delta, as.character)
  log_file_current <- if (is.null(LOG_FILE)) {
    file.path(out_dir, "DELTA_framework", "log.txt")
  } else {
    LOG_FILE
  }

  stat_mode <- "srlr"
  paired_by <- NULL

  res <- phiper::compute_delta(
    x = data_frameworks,
    exist_col = "exist",
    rank_cols = "is_flagellin",
    group_cols = "group_char",
    peptide_library = peplib_delta,
    B_permutations = 150000L,
    weight_mode = "equal",
    stat_mode = stat_mode,
    winsor_z = Inf,
    rank_feature_keep = c(is_flagellin = "1"),
    log = LOG,
    log_file = log_file_current,
    paired_by = paired_by
  )
  res <- as.data.frame(res)

  write.csv(res, file = file.path(out_dir, "DELTA_framework", "delta_table.csv"))
  tax_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  std_idx <- res$rank %in% tax_ranks
  res$feature[!std_idx] <- as.character(res$rank[!std_idx])
  res$rank <- "all"

  CairoSVG(file.path(out_dir, "DELTA_framework", "uncorrected_significant_static_forestplot.svg"),
    dpi = 300,
    height = 30, width = 30, unit = "cm", bg = "white"
  )
  p_forest_unc <- phiper::forestplot(
    results_tbl = res,
    rank_of_interest = "all",
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", df_rank$group1[1]),
    right_label = paste0("More in ", df_rank$group2[1]),
    label_vjust = -0.9,
    y_pad = 0.3,
    label_x_gap_frac = -0.3,
    statistic_to_plot = "T_stand"
  )
  print(p_forest_unc)
  dev.off()

  p_inter <- phiper::forestplot_interactive(
    results_tbl = res,
    rank_of_interest = "all",
    statistic_to_plot = "T_stand",
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", df_rank$group1[1]),
    right_label = paste0("More in ", df_rank$group2[1]),
    arrow_length_frac = 0.35,
    label_x_gap_frac = -0.3,
    label_y_offset = -0.9
  )$plot
  htmlwidgets::saveWidget(
    p_inter,
    file = file.path(out_dir, "DELTA_framework", "uncorrected_significant_interactive_forestplot.html"),
    selfcontained = FALSE
  )

  # ------------------ plot interesting features ---------
  # Always include these features in the "interesting" set (even if not significant)
  always_keep <- c("")

  res_filtered <- res %>%
    dplyr::mutate(
      .force_keep = (.data$feature %in% always_keep) | (.data$rank %in% always_keep)
    ) %>%
    dplyr::filter(
      .force_keep | (.data$p_perm < 0.05)
    ) %>%
    dplyr::arrange(
      dplyr::desc(.force_keep), # optional: keep forced ones at the top
      dplyr::desc(.data$T_obs_stand)
    ) %>%
    dplyr::select(-.force_keep)

  write.csv(
    res_filtered,
    file      = file.path(out_dir, "DELTA_framework", "delta_table_interesting.csv"),
    row.names = FALSE
  )

  tax_cols <- intersect(tax_ranks, names(peplib))

  special_features <- c(
    "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
    "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
    "is_probio", "is_IgA", "is_flagellum", "is_flagellin", "signalp6_slow",
    "is_topgraph_new", "is_allergens",
    "anno_is_fungi", "anno_is_food", "anno_is_homo_sapiens", "anno_is_lacto_phage"
  )

  get_binary_and_ids <- function(feature, peplib, tax_cols,
                                 peptide_col = "peptide_id") {
    if (feature %in% special_features) {
      vals <- peplib[[feature]]
      present <- !is.na(vals) & as.logical(vals)
    } else if (feature %in% names(peplib)) {
      vals <- peplib[[feature]]
      present <- as.logical(vals)
      present[is.na(present)] <- FALSE
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
    list(present = present, peptide_ids = peptide_ids)
  }

  res_with_pep <- res_filtered %>%
    mutate(
      match_info = map(feature, ~ get_binary_and_ids(.x, peplib, tax_cols)),
      binary_in_peplib = map(match_info, "present"),
      peptide_ids = map(match_info, "peptide_ids"),
      pep_tbl_subset = map(peptide_ids, ~ pep_tbl %>% filter(feature %in% .x))
    ) %>%
    select(-match_info)

  dir.create(file.path(out_dir, "DELTA_framework", "interesting_features"), recursive = TRUE, showWarnings = FALSE)

  add_background_static <- function(p, bg,
                                    size = 0.8,
                                    alpha = 0.12,
                                    color = "#808080") {
    if (is.null(bg) || !nrow(bg)) {
      return(p)
    }

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

  plot_feature_all <- function(feature_name, group1, group2, feature_data, out_dir) {
    if (is.null(feature_data) || nrow(feature_data) == 0L) {
      message("Skipping ", feature_name, " (no peptide data)")
      return(invisible(NULL))
    }
    safe_name <- gsub("[^A-Za-z0-9_-]+", "_", as.character(feature_name))
    base_dir <- file.path(out_dir, "DELTA_framework", "interesting_features")
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    scatter_dir <- file.path(base_dir, "scatter")
    dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
    file_prefix <- file.path(base_dir, safe_name)
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

      # keep only unique (percent1, percent2) combos; drop the rest at random
      if (all(c("percent1", "percent2") %in% names(bg_df))) {
        set.seed(BG_SEED) # makes the random pick reproducible
        bg_df <- bg_df %>%
          dplyr::slice_sample(prop = 1) %>% # shuffle rows
          dplyr::distinct(percent1, percent2, .keep_all = TRUE)
      }

      if (nrow(bg_df) > BG_MAX_N_INTERACTIVE) {
        set.seed(BG_SEED)
        bg_df <- bg_df %>% dplyr::slice_sample(n = BG_MAX_N_INTERACTIVE)
      }
    }

    feature_data_static <- downsample_for_static(feature_data, prop = 1, seed = 1L)
    bg_df_static <- downsample_for_static(bg_df, prop = 1, seed = BG_SEED)

    ## ---------------- SCATTER STATIC ----------------
    CairoSVG(paste0(scatter_prefix, "_scatter_static.svg"),
      dpi = 300,
      height = 30, width = 30, unit = "cm", bg = "white"
    )
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
        text        = ggplot2::element_text(family = "Montserrat")
      ) +
      ggplot2::ggtitle(feature_name)

    p_scatter <- add_background_static(
      p_scatter,
      bg = bg_df_static,
      size = 1,
      alpha = 0.35
    )

    print(p_scatter)
    dev.off()

    ## ---------------- SCATTER INTERACTIVE ----------------
    # Make sure background cannot “inherit” any categories
    if (!is.null(bg_df) && nrow(bg_df)) {
      bg_df <- bg_df %>%
        dplyr::select(-dplyr::any_of(c("category"))) %>% # safety
        dplyr::mutate(category = "all peptides") # force separate legend/trace
    }

    cat_cols <- c(
      "significant (wBH, per rank)" = "#FF1744",
      "nominal only"                = "#00E676",
      "not significant"             = "#2979FF",
      "all peptides"                = "#7A7A7A"
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
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )

    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(file_prefix, "_scatter_interactive.html"),
      selfcontained = FALSE
    )

    ## ---------------- DELTA PLOT: conditional smooth ----------------
    use_smooth <- nrow(feature_data) >= 7
    smooth_k <- if (use_smooth) 3L else 1L

    ## ---------------- DELTA PLOT STATIC ----------------
    CairoSVG(paste0(file_prefix, "_deltaplot_static.svg"),
      dpi = 300,
      height = 30, width = 30, unit = "cm", bg = "white"
    )
    p_delta_static <- tryCatch(
      {
        deltaplot(
          prev_tbl              = feature_data,
          group_pair_values     = c(group1, group2),
          group_labels          = c(group1, group2),
          point_jitter_width    = 0.01,
          point_jitter_height   = 0.01,
          point_alpha           = 0.6,
          point_size            = 6,
          add_smooth            = use_smooth,
          smooth_k              = smooth_k,
          arrow_head_length_mm  = 4
        ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12))
      },
      error = function(e) {
        message(
          "Delta static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No delta plot for ", feature_name, "\n(", group1, " vs ", group2, ")")
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat"))
      }
    )
    print(p_delta_static)
    dev.off()

    ## ---------------- DELTA PLOT INTERACTIVE ----------------
    # p_delta <- tryCatch(
    #   deltaplot_interactive(
    #     prev_tbl            = feature_data,
    #     group_pair_values   = c(group1, group2),
    #     group_labels        = c(group1, group2),
    #     point_alpha         = 0.6,
    #     point_size          = 6,
    #     add_smooth          = use_smooth,
    #     smooth_k            = smooth_k,
    #     arrow_length_frac   = 0.35, # old arrow_frac_h
    #     point_jitter_width  = 0.01,
    #     point_jitter_height = 0.01
    #   ),
    #   error = function(e) {
    #     message(
    #       "Delta interactive plot failed for ", feature_name,
    #       " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #     )
    #     NULL
    #   }
    # )
    # if (!is.null(p_delta)) {
    #   htmlwidgets::saveWidget(
    #     p_delta,
    #     file = paste0(file_prefix, "_deltaplot_interactive.html"),
    #     selfcontained = FALSE
    #   )
    # }

    ## ---------------- ECDF STATIC ----------------
    CairoSVG(paste0(file_prefix, "_ecdfplot_static.svg"),
      dpi = 300,
      height = 30, width = 30, unit = "cm", bg = "white"
    )
    p_ecdf_static <- tryCatch(
      {
        ecdf_plot(
          prev_tbl            = feature_data,
          group_pair_values   = c(group1, group2),
          group_labels        = c(group1, group2),
          line_width_pt       = 1,
          line_alpha          = 1,
          show_median_lines   = TRUE,
          show_ks_test        = TRUE
        ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12))
      },
      error = function(e) {
        message(
          "ECDF static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No ECDF plot for ", feature_name, "\n(", group1, " vs ", group2, ")")
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat"))
      }
    )
    print(p_ecdf_static)
    dev.off()

    ## ---------------- ECDF INTERACTIVE ----------------
    # p_ecdf <- tryCatch(
    #   ecdf_plot_interactive(
    #     prev_tbl            = feature_data,
    #     group_pair_values   = c(group1, group2),
    #     group_labels        = c(group1, group2),
    #     line_width_px       = 2,
    #     line_alpha          = 1,
    #     show_median_lines   = TRUE,
    #     show_ks_test        = TRUE
    #   ),
    #   error = function(e) {
    #     message(
    #       "ECDF interactive plot failed for ", feature_name,
    #       " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #     )
    #     NULL
    #   }
    # )
    # if (!is.null(p_ecdf)) {
    #   htmlwidgets::saveWidget(
    #     p_ecdf,
    #     file = paste0(file_prefix, "_ecdfplot_interactive.html"),
    #     selfcontained = FALSE
    #   )
    # }

    invisible(NULL)
  }

  n_features <- nrow(res_with_pep)
  for (i in seq_len(n_features)) {
    feature_name <- res_with_pep$feature[i]
    group1 <- res_with_pep$pep_tbl_subset[[i]]$group1[1]
    group2 <- res_with_pep$pep_tbl_subset[[i]]$group2[1]
    feature_data <- res_with_pep$pep_tbl_subset[[i]]
    message("Plotting [", i, "/", n_features, "]: ", feature_name)
    plot_feature_all(feature_name, group1, group2, feature_data, out_dir)
  }
}

# restore original future plan
future::plan(original_plan)
