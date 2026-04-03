# ==============================================================================
# 27-richness.R
# exports: alpha diversity tables + enrichment plots
#
# Outputs saved to results/richness/:
#   time_alpha__big_group_timepoint_recoded.{csv,rds}
#   time_alpha__kid__<var>_*_timepoint_recoded.{csv,rds}   (one per meta_var)
#   enrichment__big_group_timepoint_recoded.png
#   enrichment__kid__<var>_*_timepoint_recoded.png          (one per meta_var)
#
# After running: manually copy results/richness/ to shiny_app/data/richness/
# ==============================================================================

suppressPackageStartupMessages({
  library(phiper)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(withr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- output directory ---------------------------------------------------------
out_dir <- "results/richness"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ------------------------------------------------------------------
.save_tbl <- function(df, stem) {
  readr::write_csv(df, file.path(out_dir, paste0(stem, ".csv")))
  saveRDS(df,         file.path(out_dir, paste0(stem, ".rds")))
}

.save_png <- function(plot_obj, file, width = 2200, height = 1400, dpi = 300, units = "px") {
  ggplot2::ggsave(
    filename = file, plot = plot_obj,
    width = width, height = height, dpi = dpi, units = units, bg = "white"
  )
}

.clean_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ==============================================================================
# Load data
# ==============================================================================
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path    = "data/phiper_validated_babies_full.parquet",
    sample_id         = "sampleID",
    peptide_id        = "peptideID",
    exist             = "exist",
    fold_change       = "fold_change",
    counts_input      = "input",
    counts_hit        = "count",
    peptide_library   = TRUE,
    materialise_table = TRUE,
    auto_expand       = FALSE,
    n_cores           = 10
  )
})

# --- join categorical metadata (required for kid subgroup enrichment plots) ---
# metas_clean.rds is produced by 01-data-prepare.R; columns used here:
#   delivery_mode, infant_sex, ever_breastfed, feedmode_m3_bin,
#   smoking, delivery_place, siblings, preterm, risk_CD_bin
metas <- readRDS("data/metas_clean.rds") %>%
  dplyr::select(
    id_infant,
    delivery_mode, infant_sex, ever_breastfed, feedmode_m3_bin,
    smoking, delivery_place, siblings, preterm, risk_CD_bin
  )

ps <- left_join(
  ps, metas,
  by    = dplyr::join_by(subject_id == id_infant),
  copy  = TRUE
)

# ==============================================================================
# A) Alpha-diversity tables — big_group × timepoint_recoded
# ==============================================================================
time_alpha <- compute_alpha_diversity(
  ps,
  ranks             = c("peptide_id", "species", "genus", "family", "order", "class", "phylum"),
  group_cols        = c("big_group", "timepoint_recoded"),
  group_interaction = TRUE,
  carry_cols        = "subject_id",
  interaction_only  = TRUE
)

purrr::iwalk(time_alpha, function(df, nm) {
  .save_tbl(df, stem = paste0("time_alpha__", .clean_name(nm)))
})

# ==============================================================================
# B) Enrichment plot — big_group × timepoint_recoded (full dataset)
# ==============================================================================
.group_map <- c(kid_serum = "I", mom_serum = "M", mom_milk = "MB")
.time_map  <- c(
  T0 = "P12", T1 = "P28", T2 = "B",  T3 = "W2",
  T4 = "M1",  T5 = "M2",  T6 = "M3", T7 = "M6", T8 = "M12"
)

cohort_code <- function(x) {
  parts <- strsplit(as.character(x), "\\*")
  bg <- trimws(vapply(parts, `[`, "", 1))
  tp <- trimws(vapply(parts, `[`, "", 2))
  g  <- unname(.group_map[bg]); t <- unname(.time_map[tp])
  g[is.na(g)] <- bg[is.na(g)]
  t[is.na(t)] <- tp[is.na(t)]
  paste0(g, "_", t)
}

combo_pal <- c(
  mom_serum__T0 = "#F6A1C7", mom_serum__T1 = "#E32877", mom_serum__T2 = "#8C0303",
  kid_serum__T2 = "#B7D68F", kid_serum__T6 = "#91BE8C", kid_serum__T8 = "#0F9C69",
  mom_milk__T3  = "#73ACF1", mom_milk__T4  = "#6CA5EB", mom_milk__T5  = "#5C97DE",
  mom_milk__T6  = "#4C8AD1", mom_milk__T7  = "#0063AB", mom_milk__T8  = "#0050A1"
)
pal_by_cohort <- setNames(unname(combo_pal), sub("__", " * ", names(combo_pal)))

plots_main <- plot_enrichment_counts(
  ps,
  group_cols        = c("big_group", "timepoint_recoded"),
  group_interaction = TRUE,
  annotation_size   = 7.5
)

p <- plots_main[["big_group * timepoint_recoded"]]
p$facet$params$labeller   <- ggplot2::labeller(Cohort = cohort_code)
p$facet$params$multi_line <- FALSE
p$data$fill_col            <- pal_by_cohort[as.character(p$data$Cohort)]
p <- p + ggplot2::theme(
  text       = ggplot2::element_text(size = 22),
  plot.title = ggplot2::element_blank()
)

.save_png(
  p,
  file.path(out_dir, "enrichment__big_group_timepoint_recoded.png"),
  height = 3500, width = 2200, dpi = 300
)

# ==============================================================================
# C) Kid-serum subgroup: metadata × timepoint_recoded
# ==============================================================================
meta_vars <- c(
  "delivery_mode", "infant_sex", "ever_breastfed", "feedmode_m3_bin",
  "smoking", "delivery_place", "siblings", "preterm", "risk_CD_bin"
)

dat_kid <- ps %>% dplyr::filter(.data$big_group == "kid_serum")

.time_map_kid   <- c(T2 = "B",  T6 = "M3", T8 = "M12")
.time_shade_kid <- c(B  = 0.60, M3 = 0.80, M12 = 1.00)
.base_blue      <- "#0050A1"
.base_red       <- "#E32877"

shade_col <- function(hex, amount = 1) {
  rgb <- grDevices::col2rgb(hex) / 255
  mix <- (1 - amount) * 1 + amount * rgb
  grDevices::rgb(mix[1], mix[2], mix[3])
}

recolor_plots_meta <- function(df) {
  stopifnot(all(c("Cohort", "fill_col") %in% names(df)))
  parts   <- strsplit(as.character(df$Cohort), "\\*")
  covar   <- trimws(vapply(parts, `[`, "", 1))
  time    <- trimws(vapply(parts, `[`, "", 2))
  levs    <- sort(unique(covar))
  base_map <- setNames(c(.base_blue, .base_red)[seq_along(levs)], levs)
  t_short  <- .time_map_kid[time];  t_short[is.na(t_short)] <- time
  shade    <- .time_shade_kid[t_short]; shade[is.na(shade)] <- 1
  new_fill  <- mapply(
    function(g, s) shade_col(base_map[[g]], s),
    covar, shade, USE.NAMES = FALSE
  )
  new_label <- paste0(covar, "_", t_short)
  df$fill_col <- new_fill
  attr(df, "._facet_labels_") <- setNames(new_label, df$Cohort)
  df
}

for (v in meta_vars) {
  dat_star <- dat_kid %>%
    dplyr::filter(!is.na(.data[[v]]), !is.na(.data$timepoint_recoded))

  nm_int <- paste0(v, " * timepoint_recoded")

  # --- enrichment plot --------------------------------------------------------
  plots_meta <- plot_enrichment_counts(
    dat_star,
    group_cols        = c(v, "timepoint_recoded"),
    group_interaction = TRUE,
    annotation_size   = 7.5
  )

  p_meta <- plots_meta[[nm_int]]
  p_meta$data <- recolor_plots_meta(p_meta$data)

  cohort_to_label <- function(x) {
    lab_map <- attr(p_meta$data, "._facet_labels_")
    out <- unname(lab_map[as.character(x)])
    out[is.na(out)] <- as.character(x)
    out
  }
  p_meta$facet$params$labeller   <- ggplot2::labeller(Cohort = cohort_to_label)
  p_meta$facet$params$multi_line <- FALSE
  p_meta <- p_meta + ggplot2::theme(
    text       = ggplot2::element_text(size = 22),
    plot.title = ggplot2::element_blank()
  )

  .save_png(
    p_meta,
    file.path(out_dir, paste0("enrichment__kid__", .clean_name(nm_int), ".png")),
    height = 1750, width = 2200
  )

  # --- alpha-diversity table --------------------------------------------------
  time_alpha_meta <- compute_alpha_diversity(
    dat_star,
    ranks             = c("peptide_id", "species", "genus", "family", "order", "class", "phylum"),
    group_cols        = c(v, "timepoint_recoded"),
    group_interaction = TRUE,
    carry_cols        = "subject_id",
    interaction_only  = TRUE
  )

  if (nm_int %in% names(time_alpha_meta)) {
    df_int <- time_alpha_meta[[nm_int]] %>%
      dplyr::filter(stringr::str_detect(.data$phip_interaction, "\\*"))
    if (nrow(df_int)) {
      .save_tbl(df_int, stem = paste0("time_alpha__kid__", .clean_name(nm_int)))
    } else {
      warning(sprintf("no interaction rows to save for '%s'", nm_int))
    }
  } else {
    warning(sprintf("interaction '%s' not found in alpha list", nm_int))
  }
}

message("Done. Results saved to: ", normalizePath(out_dir))
