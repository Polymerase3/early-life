# =============================================================================
# setup
# =============================================================================
# load required packages
library(phiper)
library(dplyr)
library(knitr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# =============================================================================
# 1) data import
# =============================================================================
# read parquet and build phip_data object
# create PHIP data object from input data with reproducible seed
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

# extract peptide library for interactive plot tooltips
peplib <- phiper::get_peptide_library(ps) %>% dplyr::collect() %>% as.data.frame()

# =============================================================================
# 2) case table
# =============================================================================
# compute unique samples per run/plate/group/timepoint
# count unique samples per group/timepoint/run/plate
cases_tbl <- ps$data_long %>%
  dplyr::distinct(
    .data$sample_id,
    .data$big_group,
    .data$timepoint_recoded,
    .data$run_id,
    .data$plate_id
  ) %>%
  dplyr::count(
    .data$big_group,
    .data$timepoint_recoded,
    .data$run_id,
    .data$plate_id,
    name = "n_samples"
  ) %>%
  dplyr::arrange(
    .data$big_group,
    .data$timepoint_recoded,
    .data$run_id,
    .data$plate_id
  ) %>%
  dplyr::collect()

# reshape to wide for excel-friendly output
cases_wide <- cases_tbl %>%
  dplyr::mutate(
    group_tp = paste(.data$big_group, .data$timepoint_recoded, sep = "*")
  ) %>%
  dplyr::select("run_id", "plate_id", "group_tp", "n_samples") %>%
  tidyr::pivot_wider(
    names_from = "group_tp",
    values_from = "n_samples",
    values_fill = 0
  ) %>%
  dplyr::arrange(.data$run_id, .data$plate_id)

# write crosstable to disk
out_dir <- "results/run-plate-diagnostics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cases_wide, file = file.path(out_dir, "run_plate_crosstable.csv"), row.names = FALSE)

# render a quick table preview
knitr::kable(
  cases_wide,
  caption = "Unique samples (n_samples) by run/plate and group*timepoint"
)

# =============================================================================
# 3) alpha diversity
# =============================================================================
# compute alpha diversity per sample (interaction-only)
alpha_div <- phiper::compute_alpha_diversity(
  ps,
  group_cols = c("run_id", "plate_id", "big_group", "timepoint_recoded"),
  group_interaction = TRUE,
  interaction_only = TRUE,
  carry_cols = c("run_id", "plate_id", "big_group", "timepoint_recoded")
)

# format output and recode timepoints
alpha_tbl <- alpha_div[[1]] %>%
  dplyr::filter(.data$rank == "peptide_id") %>%
  dplyr::mutate(
    timepoint_label = dplyr::recode(
      .data$timepoint_recoded,
      "T0" = "P12",
      "T1" = "P28",
      "T2" = "B",
      "T3" = "W2",
      "T4" = "M1",
      "T5" = "M2",
      "T6" = "M3",
      "T7" = "M6",
      "T8" = "M12"
    ),
    timepoint_label = factor(
      .data$timepoint_label,
      levels = c("P12", "P28", "B", "W2", "M1", "M2", "M3", "M6", "M12")
    )
  )

# enable montserrat font for plots
phiper::phip_use_montserrat()

# apply richness stability filter
alpha_tbl_filt <- alpha_tbl %>% dplyr::filter(.data$richness < 3000)

# -----------------------------------------------------------------------------
# 3a) marginal run analysis (most general)
# -----------------------------------------------------------------------------
# run_id-only marginal plots
run_levels <- sort(unique(alpha_tbl_filt$run_id))
run_richness_stats <- alpha_tbl_filt %>%
  rstatix::pairwise_t_test(richness ~ run_id, p.adjust.method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "run_id")

run_shannon_stats <- alpha_tbl_filt %>%
  rstatix::pairwise_t_test(shannon_diversity ~ run_id, p.adjust.method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "run_id")

run_richness_plot <- ggplot(
  alpha_tbl_filt,
  aes(x = run_id, y = richness, fill = run_id)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_violin(fill = NA, linewidth = 0.6, alpha = 1, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 3, color = "black", shape = 16) +
  ggpubr::stat_pvalue_manual(
    run_richness_stats,
    label = "p.adj.signif",
    hide.ns = FALSE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::scale_fill_phip() +
  phiper::scale_color_phip() +
  labs(x = "Run", y = "Richness", fill = "Run", color = "Run") +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

run_shannon_plot <- ggplot(
  alpha_tbl_filt,
  aes(x = run_id, y = shannon_diversity, fill = run_id)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_violin(fill = NA, linewidth = 0.6, alpha = 1, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 3, color = "black", shape = 16) +
  ggpubr::stat_pvalue_manual(
    run_shannon_stats,
    label = "p.adj.signif",
    hide.ns = FALSE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::scale_fill_phip() +
  phiper::scale_color_phip() +
  labs(x = "Run", y = "Shannon diversity", fill = "Run", color = "Run") +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

ggsave(
  filename = file.path(out_dir, "alpha_richness_by_run.png"),
  plot = run_richness_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "alpha_shannon_by_run.png"),
  plot = run_shannon_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 3b) marginal run * plate analysis (more specific)
# -----------------------------------------------------------------------------
# run * plate marginal plots
run_plate_counts <- alpha_tbl_filt %>%
  dplyr::count(.data$run_id, .data$plate_id, name = "n") %>%
  dplyr::filter(.data$n >= 2)

alpha_run_plate <- alpha_tbl_filt %>%
  dplyr::inner_join(run_plate_counts, by = c("run_id", "plate_id"))

run_plate_richness_stats <- alpha_run_plate %>%
  dplyr::group_by(.data$run_id) %>%
  rstatix::pairwise_t_test(richness ~ plate_id, p.adjust.method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "run_id", group = "plate_id", dodge = 0.8)

run_plate_shannon_stats <- alpha_run_plate %>%
  dplyr::group_by(.data$run_id) %>%
  rstatix::pairwise_t_test(shannon_diversity ~ plate_id, p.adjust.method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "run_id", group = "plate_id", dodge = 0.8)

run_plate_richness_plot <- ggplot(
  alpha_run_plate,
  aes(x = run_id, y = richness, fill = plate_id)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    alpha = 0.5,
    size = 3,
    color = "black",
    shape = 16,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)
  ) +
  ggpubr::stat_pvalue_manual(
    run_plate_richness_stats,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::scale_fill_phip() +
  phiper::scale_color_phip() +
  labs(x = "Run", y = "Richness", fill = "Plate", color = "Plate") +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

run_plate_shannon_plot <- ggplot(
  alpha_run_plate,
  aes(x = run_id, y = shannon_diversity, fill = plate_id)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    alpha = 0.5,
    size = 3,
    color = "black",
    shape = 16,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)
  ) +
  ggpubr::stat_pvalue_manual(
    run_plate_shannon_stats,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::scale_fill_phip() +
  phiper::scale_color_phip() +
  labs(x = "Run", y = "Shannon diversity", fill = "Plate", color = "Plate") +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

ggsave(
  filename = file.path(out_dir, "alpha_richness_by_run_plate.png"),
  plot = run_plate_richness_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "alpha_shannon_by_run_plate.png"),
  plot = run_plate_shannon_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 3c) run * plate * group * timepoint analysis (most specific)
# -----------------------------------------------------------------------------
# cache run order and group ids for dodged plotting
run_levels <- sort(unique(alpha_tbl_filt$run_id))
run_comparisons <- if (length(run_levels) > 1) {
  utils::combn(run_levels, 2, simplify = FALSE)
} else {
  list()
}

alpha_tbl_filt <- alpha_tbl_filt %>%
  dplyr::mutate(
    run_id = factor(.data$run_id, levels = run_levels),
    group_id = interaction(.data$timepoint_label, .data$run_id, sep = "__", drop = TRUE)
  )

# keep only groups with >= 2 runs for stats
alpha_tbl_test <- alpha_tbl_filt %>%
  dplyr::group_by(.data$big_group, .data$timepoint_label) %>%
  dplyr::filter(dplyr::n_distinct(.data$run_id) >= 2) %>%
  dplyr::ungroup()

# helper to compute pairwise tests and map bracket positions
build_stats <- function(df, y_col, p_base) {
  test_results <- df %>%
    dplyr::group_by(.data$big_group, .data$timepoint_label) %>%
    rstatix::pairwise_wilcox_test(
      stats::as.formula(paste0(y_col, " ~ run_id")),
      p.adjust.method = "holm",
      exact = FALSE
    ) %>%
    rstatix::add_significance("p.adj")

  if (nrow(test_results) == 0) {
    return(test_results)
  }

  p_build <- ggplot2::ggplot_build(p_base)
  panel_map <- p_build$layout$layout[, c("PANEL", "big_group"), drop = FALSE]
  box_map <- p_build$data[[1]] %>%
    dplyr::select(.data$PANEL, .data$group, .data$x) %>%
    dplyr::distinct() %>%
    dplyr::left_join(panel_map, by = "PANEL")

  group_levels <- levels(df$group_id)
  group_lookup <- tibble::tibble(
    group = seq_along(group_levels),
    group_id = group_levels
  ) %>%
    tidyr::separate(.data$group_id, into = c("timepoint_label", "run_id"), sep = "__") %>%
    dplyr::mutate(run_id = factor(.data$run_id, levels = run_levels))

  box_map <- box_map %>%
    dplyr::left_join(group_lookup, by = "group")

  y_max_tbl <- df %>%
    dplyr::group_by(.data$big_group, .data$timepoint_label) %>%
    dplyr::summarise(y_max = max(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(step = pmax(0.05 * .data$y_max, 0.1))

  test_results <- test_results %>%
    dplyr::left_join(y_max_tbl, by = c("big_group", "timepoint_label")) %>%
    dplyr::group_by(.data$big_group, .data$timepoint_label) %>%
    dplyr::arrange(.data$p.adj, .by_group = TRUE) %>%
    dplyr::mutate(y.position = .data$y_max + .data$step * dplyr::row_number()) %>%
    dplyr::ungroup()

  test_results <- test_results %>%
    dplyr::left_join(
      box_map %>%
        dplyr::select(.data$big_group, .data$timepoint_label, .data$run_id, xmin = .data$x),
      by = c("big_group", "timepoint_label", "group1" = "run_id")
    ) %>%
    dplyr::left_join(
      box_map %>%
        dplyr::select(.data$big_group, .data$timepoint_label, .data$run_id, xmax = .data$x),
      by = c("big_group", "timepoint_label", "group2" = "run_id")
    )

  test_results
}

# base plot for richness
richness_base <- ggplot(
  alpha_tbl_filt,
  aes(x = timepoint_label, y = richness, fill = run_id, group = group_id)
) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.6) +
  phiper::scale_fill_phip() +
  facet_wrap(~big_group, scales = "free_x") +
  labs(x = "Timepoint", y = "Richness", fill = "Run") +
  phiper::theme_phip()

# compute bracket positions and add to plot
richness_stats <- build_stats(alpha_tbl_test, "richness", richness_base)

richness_plot <- richness_base +
  ggpubr::stat_pvalue_manual(
    richness_stats,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::scale_fill_phip() +
  facet_wrap(~big_group, scales = "free_x") +
  labs(x = "Timepoint", y = "Richness", fill = "Run") +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

# base plot for shannon diversity
shannon_base <- ggplot(
  alpha_tbl_filt,
  aes(x = timepoint_label, y = shannon_diversity, fill = run_id, group = group_id)
) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.6) +
  phiper::scale_fill_phip() +
  facet_wrap(~big_group, scales = "free_x") +
  labs(x = "Timepoint", y = "Shannon diversity", fill = "Run") +
  phiper::theme_phip()

# compute bracket positions and add to plot
shannon_stats <- build_stats(alpha_tbl_test, "shannon_diversity", shannon_base)

shannon_plot <- shannon_base +
  ggpubr::stat_pvalue_manual(
    shannon_stats,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 18
  ) +
  phiper::theme_phip() +
  theme(text = element_text(size = 60),
        legend.key.size = grid::unit(1.2, "cm"))

# save plots
ggsave(
  filename = file.path(out_dir, "alpha_richness_by_timepoint_run.png"),
  plot = richness_plot,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "alpha_shannon_by_timepoint_run.png"),
  plot = shannon_plot,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 3d) richness summary table (median)
# -----------------------------------------------------------------------------
# summarize richness per run/plate/group/timepoint
richness_summary <- alpha_tbl_filt %>%
  dplyr::group_by(.data$run_id, .data$plate_id, .data$big_group, .data$timepoint_recoded) %>%
  dplyr::summarise(
    med_richness = stats::median(.data$richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    med = round(.data$med_richness, 1),
    group_tp = paste(.data$big_group, .data$timepoint_recoded, sep = "*")
  ) %>%
  dplyr::select(.data$run_id, .data$plate_id, .data$group_tp, .data$med)

richness_wide <- richness_summary %>%
  tidyr::pivot_wider(
    names_from = "group_tp",
    values_from = "med"
  ) %>%
  dplyr::arrange(.data$run_id, .data$plate_id)

write.csv(
  richness_wide,
  file = file.path(out_dir, "richness_median_crosstable.csv"),
  row.names = FALSE
)

knitr::kable(
  richness_wide,
  caption = "Richness median by run/plate and group*timepoint"
)

# =============================================================================
# 4) POP run comparisons (static scatterplots)
# =============================================================================
# helper to extract a plain tibble from phiper outputs
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

pop_dir <- file.path(out_dir, "POP_run_compare")
dir.create(pop_dir, recursive = TRUE, showWarnings = FALSE)

run_pairs <- list(
  c("R25", "R26"),
  c("R25", "R27"),
  c("R26", "R27")
)

for (pair in run_pairs) {
  pair_label <- paste(pair, collapse = "_vs_")

  # stratified sampling by pooled prevalence bin (SQL-side)
  base_tbl <- ps$data_long %>%
    dplyr::filter(.data$run_id %in% pair)

  prev_tbl <- base_tbl %>%
    dplyr::group_by(.data$peptide_id) %>%
    dplyr::summarise(prev_pct = 100 * mean(.data$exist), .groups = "drop") %>%
    dplyr::mutate(
      prev_bin = dplyr::case_when(
        .data$prev_pct < 10 ~ "0-10",
        .data$prev_pct < 20 ~ "10-20",
        .data$prev_pct < 30 ~ "20-30",
        .data$prev_pct < 35 ~ "30-35",
        .data$prev_pct < 40 ~ "35-40",
        .data$prev_pct < 45 ~ "40-45",
        .data$prev_pct < 50 ~ "45-50",
        .data$prev_pct < 55 ~ "50-55",
        .data$prev_pct < 60 ~ "55-60",
        .data$prev_pct < 65 ~ "60-65",
        .data$prev_pct < 70 ~ "65-70",
        .data$prev_pct < 75 ~ "70-75",
        .data$prev_pct < 80 ~ "75-80",
        .data$prev_pct < 85 ~ "80-85",
        .data$prev_pct < 90 ~ "85-90",
        .data$prev_pct < 95 ~ "90-95",
        .data$prev_pct < 100 ~ "95-100",
        TRUE ~ "100"
      )
    )

  sampled_peptides <- prev_tbl %>%
    dplyr::mutate(
      .rand = dbplyr::sql("random()"),
      sample_rate = 1
    ) %>%
    dplyr::group_by(.data$prev_bin) %>%
    dplyr::filter(.data$.rand < .data$sample_rate) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$peptide_id)

  # collect only the columns needed for prevalence comparison
  data_frameworks <- base_tbl %>%
    dplyr::semi_join(sampled_peptides, by = "peptide_id") %>%
    dplyr::select("sample_id", "peptide_id", "run_id", "exist") %>%
    dplyr::collect() %>%
    dplyr::mutate(group_char = dplyr::if_else(.data$run_id == pair[1], pair[1], pair[2]))

  prev_res_pep <- phiper::ph_prevalence_compare(
    x                 = data_frameworks,
    group_cols        = "group_char",
    rank_cols         = "peptide_id",
    compute_ratios_db = TRUE,
    paired            = FALSE,
    parallel          = TRUE,
    collect           = TRUE
  )
  rm(data_frameworks)

  pep_tbl <- extract_tbl(prev_res_pep)
  rm(prev_res_pep)

  p_static <- phiper::scatter_static(
    df = pep_tbl,
    rank = "peptide_id",
    xlab = pep_tbl$group1[1],
    ylab = pep_tbl$group2[1],
    point_size = 3,
    jitter_width_pp = 0.15,
    jitter_height_pp = 0.15,
    point_alpha = 0.85,
    font_family = "Montserrat",
    font_size = 23
  ) +
    ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
    phiper::theme_phip() +
    ggplot2::theme(
      plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
      text        = ggplot2::element_text(size = 23, family = "Montserrat"),
      legend.key.size = grid::unit(1.2, "cm")
    )

  ggsave(
    filename = file.path(pop_dir, paste0(pair_label, "_static.svg")),
    plot = p_static,
    dpi = 300,
    height = 30,
    width = 30,
    units = "cm",
    bg = "white"
  )
  rm(p_static)

  p_inter <- phiper::scatter_interactive(
    df = pep_tbl,
    rank = "peptide_id",
    xlab = pep_tbl$group1[1],
    ylab = pep_tbl$group2[1],
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
    file = file.path(pop_dir, paste0(pair_label, "_interactive.html")),
    selfcontained = TRUE
  )
  rm(p_inter, pep_tbl)
  gc()
}
