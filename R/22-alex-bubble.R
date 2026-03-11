#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scico)
})

args <- commandArgs(trailingOnly = TRUE)

summary_dir_default <- file.path("results", "alex_flags", "summary")
in_path_default <- file.path(summary_dir_default, "delta_table_all_microbiome_recoded.csv")

in_path <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else in_path_default
out_dir <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else summary_dir_default

nominal_threshold <- 0.05
bh_threshold <- 0.10
timepoint_levels <- c("m6", "m9", "m12", "combined")
plot_base_size <- 5
plot_font_family <- "Nimbus Sans"

if (!file.exists(in_path)) {
  stop("Missing input table: ", in_path)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")

required_cols <- c("microbiome_tested", "timepoint", "p_perm", "p_adj_bh", "D_standardized", "max_delta")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "))
}

d <- d %>%
  mutate(
    microbiome_tested = as.character(.data$microbiome_tested),
    timepoint = as.character(.data$timepoint),
    p_perm = as.numeric(.data$p_perm),
    p_adj_bh = as.numeric(.data$p_adj_bh),
    D_standardized = as.numeric(.data$D_standardized),
    max_delta = as.numeric(.data$max_delta)
  ) %>%
  filter(.data$timepoint %in% timepoint_levels) %>%
  group_by(.data$microbiome_tested, .data$timepoint) %>%
  arrange(.data$p_perm, desc(abs(.data$D_standardized)), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

sig_microbes <- d %>%
  group_by(.data$microbiome_tested) %>%
  summarise(has_nominal = any(!is.na(.data$p_perm) & .data$p_perm < nominal_threshold), .groups = "drop") %>%
  filter(.data$has_nominal) %>%
  pull(.data$microbiome_tested)

if (length(sig_microbes) == 0L) {
  stop("No microbiome_tested has nominal significance at p_perm < ", nominal_threshold)
}

bubble_df <- d %>%
  filter(.data$microbiome_tested %in% sig_microbes) %>%
  mutate(timepoint = factor(.data$timepoint, levels = timepoint_levels)) %>%
  complete(
    microbiome_tested,
    timepoint = factor(timepoint_levels, levels = timepoint_levels),
    fill = list(
      p_perm = NA_real_,
      p_adj_bh = NA_real_,
      D_standardized = NA_real_,
      max_delta = NA_real_
    )
  )

preferred_row_levels <- c(
  "Blautia_faecis.t__SGB4820",
  "Eubacterium_rectale.t__SGB4933_group",
  "Lacrimispora_amygdalina.t__SGB4716",
  "Roseburia_faecis.t__SGB4925",
  "Anaerostipes_hadrus.t__SGB4540_group",
  "Enterocloster_aldensis.t__SGB4762",
  "Blautia_wexlerae.t__SGB4837_group",
  "Fusicatenibacter_saccharivorans.t__SGB4874",
  "Clostridiaceae_bacterium.t__SGB4269",
  "Clostridium_sp_AF34_10BH.t__SGB4914",
  "Clostridium_symbiosum.t__SGB4699",
  "Clostridium_neonatale.t__SGB6169",
  "Agathobaculum_butyriciproducens.t__SGB14993_group",
  "Faecalibacterium_prausnitzii.t__SGB15316_group",
  "Faecalibacterium_prausnitzii.t__SGB15342",
  "Bifidobacterium_longum.t__SGB17248",
  "Bifidobacterium_bifidum.t__SGB17256",
  "Erysipelatoclostridium_ramosum.t__SGB6744",
  "Eggerthella_lenta.t__SGB14809",
  "Veillonella_atypica.t__SGB6936",
  "Veillonella_dispar.t__SGB6952",
  "Klebsiella_oxytoca.t__SGB10118",
  "Phocaeicola_dorei.t__SGB1815",
  "Bacteroides_ovatus.t__SGB1871",
  "Lacticaseibacillus_rhamnosus.t__SGB7144",
  "Lacticaseibacillus_paracasei.t__SGB7142"
)

preferred_row_levels <- unique(preferred_row_levels)
present_microbes <- bubble_df %>%
  distinct(.data$microbiome_tested) %>%
  pull(.data$microbiome_tested)

missing_in_preferred <- setdiff(present_microbes, preferred_row_levels)
extra_in_preferred <- setdiff(preferred_row_levels, present_microbes)

if (length(extra_in_preferred) > 0L) {
  warning("Some preferred microbiome_tested entries are not present and will be ignored.")
}
if (length(missing_in_preferred) > 0L) {
  warning("Some microbiome_tested entries are missing from preferred order and will be appended.")
}

row_levels <- c(
  preferred_row_levels[preferred_row_levels %in% present_microbes],
  sort(missing_in_preferred)
)

bubble_df <- bubble_df %>%
  mutate(
    microbiome_row = factor(.data$microbiome_tested, levels = rev(row_levels)),
    is_nominal = !is.na(.data$p_perm) & .data$p_perm < nominal_threshold,
    is_bh = !is.na(.data$p_adj_bh) & .data$p_adj_bh <= bh_threshold,
    max_delta_pct = 100 * .data$max_delta,
    max_delta_bin = case_when(
      is.na(.data$max_delta_pct) ~ "0-5%",
      .data$max_delta_pct <= 5 ~ "0-5%",
      .data$max_delta_pct <= 10 ~ "5-10%",
      .data$max_delta_pct <= 20 ~ "10-20%",
      .data$max_delta_pct <= 30 ~ "20-30%",
      TRUE ~ ">30%"
    ),
    max_delta_bin = factor(.data$max_delta_bin, levels = c("0-5%", "5-10%", "10-20%", "20-30%", ">30%")),
    x_num = as.numeric(.data$timepoint)
  )

line_gap <- 0.12
line_df <- bubble_df %>%
  distinct(.data$microbiome_row) %>%
  crossing(segment = c("m6_m9", "m9_m12", "m12_combined")) %>%
  mutate(
    x = case_when(
      .data$segment == "m6_m9" ~ 1 + line_gap,
      .data$segment == "m9_m12" ~ 2 + line_gap,
      TRUE ~ 3 + line_gap
    ),
    xend = case_when(
      .data$segment == "m6_m9" ~ 2 - line_gap,
      .data$segment == "m9_m12" ~ 3 - line_gap,
      TRUE ~ 4 - line_gap
    ),
    y = .data$microbiome_row,
    yend = .data$microbiome_row
  ) %>%
  select("x", "xend", "y", "yend")

# Match the scaling approach used in R/14-sankey-summary.R
plot_height_mm <- max(120, 2.6 * length(row_levels) + 40)
old_height_cm <- 2.5 * max(12, 0.35 * length(row_levels))
size_scale <- plot_height_mm / (old_height_cm * 10)

p <- ggplot(bubble_df, aes(x = .data$x_num, y = .data$microbiome_row)) +
  geom_segment(
    data = line_df,
    aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
    inherit.aes = FALSE,
    linewidth = 0.14,
    color = "grey45",
    alpha = 0.75
  ) +
  geom_point(
    aes(size = .data$max_delta_bin, fill = .data$D_standardized),
    shape = 21,
    color = "grey20",
    stroke = 0.2 * size_scale,
    na.rm = TRUE
  ) +
  geom_point(
    data = dplyr::filter(bubble_df, .data$is_nominal),
    aes(x = .data$x_num, y = .data$microbiome_row),
    inherit.aes = FALSE,
    shape = 16,
    size = 1.35 * size_scale,
    color = "black"
  ) +
  geom_point(
    data = dplyr::filter(bubble_df, .data$is_bh),
    aes(x = .data$x_num, y = .data$microbiome_row),
    inherit.aes = FALSE,
    shape = 8,
    size = 2.65 * size_scale,
    color = "black",
    stroke = 0.45 * size_scale
  ) +
  scico::scale_fill_scico(
    palette = "roma",
    begin = 0.02,
    end = 0.85,
    direction = 1,
    name = "D standardized",
    limits = c(-3, NA),
    oob = scales::squish
  ) +
  scale_size_manual(
    values = c(
      "0-5%" = 4.97,
      "5-10%" = 6.46,
      "10-20%" = 7.95,
      "20-30%" = 9.44,
      ">30%" = 10.93
    ) * size_scale,
    name = "Max delta",
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = timepoint_levels,
    limits = c(0.8, 4.2)
  ) +
  labs(
    x = "Timepoint",
    y = "microbiome_tested"
  ) +
  theme_bw(base_size = plot_base_size, base_family = plot_font_family) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = plot_base_size),
    axis.text.x = element_text(size = plot_base_size, face = "bold"),
    axis.title = element_text(size = plot_base_size),
    legend.text = element_text(size = plot_base_size),
    legend.title = element_text(size = plot_base_size),
    legend.position = "right"
  )

out_long <- file.path(out_dir, "alex_bubble_long.csv")
out_plot <- file.path(out_dir, "alex_bubbleplot.svg")

readr::write_csv(bubble_df, out_long)

ggplot2::ggsave(
  filename = out_plot,
  plot = p,
  width = 120,
  height = plot_height_mm,
  units = "mm",
  dpi = 300,
  device = "svg",
  bg = "white"
)

cat("Input:", in_path, "\n")
cat("Output long table:", out_long, "\n")
cat("Output plot:", out_plot, "\n")
cat("Microbiome features retained (>=1 nominal timepoint):", length(sig_microbes), "\n")
