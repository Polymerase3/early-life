#!/usr/bin/env Rscript

# Principal Coordinates Analysis (PCoA) on the Kulczynski distance matrix
# reused from the t-SNE workflow (R/03-tsne-multiple-groups.R).

suppressPackageStartupMessages({
  library(phiper)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

set.seed(15092025)

# ------------------------------------------------------------------
# Load PHIP object and metadata (same as t-SNE script)
# ------------------------------------------------------------------
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

source("R/99-utils.R")

ps_cmp <- ps
ps_cmp$data_long <- ps_cmp$data_long %>%
  filter(
    !is.na(subject_id),
    !is.na(big_group),
    !is.na(timepoint_recoded),
    !is.na(dyade_recoded)
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    timepoint_factor = timepoint_recoded,
    dyade = dyade_recoded
  )

sample_meta <- ps_cmp$data_long %>%
  select(
    sample_id,
    subject_id,
    big_group,
    timepoint_recoded,
    dyade_recoded,
    run_id,
    plate_id,
    months_since_birth
  ) %>%
  distinct() %>%
  collect()

# Descriptives: months_since_birth for kid serum T6 and T8
msb_desc <- sample_meta %>%
  filter(big_group == "kid_serum", timepoint_recoded %in% c("T6", "T8")) %>%
  group_by(timepoint_recoded) %>%
  summarise(
    n = sum(!is.na(months_since_birth)),
    mean = mean(months_since_birth, na.rm = TRUE),
    sd = sd(months_since_birth, na.rm = TRUE),
    min = min(months_since_birth, na.rm = TRUE),
    max = max(months_since_birth, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nMonths since birth (kid_serum):\n")
print(msb_desc)


# ------------------------------------------------------------------
# Distance matrix (reuse cached version from t-SNE)
# ------------------------------------------------------------------
dir.create("results/tsne/data", recursive = TRUE, showWarnings = FALSE)
dist_cache_path <- file.path("results", "tsne", "data", "dist_bc.rds")

if (file.exists(dist_cache_path)) {
  dist_bc <- readRDS(dist_cache_path)
} else {
  dist_bc <- phiper:::compute_distance(
    ps_cmp,
    value_col = "exist",
    method_normalization = "none",
    distance = "kulczynski",
    n_threads = 10
  )
  saveRDS(dist_bc, dist_cache_path)
}

D <- as.matrix(dist_bc)
stopifnot(identical(rownames(D), colnames(D)))

# ------------------------------------------------------------------
# Align metadata to distance matrix order
# ------------------------------------------------------------------
nms <- rownames(D)
parsed <- sample_meta %>%
  filter(sample_id %in% nms) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  arrange(match(sample_id, nms)) %>%
  rename(id_full = sample_id)
stopifnot(identical(parsed$id_full, nms))

# ------------------------------------------------------------------
# Compute PCoA
# ------------------------------------------------------------------
pcoa_res <- phiper::compute_pcoa(dist_bc)

# ------------------------------------------------------------------
# Prepare coordinates + metadata and plot PC1 vs PC2
# ------------------------------------------------------------------
pcoa_pal <- c(
  mom_serum_T0 = "#F6A1C7",
  mom_serum_T1 = "#E32877",
  mom_serum_T2 = "#8C0303",
  kid_serum_T2 = "#B7D68F",
  kid_serum_T6 = "#91BE8C",
  kid_serum_T8 = "#0F9C69",
  mom_milk_T3 = "#73ACF1",
  mom_milk_T4 = "#6CA5EB",
  mom_milk_T5 = "#5C97DE",
  mom_milk_T6 = "#4C8AD1",
  mom_milk_T7 = "#0063AB",
  mom_milk_T8 = "#0050A1"
)

pcoa_df <- pcoa_res$sample_coords %>%
  rename(id_full = sample_id) %>%
  left_join(parsed, by = "id_full") %>%
  mutate(
    group_time = paste(big_group, timepoint_recoded, sep = "_")
  )

pcoa_centroids <- pcoa_df %>%
  group_by(group_time) %>%
  summarise(
    PCoA1 = mean(PCoA1, na.rm = TRUE),
    PCoA2 = mean(PCoA2, na.rm = TRUE),
    .groups = "drop"
  )

kid_order <- c("kid_serum_T2", "kid_serum_T6", "kid_serum_T8")
kid_path <- pcoa_centroids %>%
  filter(group_time %in% kid_order) %>%
  mutate(group_time = factor(group_time, levels = kid_order)) %>%
  arrange(group_time)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group_time)) +
  geom_point(size = 2.2, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.8, alpha = 0.4, show.legend = FALSE) +
  geom_point(
    data = pcoa_centroids %>% filter(!grepl("^kid_serum", group_time)),
    aes(x = PCoA1, y = PCoA2, color = group_time),
    size = 12
  ) +
  geom_point(
    data = pcoa_centroids %>% filter(grepl("^kid_serum", group_time)),
    aes(x = PCoA1, y = PCoA2, color = group_time),
    size = 12
  ) +
  labs(
    x = "PCoA1",
    y = "PCoA2",
    color = "Group_Time"
  ) +
  scale_color_manual(values = pcoa_pal, drop = FALSE) +
  theme_phip() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  coord_equal()

p_pcoa

dir.create("results/other_plots", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "results/other_plots/pcoa_pc1_pc2.svg",
  plot = p_pcoa,
  width = 12, height = 12, units = "in", dpi = 300, bg = "white"
)

cat("Saved PCoA plot to results/pcoa/pcoa_pc1_pc2.svg\n")

# ------------------------------------------------------------------
# Compute NMDS on the same distance matrix and plot NMDS1 vs NMDS2
# ------------------------------------------------------------------
if (!requireNamespace("vegan", quietly = TRUE)) {
  warning("Package 'vegan' is not installed; skipping NMDS.")
} else {
  nmds_res <- tryCatch(
    vegan::metaMDS(
      dist_bc,
      k = 2,
      trymax = 100,
      autotransform = FALSE,
      trace = FALSE
    ),
    error = function(e) {
      warning("NMDS computation failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(nmds_res)) {
    nmds_df <- as.data.frame(vegan::scores(nmds_res, display = "sites")) %>%
      tibble::rownames_to_column(var = "id_full") %>%
      left_join(parsed, by = "id_full") %>%
      mutate(group_time = paste(big_group, timepoint_recoded, sep = "_"))

    nmds_centroids <- nmds_df %>%
      group_by(group_time) %>%
      summarise(
        NMDS1 = mean(NMDS1, na.rm = TRUE),
        NMDS2 = mean(NMDS2, na.rm = TRUE),
        .groups = "drop"
      )

    p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = group_time)) +
      geom_point(size = 2.2, alpha = 0.8) +
      stat_ellipse(type = "norm", level = 0.68, linewidth = 0.8, alpha = 0.4, show.legend = FALSE) +
      geom_point(
        data = nmds_centroids %>% filter(!grepl("^kid_serum", group_time)),
        aes(x = NMDS1, y = NMDS2, color = group_time),
        size = 12
      ) +
      geom_point(
        data = nmds_centroids %>% filter(grepl("^kid_serum", group_time)),
        aes(x = NMDS1, y = NMDS2, color = group_time),
        size = 12
      ) +
      labs(
        x = "NMDS1",
        y = "NMDS2",
        color = "Group_Time",
        subtitle = paste0("Stress = ", signif(nmds_res$stress, 3))
      ) +
      scale_color_manual(values = pcoa_pal, drop = FALSE) +
      theme_phip() +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      ) +
      coord_equal()

    ggsave(
      filename = "results/other_plots/nmds1_nmds2.svg",
      plot = p_nmds,
      width = 12, height = 12, units = "in", dpi = 300, bg = "white"
    )

    cat("Saved NMDS plot to results/other_plots/nmds1_nmds2.svg\n")
  }
}
