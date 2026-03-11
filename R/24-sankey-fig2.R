#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
  library(Rtsne)
})

# ── Config ────────────────────────────────────────────────────────────────────

comparisons <- c(
  "kid_serum_T2_vs_kid_serum_T6",
  "kid_serum_T6_vs_mom_milk_T6"
)
comp_labels <- c(
  "kid_serum_T2_vs_kid_serum_T6" = "T2 vs T6\n(kid serum)",
  "kid_serum_T6_vs_mom_milk_T6"  = "T6 kid serum\nvs T6 mom milk"
)
sig_thr   <- 0.05
base_dir  <- file.path("results", "results")
out_dir   <- file.path("results", "sankey_fig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load delta_table_fdr.csv ──────────────────────────────────────────────────

load_fdr <- function(cmp_id) {
  path <- file.path(base_dir, cmp_id, "DELTA_framework", "delta_table_fdr.csv")
  if (!file.exists(path)) stop("Missing file: ", path)
  d <- readr::read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet")
  d <- d %>% dplyr::select(-dplyr::any_of("...1"))

  required <- c("feature", "p_perm", "T_obs_stand")
  missing  <- setdiff(required, names(d))
  if (length(missing) > 0) stop("Missing columns in ", cmp_id, ": ", paste(missing, collapse = ", "))

  d %>%
    dplyr::mutate(
      feature     = as.character(feature),
      p_perm      = as.numeric(p_perm),
      T_obs_stand = as.numeric(T_obs_stand)
    ) %>%
    # keep one row per feature (best = lowest p_perm)
    dplyr::group_by(feature) %>%
    dplyr::slice_min(order_by = p_perm, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(feature, p_perm, T_obs_stand)
}

dlist <- lapply(comparisons, load_fdr)
names(dlist) <- comparisons

# ── Build binary significance matrix ─────────────────────────────────────────

all_features <- sort(unique(unlist(lapply(dlist, `[[`, "feature"))))
message("Total unique features across both comparisons: ", length(all_features))

sig_mat <- tibble::tibble(feature = all_features)
for (cmp in comparisons) {
  sig_mat <- sig_mat %>%
    dplyr::left_join(
      dlist[[cmp]] %>% dplyr::rename(p = p_perm, z = T_obs_stand),
      by = "feature"
    ) %>%
    dplyr::mutate(
      !!cmp := dplyr::if_else(!is.na(p) & p < sig_thr, "Significant", "Not significant"),
      !!paste0(cmp, "_signed") := dplyr::case_when(
        !is.na(p) & p < sig_thr & z > 0  ~ "Sig +",
        !is.na(p) & p < sig_thr & z < 0  ~ "Sig -",
        !is.na(p) & p < sig_thr & z == 0 ~ "Sig +",  # edge case: treat zero as positive
        TRUE ~ "Not significant"
      )
    ) %>%
    dplyr::select(-p, -z)
}

readr::write_csv(sig_mat, file.path(out_dir, "binary_significance_matrix.csv"))
message("Saved: binary_significance_matrix.csv (", nrow(sig_mat), " features)")

# ── Confusion matrix counts ───────────────────────────────────────────────────

cmp1 <- comparisons[1]
cmp2 <- comparisons[2]

confusion <- sig_mat %>%
  dplyr::count(.data[[cmp1]], .data[[cmp2]], name = "n_features") %>%
  dplyr::rename(status_cmp1 = all_of(cmp1), status_cmp2 = all_of(cmp2))

readr::write_csv(confusion, file.path(out_dir, "confusion_counts.csv"))
message("Saved: confusion_counts.csv")
print(confusion)

# ── Jaccard and Kulczynski distances ─────────────────────────────────────────

sig_A <- sig_mat[["feature"]][sig_mat[[cmp1]] == "Significant"]
sig_B <- sig_mat[["feature"]][sig_mat[[cmp2]] == "Significant"]

n_A       <- length(sig_A)
n_B       <- length(sig_B)
n_inter   <- length(intersect(sig_A, sig_B))
n_union   <- length(union(sig_A, sig_B))

jaccard_sim   <- if (n_union == 0) NA_real_ else n_inter / n_union
jaccard_dist  <- if (is.na(jaccard_sim)) NA_real_ else 1 - jaccard_sim

kulc_sim  <- if (n_A == 0 || n_B == 0) NA_real_ else
  0.5 * (n_inter / n_A + n_inter / n_B)
kulc_dist <- if (is.na(kulc_sim)) NA_real_ else 1 - kulc_sim

dist_tbl <- tibble::tibble(
  metric         = c("Jaccard similarity", "Jaccard distance", "Kulczynski similarity", "Kulczynski distance"),
  value          = c(jaccard_sim, jaccard_dist, kulc_sim, kulc_dist),
  n_sig_cmp1     = n_A,
  n_sig_cmp2     = n_B,
  n_intersection = n_inter,
  n_union        = n_union
)

readr::write_csv(dist_tbl, file.path(out_dir, "jaccard_kulczynski_distances.csv"))
message("Saved: jaccard_kulczynski_distances.csv")
print(dist_tbl)

# ── Sankey / alluvial plot ────────────────────────────────────────────────────

sig_levels <- c("Significant", "Not significant")

flow_df <- confusion %>%
  dplyr::mutate(
    status_cmp1 = factor(status_cmp1, levels = sig_levels),
    status_cmp2 = factor(status_cmp2, levels = sig_levels)
  )

lodes_df <- ggalluvial::to_lodes_form(
  data  = flow_df,
  axes  = c("status_cmp1", "status_cmp2"),
  key   = "comparison",
  value = "status",
  id    = "alluvium"
)

lodes_df <- lodes_df %>%
  dplyr::mutate(
    comparison_label = dplyr::recode(
      comparison,
      status_cmp1 = comp_labels[[cmp1]],
      status_cmp2 = comp_labels[[cmp2]]
    )
  )

# ordered factor so the x axis follows cmp1 → cmp2
lodes_df$comparison_label <- factor(
  lodes_df$comparison_label,
  levels = unname(comp_labels[comparisons])
)

sig_colors <- c(
  "Significant"     = "#2E7D32",
  "Not significant" = "#B0BEC5"
)

p_sankey <- ggplot(
  lodes_df,
  aes(
    x        = comparison_label,
    stratum  = status,
    alluvium = alluvium,
    y        = n_features
  )
) +
  ggalluvial::geom_flow(aes(fill = status), alpha = 0.55) +
  ggalluvial::geom_stratum(
    aes(fill = status),
    width = 1 / 6,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3.5
  ) +
  scale_fill_manual(values = sig_colors, breaks = sig_levels, drop = FALSE) +
  scale_x_discrete(expand = c(0.12, 0.12)) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Feature significance overlap between comparisons",
    subtitle = paste0(
      "p_perm < ", sig_thr, "; total features = ", length(all_features),
      "; sig in T2vsT6 = ", n_A, "; sig in T6vsMilk = ", n_B,
      "; intersection = ", n_inter
    ),
    x    = "Comparison",
    y    = "Number of features",
    fill = "Significance status"
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

svg_path <- file.path(out_dir, "sankey_fig2.svg")
ggsave(svg_path, p_sankey, width = 28, height = 18, units = "cm", dpi = 300, bg = "white")
message("Saved: sankey_fig2.svg")

# ── Sankey signed by T_obs_stand ─────────────────────────────────────────────

signed_levels <- c("Sig +", "Sig -", "Not significant")

signed_col1 <- paste0(cmp1, "_signed")
signed_col2 <- paste0(cmp2, "_signed")

signed_confusion <- sig_mat %>%
  dplyr::count(.data[[signed_col1]], .data[[signed_col2]], name = "n_features") %>%
  dplyr::rename(status_cmp1 = all_of(signed_col1), status_cmp2 = all_of(signed_col2)) %>%
  dplyr::mutate(
    status_cmp1 = factor(status_cmp1, levels = signed_levels),
    status_cmp2 = factor(status_cmp2, levels = signed_levels)
  )

readr::write_csv(signed_confusion, file.path(out_dir, "confusion_counts_signed.csv"))
message("Saved: confusion_counts_signed.csv")
print(signed_confusion)

signed_lodes <- ggalluvial::to_lodes_form(
  data  = signed_confusion,
  axes  = c("status_cmp1", "status_cmp2"),
  key   = "comparison",
  value = "status",
  id    = "alluvium"
) %>%
  dplyr::mutate(
    comparison_label = factor(
      dplyr::recode(
        comparison,
        status_cmp1 = comp_labels[[cmp1]],
        status_cmp2 = comp_labels[[cmp2]]
      ),
      levels = unname(comp_labels[comparisons])
    ),
    status = factor(status, levels = signed_levels)
  )

signed_colors <- c(
  "Sig +"          = "#2E7D32",
  "Sig -"          = "#C62828",
  "Not significant" = "#B0BEC5"
)

p_sankey_signed <- ggplot(
  signed_lodes,
  aes(
    x        = comparison_label,
    stratum  = status,
    alluvium = alluvium,
    y        = n_features
  )
) +
  ggalluvial::geom_flow(aes(fill = status), alpha = 0.55) +
  ggalluvial::geom_stratum(
    aes(fill = status),
    width = 1 / 6,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3.5
  ) +
  scale_fill_manual(values = signed_colors, breaks = signed_levels, drop = FALSE) +
  scale_x_discrete(expand = c(0.12, 0.12)) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Feature significance overlap between comparisons (signed by T_obs_stand)",
    subtitle = paste0(
      "p_perm < ", sig_thr, "; total features = ", length(all_features),
      "; sig in T2vsT6 = ", n_A, "; sig in T6vsMilk = ", n_B,
      "; intersection = ", n_inter
    ),
    x    = "Comparison",
    y    = "Number of features",
    fill = "Significance status"
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

svg_signed_path <- file.path(out_dir, "sankey_fig2_signed.svg")
ggsave(svg_signed_path, p_sankey_signed, width = 28, height = 18, units = "cm", dpi = 300, bg = "white")
message("Saved: sankey_fig2_signed.svg")

# ── Heatmap: T_obs_stand (kid serum | mom milk reversed | difference) ─────────
# Reverse sign for mom_milk comparison so both comparisons express the same
# directional meaning (higher = more abundant in kid serum T6 relative to
# the other group).  Difference = z_kid − z_milk_rev.

hm_df <- tibble::tibble(feature = all_features) %>%
  dplyr::left_join(
    dlist[[cmp1]] %>% dplyr::select(feature, p_perm, T_obs_stand) %>%
      dplyr::rename(p_cmp1 = p_perm, z_kid = T_obs_stand),
    by = "feature"
  ) %>%
  dplyr::left_join(
    dlist[[cmp2]] %>% dplyr::select(feature, p_perm, T_obs_stand) %>%
      dplyr::rename(p_cmp2 = p_perm, z_milk_rev = T_obs_stand),
    by = "feature"
  ) %>%
  # reverse sign for mom_milk comparison
  dplyr::mutate(z_milk_rev = -z_milk_rev) %>%
  # keep only union-significant features
  dplyr::filter(
    (!is.na(p_cmp1) & p_cmp1 < sig_thr) |
    (!is.na(p_cmp2) & p_cmp2 < sig_thr)
  ) %>%
  dplyr::mutate(diff = z_kid - z_milk_rev) %>%
  dplyr::select(feature, z_kid, z_milk_rev, diff)

message("Significant features for heatmap: ", nrow(hm_df))

# Ward clustering on the two T_obs_stand columns (NAs → 0 for distance)
hm_mat <- as.matrix(hm_df[, c("z_kid", "z_milk_rev")])
hm_mat[is.na(hm_mat)] <- 0
hc <- hclust(dist(hm_mat), method = "ward.D2")
feature_order <- hm_df$feature[hc$order]

panel_labels <- c(
  z_kid      = "T2 vs T6\n(kid serum)",
  z_milk_rev = "T6 kid serum vs T6 mom milk\n(sign reversed)",
  diff       = "Difference\n(kid − milk reversed)"
)

hm_long <- hm_df %>%
  tidyr::pivot_longer(
    cols      = c(z_kid, z_milk_rev, diff),
    names_to  = "panel",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    feature = factor(feature, levels = feature_order),
    panel   = factor(panel, levels = c("z_kid", "z_milk_rev", "diff"),
                     labels = unname(panel_labels))
  )

# symmetric colour limit (95th percentile of absolute values, excluding diff)
diff_label <- panel_labels[["diff"]]
z_abs_max <- quantile(
  abs(hm_long$value[hm_long$panel != diff_label]),
  na.rm = TRUE, probs = 0.95
)

p_heatmap <- ggplot(hm_long, aes(x = panel, y = feature, fill = value)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low      = "#053061",
    mid      = "white",
    high     = "#67001F",
    midpoint = 0,
    limits   = c(-z_abs_max, z_abs_max),
    oob      = scales::squish,
    name     = "T_obs_stand\n(or difference)"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_text(size = 8, face = "bold"),
    axis.title       = element_blank(),
    plot.title       = element_text(face = "bold"),
    legend.position  = "right",
    panel.border     = element_rect(color = "grey40", fill = NA)
  ) +
  labs(
    title    = "T_obs_stand heatmap — kid serum vs mom milk (Ward clustering)",
    subtitle = paste0(
      "Union significant features (p_perm < ", sig_thr, ", n = ", nrow(hm_df), "); ",
      "mom milk sign reversed; features ordered by ward.D2"
    )
  )

# height scales with number of features (min 15 cm)
hm_height_cm <- max(15, nrow(hm_df) * 0.06)

svg_hm_path <- file.path(out_dir, "heatmap_ward.svg")
ggsave(svg_hm_path, p_heatmap, width = 22, height = hm_height_cm, units = "cm", dpi = 300, bg = "white")
message("Saved: heatmap_ward.svg")

readr::write_csv(hm_df %>% dplyr::mutate(ward_order = match(feature, feature_order)),
                 file.path(out_dir, "heatmap_data.csv"))
message("Saved: heatmap_data.csv")

# ── Scatter plot: z_kid vs z_milk_rev ────────────────────────────────────────

# Label threshold: both axes exceed this absolute value AND in opposite quadrants
label_thr <- 1.0

scatter_df <- hm_df %>%
  dplyr::mutate(
    quadrant = dplyr::case_when(
      z_kid >= 0 & z_milk_rev >= 0 ~ "++ (concordant high)",
      z_kid <  0 & z_milk_rev <  0 ~ "-- (concordant low)",
      z_kid >= 0 & z_milk_rev <  0 ~ "+- (discordant)",
      z_kid <  0 & z_milk_rev >= 0 ~ "-+ (discordant)",
      TRUE ~ NA_character_
    ),
    discordant = grepl("discordant", quadrant) &
                   abs(z_kid) > label_thr & abs(z_milk_rev) > label_thr,
    label_conc_low = quadrant == "-- (concordant low)" &
                   abs(z_kid) > label_thr & abs(z_milk_rev) > label_thr
  )

quad_colors <- c(
  "++ (concordant high)" = "#B2182B",
  "-- (concordant low)"  = "#1F78B4",
  "+- (discordant)"      = "#6A3D9A",
  "-+ (discordant)"      = "#E6AB02"
)

r_pearson  <- cor(scatter_df$z_kid, scatter_df$z_milk_rev, use = "complete.obs", method = "pearson")
r_spearman <- cor(scatter_df$z_kid, scatter_df$z_milk_rev, use = "complete.obs", method = "spearman")

n_discord  <- sum(grepl("discordant", scatter_df$quadrant), na.rm = TRUE)
n_by_quad  <- table(scatter_df$quadrant)

scatter_notes <- c(
  paste0("n features: ", nrow(scatter_df)),
  paste0("Pearson r: ", round(r_pearson, 4)),
  paste0("Spearman rho: ", round(r_spearman, 4)),
  paste0("Discordant (|NES| > ", label_thr, " on both axes): ", n_discord),
  "",
  "Quadrant counts:",
  paste0("  ", names(n_by_quad), ": ", as.integer(n_by_quad)),
  "",
  "Axis definitions:",
  "  x: NES from kid_serum_T2_vs_kid_serum_T6  (group1=B, group2=M3; positive = higher at M3)",
  "  y: NES from kid_serum_T6_vs_mom_milk_T6, sign reversed  (original positive = higher in mom milk M3)"
)
writeLines(scatter_notes, file.path(out_dir, "scatter_concordance_notes.txt"))
message("Saved: scatter_concordance_notes.txt")

p_scatter <- ggplot(scatter_df, aes(x = z_kid, y = z_milk_rev)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "solid", color = "black") +
  # concordant low (--) drawn first so discordant points sit on top
  geom_point(
    data = dplyr::filter(scatter_df, quadrant == "++ (concordant high)"),
    aes(color = quadrant), size = 0.4, alpha = 0.65
  ) +
  geom_point(
    data = dplyr::filter(scatter_df, quadrant == "-- (concordant low)"),
    aes(color = quadrant), size = 0.4, alpha = 0.65
  ) +
  geom_point(
    data = dplyr::filter(scatter_df, grepl("discordant", quadrant)),
    aes(color = quadrant), size = 0.6, alpha = 0.90
  ) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(scatter_df, feature %in% c(
      "Bubalus bubalis", "milk", "Ovis aries", "Bos taurus", "Capra hircus",
      "Onchocercidae", "Brugia malayi", "Nicotiana tabacum", "Bacteroides nordii",
      "Phocaeicola coprocola", "Pecentumvirus LP0832", "Porcine epidemic diarrhea virus",
      "Pisum sativum", "Lupinus angustifolius", "Euroglyphus maynei", "Nectriaceae",
      "Fusarium", "Aggregatibacter segnis", "Prevotella sp. CAG:873", "ragweed",
      "Arthroderma benhamiae", "Oreochromis niloticus", "Bacteroides acidifaciens",
      "Rangifer tarandus", "Cervidae", "Azospirillum sp. 51_20"
    )),
    aes(label    = feature, color = quadrant),
    size         = 1.2,
    max.overlaps = Inf,
    segment.size = 0.15,
    show.legend  = FALSE,
    family       = "Montserrat"
  ) +
  scale_color_manual(values = quad_colors, name = NULL, drop = FALSE) +
  phiper::theme_phip(base_size = 5) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    axis.line        = element_line(colour = "black", linewidth = 0.3),
    legend.key.size  = unit(0.25, "cm")
  ) +
  labs(x = NULL, y = NULL)

svg_scatter_path <- file.path(out_dir, "scatter_concordance.svg")
ggsave(svg_scatter_path, p_scatter, width = 90, height = 65, units = "mm", dpi = 300, bg = "white")
message("Saved: scatter_concordance.svg")

feature_table <- scatter_df %>%
  dplyr::filter(quadrant %in% c("-- (concordant low)", "-+ (discordant)")) %>%
  dplyr::arrange(quadrant, z_kid) %>%
  dplyr::select(feature, quadrant, z_kid, z_milk_rev, diff)

readr::write_csv(feature_table, file.path(out_dir, "features_concordant_low_and_discordant.csv"))
message("Saved: features_concordant_low_and_discordant.csv (", nrow(feature_table), " features)")

message("\nAll outputs written to: ", out_dir)
message(" - binary_significance_matrix.csv")
message(" - confusion_counts.csv")
message(" - confusion_counts_signed.csv")
message(" - jaccard_kulczynski_distances.csv")
message(" - sankey_fig2.svg")
message(" - sankey_fig2_signed.svg")
message(" - heatmap_ward.svg")
message(" - heatmap_data.csv")
message(" - scatter_concordance.svg")
message(" - scatter_concordance_notes.txt")
message(" - features_concordant_low_and_discordant.csv")

# ── t-SNE on raw z-scores (z_kid, z_milk_rev) ────────────────────────────────

tsne_mat <- as.matrix(hm_df[, c("z_kid", "z_milk_rev")])
tsne_mat[is.na(tsne_mat)] <- 0

n_tsne        <- nrow(tsne_mat)
max_perplexity <- floor((n_tsne - 1) / 3)
perplexity     <- min(30, max_perplexity)

message("Running t-SNE on ", n_tsne, " features with perplexity=", perplexity)
set.seed(42)
tsne_fit <- Rtsne::Rtsne(
  tsne_mat,
  dims             = 2,
  perplexity       = perplexity,
  check_duplicates = FALSE,
  pca              = FALSE,
  verbose          = FALSE
)

tsne_df <- hm_df %>%
  dplyr::mutate(
    tsne1    = tsne_fit$Y[, 1],
    tsne2    = tsne_fit$Y[, 2],
    quadrant = dplyr::case_when(
      z_kid >= 0 & z_milk_rev >= 0 ~ "++ (concordant high)",
      z_kid <  0 & z_milk_rev <  0 ~ "-- (concordant low)",
      z_kid >= 0 & z_milk_rev <  0 ~ "+- (discordant)",
      z_kid <  0 & z_milk_rev >= 0 ~ "-+ (discordant)",
      TRUE ~ NA_character_
    ),
    discordant = grepl("discordant", quadrant) &
                   abs(z_kid) > label_thr & abs(z_milk_rev) > label_thr
  )

readr::write_csv(tsne_df, file.path(out_dir, "tsne_data.csv"))
message("Saved: tsne_data.csv")

p_tsne <- ggplot(tsne_df, aes(x = tsne1, y = tsne2)) +
  geom_point(aes(color = quadrant), size = 1.6, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data         = dplyr::filter(tsne_df, discordant),
    aes(label    = feature, color = quadrant),
    size         = 2.8,
    max.overlaps = 30,
    segment.size = 0.3,
    show.legend  = FALSE
  ) +
  scale_color_manual(values = quad_colors, name = "Quadrant") +
  theme_bw(base_size = 12) +
  labs(
    title    = "t-SNE of union-significant features (input: z_kid, z_milk_rev)",
    subtitle = paste0(
      "n = ", n_tsne, " features; perplexity = ", perplexity,
      "; NAs → 0; quadrant based on sign of z-scores"
    ),
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )

svg_tsne_path <- file.path(out_dir, "tsne_concordance.svg")
ggsave(svg_tsne_path, p_tsne, width = 28, height = 22, units = "cm", dpi = 300, bg = "white")
message("Saved: tsne_concordance.svg")
message(" - tsne_data.csv")
message(" - tsne_concordance.svg")

# ── Scatter: mean-centred Z-scores, ALL features (no p_perm filter) ───────────
# For each comparison independently, subtract the mean T_obs_stand across ALL
# features (before any significance filtering). The sign of cmp2 is reversed
# to keep the same directional convention as the main scatter above.

full_df <- tibble::tibble(feature = all_features) %>%
  dplyr::left_join(
    dlist[[cmp1]] %>% dplyr::select(feature, p_perm, T_obs_stand) %>%
      dplyr::rename(p_cmp1 = p_perm, z_kid_raw = T_obs_stand),
    by = "feature"
  ) %>%
  dplyr::left_join(
    dlist[[cmp2]] %>% dplyr::select(feature, p_perm, T_obs_stand) %>%
      dplyr::rename(p_cmp2 = p_perm, z_milk_raw = T_obs_stand),
    by = "feature"
  ) %>%
  # centre on ALL features, then filter to union-significant for plotting
  dplyr::mutate(
    z_kid      = z_kid_raw  - mean(z_kid_raw,  na.rm = TRUE),
    z_milk_rev = -(z_milk_raw - mean(z_milk_raw, na.rm = TRUE))
  ) %>%
  dplyr::filter(
    (!is.na(p_cmp1) & p_cmp1 < sig_thr) |
    (!is.na(p_cmp2) & p_cmp2 < sig_thr)
  ) %>%
  dplyr::select(feature, z_kid, z_milk_rev)

scatter_df_cent <- full_df %>%
  dplyr::mutate(
    quadrant = dplyr::case_when(
      z_kid >= 0 & z_milk_rev >= 0 ~ "++ (concordant high)",
      z_kid <  0 & z_milk_rev <  0 ~ "-- (concordant low)",
      z_kid >= 0 & z_milk_rev <  0 ~ "+- (discordant)",
      z_kid <  0 & z_milk_rev >= 0 ~ "-+ (discordant)",
      TRUE ~ NA_character_
    ),
    discordant = grepl("discordant", quadrant) &
                   abs(z_kid) > label_thr & abs(z_milk_rev) > label_thr,
    label_conc_low = quadrant == "-- (concordant low)" &
                   abs(z_kid) > label_thr & abs(z_milk_rev) > label_thr
  )

r_pearson_cent  <- cor(scatter_df_cent$z_kid, scatter_df_cent$z_milk_rev,
                       use = "complete.obs", method = "pearson")
r_spearman_cent <- cor(scatter_df_cent$z_kid, scatter_df_cent$z_milk_rev,
                       use = "complete.obs", method = "spearman")

n_discord_cent <- sum(grepl("discordant", scatter_df_cent$quadrant), na.rm = TRUE)
n_by_quad_cent <- table(scatter_df_cent$quadrant)

scatter_cent_notes <- c(
  "Mean-centred Z-scores; ALL features included (no p_perm filter)",
  paste0("n features: ", nrow(scatter_df_cent)),
  paste0("Mean subtracted per comparison before sign reversal of cmp2"),
  paste0("Pearson r: ", round(r_pearson_cent, 4)),
  paste0("Spearman rho: ", round(r_spearman_cent, 4)),
  paste0("Discordant (|NES| > ", label_thr, " on both axes): ", n_discord_cent),
  "",
  "Quadrant counts:",
  paste0("  ", names(n_by_quad_cent), ": ", as.integer(n_by_quad_cent)),
  "",
  "Axis definitions:",
  "  x: mean-centred NES from kid_serum_T2_vs_kid_serum_T6",
  "  y: mean-centred NES from kid_serum_T6_vs_mom_milk_T6, sign reversed"
)
writeLines(scatter_cent_notes, file.path(out_dir, "scatter_centred_notes.txt"))
readr::write_csv(scatter_df_cent, file.path(out_dir, "scatter_centred_data.csv"))

p_scatter_cent <- ggplot(scatter_df_cent, aes(x = z_kid, y = z_milk_rev)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "solid", color = "black") +
  geom_point(
    data = dplyr::filter(scatter_df_cent, quadrant == "++ (concordant high)"),
    aes(color = quadrant), size = 0.4, alpha = 0.65
  ) +
  geom_point(
    data = dplyr::filter(scatter_df_cent, quadrant == "-- (concordant low)"),
    aes(color = quadrant), size = 0.4, alpha = 0.65
  ) +
  geom_point(
    data = dplyr::filter(scatter_df_cent, grepl("discordant", quadrant)),
    aes(color = quadrant), size = 0.6, alpha = 0.90
  ) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(scatter_df_cent, feature %in% c(
      "Bubalus bubalis", "milk", "Ovis aries", "Bos taurus", "Capra hircus",
      "Onchocercidae", "Brugia malayi", "Nicotiana tabacum", "Bacteroides nordii",
      "Phocaeicola coprocola", "Pecentumvirus LP0832", "Porcine epidemic diarrhea virus",
      "Pisum sativum", "Lupinus angustifolius", "Euroglyphus maynei", "Nectriaceae",
      "Fusarium", "Aggregatibacter segnis", "Prevotella sp. CAG:873", "ragweed",
      "Arthroderma benhamiae", "Oreochromis niloticus", "Bacteroides acidifaciens",
      "Rangifer tarandus", "Cervidae", "Azospirillum sp. 51_20"
    )),
    aes(label = feature, color = quadrant),
    size         = 1.2,
    max.overlaps = Inf,
    segment.size = 0.15,
    show.legend  = FALSE,
    family       = "Montserrat"
  ) +
  scale_color_manual(values = quad_colors, name = NULL, drop = FALSE) +
  phiper::theme_phip(base_size = 5) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    axis.line        = element_line(colour = "black", linewidth = 0.3),
    legend.key.size  = unit(0.25, "cm")
  ) +
  labs(x = NULL, y = NULL)

svg_scatter_cent_path <- file.path(out_dir, "scatter_centred.svg")
ggsave(svg_scatter_cent_path, p_scatter_cent,
       width = 90, height = 65, units = "mm", dpi = 300, bg = "white")
message("Saved: scatter_centred.svg")
message(" - scatter_centred_data.csv")
message(" - scatter_centred_notes.txt")
