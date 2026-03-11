################################################################################
########################### RICHNESS PLOT FIGURE 1 #############################
################################################################################
source("R/99-utils.R")

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

## computing alpha diversity
## this is the main workhorse for computing the alpha diversity in phiper - it
## does everything, calculate the coefs and absolute richness by group/time and
## combined groups and carries the columns important for later plotting
ranks <- c("peptide_id", "species", "genus", "order")
alpha_df <- compute_alpha_diversity(
  ps,
  group_col   = c("big_group", "timepoint_recoded"),
  group_interaction = TRUE,
  ranks       = ranks,
  carry_cols  = c("months_since_birth", "subject_id")
)

# define group-specific knots (months) and colors
group_time_palettes <- list(
  mom_serum = list(
    knots = c(-7, -3, 0),
    cols = c("#F6A1C7", "#E32877", "#8C0303")
  ),
  kid_serum = list(
    knots = c(0, 6, 12),
    cols = c("#B7D68F", "#68A689", "#0f9c69")
  ), # "#06402B"
  mom_milk = list(
    knots = c(0, 6, 12, 16),
    cols = c("#7bb3f8", "#0063ab", "#0050a1", "#002765")
  ) # "#372569"
)

# returns df with a new column `point_color` (hex per row)
make_group_time_colors <- function(df,
                                   group_col = "group", # can be NULL now
                                   time_col = "months_since_birth",
                                   palettes = group_time_palettes,
                                   clamp = TRUE,
                                   # choose a palette when group_col is NULL:
                                   palette_if_null = c("auto", "mom_serum", "kid_serum", "mom_milk"),
                                   custom_if_null = NULL # vector of hex OR list(knots, cols)
) {
  palette_if_null <- match.arg(palette_if_null)

  # extract time, robust-ish
  t_raw <- df[[time_col]]
  t_num <- if (is.numeric(t_raw)) t_raw else suppressWarnings(as.numeric(t_raw))

  # ---------- helpers ----------
  build_grad <- function(pal, tt) {
    if (is.null(pal) || is.na(tt)) {
      return(NA_character_)
    }
    rng <- range(pal$knots, finite = TRUE)
    vals <- rescale(pal$knots, to = c(0, 1), from = rng)
    grad <- gradient_n_pal(pal$cols, values = vals)
    ts <- rescale(tt, to = c(0, 1), from = rng)
    if (clamp) ts <- pmin(pmax(ts, 0), 1)
    grad(ts)
  }

  auto_palette_from_range <- function(tt) {
    r <- range(tt, finite = TRUE)
    if (!all(is.finite(r))) r <- c(0, 1)
    if (r[2] <= 0) {
      return(palettes$mom_serum)
    }
    if (r[1] >= 0 && r[2] <= 12) {
      return(palettes$kid_serum)
    }
    if (r[1] >= 0 && r[2] > 12) {
      return(palettes$mom_milk)
    }
    # spans negative->positive: stitch through 0
    if (r[2] <= 12) {
      list(
        knots = c(palettes$mom_serum$knots, palettes$kid_serum$knots[-1]),
        cols  = c(palettes$mom_serum$cols, palettes$kid_serum$cols[-1])
      )
    } else {
      list(
        knots = c(palettes$mom_serum$knots, palettes$mom_milk$knots[-1]),
        cols  = c(palettes$mom_serum$cols, palettes$mom_milk$cols[-1])
      )
    }
  }

  make_from_hex_vec <- function(tt, cols) {
    r <- range(tt, finite = TRUE)
    if (!all(is.finite(r))) r <- c(0, 1)
    list(knots = seq(r[1], r[2], length.out = length(cols)), cols = cols)
  }

  # ---------- main logic ----------
  if (is.null(group_col)) {
    # pick the single palette
    pal_single <-
      if (!is.null(custom_if_null)) {
        if (is.list(custom_if_null)) custom_if_null else make_from_hex_vec(t_num, custom_if_null)
      } else if (palette_if_null == "auto") {
        auto_palette_from_range(t_num)
      } else {
        palettes[[palette_if_null]]
      }

    df$point_color <- vapply(t_num, function(tt) build_grad(pal_single, tt), character(1))
    return(df)
  }

  # with groups: use the per-group palettes as before
  g <- df[[group_col]]
  df$point_color <- mapply(function(grp, tt) build_grad(palettes[[as.character(grp)]], tt),
                           g, t_num,
                           USE.NAMES = FALSE
  )
  df
}

################################################################################
######################### PLOTTING #############################################
################################################################################
# i plot it using the phiper function plot_alpha_diversity; it is not finished
# and is a dev-version, actually as the whole package, cause i didnt have time
# to finish it...; but nvm it actually works for now, the user interface is only
# not the best

set.seed(14082025)

# richness + GAM smooth
p1 <- alpha_df %>%
  plot_alpha_diversity(
    metric          = "richness",
    group_col       = "big_group",
    rank_col        = "rank",
    filter_ranks    = "peptide_id",
    time_col        = "months_since_birth",
    continuous_mode = "gam",
    gam_k           = 9,
    ci_method       = "bootstrap",
    ci_level        = 0.95,
    boot_R          = 5000,
    boot_seed       = 2137,
    boot_progress   = TRUE,
    ci_fill         = "grey10",
    ci_alpha        = 0.30,
    point_alpha     = 0.0
  ) +
  theme(
    text = element_text(family = "monte", size = 20),
    legend.position = "top", # move it above the plot
    legend.box = "horizontal", # keep items in a row
    legend.margin = margin(b = -5), # fine-tune spacing) +
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    axis.title.y = element_text(lineheight = 0.33, hjust = 0.5),
    plot.margin = margin(t = 7, r = 3, b = 8, l = 9)
  ) +
  scale_colour_manual(
    name = "Group",
    labels = c("Kid serum", "Mom milk", "Mom serum"),
    values = c("#68a689", "#07338c", "#8c0303")
  ) +
  ylim(0, 2260) +
  scale_x_continuous(breaks = seq(-8, 16, 2)) +
  ylab("Number of enriched\n peptides per subject") +
  xlab("Months since birth") +
  geom_vline(xintercept = 4, color = "red", linewidth = 1, linetype = "dashed") +
  annotate(
    geom = "text", x = -5.8, y = 1900, label = "M_P12", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = -2.2, y = 2100, label = "M_P28", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 1.4, y = 550, label = "M_B/I_B", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 0.9, y = 2150, label = "BM_M1", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 3, y = 2225, label = "BM_M3", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 6, y = 2250, label = "BM_M6", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 12, y = 2050, label = "BM_M12", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 5.1, y = 600, label = "I_M3", size = 4,
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 14.1, y = 850, label = "I_M12", size = 4,
    family = "monte"
  )

# make colored points
pts <- alpha_df$big_group %>%
  dplyr::filter(rank == "peptide_id") %>%
  make_group_time_colors(
    group_col = "big_group",
    time_col = "months_since_birth"
  )

# overlay with a separate color scale just for points
p_colored <- p1 +
  ggnewscale::new_scale_color() +
  geom_point(
    data = pts,
    aes(x = months_since_birth, y = richness, color = point_color),
    size = 1.6, alpha = 0.7, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_identity()

ggsave("results/other_plots/panel1d_richness_final.svg",
       height = 10,
       width = 28, units = "cm", dpi = 300, bg = "white"
)

################################################################################
################################ STABPLOT ######################################
################################################################################
source("099-old_phiper_utils.R") # i used old functions which were not in the
# phiper package anymore

set.seed(14082025)

## create phip_data PREDICTS; MAY BLOW UP YOUR RAM IF YOU HAVE LOWER THAN 64GB
withr::with_preserve_seed({
  ps1 <- phip_convert(
    data_long_path = "peptides_joined.parquet",
    backend = "duckdb",
    peptide_library = TRUE,
    subject_id = "SUBJECTID",
    sample_id = "sample_id",
    peptide_id = "peptide_name",
    fold_change = "fold_change",
    timepoint = "timepoint",
    materialise_table = TRUE,
    auto_expand = TRUE,
    n_cores = 12
  )
})

# add binary exist columns + filter only to the healthy adults
ps1 <- add_exist(ps1)
ps1 %<>% filter(GROUP == "Healthy controls")

# 1) compute kulczynski for babies data
stab_df <- compute_repertoire_stability(
  x = ps_merged_box_bin,
  group_col = "big_group",
  time_col = "months_since_birth",
  similarity = "kulczynski",
  drop_self = TRUE
)

# 1) compute kulczynski for PREDICTS data
stab_df2 <- compute_repertoire_stability(
  x = ps1,
  group_col = NULL,
  time_col = "timepoint",
  similarity = "kulczynski",
  drop_self = TRUE
)

# get the marker tables
df <- tibble::as_tibble(stab_df)
marker_tbl <- attr(stab_df, "baseline_marker", exact = TRUE)
drop_self <- isTRUE(attr(stab_df, "drop_self", exact = TRUE))

# block the group levels
grp_levels <- unique(df$group)
df$group <- factor(df$group, levels = grp_levels)
if (is.data.frame(marker_tbl) && nrow(marker_tbl)) {
  marker_tbl$group <- factor(marker_tbl$group, levels = grp_levels)
}

# --- bootstrap GAM per group: k=9, R=500, level=0.95, seed=2137 ---------------
## helper
preds <- split(df, df$group) |>
  purrr::map_dfr(function(d) {
    out <- .bootstrap_gam_one(
      as.data.frame(d),
      xvar = "time",
      yvar = "similarity",
      k_req = 9,
      R = 5000,
      level = 0.95,
      seed = 2137,
      progress = TRUE,
      nonneg = TRUE
    )
    out$group <- unique(d$group)
    out
  })
preds$group <- factor(preds$group, levels = grp_levels)
preds$.grp <- preds$group

# marker fill colors (unnamed -> map to levels)
# its not important anymore, as sasha wanted the markers removed
marker_fill_colors <- c("#07338c", "#8c0303", "#68a689")
if (length(marker_fill_colors) == length(grp_levels) &&
    is.null(names(marker_fill_colors))) {
  names(marker_fill_colors) <- grp_levels
}

# --- final ggplot --------------------------------------------------------------
p_stab <- ggplot() +
  # in this setting we do not draw points (point_alpha = 0)
  geom_ribbon(
    data = preds,
    aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
    fill = "grey40", alpha = 0.30, colour = NA, inherit.aes = FALSE
  ) +
  geom_line(
    data = preds,
    aes(x = .data$.x, y = .data$.y, color = .data$group, group = .data$.grp),
    linewidth = 1
  ) +
  labs(color = "Group") +
  theme_phip() +
  # baseline markers with custom fills
  (if (isTRUE(drop_self) && is.data.frame(marker_tbl) && nrow(marker_tbl)) ggnewscale::new_scale_fill() else NULL) +
  (if (isTRUE(drop_self) && is.data.frame(marker_tbl) && nrow(marker_tbl) && FALSE) {
    geom_point(
      data = marker_tbl,
      aes(x = .data$x_marker, y = .data$y_marker, fill = .data$group),
      shape = 21, size = 2.5, stroke = 0.6, colour = "black",
      inherit.aes = FALSE, show.legend = FALSE
    )
  } else {
    NULL
  }) +
  (if (isTRUE(drop_self) && is.data.frame(marker_tbl) && nrow(marker_tbl)) {
    scale_fill_manual(values = marker_fill_colors, guide = "none")
  } else {
    NULL
  }) +
  coord_cartesian(ylim = c(0, 1)) +
  # scale_x_continuous(breaks = seq(-8, 16, 2)) +
  ylab("Repertoire stability (Kulczynski)") +
  xlab("Time since baseline (months)") +
  scale_colour_manual(
    name = "Group",
    labels = c("Kid serum", "Mom milk", "Mom serum"),
    values = c(
      "mom_serum" = "#8c0303",
      "mom_milk" = "#07338c",
      "kid_serum" = "#68a689"
    )
  )

# 1) years -> months
df2 <- stab_df2 %>%
  transmute(
    subject_id,
    group = "all",
    time = (as.numeric(time) * 12) + 16,
    similarity
  )

# 2) Bootstrap GAM df2
preds2 <- .bootstrap_gam_one(
  as.data.frame(df2),
  xvar = "time",
  yvar = "similarity",
  k_req = 9,
  R = 5000,
  level = 0.95,
  seed = 2137,
  progress = TRUE,
  nonneg = TRUE
)

# 3) overlay on the same plot: gray ribbon + black dashed line
p_stab2 <- p_stab +
  ggplot2::geom_ribbon(
    data = preds2,
    ggplot2::aes(x = .x, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "grey40", alpha = 0.20
  ) +
  ggplot2::geom_line(
    data = preds2,
    ggplot2::aes(x = .x, y = .y, linetype = "Overall (df2)"),
    inherit.aes = FALSE,
    linewidth = 1.2,
    color = "#9a37de",
    alpha = 0.75
  ) +
  ggplot2::scale_linetype_manual(
    name = NULL,
    values = c("Overall (df2)" = "longdash"),
    labels = "Healthy adults"
  ) #+
# scale_x_break(c(-9.5 * 12 + 16, -4.5 * 12 + 16))

p_stab3 <- p_stab2 +
  ggplot2::geom_vline(
    xintercept = -104,
    linetype = "dotted",
    linewidth = 1,
    alpha = 0.6,
    color = "#9a37de"
  ) +
  ggplot2::geom_vline(
    xintercept = -6,
    linetype = "dotted",
    linewidth = 1,
    alpha = 0.6,
    color = "#8c0303"
  ) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    linewidth = 1,
    alpha = 0.6,
    color = "#07338c"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 1,
    alpha = 0.6,
    color = "#68a689"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-106, 17),
    breaks = c(-104, -32, seq(-10, 16, 2)),
    minor_breaks = NULL
  ) +
  ggbreak::scale_x_break(
    c(-103, -33, -31, -9),
    scales = "fixed",
    space  = 0.02
  ) +
  ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())

# --- 1) POINTS for stab_df (3 groups) -----------------------------------------
pts1 <- stab_df %>%
  make_group_time_colors(group_col = "group", time_col = "time")

p_stab4 <- p_stab3 +
  ggnewscale::new_scale_color() +
  geom_point(
    data = pts1,
    aes(x = time, y = similarity, color = point_color),
    size = 1.4, alpha = 0.5, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_identity()

pts2 <- stab_df2 %>%
  mutate(time2 = time * 12 + 16) %>%
  make_group_time_colors(
    group_col = NULL,
    time_col = "time2",
    custom_if_null = c("#e0087f", "#e0087f", "#d037de", "#9a37de")
  )

p_stab5 <- p_stab4 +
  ggnewscale::new_scale_color() +
  geom_jitter(
    data = pts2, width = 0.30,
    aes(x = time2, y = similarity, color = point_color),
    size = 1.4, alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_identity() +
  guides(
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      keywidth = unit(0.1, "cm"),
      keyheight = unit(0.9, "cm"),
      override.aes = list(linewidth = 3)
    )
  ) +
  annotate(
    geom = "text", label = "Adults\n baseline", x = -104, y = 0.1,
    size = 12, color = "#9a37de", family = "monte", fontface = "bold",
    lineheight = 0.35
  ) +
  annotate(
    geom = "text", label = "Mom serum\n baseline", x = -7.5, y = 0.1,
    size = 12, color = "#8c0303", family = "monte", fontface = "bold",
    lineheight = 0.35
  ) +
  annotate(
    geom = "text", label = "Kid serum\n baseline", x = -1.5, y = 0.1,
    size = 12, color = "#68a689", family = "monte", fontface = "bold",
    lineheight = 0.35
  ) +
  annotate(
    geom = "text", label = "Mom milk\n baseline", x = 1.8, y = 0.1,
    size = 12, color = "#07338c", family = "monte", fontface = "bold",
    lineheight = 0.35
  ) +
  annotate(
    geom = "text", x = -32, y = 0.935, label = "Adults after\n 6 years",
    size = 10, lineheight = 0.35, fontface = "bold",
    family = "monte"
  ) +
  annotate(
    geom = "text", x = -8, y = 0.92, label = "Adults after\n 8 years",
    size = 10, lineheight = 0.35, fontface = "bold",
    family = "monte"
  ) +
  annotate(
    geom = "text", x = 16, y = 0.835, label = "Adults after\n 10 years",
    size = 10, lineheight = 0.35, fontface = "bold",
    family = "monte"
  ) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    legend.title = element_text(size = 55),
    legend.text = element_text(size = 50),
    legend.spacing.x = grid::unit(0.1, "cm"),
    legend.spacing.y = grid::unit(0.06, "cm"),
    legend.key.width = grid::unit(1.1, "cm"),
    legend.key.height = grid::unit(0.6, "cm"),
    legend.box.spacing = grid::unit(0.08, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.position = "top",
    text = element_text(family = "monte", size = 55),
    axis.title.y = element_text(
      margin = margin(r = -1.7, unit = "lines"),
      lineheight = 0.01,
      hjust = 0.33,
      vjust = 3.2
    ),
    axis.title.x = element_text(
      margin = margin(
        t = -1.3,
        b = 0.2,
        unit = "lines"
      ),
      vjust = 0.12,
      hjust = 0.57
    )
  )

ggsave("results/stab-rich_plots/stability_natmed_final.png",
       height = 20, width = 40, units = "cm", dpi = 300, bg = "white"
)
ggsave("results/stab-rich_plots/stability_natmed_legend.png",
       height = 15, width = 60, units = "cm", dpi = 300, bg = "white"
)

###############################################################################
############################ BETA DIVERSITY TO SHARE ###########################
################################################################################
