#!/usr/bin/env Rscript

# Minimal scatterplot for moms: birth (T0) vs week P12 (T1)
# Uses existing POP/DELTA outputs for the comparison "mom_serum_T0_vs_mom_serum_T1".
# Generates one static PNG and one interactive HTML scatter of percent positives.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(phiper)
  library(RColorBrewer)
})

# ------------------------------------------------------------------
# Paths and comparison
# ------------------------------------------------------------------
base_dir <- file.path("results", "results")
# T0 corresponds to P12; T2 corresponds to birth (B). We want P12 vs birth.
comparison <- "mom_serum_T0_vs_mom_serum_T2"
single_pep_path <- file.path(base_dir, comparison, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path)) {
  stop("Missing POP single_peptide.csv for ", comparison, " at ", single_pep_path)
}

# Output locations
out_dir <- file.path("results", "other_plots")

static_path <- file.path(out_dir, "mom_P12_vs_birth.svg")
interactive_path <- file.path(out_dir, "mom_P12_vs_birth.html")

# Peptide metadata and shared color map for category flags
peplib <- phiper::get_peptide_meta() %>% dplyr::collect()
highlight_points <- tibble::tribble(
  ~rank, ~feature,
  "anno_food_item", "milk",
  "is_flagellum", "TRUE",
  "anno_food_item", "wheat",
  "anno_viruses_bacteriophage", "EnterovirusA",
  "anno_viruses_bacteriophage", "HHV4",
  "anno_viruses_bacteriophage", "HHV-3",
  "genus", "Prevotella",
  "species", "Barnesiella intestinihominis",
  "genus", "Streptococcus",
  "anno_viruses_bacteriophage", "RhinovirusA",
  "anno_viruses_bacteriophage", "RhinovirusB",
  "species", "Parabacteroides distasonis",
  "genus", "Lactobacillus"
) %>%
  mutate(cat_label = feature)

# custom palette for highlight categories (13 distinct colors)
highlight_labels <- unique(highlight_points$cat_label)
palette_manual <- c(
  "#1B9E77", # teal
  "#B2182B", # dark red  (was orange)
  "#7570B3", # purple
  "#E7298A", # magenta
  "#c96b6b", # light red (was green)
  "#E6AB02", # mustard
  "#A6761D", # brown
  "#1F78B4", # blue
  "#B2DF8A", # light green
  "#FB9A99", # light red
  "#CAB2D6", # light purple
  "#6A3D9A", # dark violet
  "#FF7F00"  # dark orange
)

pal_fun <- colorRampPalette(palette_manual)
highlight_cols <- setNames(
  pal_fun(max(3, length(highlight_labels)))[seq_along(highlight_labels)],
  highlight_labels
)
build_highlight_peptides <- function(peplib_df, hp_tbl) {
  rows <- lapply(seq_len(nrow(hp_tbl)), function(i) {
    r <- hp_tbl$rank[i]
    f <- hp_tbl$feature[i]
    f_chr <- as.character(f)
    lbl <- hp_tbl$cat_label[i]
    if (!r %in% names(peplib_df)) return(NULL)
    hits <- peplib_df %>%
      filter(!is.na(.data[[r]])) %>%
      filter(as.character(.data[[r]]) == f_chr) %>%
      select(peptide_id) %>%
      mutate(cat_label = lbl)
    if (nrow(hits) == 0) return(NULL)
    hits
  })
  bind_rows(rows) %>% distinct(peptide_id, cat_label)
}
highlight_peptides <- build_highlight_peptides(peplib, highlight_points)
highlight_peptides_map <- highlight_peptides %>% rename(feature = peptide_id)

build_highlight_df <- function(df) {
  df %>%
    inner_join(highlight_peptides_map, by = "feature") %>%
    mutate(cat_label = factor(cat_label, levels = names(highlight_cols))) %>%
    tidyr::drop_na(cat_label, percent1, percent2, x_jit, y_jit)
}
add_categories <- function(df) {
  df %>%
    left_join(peplib, by = c("feature" = "peptide_id")) %>%
    mutate(
      is_milk          = anno_food_item == "milk",
      is_enterovirus   = genus == "Enterovirus",
      is_streptococcus = genus == "Streptococcus",
      is_rhinovirus    = !is.na(anno_viruses_bacteriophage) &
        grepl("Rhinovirus", anno_viruses_bacteriophage, ignore.case = TRUE),
      is_hhv3          = anno_viruses_bacteriophage %in% c("HHV-3", "HHV3"),
      is_hhv4          = anno_viruses_bacteriophage %in% c("HHV-4", "HHV4"),
      category = dplyr::case_when(
        is_milk ~ "milk",
        is_rhinovirus ~ "rhinovirus",
        is_hhv3 ~ "HHV-3",
        is_hhv4 ~ "HHV-4",
        is_enterovirus ~ "enterovirus",
        is_streptococcus ~ "streptococcus",
        TRUE ~ "other"
      )
    )
}
# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
pep <- readr::read_csv(single_pep_path, show_col_types = FALSE)

# Keep only the needed columns and drop rows without percentages
pep_small <- pep %>%
  select(feature, rank, percent1, percent2, category_rank_bh) %>%
  filter(!is.na(percent1) & !is.na(percent2))
highlight_p1 <- pep_small %>%
  mutate(x_jit = percent1, y_jit = percent2) %>%
  build_highlight_df()

# ------------------------------------------------------------------
# Static scatter
# ------------------------------------------------------------------
p <- ggplot(pep_small, aes(x = percent1, y = percent2)) +
  geom_point(alpha = 0.7, size = 1.8, color = "grey70") +
  geom_point(
    data = highlight_p1,
    inherit.aes = FALSE,
    aes(x = x_jit, y = y_jit, colour = cat_label),
    size = 3.6, alpha = 0.8, na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Mothers: P12 (T0) vs Birth (T2)",
    x = "Percent positive at P12 (T0)",
    y = "Percent positive at Birth (T2)"
  ) +
  scale_color_manual(values = highlight_cols, name = "Highlight", na.value = "grey50", drop = FALSE) +
  theme_phip() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

ggsave(static_path, p, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

# ------------------------------------------------------------------
# Interactive scatter with hover info
# ------------------------------------------------------------------
p_inter <- plot_ly(
  pep_small,
  x = ~percent1,
  y = ~percent2,
  type = "scatter",
  mode = "markers",
  marker = list(size = 6, color = "grey75", opacity = 0.35),
  text = ~paste0(
    "Feature: ", feature,
    "<br>Rank: ", rank,
    "<br>P12 (T0) %+: ", round(percent1, 2),
    "<br>Birth (T2) %+: ", round(percent2, 2),
    "<br>BH category: ", category_rank_bh
  ),
  hoverinfo = "text"
) %>%
  add_markers(
    data = highlight_p1,
    x = ~percent1, y = ~percent2,
    type = "scatter", mode = "markers",
    color = ~cat_label, colors = highlight_cols,
    marker = list(
      size = 11,
      symbol = "circle",
      opacity = 0.8,
      line = list(width = 0, color = "rgba(0,0,0,0)")
    ),
    text = ~paste0(
      "Feature: ", feature,
      "<br>Rank: ", rank,
      "<br>P12 (T0) %+: ", round(percent1, 2),
      "<br>Birth (T2) %+: ", round(percent2, 2)
    ),
    hoverinfo = "text",
    showlegend = TRUE
  ) %>%
  layout(
    title = "Mothers: P12 (T0) vs Birth (T2)",
    xaxis = list(title = "Percent positive at P12 (T0)", range = c(0, 100)),
    yaxis = list(title = "Percent positive at Birth (T2)", range = c(0, 100)),
    shapes = list(
      list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
           line = list(dash = "dash", color = "grey"))
    )
  )

htmlwidgets::saveWidget(p_inter, interactive_path, selfcontained = FALSE)

cat("Saved static scatter to:", static_path, "\n")
cat("Saved interactive scatter to:", interactive_path, "\n")

# ------------------------------------------------------------------
# Minimal scatterplot for moms vs kids: birth (T2) vs birth (T2)
# Uses existing POP/DELTA outputs for the comparison "mom_serum_T2_vs_kid_serum_T2".
# Generates one static PNG and one interactive HTML scatter of percent positives.

comparison2 <- "mom_serum_T2_vs_kid_serum_T2"
single_pep_path2 <- file.path(base_dir, comparison2, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path2)) {
  stop("Missing POP single_peptide.csv for ", comparison2, " at ", single_pep_path2)
}

static_path2 <- file.path(out_dir, "mom_birth_vs_kid_birth.svg")
interactive_path2 <- file.path(out_dir, "mom_birth_vs_kid_birth.html")

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
pep2 <- readr::read_csv(single_pep_path2, show_col_types = FALSE)

pep_small2 <- pep2 %>%
  select(feature, rank, percent1, percent2, category_rank_bh) %>%
  filter(!is.na(percent1) & !is.na(percent2)) %>%
  add_categories() %>%
  mutate(
    x_jit = jitter(percent1, amount = 0.5),
    y_jit = jitter(percent2, amount = 0.5),
    x_jit = pmin(pmax(x_jit, 0), 100),
    y_jit = pmin(pmax(y_jit, 0), 100)
  )
highlight_p2 <- pep_small2 %>% build_highlight_df()

# ------------------------------------------------------------------
# Static scatter
# ------------------------------------------------------------------
p2 <- ggplot(pep_small2, aes(x = x_jit, y = y_jit)) +
  geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
  geom_point(
    data = highlight_p2,
    inherit.aes = FALSE,
    aes(x = x_jit, y = y_jit, colour = cat_label),
    size = 3.6
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Mothers vs Kids: Birth (T2)",
    x = "Percent positive in Mothers at Birth (T2)",
    y = "Percent positive in Kids at Birth (T2)"
  ) +
  scale_color_manual(values = highlight_cols, name = "Highlight",
                    na.value = "grey50", drop = FALSE) +
  theme_phip() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank()
  )
p2
ggsave(static_path2, p2, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

# ------------------------------------------------------------------
# Interactive scatter with hover info
# ------------------------------------------------------------------
p_inter2 <- plot_ly() %>%
  add_markers(
    data = pep_small2,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    marker = list(size = 5, color = "grey80", opacity = 0.3),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  add_markers(
    data = highlight_p2,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    color = ~cat_label, colors = highlight_cols,
    marker = list(
      size = 12,
      symbol = "circle",
      opacity = 0.8,
      line = list(width = 0, color = "rgba(0,0,0,0)")
    ),
    text = ~paste0(
      "Feature: ", feature,
      "<br>Rank: ", rank,
      "<br>Mothers (T2) %+: ", round(percent1, 2),
      "<br>Kids (T2) %+: ", round(percent2, 2),
      "<br>BH category: ", category_rank_bh
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
  layout(
    title = "Mothers vs Kids: Birth (T2)",
    xaxis = list(title = "Percent positive in Mothers at Birth (T2)", range = c(0, 100)),
    yaxis = list(title = "Percent positive in Kids at Birth (T2)", range = c(0, 100)),
    shapes = list(
      list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
           line = list(dash = "dash", color = "grey"))
    ),
    legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
  )

htmlwidgets::saveWidget(p_inter2, interactive_path2, selfcontained = FALSE)

cat("Saved static scatter to:", static_path2, "\n")
cat("Saved interactive scatter to:", interactive_path2, "\n")

# ------------------------------------------------------------------
# Kids: birth (T2) vs follow-up (T6)
# ------------------------------------------------------------------
comparison3 <- "kid_serum_T2_vs_kid_serum_T6"
single_pep_path3 <- file.path(base_dir, comparison3, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path3)) {
  stop("Missing POP single_peptide.csv for ", comparison3, " at ", single_pep_path3)
}

static_path3 <- file.path(out_dir, "kid_T2_vs_T6.svg")
interactive_path3 <- file.path(out_dir, "kid_T2_vs_T6.html")

pep3 <- readr::read_csv(single_pep_path3, show_col_types = FALSE)

pep_small3 <- pep3 %>%
  select(feature, rank, percent1, percent2, category_rank_bh) %>%
  filter(!is.na(percent1) & !is.na(percent2)) %>%
  add_categories()

# jitter positions for background and markers
set.seed(123)
pep_small3 <- pep_small3 %>%
  mutate(
    x_jit = jitter(percent1, amount = 0.5),
    y_jit = jitter(percent2, amount = 0.5),
    x_jit = pmin(pmax(x_jit, 0), 100),
    y_jit = pmin(pmax(y_jit, 0), 100)
  )
highlight_p3 <- pep_small3 %>% build_highlight_df()

p3 <- ggplot(pep_small3, aes(x = x_jit, y = y_jit)) +
  geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
  geom_point(
    data = highlight_p3,
    inherit.aes = FALSE,
    aes(x = x_jit, y = y_jit, colour = cat_label),
    size = 3.6, alpha = 0.8, na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Kids: Birth (T2) vs T6",
    x = "Percent positive at Birth (T2)",
    y = "Percent positive at T6"
  ) +
  scale_color_manual(values = highlight_cols, name = "Highlight", na.value = "grey50", drop = FALSE) +
  theme_phip() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank()
  )
p3
ggsave(static_path3, p3, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

p_inter3 <- plot_ly() %>%
  add_markers(
    data = pep_small3,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    marker = list(size = 5, color = "grey80", opacity = 0.3),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  add_markers(
    data = highlight_p3,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    color = ~cat_label, colors = highlight_cols,
    marker = list(
      size = 12,
      symbol = "circle",
      opacity = 0.8,
      line = list(width = 0, color = "rgba(0,0,0,0)")
    ),
    text = ~paste0(
      "Feature: ", feature,
      "<br>Rank: ", rank,
      "<br>Birth (T2) %+: ", round(percent1, 2),
      "<br>T6 %+: ", round(percent2, 2),
      "<br>BH category: ", category_rank_bh
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
  layout(
    title = "Kids: Birth (T2) vs T6",
    xaxis = list(title = "Percent positive at Birth (T2)", range = c(0, 100)),
    yaxis = list(title = "Percent positive at T6", range = c(0, 100)),
    shapes = list(
      list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
           line = list(dash = "dash", color = "grey"))
    ),
    legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
  )

htmlwidgets::saveWidget(p_inter3, interactive_path3, selfcontained = FALSE)

cat("Saved static scatter to:", static_path3, "\n")
cat("Saved interactive scatter to:", interactive_path3, "\n")

# ------------------------------------------------------------------
# Kids: birth (T2) vs follow-up (T8)
# ------------------------------------------------------------------
comparison4 <- "kid_serum_T2_vs_kid_serum_T8"
single_pep_path4 <- file.path(base_dir, comparison4, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path4)) {
  stop("Missing POP single_peptide.csv for ", comparison4, " at ", single_pep_path4)
}

static_path4 <- file.path(out_dir, "kid_T2_vs_T8.svg")
interactive_path4 <- file.path(out_dir, "kid_T2_vs_T8.html")

pep4 <- readr::read_csv(single_pep_path4, show_col_types = FALSE)

pep_small4 <- pep4 %>%
  select(feature, rank, percent1, percent2, category_rank_bh) %>%
  filter(!is.na(percent1) & !is.na(percent2))

# reuse peplib and categorization logic
pep_small4 <- pep_small4 %>%
  add_categories() %>%
  mutate(
    x_jit = jitter(percent1, amount = 0.5),
    y_jit = jitter(percent2, amount = 0.5),
    x_jit = pmin(pmax(x_jit, 0), 100),
    y_jit = pmin(pmax(y_jit, 0), 100)
  )
highlight_p4 <- pep_small4 %>% build_highlight_df()

p4 <- ggplot(pep_small4, aes(x = x_jit, y = y_jit)) +
  geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
  geom_point(
    data = highlight_p4,
    inherit.aes = FALSE,
    aes(x = x_jit, y = y_jit, colour = cat_label),
    size = 3.6, alpha = 0.8, na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Kids: Birth (T2) vs T8",
    x = "Percent positive at Birth (T2)",
    y = "Percent positive at T8"
  ) +
  scale_color_manual(values = highlight_cols, name = "Highlight", na.value = "grey50", drop = FALSE) +
  theme_phip() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank()
  )
p4
ggsave(static_path4, p4, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

p_inter4 <- plot_ly() %>%
  add_markers(
    data = pep_small4,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    marker = list(size = 5, color = "grey80", opacity = 0.3),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  add_markers(
    data = highlight_p4,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    color = ~cat_label, colors = highlight_cols,
    marker = list(
      size = 12,
      symbol = "circle",
      opacity = 0.8,
      line = list(color = "rgba(0,0,0,0)", width = 0)
    ),
    text = ~paste0(
      "Feature: ", feature,
      "<br>Rank: ", rank,
      "<br>Birth (T2) %+: ", round(percent1, 2),
      "<br>T8 %+: ", round(percent2, 2),
      "<br>BH category: ", category_rank_bh
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
  layout(
    title = "Kids: Birth (T2) vs T8",
    xaxis = list(title = "Percent positive at Birth (T2)", range = c(0, 100)),
    yaxis = list(title = "Percent positive at T8", range = c(0, 100)),
    shapes = list(
      list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
           line = list(dash = "dash", color = "grey"))
    ),
    legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
  )

htmlwidgets::saveWidget(p_inter4, interactive_path4, selfcontained = FALSE)

cat("Saved static scatter to:", static_path4, "\n")
cat("Saved interactive scatter to:", interactive_path4, "\n")

# ------------------------------------------------------------------
# Moms: serum birth (T2) vs milk follow-up (T6)
# ------------------------------------------------------------------
comparison5 <- "mom_serum_T2_vs_mom_milk_T7"
single_pep_path5 <- file.path(base_dir, comparison5, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path5)) {
  stop("Missing POP single_peptide.csv for ", comparison5, " at ", single_pep_path5)
}

static_path5 <- file.path(out_dir, "mom_serum_T2_vs_mom_milk_T6.svg")
interactive_path5 <- file.path(out_dir, "mom_serum_T2_vs_mom_milk_T6.html")

pep5 <- readr::read_csv(single_pep_path5, show_col_types = FALSE)

highlight_points5 <- tibble::tribble(
  ~rank, ~feature,
  "class", "Mammalia",
  #"species", "Salmonella enterica",
  #"genus", "Parabacteroides",
  "species", "Streptococcus pneumoniae",
  "species", "Neisseria meningitidis",
  "species", "Staphylococcus aureus"
) %>%
  mutate(cat_label = feature)

highlight_labels5 <- unique(highlight_points5$cat_label)
highlight_cols5 <- setNames(
  pal_fun(max(3, length(highlight_labels5)))[seq_along(highlight_labels5)],
  highlight_labels5
)
highlight_peptides5 <- build_highlight_peptides(peplib, highlight_points5)
highlight_peptides_map5 <- highlight_peptides5 %>% rename(feature = peptide_id)

pep_small5 <- pep5 %>%
  select(feature, rank, percent1, percent2, category_rank_bh) %>%
  filter(!is.na(percent1) & !is.na(percent2)) %>%
  add_categories() %>%
  mutate(
    x_jit = jitter(percent1, amount = 0.5),
    y_jit = jitter(percent2, amount = 0.5),
    x_jit = pmin(pmax(x_jit, 0), 100),
    y_jit = pmin(pmax(y_jit, 0), 100)
  )
highlight_p5 <- pep_small5 %>%
  inner_join(highlight_peptides_map5, by = "feature") %>%
  mutate(cat_label = factor(cat_label, levels = names(highlight_cols5))) %>%
  tidyr::drop_na(cat_label, percent1, percent2, x_jit, y_jit)

p5 <- ggplot(pep_small5, aes(x = x_jit, y = y_jit)) +
  geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
  geom_point(
    data = highlight_p5,
    inherit.aes = FALSE,
    aes(x = x_jit, y = y_jit, colour = cat_label),
    size = 3.6, alpha = 0.8, na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Mothers: Serum Birth (T2) vs Milk T6",
    x = "Percent positive in Maternal Serum at Birth (T2)",
    y = "Percent positive in Maternal Milk at T6"
  ) +
  scale_color_manual(values = highlight_cols5, name = "Highlight", na.value = "grey50", drop = FALSE) +
  theme_phip() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank()
  )
p5
ggsave(static_path5, p5, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

p_inter5 <- plot_ly() %>%
  add_markers(
    data = pep_small5,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    marker = list(size = 5, color = "grey80", opacity = 0.3),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  add_markers(
    data = highlight_p5,
    x = ~x_jit, y = ~y_jit,
    type = "scatter", mode = "markers",
    color = ~cat_label, colors = highlight_cols5,
    marker = list(
      size = 12,
      symbol = "circle",
      opacity = 0.8,
      line = list(color = "rgba(0,0,0,0)", width = 0)
    ),
    text = ~paste0(
      "Feature: ", feature,
      "<br>Rank: ", rank,
      "<br>Maternal Serum Birth (T2) %+: ", round(percent1, 2),
      "<br>Maternal Milk T6 %+: ", round(percent2, 2),
      "<br>BH category: ", category_rank_bh
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
  layout(
    title = "Mothers: Serum Birth (T2) vs Milk T6",
    xaxis = list(title = "Percent positive in Maternal Serum at Birth (T2)", range = c(0, 100)),
    yaxis = list(title = "Percent positive in Maternal Milk at T6", range = c(0, 100)),
    shapes = list(
      list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
           line = list(dash = "dash", color = "grey"))
    ),
    legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
  )

htmlwidgets::saveWidget(p_inter5, interactive_path5, selfcontained = FALSE)

cat("Saved static scatter to:", static_path5, "\n")
cat("Saved interactive scatter to:", interactive_path5, "\n")

# ==================================================================
# ==================== BREAST MILK FEATURE SCATTERS =================
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("BREAST MILK FEATURE SCATTERS\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

base_dir_milk <- file.path("results", "results_milk")

milk_targets <- list(
  list(
    comparison = "mom_milk_T4_smoking_yes_vs_mom_milk_T4_smoking_no",
    highlights = tibble::tribble(
      ~rank, ~feature,
      "species", "Mesomycoplasma hyopneumoniae",
      "species", "Periplaneta americana",
      "family", "Blattidae"
    )
  ),
  list(
    comparison = "mom_milk_T6_sex_male_vs_mom_milk_T6_sex_female",
    highlights = tibble::tribble(
      ~rank, ~feature,
      "genus", "Cladosporium"
    )
  )
)

for (target in milk_targets) {
  comparison_milk <- target$comparison
  single_pep_path_milk <- file.path(
    base_dir_milk, comparison_milk, "POP_framework", "single_peptide.csv"
  )

  if (!file.exists(single_pep_path_milk)) {
    warning("Missing POP single_peptide.csv for ", comparison_milk, " at ", single_pep_path_milk)
    next
  }

  comparison_safe <- gsub("[^A-Za-z0-9_-]+", "_", comparison_milk)
  static_path_milk <- file.path(out_dir, paste0(comparison_safe, "_selected_features.svg"))
  interactive_path_milk <- file.path(out_dir, paste0(comparison_safe, "_selected_features.html"))

  groups <- strsplit(comparison_milk, "_vs_", fixed = TRUE)[[1]]
  group1_lbl <- gsub("_", " ", groups[1])
  group2_lbl <- gsub("_", " ", groups[2])

  hp_milk <- target$highlights %>%
    mutate(cat_label = feature)
  hlabels_milk <- unique(hp_milk$cat_label)
  hcols_milk <- setNames(
    pal_fun(max(3, length(hlabels_milk)))[seq_along(hlabels_milk)],
    hlabels_milk
  )
  hpep_map_milk <- build_highlight_peptides(peplib, hp_milk) %>%
    rename(feature = peptide_id)

  pep_milk <- readr::read_csv(single_pep_path_milk, show_col_types = FALSE)
  pep_small_milk <- pep_milk %>%
    select(feature, rank, percent1, percent2, category_rank_bh) %>%
    filter(!is.na(percent1) & !is.na(percent2)) %>%
    mutate(
      x_jit = jitter(percent1, amount = 0.5),
      y_jit = jitter(percent2, amount = 0.5),
      x_jit = pmin(pmax(x_jit, 0), 100),
      y_jit = pmin(pmax(y_jit, 0), 100)
    )

  highlight_milk <- pep_small_milk %>%
    inner_join(hpep_map_milk, by = "feature") %>%
    mutate(cat_label = factor(cat_label, levels = names(hcols_milk))) %>%
    tidyr::drop_na(cat_label, percent1, percent2, x_jit, y_jit)

  p_milk <- ggplot(pep_small_milk, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_milk,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = paste0("Breast milk: ", group1_lbl, " vs ", group2_lbl),
      x = paste0("Percent positive in ", group1_lbl),
      y = paste0("Percent positive in ", group2_lbl)
    ) +
    scale_color_manual(values = hcols_milk, name = "Selected features", na.value = "grey50", drop = FALSE) +
    theme_phip() +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )
  ggsave(static_path_milk, p_milk, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

  p_inter_milk <- plot_ly() %>%
    add_markers(
      data = pep_small_milk,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data = highlight_milk,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      color = ~cat_label, colors = hcols_milk,
      marker = list(
        size = 12,
        symbol = "circle",
        opacity = 0.85,
        line = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>", group1_lbl, " %+: ", round(percent1, 2),
        "<br>", group2_lbl, " %+: ", round(percent2, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title = paste0("Breast milk: ", group1_lbl, " vs ", group2_lbl),
      xaxis = list(title = paste0("Percent positive in ", group1_lbl), range = c(0, 100)),
      yaxis = list(title = paste0("Percent positive in ", group2_lbl), range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_milk, interactive_path_milk, selfcontained = FALSE)
  cat("Saved breast-milk static scatter to:", static_path_milk, "\n")
  cat("Saved breast-milk interactive scatter to:", interactive_path_milk, "\n")
}

# ==================================================================
# ================== KID DELMODE VG VS CS SCATTERS =================
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("KID DELMODE VG VS CS SCATTERS\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

normalize_feature_key <- function(x) {
  tolower(gsub("[^a-z0-9]", "", as.character(x)))
}

build_highlight_peptides_from_labels <- function(peplib_df, labels) {
  target_tbl <- tibble::tibble(
    cat_label = labels,
    feature_key = normalize_feature_key(labels)
  )

  tax_cols <- intersect(
    c("species", "genus", "family", "order", "class", "phylum", "anno_viruses_bacteriophage"),
    names(peplib_df)
  )

  peplib_df %>%
    dplyr::select(peptide_id, dplyr::all_of(tax_cols)) %>%
    tidyr::pivot_longer(
      cols = -peptide_id,
      names_to = "rank_source",
      values_to = "feature_value"
    ) %>%
    dplyr::filter(!is.na(.data$feature_value)) %>%
    dplyr::mutate(feature_key = normalize_feature_key(.data$feature_value)) %>%
    dplyr::inner_join(target_tbl, by = "feature_key") %>%
    dplyr::select(peptide_id, cat_label) %>%
    dplyr::distinct()
}

# Shared metadata color map across breastfeeding/delmode/siblings/smoking plots.
# Uses the user-defined manual palette first, then extras, and assigns colors in an
# interleaved order so each single plot gets well-separated colors.
metadata_group_labels <- list(
  delmode = c(
    "Rhinovirus B",
    "Poales",
    "Streptococcus pneumoniae",
    "Varicellovirus humanalpha3",
    "Simplexvirus humanalpha2"
  ),
  feeding = c(
    "Milk proteins",
    "Orthoherpesviridae",
    "Simplexvirus",
    "Mycoplasmoides pneumoniae",
    "Pseudomonas aeruginosa"
  ),
  siblings = c(
    "Streptococcus pneumoniae",
    "Rhinovirus A",
    "Clostridium perfringens",
    "Staphylococcus aureus",
    "Pisum sativum"
  ),
  smoking = c(
    "Orthopneumovirus hominis",
    "Cytomegalovirus humanbeta5"
  )
)

metadata_feature_levels <- unique(unlist(metadata_group_labels, use.names = FALSE))
metadata_extra_palette <- c(
  "#00A6D6", # cyan-blue
  "#D55E00", # vermillion
  "#66A61E", # olive green
  "#332288", # indigo
  "#F0E442"  # yellow
)
metadata_palette <- unique(c(palette_manual, metadata_extra_palette))

if (length(metadata_palette) < length(metadata_feature_levels)) {
  stop("Not enough colors for metadata categories.")
}

metadata_cols <- c(
  "Rhinovirus B" = "#1B9E77",
  "Poales" = "#E7298A",
  "Streptococcus pneumoniae" = "#A6761D",
  "Varicellovirus humanalpha3" = "#00A6D6",
  "Simplexvirus humanalpha2" = "#B2DF8A",
  "Milk proteins" = "#B2182B",
  "Orthoherpesviridae" = "#1F78B4",
  "Simplexvirus" = "#E6AB02",
  "Mycoplasmoides pneumoniae" = "#6A3D9A",
  "Pseudomonas aeruginosa" = "#66A61E",
  "Rhinovirus A" = "#7570B3",
  "Clostridium perfringens" = "#FF7F00",
  "Staphylococcus aureus" = "#FB9A99",
  "Pisum sativum" = "#F0E442",
  "Orthopneumovirus hominis" = "#332288",
  "Cytomegalovirus humanbeta5" = "#D55E00"
)

if (!all(metadata_feature_levels %in% names(metadata_cols))) {
  missing_metadata_cols <- setdiff(metadata_feature_levels, names(metadata_cols))
  stop("Missing metadata color assignments for: ", paste(missing_metadata_cols, collapse = ", "))
}
metadata_cols <- metadata_cols[metadata_feature_levels]

delmode_feature_groups <- metadata_group_labels$delmode

delmode_cols <- metadata_cols[delmode_feature_groups]

highlight_peptides_delmode <- build_highlight_peptides_from_labels(peplib, delmode_feature_groups)
highlight_peptides_delmode_map <- highlight_peptides_delmode %>% dplyr::rename(feature = peptide_id)

delmode_comparisons <- list(
  list(
    comparison_vg_first = "kid_serum_T6_delmode_VG_vs_kid_serum_T6_delmode_CS",
    title = "Kids T6: Delmode VG vs Delmode CS",
    xlab = "Percent positive in kid serum T6 delmode VG",
    ylab = "Percent positive in kid serum T6 delmode CS"
  ),
  list(
    comparison_vg_first = "kid_serum_T8_delmode_VG_vs_kid_serum_T8_delmode_CS",
    title = "Kids T8: Delmode VG vs Delmode CS",
    xlab = "Percent positive in kid serum T8 delmode VG",
    ylab = "Percent positive in kid serum T8 delmode CS"
  )
)

for (cfg in delmode_comparisons) {
  cmp_vg_first <- cfg$comparison_vg_first
  groups <- strsplit(cmp_vg_first, "_vs_", fixed = TRUE)[[1]]
  cmp_cs_first <- paste0(groups[2], "_vs_", groups[1])

  path_vg_first <- file.path(base_dir, cmp_vg_first, "POP_framework", "single_peptide.csv")
  path_cs_first <- file.path(base_dir, cmp_cs_first, "POP_framework", "single_peptide.csv")

  swap_axes <- FALSE
  single_pep_path_delmode <- path_vg_first

  if (!file.exists(path_vg_first)) {
    if (file.exists(path_cs_first)) {
      swap_axes <- TRUE
      single_pep_path_delmode <- path_cs_first
      message("Using ", cmp_cs_first, " and swapping axes to keep VG on x-axis.")
    } else {
      warning("Missing POP single_peptide.csv for both: ", cmp_vg_first, " and ", cmp_cs_first)
      next
    }
  }

  cmp_safe <- gsub("[^A-Za-z0-9_-]+", "_", cmp_vg_first)
  static_path_delmode <- file.path(out_dir, paste0(cmp_safe, "_selected_groups.svg"))
  interactive_path_delmode <- file.path(out_dir, paste0(cmp_safe, "_selected_groups.html"))

  pep_delmode <- readr::read_csv(
    single_pep_path_delmode,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )

  group1_label <- NA_character_
  group2_label <- NA_character_
  if (all(c("group1", "group2") %in% names(pep_delmode))) {
    g1_vals <- unique(stats::na.omit(as.character(pep_delmode$group1)))
    g2_vals <- unique(stats::na.omit(as.character(pep_delmode$group2)))
    if (length(g1_vals) >= 1L) group1_label <- g1_vals[[1]]
    if (length(g2_vals) >= 1L) group2_label <- g2_vals[[1]]

    g1_is_vg <- isTRUE(grepl("delmode_VG", group1_label, fixed = TRUE))
    g2_is_vg <- isTRUE(grepl("delmode_VG", group2_label, fixed = TRUE))
    if (xor(g1_is_vg, g2_is_vg)) {
      swap_axes <- g2_is_vg
    } else {
      warning(
        "Could not infer unique VG side from group1/group2 for ", cmp_vg_first,
        "; using file-order fallback."
      )
    }
  } else {
    warning(
      "Columns group1/group2 missing in ", basename(single_pep_path_delmode),
      "; using file-order fallback."
    )
  }

  x_group_label <- if (swap_axes) group2_label else group1_label
  y_group_label <- if (swap_axes) group1_label else group2_label
  xlab_use <- if (!is.na(x_group_label) && nzchar(x_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", x_group_label))
  } else {
    cfg$xlab
  }
  ylab_use <- if (!is.na(y_group_label) && nzchar(y_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", y_group_label))
  } else {
    cfg$ylab
  }
  x_group_hover <- if (!is.na(x_group_label) && nzchar(x_group_label)) x_group_label else "group_x"
  y_group_hover <- if (!is.na(y_group_label) && nzchar(y_group_label)) y_group_label else "group_y"

  pep_small_delmode <- pep_delmode %>%
    dplyr::select(feature, rank, percent1, percent2, category_rank_bh) %>%
    dplyr::filter(!is.na(percent1) & !is.na(percent2)) %>%
    dplyr::mutate(
      x_pct = if (swap_axes) .data$percent2 else .data$percent1,
      y_pct = if (swap_axes) .data$percent1 else .data$percent2,
      x_jit = jitter(.data$x_pct, amount = 0.5),
      y_jit = jitter(.data$y_pct, amount = 0.5),
      x_jit = pmin(pmax(.data$x_jit, 0), 100),
      y_jit = pmin(pmax(.data$y_jit, 0), 100)
    )

  highlight_delmode <- pep_small_delmode %>%
    dplyr::inner_join(highlight_peptides_delmode_map, by = "feature") %>%
    dplyr::mutate(cat_label = factor(cat_label, levels = names(delmode_cols))) %>%
    dplyr::arrange(.data$cat_label) %>%
    tidyr::drop_na(cat_label, x_jit, y_jit)

  p_delmode <- ggplot(pep_small_delmode, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_delmode,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = cfg$title,
      x = xlab_use,
      y = ylab_use
    ) +
    scale_color_manual(values = delmode_cols, name = "Selected features", na.value = "grey50", drop = TRUE) +
    theme_phip() +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )

  ggsave(static_path_delmode, p_delmode, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

  p_inter_delmode <- plot_ly() %>%
    add_markers(
      data = pep_small_delmode,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data = highlight_delmode,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      color = ~cat_label, colors = delmode_cols,
      marker = list(
        size = 12,
        symbol = "circle",
        opacity = 0.85,
        line = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>", x_group_hover, " %+: ", round(x_pct, 2),
        "<br>", y_group_hover, " %+: ", round(y_pct, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title = cfg$title,
      xaxis = list(title = xlab_use, range = c(0, 100)),
      yaxis = list(title = ylab_use, range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_delmode, interactive_path_delmode, selfcontained = FALSE)
  cat("Saved delmode static scatter to:", static_path_delmode, "\n")
  cat("Saved delmode interactive scatter to:", interactive_path_delmode, "\n")
}

# ==================================================================
# ================== KID FEEDING BF VS MF SCATTERS =================
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("KID FEEDING BF VS MF SCATTERS\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

build_highlight_peptides_from_map <- function(peplib_df, feature_map_tbl) {
  target_tbl <- feature_map_tbl %>%
    dplyr::mutate(feature_key = normalize_feature_key(.data$feature_label))

  tax_cols <- intersect(
    c("species", "genus", "family", "order", "class", "phylum", "anno_viruses_bacteriophage"),
    names(peplib_df)
  )

  peplib_df %>%
    dplyr::select(peptide_id, dplyr::all_of(tax_cols)) %>%
    tidyr::pivot_longer(
      cols = -peptide_id,
      names_to = "rank_source",
      values_to = "feature_value"
    ) %>%
    dplyr::filter(!is.na(.data$feature_value)) %>%
    dplyr::mutate(feature_key = normalize_feature_key(.data$feature_value)) %>%
    dplyr::inner_join(target_tbl, by = "feature_key") %>%
    dplyr::select(peptide_id, cat_label) %>%
    dplyr::distinct()
}

feeding_feature_map <- tibble::tribble(
  ~feature_label, ~cat_label,
  "Bos taurus", "Milk proteins",
  "Bubalus bubalis", "Milk proteins",
  "Capra hircus", "Milk proteins",
  "Ovis aries", "Milk proteins",
  "Rangifer tarandus", "Milk proteins",
  "Sus scrofa", "Milk proteins",
  "Orthoherpesviridae", "Orthoherpesviridae",
  "Simplexvirus", "Simplexvirus",
  "Mycoplasmoides pneumoniae", "Mycoplasmoides pneumoniae",
  "Pseudomonas aeruginosa", "Pseudomonas aeruginosa"
)

feeding_cat_levels <- metadata_group_labels$feeding

feeding_cols <- metadata_cols[feeding_cat_levels]

highlight_peptides_feeding <- build_highlight_peptides_from_map(peplib, feeding_feature_map)
highlight_peptides_feeding_map <- highlight_peptides_feeding %>% dplyr::rename(feature = peptide_id)

feeding_comparisons <- list(
  list(
    comparison_bf_first = "kid_serum_T6_feeding_BF_vs_kid_serum_T6_feeding_MF",
    title = "Kids T6: Feeding BF vs Feeding MF",
    xlab = "Percent positive in kid serum T6 feeding BF",
    ylab = "Percent positive in kid serum T6 feeding MF"
  ),
  list(
    comparison_bf_first = "kid_serum_T8_feeding_BF_vs_kid_serum_T8_feeding_MF",
    title = "Kids T8: Feeding BF vs Feeding MF",
    xlab = "Percent positive in kid serum T8 feeding BF",
    ylab = "Percent positive in kid serum T8 feeding MF"
  )
)

for (cfg in feeding_comparisons) {
  cmp_bf_first <- cfg$comparison_bf_first
  groups <- strsplit(cmp_bf_first, "_vs_", fixed = TRUE)[[1]]
  cmp_mf_first <- paste0(groups[2], "_vs_", groups[1])

  path_bf_first <- file.path(base_dir, cmp_bf_first, "POP_framework", "single_peptide.csv")
  path_mf_first <- file.path(base_dir, cmp_mf_first, "POP_framework", "single_peptide.csv")

  swap_axes <- FALSE
  single_pep_path_feeding <- path_bf_first

  if (!file.exists(path_bf_first)) {
    if (file.exists(path_mf_first)) {
      swap_axes <- TRUE
      single_pep_path_feeding <- path_mf_first
      message("Using ", cmp_mf_first, " and swapping axes to keep BF on x-axis.")
    } else {
      warning("Missing POP single_peptide.csv for both: ", cmp_bf_first, " and ", cmp_mf_first)
      next
    }
  }

  cmp_safe <- gsub("[^A-Za-z0-9_-]+", "_", cmp_bf_first)
  static_path_feeding <- file.path(out_dir, paste0(cmp_safe, "_selected_groups.svg"))
  interactive_path_feeding <- file.path(out_dir, paste0(cmp_safe, "_selected_groups.html"))

  pep_feeding <- readr::read_csv(
    single_pep_path_feeding,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )

  group1_label <- NA_character_
  group2_label <- NA_character_
  if (all(c("group1", "group2") %in% names(pep_feeding))) {
    g1_vals <- unique(stats::na.omit(as.character(pep_feeding$group1)))
    g2_vals <- unique(stats::na.omit(as.character(pep_feeding$group2)))
    if (length(g1_vals) >= 1L) group1_label <- g1_vals[[1]]
    if (length(g2_vals) >= 1L) group2_label <- g2_vals[[1]]

    g1_is_bf <- isTRUE(grepl("feeding_BF", group1_label, fixed = TRUE))
    g2_is_bf <- isTRUE(grepl("feeding_BF", group2_label, fixed = TRUE))
    if (xor(g1_is_bf, g2_is_bf)) {
      swap_axes <- g2_is_bf
    } else {
      warning(
        "Could not infer unique BF side from group1/group2 for ", cmp_bf_first,
        "; using file-order fallback."
      )
    }
  } else {
    warning(
      "Columns group1/group2 missing in ", basename(single_pep_path_feeding),
      "; using file-order fallback."
    )
  }

  x_group_label <- if (swap_axes) group2_label else group1_label
  y_group_label <- if (swap_axes) group1_label else group2_label
  xlab_use <- if (!is.na(x_group_label) && nzchar(x_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", x_group_label))
  } else {
    cfg$xlab
  }
  ylab_use <- if (!is.na(y_group_label) && nzchar(y_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", y_group_label))
  } else {
    cfg$ylab
  }
  x_group_hover <- if (!is.na(x_group_label) && nzchar(x_group_label)) x_group_label else "group_x"
  y_group_hover <- if (!is.na(y_group_label) && nzchar(y_group_label)) y_group_label else "group_y"

  pep_small_feeding <- pep_feeding %>%
    dplyr::select(feature, rank, percent1, percent2, category_rank_bh) %>%
    dplyr::filter(!is.na(percent1) & !is.na(percent2)) %>%
    dplyr::mutate(
      x_pct = if (swap_axes) .data$percent2 else .data$percent1,
      y_pct = if (swap_axes) .data$percent1 else .data$percent2,
      x_jit = jitter(.data$x_pct, amount = 0.5),
      y_jit = jitter(.data$y_pct, amount = 0.5),
      x_jit = pmin(pmax(.data$x_jit, 0), 100),
      y_jit = pmin(pmax(.data$y_jit, 0), 100)
    )

  highlight_feeding <- pep_small_feeding %>%
    dplyr::inner_join(highlight_peptides_feeding_map, by = "feature") %>%
    dplyr::mutate(cat_label = factor(cat_label, levels = feeding_cat_levels)) %>%
    dplyr::arrange(.data$cat_label) %>%
    tidyr::drop_na(cat_label, x_jit, y_jit)

  p_feeding <- ggplot(pep_small_feeding, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_feeding,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = cfg$title,
      x = xlab_use,
      y = ylab_use
    ) +
    scale_color_manual(values = feeding_cols, name = "Selected features", na.value = "grey50", drop = TRUE) +
    theme_phip() +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )

  ggsave(static_path_feeding, p_feeding, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

  p_inter_feeding <- plot_ly() %>%
    add_markers(
      data = pep_small_feeding,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data = highlight_feeding,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      color = ~cat_label, colors = feeding_cols,
      marker = list(
        size = 12,
        symbol = "circle",
        opacity = 0.85,
        line = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>", x_group_hover, " %+: ", round(x_pct, 2),
        "<br>", y_group_hover, " %+: ", round(y_pct, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title = cfg$title,
      xaxis = list(title = xlab_use, range = c(0, 100)),
      yaxis = list(title = ylab_use, range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_feeding, interactive_path_feeding, selfcontained = FALSE)
  cat("Saved feeding static scatter to:", static_path_feeding, "\n")
  cat("Saved feeding interactive scatter to:", interactive_path_feeding, "\n")
}

# ==================================================================
# ================= KID T8 SIBLINGS VS NO-SIBLINGS ================
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("KID T8 SIBLINGS VS NO-SIBLINGS\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

siblings_feature_groups <- metadata_group_labels$siblings

siblings_cols <- metadata_cols[siblings_feature_groups]

highlight_peptides_siblings <- build_highlight_peptides_from_labels(peplib, siblings_feature_groups)
highlight_peptides_siblings_map <- highlight_peptides_siblings %>% dplyr::rename(feature = peptide_id)

comparison_siblings <- "kid_serum_T8_siblings_vs_kid_serum_T8_no_siblings"
single_pep_path_siblings <- file.path(base_dir, comparison_siblings, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path_siblings)) {
  warning("Missing POP single_peptide.csv for ", comparison_siblings, " at ", single_pep_path_siblings)
} else {
  static_path_siblings <- file.path(out_dir, "kid_serum_T8_siblings_vs_kid_serum_T8_no_siblings_selected_groups.svg")
  interactive_path_siblings <- file.path(out_dir, "kid_serum_T8_siblings_vs_kid_serum_T8_no_siblings_selected_groups.html")

  pep_siblings <- readr::read_csv(
    single_pep_path_siblings,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )

  group1_label <- NA_character_
  group2_label <- NA_character_
  swap_axes <- FALSE
  if (all(c("group1", "group2") %in% names(pep_siblings))) {
    g1_vals <- unique(stats::na.omit(as.character(pep_siblings$group1)))
    g2_vals <- unique(stats::na.omit(as.character(pep_siblings$group2)))
    if (length(g1_vals) >= 1L) group1_label <- g1_vals[[1]]
    if (length(g2_vals) >= 1L) group2_label <- g2_vals[[1]]

    g1_is_no_siblings <- grepl("no_siblings", group1_label, fixed = TRUE)
    g2_is_no_siblings <- grepl("no_siblings", group2_label, fixed = TRUE)
    if (xor(g1_is_no_siblings, g2_is_no_siblings)) {
      swap_axes <- g2_is_no_siblings
    } else {
      warning("Could not infer unique no_siblings side from group1/group2; using file-order fallback.")
    }
  } else {
    warning("Columns group1/group2 missing in ", basename(single_pep_path_siblings), "; using file-order fallback.")
  }

  x_group_label <- if (swap_axes) group2_label else group1_label
  y_group_label <- if (swap_axes) group1_label else group2_label
  xlab_use <- if (!is.na(x_group_label) && nzchar(x_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", x_group_label))
  } else {
    "Percent positive in kid serum T8 no siblings"
  }
  ylab_use <- if (!is.na(y_group_label) && nzchar(y_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", y_group_label))
  } else {
    "Percent positive in kid serum T8 siblings"
  }
  x_group_hover <- if (!is.na(x_group_label) && nzchar(x_group_label)) x_group_label else "group_x"
  y_group_hover <- if (!is.na(y_group_label) && nzchar(y_group_label)) y_group_label else "group_y"

  pep_small_siblings <- pep_siblings %>%
    dplyr::select(feature, rank, percent1, percent2, category_rank_bh) %>%
    dplyr::filter(!is.na(percent1) & !is.na(percent2)) %>%
    dplyr::mutate(
      x_pct = if (swap_axes) .data$percent2 else .data$percent1,
      y_pct = if (swap_axes) .data$percent1 else .data$percent2,
      x_jit = jitter(.data$x_pct, amount = 0.5),
      y_jit = jitter(.data$y_pct, amount = 0.5),
      x_jit = pmin(pmax(.data$x_jit, 0), 100),
      y_jit = pmin(pmax(.data$y_jit, 0), 100)
    )

  highlight_siblings <- pep_small_siblings %>%
    dplyr::inner_join(highlight_peptides_siblings_map, by = "feature") %>%
    dplyr::mutate(cat_label = factor(cat_label, levels = siblings_feature_groups)) %>%
    dplyr::arrange(.data$cat_label) %>%
    tidyr::drop_na(cat_label, x_jit, y_jit)

  p_siblings <- ggplot(pep_small_siblings, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_siblings,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = "Kids T8: Siblings vs No siblings",
      x = xlab_use,
      y = ylab_use
    ) +
    scale_color_manual(values = siblings_cols, name = "Selected features", na.value = "grey50", drop = TRUE) +
    theme_phip() +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )

  ggsave(static_path_siblings, p_siblings, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

  p_inter_siblings <- plot_ly() %>%
    add_markers(
      data = pep_small_siblings,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data = highlight_siblings,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      color = ~cat_label, colors = siblings_cols,
      marker = list(
        size = 12,
        symbol = "circle",
        opacity = 0.85,
        line = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>", x_group_hover, " %+: ", round(x_pct, 2),
        "<br>", y_group_hover, " %+: ", round(y_pct, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title = "Kids T8: Siblings vs No siblings",
      xaxis = list(title = xlab_use, range = c(0, 100)),
      yaxis = list(title = ylab_use, range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_siblings, interactive_path_siblings, selfcontained = FALSE)
  cat("Saved siblings static scatter to:", static_path_siblings, "\n")
  cat("Saved siblings interactive scatter to:", interactive_path_siblings, "\n")
}

# ==================================================================
# ================= KID T2 SMOKING YES VS NO =======================
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("KID T2 SMOKING YES VS NO\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

smoking_feature_map <- tibble::tribble(
  ~feature_label, ~cat_label,
  "Orthopneumovirus hominis", "Orthopneumovirus hominis",
  "Cytomeglovirus humanbeta5", "Cytomegalovirus humanbeta5",
  "Cytomegalovirus humanbeta5", "Cytomegalovirus humanbeta5"
)

smoking_cat_levels <- metadata_group_labels$smoking

smoking_cols <- metadata_cols[smoking_cat_levels]

highlight_peptides_smoking <- build_highlight_peptides_from_map(peplib, smoking_feature_map)
highlight_peptides_smoking_map <- highlight_peptides_smoking %>% dplyr::rename(feature = peptide_id)

comparison_smoking <- "kid_serum_T2_smoking_yes_vs_kid_serum_T2_smoking_no"
single_pep_path_smoking <- file.path(base_dir, comparison_smoking, "POP_framework", "single_peptide.csv")

if (!file.exists(single_pep_path_smoking)) {
  warning("Missing POP single_peptide.csv for ", comparison_smoking, " at ", single_pep_path_smoking)
} else {
  static_path_smoking <- file.path(out_dir, "kid_serum_T2_smoking_yes_vs_kid_serum_T2_smoking_no_selected_groups.svg")
  interactive_path_smoking <- file.path(out_dir, "kid_serum_T2_smoking_yes_vs_kid_serum_T2_smoking_no_selected_groups.html")

  pep_smoking <- readr::read_csv(
    single_pep_path_smoking,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )

  group1_label <- NA_character_
  group2_label <- NA_character_
  swap_axes <- FALSE
  if (all(c("group1", "group2") %in% names(pep_smoking))) {
    g1_vals <- unique(stats::na.omit(as.character(pep_smoking$group1)))
    g2_vals <- unique(stats::na.omit(as.character(pep_smoking$group2)))
    if (length(g1_vals) >= 1L) group1_label <- g1_vals[[1]]
    if (length(g2_vals) >= 1L) group2_label <- g2_vals[[1]]

    g1_is_no <- grepl("smoking_no", group1_label, fixed = TRUE)
    g2_is_no <- grepl("smoking_no", group2_label, fixed = TRUE)
    if (xor(g1_is_no, g2_is_no)) {
      swap_axes <- g2_is_no
    } else {
      warning("Could not infer unique smoking_no side from group1/group2; using file-order fallback.")
    }
  } else {
    warning("Columns group1/group2 missing in ", basename(single_pep_path_smoking), "; using file-order fallback.")
  }

  x_group_label <- if (swap_axes) group2_label else group1_label
  y_group_label <- if (swap_axes) group1_label else group2_label
  xlab_use <- if (!is.na(x_group_label) && nzchar(x_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", x_group_label))
  } else {
    "Percent positive in kid serum T2 smoking no"
  }
  ylab_use <- if (!is.na(y_group_label) && nzchar(y_group_label)) {
    paste0("Percent positive in ", gsub("_", " ", y_group_label))
  } else {
    "Percent positive in kid serum T2 smoking yes"
  }
  x_group_hover <- if (!is.na(x_group_label) && nzchar(x_group_label)) x_group_label else "group_x"
  y_group_hover <- if (!is.na(y_group_label) && nzchar(y_group_label)) y_group_label else "group_y"

  pep_small_smoking <- pep_smoking %>%
    dplyr::select(feature, rank, percent1, percent2, category_rank_bh) %>%
    dplyr::filter(!is.na(percent1) & !is.na(percent2)) %>%
    dplyr::mutate(
      x_pct = if (swap_axes) .data$percent2 else .data$percent1,
      y_pct = if (swap_axes) .data$percent1 else .data$percent2,
      x_jit = jitter(.data$x_pct, amount = 0.5),
      y_jit = jitter(.data$y_pct, amount = 0.5),
      x_jit = pmin(pmax(.data$x_jit, 0), 100),
      y_jit = pmin(pmax(.data$y_jit, 0), 100)
    )

  highlight_smoking <- pep_small_smoking %>%
    dplyr::inner_join(highlight_peptides_smoking_map, by = "feature") %>%
    dplyr::mutate(cat_label = factor(cat_label, levels = smoking_cat_levels)) %>%
    dplyr::arrange(.data$cat_label) %>%
    tidyr::drop_na(cat_label, x_jit, y_jit)

  p_smoking <- ggplot(pep_small_smoking, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_smoking,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = "Kids T2: Smoking yes vs no",
      x = xlab_use,
      y = ylab_use
    ) +
    scale_color_manual(values = smoking_cols, name = "Selected features", na.value = "grey50", drop = TRUE) +
    theme_phip() +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )

  ggsave(static_path_smoking, p_smoking, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")

  p_inter_smoking <- plot_ly() %>%
    add_markers(
      data = pep_small_smoking,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data = highlight_smoking,
      x = ~x_jit, y = ~y_jit,
      type = "scatter", mode = "markers",
      color = ~cat_label, colors = smoking_cols,
      marker = list(
        size = 12,
        symbol = "circle",
        opacity = 0.85,
        line = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>", x_group_hover, " %+: ", round(x_pct, 2),
        "<br>", y_group_hover, " %+: ", round(y_pct, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title = "Kids T2: Smoking yes vs no",
      xaxis = list(title = xlab_use, range = c(0, 100)),
      yaxis = list(title = ylab_use, range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_smoking, interactive_path_smoking, selfcontained = FALSE)
  cat("Saved smoking static scatter to:", static_path_smoking, "\n")
  cat("Saved smoking interactive scatter to:", interactive_path_smoking, "\n")
}

# ==================================================================
# ========= BREAST MILK COMPARISONS: KID/MOM CROSS-SERIES =========
# ==================================================================
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("BREAST MILK COMPARISONS: KID/MOM CROSS-SERIES\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n", sep = "")

breastmilk_cmp_species <- c(
  "Staphylococcus aureus",
  "Barnesiella intestinihominis",
  "Lactobacillales",
  "Lymphocryptovirus humangamma4",
  "Escherichia coli",
  "Alphainfluenzavirus influenzae"  # last = drawn on top
)

breastmilk_cmp_cols <- setNames(
  pal_fun(max(3, length(breastmilk_cmp_species)))[seq_along(breastmilk_cmp_species)],
  breastmilk_cmp_species
)
breastmilk_cmp_cols["Staphylococcus aureus"]          <- palette_manual[2]  # dark red
breastmilk_cmp_cols["Alphainfluenzavirus influenzae"] <- palette_manual[1]  # teal

highlight_peptides_bmc <- build_highlight_peptides_from_labels(peplib, breastmilk_cmp_species)
highlight_peptides_bmc_map <- highlight_peptides_bmc %>% dplyr::rename(feature = peptide_id)

# swap = TRUE means percent1/percent2 are swapped so mom_milk lands on x-axis
breastmilk_cmp_configs <- list(
  list(
    comparison = "kid_serum_T6_vs_mom_milk_T6",
    title      = "Mom Milk T6 vs Kid Serum T6",
    xlab       = "Percent positive in mom milk T6",
    ylab       = "Percent positive in kid serum T6",
    file_stem  = "kid_serum_T6_vs_mom_milk_T6_breastmilk_species",
    swap       = TRUE
  ),
  list(
    comparison = "mom_serum_T2_vs_mom_milk_T4",
    title      = "Mom Milk T4 vs Mom Serum T2",
    xlab       = "Percent positive in mom milk T4",
    ylab       = "Percent positive in mom serum T2",
    file_stem  = "mom_serum_T2_vs_mom_milk_T4_breastmilk_species",
    swap       = FALSE
  ),
  list(
    comparison = "mom_milk_T4_vs_mom_milk_T7",
    title      = "Mom Milk T4 vs Mom Milk T7",
    xlab       = "Percent positive in mom milk T4",
    ylab       = "Percent positive in mom milk T7",
    file_stem  = "mom_milk_T4_vs_mom_milk_T7_breastmilk_species",
    swap       = FALSE
  )
)

for (cfg in breastmilk_cmp_configs) {
  cmp <- cfg$comparison
  single_pep_path_bmc <- file.path(base_dir, cmp, "POP_framework", "single_peptide.csv")

  if (!file.exists(single_pep_path_bmc)) {
    warning("Missing POP single_peptide.csv for ", cmp, " at ", single_pep_path_bmc)
    next
  }

  static_path_bmc      <- file.path(out_dir, paste0(cfg$file_stem, ".svg"))
  interactive_path_bmc <- file.path(out_dir, paste0(cfg$file_stem, ".html"))

  pep_bmc <- readr::read_csv(single_pep_path_bmc, show_col_types = FALSE)

  pep_small_bmc <- pep_bmc %>%
    dplyr::select(feature, rank, percent1, percent2, category_rank_bh) %>%
    dplyr::filter(!is.na(percent1) & !is.na(percent2)) %>%
    dplyr::mutate(
      x_pct = if (cfg$swap) .data$percent2 else .data$percent1,
      y_pct = if (cfg$swap) .data$percent1 else .data$percent2,
      x_jit = jitter(x_pct, amount = 0.5),
      y_jit = jitter(y_pct, amount = 0.5),
      x_jit = pmin(pmax(x_jit, 0), 100),
      y_jit = pmin(pmax(y_jit, 0), 100)
    )

  highlight_bmc <- pep_small_bmc %>%
    dplyr::inner_join(highlight_peptides_bmc_map, by = "feature") %>%
    dplyr::mutate(cat_label = factor(cat_label, levels = breastmilk_cmp_species)) %>%
    dplyr::arrange(cat_label) %>%
    tidyr::drop_na(cat_label, x_jit, y_jit)

  p_bmc <- ggplot(pep_small_bmc, aes(x = x_jit, y = y_jit)) +
    geom_point(color = "grey70", alpha = 0.75, size = 1.5) +
    geom_point(
      data = highlight_bmc,
      inherit.aes = FALSE,
      aes(x = x_jit, y = y_jit, colour = cat_label),
      size = 3.8, alpha = 0.85, na.rm = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey20") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      title = cfg$title,
      x     = cfg$xlab,
      y     = cfg$ylab
    ) +
    scale_color_manual(
      values = breastmilk_cmp_cols, name = "Species",
      na.value = "grey50", drop = FALSE
    ) +
    theme_phip() +
    theme(
      legend.position      = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title         = element_blank()
    )

  ggsave(static_path_bmc, p_bmc, width = 10, height = 10, dpi = 300, device = "svg", bg = "white")
  ggsave(sub("\\.svg$", ".png", static_path_bmc), p_bmc, width = 10, height = 10, dpi = 300, device = "png", bg = "white")

  p_inter_bmc <- plot_ly() %>%
    add_markers(
      data   = pep_small_bmc,
      x      = ~x_jit, y = ~y_jit,
      type   = "scatter", mode = "markers",
      marker = list(size = 5, color = "grey80", opacity = 0.3),
      hoverinfo  = "skip",
      showlegend = FALSE
    ) %>%
    add_markers(
      data   = highlight_bmc,
      x      = ~x_jit, y = ~y_jit,
      type   = "scatter", mode = "markers",
      color  = ~cat_label, colors = breastmilk_cmp_cols,
      marker = list(
        size    = 12,
        symbol  = "circle",
        opacity = 0.85,
        line    = list(color = "rgba(0,0,0,0)", width = 0)
      ),
      text = ~paste0(
        "Feature: ", feature,
        "<br>Rank: ", rank,
        "<br>x %+: ", round(x_pct, 2),
        "<br>y %+: ", round(y_pct, 2),
        "<br>BH category: ", category_rank_bh
      ),
      hoverinfo  = "text",
      showlegend = TRUE
    ) %>%
    layout(
      title  = cfg$title,
      xaxis  = list(title = cfg$xlab, range = c(0, 100)),
      yaxis  = list(title = cfg$ylab, range = c(0, 100)),
      shapes = list(
        list(type = "line", x0 = 0, y0 = 0, x1 = 100, y1 = 100,
             line = list(dash = "dash", color = "grey"))
      ),
      legend = list(x = 0, y = 1, xanchor = "left", yanchor = "top")
    )

  htmlwidgets::saveWidget(p_inter_bmc, interactive_path_bmc, selfcontained = FALSE)
  cat("Saved breastmilk comparison scatter to:", static_path_bmc, "\n")
  cat("Saved breastmilk comparison interactive to:", interactive_path_bmc, "\n")
}
