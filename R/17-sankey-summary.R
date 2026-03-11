#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

timepoints <- c("T2", "T6", "T8")
z_cols <- paste0("z_", timepoints)
top_n <- 4L
timepoint_spacing <- 0.7
plot_width_mm <- 100
plot_height_mm <- 105
plot_font_family <- "Nimbus Sans"
plot_base_size <- 5
bh_sig_threshold <- 0.10
excluded_metadata <- c("CDrisk", "lockdown")

root_in_dir <- file.path("results", "sankey_meta")
out_dir <- file.path(root_in_dir, "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(out_dir, c("summary_heatmap_long.csv", "summary_heatmap_T_obs_stand.svg")))

metadata_specs <- list(
  siblings = c("siblings", "no_siblings"),
  delmode = c("delmode_VG", "delmode_CS"),
  delplace = c("delplace_home", "delplace_hospital"),
  feeding = c("feeding_BF", "feeding_MF"),
  preterm = c("preterm_yes", "preterm_no"),
  smoking = c("smoking_yes", "smoking_no"),
  sex = c("sex_male", "sex_female")
)

metadata_specs <- metadata_specs[setdiff(names(metadata_specs), excluded_metadata)]

forced_feeding_features <- tibble::tribble(
  ~forced_rank, ~feature, ~forced_order,
  "anno_food_item", "milk", 1L,
  "anno_food_item", "pig", 2L,
  "is_allergens", "TRUE", 3L,
  "anno_is_food", "TRUE", 4L,
  "species", "Listeria innocua", 5L
)
excluded_feeding_features <- c("Sus scrofa")

report_helpers_path <- file.path("results", "report_interesting_metadata", "report_helpers.R")
base_results_dir <- file.path("results", "results")
if (file.exists(report_helpers_path)) {
  report_env <- new.env(parent = baseenv())
  suppressWarnings(source(report_helpers_path, local = report_env))
  if (exists("root_results_dir", envir = report_env, inherits = FALSE)) {
    base_results_dir <- get("root_results_dir", envir = report_env, inherits = FALSE)
  }
}

make_comparison_ids <- function(lhs_suffix, rhs_suffix, tps = timepoints) {
  ids <- vapply(
    tps,
    function(tp) {
      paste0("kid_serum_", tp, "_", lhs_suffix, "_vs_kid_serum_", tp, "_", rhs_suffix)
    },
    FUN.VALUE = character(1)
  )
  names(ids) <- tps
  ids
}

calc_median_z <- function(z1, z2, z3) {
  stats::median(c(z1, z2, z3), na.rm = TRUE)
}

format_feature_label <- function(x) {
  dplyr::if_else(
    x == "Middle East respiratory syndrome-related coronavirus",
    "Middle East respiratory syndrome-\nrelated coronavirus",
    x
  )
}

load_max_delta_by_feature <- function(cmp_id) {
  delta_fdr_path <- file.path(base_results_dir, cmp_id, "DELTA_framework", "delta_table_fdr.csv")
  delta_path <- file.path(base_results_dir, cmp_id, "DELTA_framework", "delta_table.csv")
  if (file.exists(delta_fdr_path)) {
    in_path <- delta_fdr_path
  } else if (file.exists(delta_path)) {
    warning("Missing delta_table_fdr.csv for ", cmp_id, "; falling back to delta_table.csv")
    in_path <- delta_path
  } else {
    warning("Missing delta tables for ", cmp_id, ": ", delta_fdr_path)
    return(tibble::tibble(feature = character(), max_delta = numeric(), rank = character(), p_adj_bh = numeric(), p_perm = numeric()))
  }

  d <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")
  if (!("feature" %in% names(d))) {
    return(tibble::tibble(feature = character(), max_delta = numeric(), rank = character(), p_adj_bh = numeric(), p_perm = numeric()))
  }
  if (!("max_delta" %in% names(d))) d$max_delta <- NA_real_
  if (!("rank" %in% names(d))) d$rank <- NA_character_
  if (!("T_obs_stand" %in% names(d))) d$T_obs_stand <- 0
  if (!("p_adj_bh" %in% names(d))) d$p_adj_bh <- NA_real_
  if (!("p_perm" %in% names(d))) d$p_perm <- NA_real_

  d %>%
    dplyr::mutate(
      feature = as.character(.data$feature),
      T_obs_stand = as.numeric(.data$T_obs_stand),
      max_delta = as.numeric(.data$max_delta),
      rank = as.character(.data$rank),
      p_adj_bh = as.numeric(.data$p_adj_bh),
      p_perm = as.numeric(.data$p_perm)
    ) %>%
    dplyr::group_by(.data$feature) %>%
    dplyr::slice_max(order_by = abs(.data$T_obs_stand), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select("feature", "max_delta", "rank", "p_adj_bh", "p_perm")
}

collect_one_metadata <- function(
  meta_name,
  meta_suffixes,
  use_max_delta_filter = FALSE,
  max_delta_threshold = 0.20,
  select_n = top_n,
  deduplicate_features = TRUE
) {
  in_path <- file.path(root_in_dir, meta_name, "feature_states_zthr0.csv")
  if (!file.exists(in_path)) {
    warning("Missing metadata input file: ", in_path)
    return(NULL)
  }

  df <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")
  required_cols <- c("feature", "class_triplet", "z_T2", "z_T6", "z_T8")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    warning("Skipping ", meta_name, "; missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  comparison_ids <- make_comparison_ids(meta_suffixes[[1]], meta_suffixes[[2]])
  comparison_name <- paste0(meta_suffixes[[1]], "_vs_", meta_suffixes[[2]])

  # Use p_perm from DELTA tables; avoid name collisions with p_perm_* potentially present
  # in feature_states_zthr0.csv.
  df <- df %>%
    dplyr::select(-dplyr::any_of(paste0("p_perm_", timepoints)))

  all_candidates <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(median_z_score = calc_median_z(.data$z_T2, .data$z_T6, .data$z_T8)) %>%
    dplyr::ungroup()

  for (tp in timepoints) {
    max_delta_tbl <- load_max_delta_by_feature(comparison_ids[[tp]]) %>%
      dplyr::rename(
        !!paste0("max_delta_", tp) := "max_delta",
        !!paste0("delta_rank_", tp) := "rank",
        !!paste0("p_adj_bh_", tp) := "p_adj_bh",
        !!paste0("p_perm_", tp) := "p_perm"
      )
    all_candidates <- all_candidates %>%
      dplyr::left_join(max_delta_tbl, by = "feature")
  }

  all_candidates <- all_candidates %>%
    dplyr::mutate(
      max_delta_any = pmax(.data$max_delta_T2, .data$max_delta_T6, .data$max_delta_T8, na.rm = TRUE),
      max_delta_any = dplyr::if_else(is.infinite(.data$max_delta_any), NA_real_, .data$max_delta_any),
      max_delta_timepoint = dplyr::case_when(
        is.na(.data$max_delta_any) ~ NA_character_,
        dplyr::coalesce(.data$max_delta_T2, -Inf) >= dplyr::coalesce(.data$max_delta_T6, -Inf) &
          dplyr::coalesce(.data$max_delta_T2, -Inf) >= dplyr::coalesce(.data$max_delta_T8, -Inf) ~ "T2",
        dplyr::coalesce(.data$max_delta_T6, -Inf) >= dplyr::coalesce(.data$max_delta_T8, -Inf) ~ "T6",
        TRUE ~ "T8"
      ),
      rank = dplyr::case_when(
        .data$max_delta_timepoint == "T2" ~ .data$delta_rank_T2,
        .data$max_delta_timepoint == "T6" ~ .data$delta_rank_T6,
        .data$max_delta_timepoint == "T8" ~ .data$delta_rank_T8,
        TRUE ~ NA_character_
      )
    )

  eligible <- all_candidates %>%
    dplyr::filter(.data$class_triplet %in% c("+++", "---"))

  if (isTRUE(use_max_delta_filter)) {
    eligible <- eligible %>%
      dplyr::filter(.data$max_delta_any > max_delta_threshold)
  }

  forced_rows <- tibble::tibble()
  if (identical(meta_name, "feeding")) {
    forced_rows <- forced_feeding_features %>%
      dplyr::left_join(all_candidates, by = "feature") %>%
      dplyr::filter(!is.na(.data$median_z_score)) %>%
      dplyr::arrange(.data$forced_order) %>%
      dplyr::mutate(rank = .data$forced_rank) %>%
      dplyr::mutate(selection_side = "forced_priority")
  }

  if (nrow(eligible) == 0L && nrow(forced_rows) == 0L) {
    warning("No eligible features for metadata ", meta_name, " under current filtering.")
    return(NULL)
  }

  selected <- tibble::tibble()
  if (nrow(eligible) > 0L) {
    selected_max <- eligible %>%
      dplyr::arrange(dplyr::desc(.data$median_z_score), .data$feature) %>%
      dplyr::slice_head(n = select_n) %>%
      dplyr::mutate(selection_side = "max_median")

    selected_min <- eligible %>%
      dplyr::arrange(.data$median_z_score, .data$feature) %>%
      dplyr::slice_head(n = select_n) %>%
      dplyr::mutate(selection_side = "min_median")

    selected <- dplyr::bind_rows(selected_max, selected_min)
  }

  if (nrow(forced_rows) > 0L) {
    selected <- dplyr::bind_rows(forced_rows, selected)
  }

  if (identical(meta_name, "feeding")) {
    selected <- selected %>%
      dplyr::filter(!(.data$feature %in% excluded_feeding_features))
  }

  if (isTRUE(deduplicate_features)) {
    selected <- selected %>%
      dplyr::distinct(.data$rank, .data$feature, .keep_all = TRUE)
  }

  selected %>%
    dplyr::mutate(
      metadata = meta_name,
      comparison_name = comparison_name,
      comparison_T2 = comparison_ids[["T2"]],
      comparison_T6 = comparison_ids[["T6"]],
      comparison_T8 = comparison_ids[["T8"]]
    )
}

write_extreme_table_by_comparison <- function(selected_tbl, suffix = "", top_k = 10L) {
  if (is.null(selected_tbl) || nrow(selected_tbl) == 0L) {
    warning("No rows available for top ", top_k, " max/min table suffix ", suffix)
    return(invisible(NULL))
  }

  out_tbl <- selected_tbl %>%
    dplyr::mutate(
      selection_group = dplyr::if_else(.data$selection_side == "max_median", "max", "min"),
      selection_order = dplyr::if_else(.data$selection_group == "max", 1L, 2L),
      score_order = dplyr::if_else(.data$selection_group == "max", -.data$median_z_score, .data$median_z_score)
    ) %>%
    dplyr::group_by(.data$comparison_name, .data$selection_group) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::arrange(.data$comparison_name, .data$selection_order, .data$score_order, .data$feature) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      comparison_name = .data$comparison_name,
      metadata = .data$metadata,
      selection_group = .data$selection_group,
      rank = .data$rank,
      feature = .data$feature,
      median_z_score = .data$median_z_score,
      z_T2 = .data$z_T2,
      z_T6 = .data$z_T6,
      z_T8 = .data$z_T8,
      max_delta_any = .data$max_delta_any
    )

  out_path <- file.path(out_dir, paste0("summary_table_top", top_k, "_max_min_by_comparison", suffix, ".csv"))
  readr::write_csv(out_tbl, out_path)
  message("Saved top ", top_k, " max/min table: ", basename(out_path))
  invisible(NULL)
}

build_summary_outputs <- function(selected_tbl, suffix = "", use_inferno = FALSE) {
  if (is.null(selected_tbl) || nrow(selected_tbl) == 0L) {
    warning("No rows available for summary suffix ", suffix)
    return(invisible(NULL))
  }

  required_p_perm_cols <- paste0("p_perm_", timepoints)
  missing_p_perm_cols <- setdiff(required_p_perm_cols, names(selected_tbl))
  if (length(missing_p_perm_cols) > 0L) {
    for (nm in missing_p_perm_cols) selected_tbl[[nm]] <- NA_real_
  }

  selected_tbl <- selected_tbl %>%
    dplyr::mutate(
      feeding_priority_order = dplyr::case_when(
        .data$metadata == "feeding" & .data$rank == "anno_food_item" & .data$feature == "milk" ~ 1L,
        .data$metadata == "feeding" & .data$rank == "anno_food_item" & .data$feature == "pig" ~ 2L,
        .data$metadata == "feeding" & .data$rank == "is_allergens" & .data$feature == "TRUE" ~ 3L,
        .data$metadata == "feeding" & .data$rank == "anno_is_food" & .data$feature == "TRUE" ~ 4L,
        .data$metadata == "feeding" & .data$rank == "species" & .data$feature == "Listeria innocua" ~ 5L,
        TRUE ~ 999L
      ),
      feature_display_plot = dplyr::case_when(
        .data$metadata == "feeding" & .data$feature == "TRUE" ~ paste(.data$rank, .data$feature),
        TRUE ~ .data$feature
      )
    ) %>%
    dplyr::arrange(.data$metadata, .data$feeding_priority_order, dplyr::desc(.data$median_z_score), .data$feature)

  summary_table <- selected_tbl %>%
    dplyr::transmute(
      comparison_name = .data$comparison_name,
      metadata = .data$metadata,
      rank = .data$rank,
      feature = .data$feature,
      median_z_score = .data$median_z_score
    )
  readr::write_csv(
    summary_table,
    file.path(out_dir, paste0("summary_table_comparison_feature_median_z", suffix, ".csv"))
  )

  detailed_table <- selected_tbl %>%
    dplyr::select(
      "metadata",
      "comparison_name",
      "comparison_T2",
      "comparison_T6",
      "comparison_T8",
      "rank",
      "feature",
      "class_triplet",
      "selection_side",
      "median_z_score",
      "max_delta_any",
      "p_adj_bh_T2",
      "p_adj_bh_T6",
      "p_adj_bh_T8",
      "p_perm_T2",
      "p_perm_T6",
      "p_perm_T8",
      "z_T2",
      "z_T6",
      "z_T8",
      "max_delta_T2",
      "max_delta_T6",
      "max_delta_T8"
    )
  readr::write_csv(
    detailed_table,
    file.path(out_dir, paste0("summary_selected_features_detailed", suffix, ".csv"))
  )

  q_cols <- paste0("p_adj_bh_", timepoints)
  p_cols <- paste0("p_perm_", timepoints)
  md_cols <- paste0("max_delta_", timepoints)
  x_ib <- 2 - timepoint_spacing
  x_im3 <- 2
  x_im12 <- 2 + timepoint_spacing

  bubble_long <- selected_tbl %>%
    dplyr::transmute(
      metadata = .data$metadata,
      comparison_name = .data$comparison_name,
      rank = .data$rank,
      feature = .data$feature,
      median_z_score = .data$median_z_score,
      feature_display = format_feature_label(.data$feature_display_plot),
      feature_row = paste0(.data$metadata, " | ", .data$feature_display),
      z_T2 = .data$z_T2,
      z_T6 = .data$z_T6,
      z_T8 = .data$z_T8,
      p_adj_bh_T2 = .data$p_adj_bh_T2,
      p_adj_bh_T6 = .data$p_adj_bh_T6,
      p_adj_bh_T8 = .data$p_adj_bh_T8,
      p_perm_T2 = .data$p_perm_T2,
      p_perm_T6 = .data$p_perm_T6,
      p_perm_T8 = .data$p_perm_T8,
      max_delta_T2 = .data$max_delta_T2,
      max_delta_T6 = .data$max_delta_T6,
      max_delta_T8 = .data$max_delta_T8
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(z_cols, q_cols, p_cols, md_cols)),
      names_to = c(".value", "timepoint"),
      names_pattern = "(z|p_adj_bh|p_perm|max_delta)_(T2|T6|T8)"
    ) %>%
    dplyr::rename(NES = "z") %>%
    dplyr::mutate(
      timepoint = dplyr::recode(.data$timepoint, T2 = "I_B", T6 = "I_M3", T8 = "I_M12"),
      timepoint = factor(.data$timepoint, levels = c("I_B", "I_M3", "I_M12")),
      is_bh_significant = !is.na(.data$p_adj_bh) & .data$p_adj_bh <= bh_sig_threshold,
      is_nominal_significant = !is.na(.data$p_perm) & .data$p_perm < 0.05,
      is_nominal_only = .data$is_nominal_significant & !.data$is_bh_significant,
      max_delta_pct = 100 * .data$max_delta,
      max_delta_bin = dplyr::case_when(
        is.na(.data$max_delta_pct) ~ "0-5%",
        .data$max_delta_pct <= 5 ~ "0-5%",
        .data$max_delta_pct <= 10 ~ "5-10%",
        .data$max_delta_pct <= 20 ~ "10-20%",
        .data$max_delta_pct <= 30 ~ "20-30%",
        TRUE ~ ">40%"
      ),
      max_delta_bin = factor(.data$max_delta_bin, levels = c("0-5%", "5-10%", "10-20%", "20-30%", ">40%")),
      x_num = dplyr::case_when(
        .data$timepoint == "I_B" ~ x_ib,
        .data$timepoint == "I_M3" ~ x_im3,
        .data$timepoint == "I_M12" ~ x_im12,
        TRUE ~ NA_real_
      )
    )

  row_levels <- selected_tbl %>%
    dplyr::transmute(
      metadata = .data$metadata,
      feeding_priority_order = .data$feeding_priority_order,
      feature_row = paste0(.data$metadata, " | ", format_feature_label(.data$feature_display_plot)),
      median_z_score = .data$median_z_score
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$metadata, .data$feeding_priority_order, dplyr::desc(.data$median_z_score), .data$feature_row) %>%
    dplyr::pull(.data$feature_row)

  bubble_long <- bubble_long %>%
    dplyr::mutate(feature_row = factor(.data$feature_row, levels = rev(row_levels)))

  readr::write_csv(
    bubble_long,
    file.path(out_dir, paste0("summary_bubble_long", suffix, ".csv"))
  )

  line_gap <- min(0.16, timepoint_spacing * 0.35)
  line_df <- bubble_long %>%
    dplyr::distinct(.data$feature_row) %>%
    tidyr::crossing(segment = c("I_B_I_M3", "I_M3_I_M12")) %>%
    dplyr::mutate(
      x = dplyr::if_else(
        .data$segment == "I_B_I_M3",
        x_ib + line_gap,
        x_im3 + line_gap
      ),
      xend = dplyr::if_else(
        .data$segment == "I_B_I_M3",
        x_im3 - line_gap,
        x_im12 - line_gap
      ),
      y = .data$feature_row,
      yend = .data$feature_row
    ) %>%
    dplyr::select("x", "xend", "y", "yend")

  fill_scale <- scale_fill_distiller(
    type = "div",
    palette = "RdYlGn",
    direction = 1,
    name = "NES",
    limits = c(-3, NA),
    oob = scales::squish
  )

  old_height_cm <- 2.5 * max(12, 0.35 * length(row_levels))
  size_scale <- plot_height_mm / (old_height_cm * 10)

  p_bubble <- ggplot(bubble_long, aes(x = .data$x_num, y = .data$feature_row)) +
    geom_segment(
      data = line_df,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      inherit.aes = FALSE,
      linewidth = 0.14,
      color = "grey45",
      alpha = 0.75
    ) +
    geom_point(
      aes(size = .data$max_delta_bin, fill = .data$NES),
      shape = 21,
      color = "grey20",
      stroke = 0.2 * size_scale
    ) +
    geom_point(
      data = dplyr::filter(bubble_long, .data$is_nominal_only),
      inherit.aes = FALSE,
      aes(x = .data$x_num, y = .data$feature_row),
      shape = 16,
      size = 1.35 * size_scale,
      color = "black"
    ) +
    geom_point(
      data = dplyr::filter(bubble_long, .data$is_bh_significant),
      inherit.aes = FALSE,
      aes(x = .data$x_num, y = .data$feature_row),
      shape = 8,
      size = 2.65 * size_scale,
      color = "black",
      stroke = 0.45 * size_scale
    ) +
    fill_scale +
    scale_size_manual(
      values = c(
        "0-5%" = 4.97,
        "5-10%" = 6.46,
        "10-20%" = 7.95,
        "20-30%" = 9.44,
        ">40%" = 10.93
      ) * size_scale,
      name = "Max delta",
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(x_ib, x_im3, x_im12),
      labels = c("I_B", "I_M3", "I_M12"),
      limits = c(1, 3)
    ) +
    labs(
      x = "Timepoint",
      y = "Metadata | Feature"
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

  out_svg <- file.path(out_dir, paste0("summary_bubbleplot_NES_nominal_p", suffix, ".svg"))
  ggplot2::ggsave(
    filename = out_svg,
    plot = p_bubble,
    width = plot_width_mm,
    height = plot_height_mm,
    units = "mm",
    dpi = 300,
    device = "svg",
    bg = "white"
  )

  message("Saved summary outputs in: ", out_dir)
  message(" - ", basename(file.path(out_dir, paste0("summary_table_comparison_feature_median_z", suffix, ".csv"))))
  message(" - ", basename(file.path(out_dir, paste0("summary_selected_features_detailed", suffix, ".csv"))))
  message(" - ", basename(file.path(out_dir, paste0("summary_bubble_long", suffix, ".csv"))))
  message(" - ", basename(out_svg))

  invisible(NULL)
}

strip_rank_prefix <- function(x) {
  x <- trimws(gsub(",$", "", x))
  x <- gsub("^(species|genus|family|order|class|phylum)\\s*;?\\s+", "", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("^anno_[A-Za-z0-9_]+\\s+", "", x, perl = TRUE)
  x <- gsub("^is_[A-Za-z0-9_]+\\s+", "", x, perl = TRUE)
  trimws(x)
}

normalize_feature_key <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

match_requested_feature <- function(requested_feature, available_features) {
  if (length(available_features) == 0L) {
    return(list(feature = NA_character_, type = "missing"))
  }

  if (requested_feature %in% available_features) {
    return(list(feature = requested_feature, type = "exact"))
  }

  requested_key <- normalize_feature_key(requested_feature)
  available_keys <- normalize_feature_key(available_features)

  idx_norm <- which(available_keys == requested_key)
  if (length(idx_norm) >= 1L) {
    return(list(feature = available_features[[idx_norm[[1]]]], type = "normalized"))
  }

  d <- as.numeric(utils::adist(requested_key, available_keys))
  best <- which(d == min(d))
  if (length(best) == 1L && d[[best[[1]]]] <= 3) {
    return(list(feature = available_features[[best[[1]]]], type = "fuzzy"))
  }

  list(feature = NA_character_, type = "missing")
}

collect_metadata_with_deltas <- function(meta_name, meta_suffixes) {
  in_path <- file.path(root_in_dir, meta_name, "feature_states_zthr0.csv")
  if (!file.exists(in_path)) {
    warning("Missing metadata input file: ", in_path)
    return(tibble::tibble())
  }

  df <- readr::read_csv(in_path, show_col_types = FALSE, name_repair = "unique_quiet")
  required_cols <- c("feature", "z_T2", "z_T6", "z_T8")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    warning("Skipping ", meta_name, "; missing columns: ", paste(missing_cols, collapse = ", "))
    return(tibble::tibble())
  }

  if (!("class_triplet" %in% names(df))) df$class_triplet <- NA_character_

  # Use p_perm from DELTA tables; avoid name collisions with p_perm_* potentially present
  # in feature_states_zthr0.csv.
  df <- df %>%
    dplyr::select(-dplyr::any_of(paste0("p_perm_", timepoints)))

  comparison_ids <- make_comparison_ids(meta_suffixes[[1]], meta_suffixes[[2]])
  comparison_name <- paste0(meta_suffixes[[1]], "_vs_", meta_suffixes[[2]])

  full_tbl <- df %>%
    dplyr::mutate(
      feature = as.character(.data$feature),
      median_z_score = apply(
        dplyr::across(dplyr::all_of(z_cols)),
        1,
        function(v) stats::median(as.numeric(v), na.rm = TRUE)
      )
    )

  for (tp in timepoints) {
    max_delta_tbl <- load_max_delta_by_feature(comparison_ids[[tp]]) %>%
      dplyr::rename(
        !!paste0("max_delta_", tp) := "max_delta",
        !!paste0("delta_rank_", tp) := "rank",
        !!paste0("p_adj_bh_", tp) := "p_adj_bh",
        !!paste0("p_perm_", tp) := "p_perm"
      )
    full_tbl <- full_tbl %>%
      dplyr::left_join(max_delta_tbl, by = "feature")
  }

  required_delta_cols <- c(
    paste0("max_delta_", timepoints),
    paste0("delta_rank_", timepoints),
    paste0("p_adj_bh_", timepoints),
    paste0("p_perm_", timepoints)
  )
  missing_delta_cols <- setdiff(required_delta_cols, names(full_tbl))
  if (length(missing_delta_cols) > 0L) {
    for (nm in missing_delta_cols) {
      full_tbl[[nm]] <- if (grepl("^delta_rank_", nm)) NA_character_ else NA_real_
    }
  }

  full_tbl %>%
    dplyr::mutate(
      max_delta_any = pmax(.data$max_delta_T2, .data$max_delta_T6, .data$max_delta_T8, na.rm = TRUE),
      max_delta_any = dplyr::if_else(is.infinite(.data$max_delta_any), NA_real_, .data$max_delta_any),
      max_delta_timepoint = dplyr::case_when(
        is.na(.data$max_delta_any) ~ NA_character_,
        dplyr::coalesce(.data$max_delta_T2, -Inf) >= dplyr::coalesce(.data$max_delta_T6, -Inf) &
          dplyr::coalesce(.data$max_delta_T2, -Inf) >= dplyr::coalesce(.data$max_delta_T8, -Inf) ~ "T2",
        dplyr::coalesce(.data$max_delta_T6, -Inf) >= dplyr::coalesce(.data$max_delta_T8, -Inf) ~ "T6",
        TRUE ~ "T8"
      ),
      rank = dplyr::case_when(
        .data$max_delta_timepoint == "T2" ~ .data$delta_rank_T2,
        .data$max_delta_timepoint == "T6" ~ .data$delta_rank_T6,
        .data$max_delta_timepoint == "T8" ~ .data$delta_rank_T8,
        TRUE ~ NA_character_
      ),
      metadata = meta_name,
      comparison_name = comparison_name,
      comparison_T2 = comparison_ids[["T2"]],
      comparison_T6 = comparison_ids[["T6"]],
      comparison_T8 = comparison_ids[["T8"]]
    )
}

build_custom_feature_plot <- function(custom_feature_spec, suffix = "_custom_ranked") {
  if (is.null(custom_feature_spec) || nrow(custom_feature_spec) == 0L) {
    warning("No custom feature specification provided for suffix ", suffix)
    return(invisible(NULL))
  }

  comparison_order <- c("feeding", "delmode", "siblings", "smoking", "preterm", "delplace", "sex")

  custom_spec <- custom_feature_spec %>%
    dplyr::mutate(
      metadata = dplyr::recode(.data$metadata, cdrisk = "CDrisk", .default = .data$metadata),
      selection_side = dplyr::recode(.data$selection_side, max = "max_median", min = "min_median"),
      feature_requested = trimws(gsub(",$", "", .data$feature_requested)),
      rank_hint = dplyr::if_else(
        grepl("^(species|genus|family|order|class|phylum|anno_[A-Za-z0-9_]+|is_[A-Za-z0-9_]+)\\s+", .data$feature_requested, perl = TRUE),
        sub("^((species|genus|family|order|class|phylum|anno_[A-Za-z0-9_]+|is_[A-Za-z0-9_]+))\\s+.*$", "\\1", .data$feature_requested, perl = TRUE),
        NA_character_
      ),
      feature_clean = strip_rank_prefix(.data$feature_requested),
      comparison_order = match(.data$metadata, comparison_order),
      side_order = dplyr::if_else(.data$selection_side == "max_median", 1L, 2L),
      feature_order = dplyr::row_number()
    ) %>%
    dplyr::filter(!.data$metadata %in% excluded_metadata) %>%
    dplyr::arrange(.data$comparison_order, .data$side_order, .data$feature_order)

  missing_meta <- setdiff(unique(custom_spec$metadata), names(metadata_specs))
  if (length(missing_meta) > 0L) {
    warning("Custom spec metadata missing from metadata_specs: ", paste(missing_meta, collapse = ", "))
  }

  meta_tables <- lapply(unique(custom_spec$metadata), function(meta_name) {
    if (!(meta_name %in% names(metadata_specs))) return(tibble::tibble())
    collect_metadata_with_deltas(meta_name, metadata_specs[[meta_name]])
  })
  names(meta_tables) <- unique(custom_spec$metadata)

  matched_rows <- vector("list", nrow(custom_spec))
  missing_lines <- character()

  for (i in seq_len(nrow(custom_spec))) {
    spec_row <- custom_spec[i, , drop = FALSE]
    meta_name <- spec_row$metadata[[1]]
    available_tbl <- meta_tables[[meta_name]]

    if (is.null(available_tbl) || nrow(available_tbl) == 0L) {
      missing_lines <- c(
        missing_lines,
        paste0(meta_name, " | ", spec_row$selection_side[[1]], " | ", spec_row$feature_requested[[1]], " (metadata table missing)")
      )
      next
    }

    rank_hint <- spec_row$rank_hint[[1]]
    candidate_tbl <- available_tbl
    if (!is.na(rank_hint)) {
      ranked_tbl <- candidate_tbl %>%
        dplyr::filter(.data$rank == rank_hint)
      if (nrow(ranked_tbl) > 0L) {
        candidate_tbl <- ranked_tbl
      }
    }

    matched <- match_requested_feature(spec_row$feature_clean[[1]], candidate_tbl$feature)
    if (is.na(matched$feature)) {
      missing_lines <- c(
        missing_lines,
        paste0(meta_name, " | ", spec_row$selection_side[[1]], " | ", spec_row$feature_requested[[1]], " (not matched)")
      )
      next
    }

    feature_candidates <- candidate_tbl %>%
      dplyr::filter(.data$feature == matched$feature) %>%
      dplyr::slice_head(n = 1)

    feature_row <- feature_candidates %>%
      dplyr::mutate(
        rank = dplyr::if_else(!is.na(rank_hint), rank_hint, .data$rank),
        selection_side = spec_row$selection_side[[1]],
        comparison_order = spec_row$comparison_order[[1]],
        side_order = spec_row$side_order[[1]],
        feature_order = spec_row$feature_order[[1]],
        feature_requested = spec_row$feature_requested[[1]],
        feature_match_type = matched$type
      )

    matched_rows[[i]] <- feature_row
  }

  custom_selected <- dplyr::bind_rows(matched_rows) %>%
    dplyr::arrange(.data$comparison_order, .data$side_order, .data$feature_order)

  if (length(missing_lines) > 0L) {
    warning("Custom feature rows not plotted:\n - ", paste(missing_lines, collapse = "\n - "))
  }

  if (nrow(custom_selected) == 0L) {
    warning("No custom features matched input data for suffix ", suffix)
    return(invisible(NULL))
  }

  required_p_perm_cols <- paste0("p_perm_", timepoints)
  missing_p_perm_cols <- setdiff(required_p_perm_cols, names(custom_selected))
  if (length(missing_p_perm_cols) > 0L) {
    for (nm in missing_p_perm_cols) custom_selected[[nm]] <- NA_real_
  }

  custom_table <- custom_selected %>%
    dplyr::transmute(
      metadata = .data$metadata,
      comparison_name = .data$comparison_name,
      selection_side = .data$selection_side,
      feature_requested = .data$feature_requested,
      feature = .data$feature,
      feature_match_type = .data$feature_match_type,
      class_triplet = .data$class_triplet,
      median_z_score = .data$median_z_score,
      rank = .data$rank,
      p_adj_bh_T2 = .data$p_adj_bh_T2,
      p_adj_bh_T6 = .data$p_adj_bh_T6,
      p_adj_bh_T8 = .data$p_adj_bh_T8,
      p_perm_T2 = .data$p_perm_T2,
      p_perm_T6 = .data$p_perm_T6,
      p_perm_T8 = .data$p_perm_T8,
      z_T2 = .data$z_T2,
      z_T6 = .data$z_T6,
      z_T8 = .data$z_T8,
      max_delta_T2 = .data$max_delta_T2,
      max_delta_T6 = .data$max_delta_T6,
      max_delta_T8 = .data$max_delta_T8,
      max_delta_any = .data$max_delta_any
    )
  readr::write_csv(custom_table, file.path(out_dir, paste0("summary_selected_features_detailed", suffix, ".csv")))

  q_cols <- paste0("p_adj_bh_", timepoints)
  p_cols <- paste0("p_perm_", timepoints)
  md_cols <- paste0("max_delta_", timepoints)

  custom_selected <- custom_selected %>%
    dplyr::mutate(
      metadata_display = dplyr::recode(.data$metadata, CDrisk = "cdrisk", .default = .data$metadata),
      side_display = dplyr::if_else(.data$selection_side == "max_median", "max", "min"),
      feature_display = dplyr::case_when(
        .data$metadata == "feeding" & .data$feature == "TRUE" ~ paste(.data$rank, .data$feature),
        TRUE ~ .data$feature
      ),
      feature_display = format_feature_label(.data$feature_display),
      feature_row = paste0(.data$metadata_display, " ", .data$side_display, " | ", .data$feature_display)
    )

  row_levels <- custom_selected %>%
    dplyr::arrange(.data$comparison_order, .data$side_order, .data$feature_order) %>%
    dplyr::pull(.data$feature_row)

  x_ib <- 2 - timepoint_spacing
  x_im3 <- 2
  x_im12 <- 2 + timepoint_spacing

  bubble_long <- custom_selected %>%
    dplyr::transmute(
      metadata = .data$metadata,
      metadata_display = .data$metadata_display,
      comparison_name = .data$comparison_name,
      selection_side = .data$selection_side,
      feature_requested = .data$feature_requested,
      rank = .data$rank,
      feature = .data$feature,
      feature_match_type = .data$feature_match_type,
      feature_row = factor(.data$feature_row, levels = rev(row_levels)),
      z_T2 = .data$z_T2,
      z_T6 = .data$z_T6,
      z_T8 = .data$z_T8,
      p_adj_bh_T2 = .data$p_adj_bh_T2,
      p_adj_bh_T6 = .data$p_adj_bh_T6,
      p_adj_bh_T8 = .data$p_adj_bh_T8,
      p_perm_T2 = .data$p_perm_T2,
      p_perm_T6 = .data$p_perm_T6,
      p_perm_T8 = .data$p_perm_T8,
      max_delta_T2 = .data$max_delta_T2,
      max_delta_T6 = .data$max_delta_T6,
      max_delta_T8 = .data$max_delta_T8
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(z_cols, q_cols, p_cols, md_cols)),
      names_to = c(".value", "timepoint"),
      names_pattern = "(z|p_adj_bh|p_perm|max_delta)_(T2|T6|T8)"
    ) %>%
    dplyr::rename(NES = "z") %>%
    dplyr::mutate(
      timepoint = dplyr::recode(.data$timepoint, T2 = "I_B", T6 = "I_M3", T8 = "I_M12"),
      timepoint = factor(.data$timepoint, levels = c("I_B", "I_M3", "I_M12")),
      is_bh_significant = !is.na(.data$p_adj_bh) & .data$p_adj_bh <= bh_sig_threshold,
      is_nominal_significant = !is.na(.data$p_perm) & .data$p_perm < 0.05,
      is_nominal_only = .data$is_nominal_significant & !.data$is_bh_significant,
      max_delta_pct = 100 * .data$max_delta,
      max_delta_bin = dplyr::case_when(
        is.na(.data$max_delta_pct) ~ "0-5%",
        .data$max_delta_pct <= 5 ~ "0-5%",
        .data$max_delta_pct <= 10 ~ "5-10%",
        .data$max_delta_pct <= 20 ~ "10-20%",
        .data$max_delta_pct <= 30 ~ "20-30%",
        TRUE ~ ">40%"
      ),
      max_delta_bin = factor(.data$max_delta_bin, levels = c("0-5%", "5-10%", "10-20%", "20-30%", ">40%")),
      x_num = dplyr::case_when(
        .data$timepoint == "I_B" ~ x_ib,
        .data$timepoint == "I_M3" ~ x_im3,
        .data$timepoint == "I_M12" ~ x_im12,
        TRUE ~ NA_real_
      )
    )

  readr::write_csv(bubble_long, file.path(out_dir, paste0("summary_bubble_long", suffix, ".csv")))

  line_gap <- min(0.16, timepoint_spacing * 0.35)
  line_df <- bubble_long %>%
    dplyr::distinct(.data$feature_row) %>%
    tidyr::crossing(segment = c("I_B_I_M3", "I_M3_I_M12")) %>%
    dplyr::mutate(
      x = dplyr::if_else(
        .data$segment == "I_B_I_M3",
        x_ib + line_gap,
        x_im3 + line_gap
      ),
      xend = dplyr::if_else(
        .data$segment == "I_B_I_M3",
        x_im3 - line_gap,
        x_im12 - line_gap
      ),
      y = .data$feature_row,
      yend = .data$feature_row
    ) %>%
    dplyr::select("x", "xend", "y", "yend")

  old_height_cm <- 2.5 * max(12, 0.35 * length(row_levels))
  size_scale <- plot_height_mm / (old_height_cm * 10)

  p_custom <- ggplot(bubble_long, aes(x = .data$x_num, y = .data$feature_row)) +
    geom_segment(
      data = line_df,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      inherit.aes = FALSE,
      linewidth = 0.14,
      color = "grey45",
      alpha = 0.75
    ) +
    geom_point(
      aes(size = .data$max_delta_bin, fill = .data$NES),
      shape = 21,
      color = "grey20",
      stroke = 0.2 * size_scale
    ) +
    geom_point(
      data = dplyr::filter(bubble_long, .data$is_nominal_only),
      inherit.aes = FALSE,
      aes(x = .data$x_num, y = .data$feature_row),
      shape = 16,
      size = 1.35 * size_scale,
      color = "black"
    ) +
    geom_point(
      data = dplyr::filter(bubble_long, .data$is_bh_significant),
      inherit.aes = FALSE,
      aes(x = .data$x_num, y = .data$feature_row),
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
      name = "NES",
      limits = c(-3, NA),
      oob = scales::squish
    ) +
    scale_size_manual(
      values = c(
        "0-5%" = 4.97,
        "5-10%" = 6.46,
        "10-20%" = 7.95,
        "20-30%" = 9.44,
        ">40%" = 10.93
      ) * size_scale,
      name = "Max delta",
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(x_ib, x_im3, x_im12),
      labels = c("I_B", "I_M3", "I_M12"),
      limits = c(1, 3)
    ) +
    labs(
      x = "Timepoint",
      y = "Comparison side | Feature"
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

  out_svg <- file.path(out_dir, paste0("summary_bubbleplot_NES_nominal_p", suffix, ".svg"))
  ggplot2::ggsave(
    filename = out_svg,
    plot = p_custom,
    width = plot_width_mm,
    height = plot_height_mm,
    units = "mm",
    dpi = 300,
    device = "svg",
    bg = "white"
  )

  message("Saved custom ranked outputs in: ", out_dir)
  message(" - ", basename(file.path(out_dir, paste0("summary_selected_features_detailed", suffix, ".csv"))))
  message(" - ", basename(file.path(out_dir, paste0("summary_bubble_long", suffix, ".csv"))))
  message(" - ", basename(out_svg))

  invisible(NULL)
}

selected_default <- dplyr::bind_rows(lapply(names(metadata_specs), function(meta_name) {
  collect_one_metadata(meta_name, metadata_specs[[meta_name]], use_max_delta_filter = FALSE)
}))

build_summary_outputs(selected_default, suffix = "", use_inferno = FALSE)

selected_maxdelta20 <- dplyr::bind_rows(lapply(names(metadata_specs), function(meta_name) {
  collect_one_metadata(
    meta_name,
    metadata_specs[[meta_name]],
    use_max_delta_filter = TRUE,
    max_delta_threshold = 0.20
  )
}))

build_summary_outputs(selected_maxdelta20, suffix = "_maxdelta20", use_inferno = TRUE)

selected_maxdelta20_top10 <- dplyr::bind_rows(lapply(names(metadata_specs), function(meta_name) {
  collect_one_metadata(
    meta_name,
    metadata_specs[[meta_name]],
    use_max_delta_filter = TRUE,
    max_delta_threshold = 0.20,
    select_n = 10L,
    deduplicate_features = FALSE
  )
}))

write_extreme_table_by_comparison(selected_maxdelta20_top10, suffix = "_maxdelta20", top_k = 10L)

custom_feature_spec <- tibble::tribble(
  ~metadata, ~selection_side, ~feature_requested,
  "feeding", "max", "anno_food_item milk",
  "feeding", "max", "anno_food_item pig",
  "feeding", "max", "is_allergens TRUE",
  "feeding", "max", "anno_is_food TRUE",
  "feeding", "max", "species Bos mutus",
  "feeding", "max", "species Porcine epidemic diarrhea virus",
  "feeding", "max", "phylum Cossaviricota",
  "feeding", "min", "species Listeria innocua",
  "feeding", "min", "family Orthoherpesviridae",
  "feeding", "min", "species Lactiplantibacillus plantarum",
  "feeding", "min", "class Betaproteobacteria",
  "feeding", "min", "genus Lachnospira",
  "delmode", "max", "anno_is_food TRUE",
  "delmode", "max", "species Human coronavirus NL63",
  "delmode", "max", "genus Triticum",
  "delmode", "max", "genus Aegilops",
  "delmode", "min", "species Bacteroides heparinolyticus",
  "delmode", "min", "species Bacteroides sp. HF-5092",
  "delmode", "min", "species Simplexvirus humanalpha2",
  "delmode", "min", "species Ashduovirus A2",
  "siblings", "max", "species Streptococcus pneumoniae",
  "siblings", "max", "species Rhinovirus A",
  "siblings", "max", "species Bacteroides eggerthii",
  "siblings", "min", "species Staphylococcus aureus",
  "siblings", "min", "species Ovis aries",
  "siblings", "min", "order Galliformes",
  "smoking", "max", "species Bifidobacterium breve",
  "smoking", "max", "species Coprococcus comes",
  "smoking", "max", "species Orthopneumovirus hominis",
  "smoking", "max", "species Lachnospiraceae bacterium 2_1_58FAA",
  "smoking", "min", "species Betainfluenzavirus influenzae",
  "smoking", "min", "species Poales",
  "smoking", "min", "genus Bacillus",
  "smoking", "min", "species Cytomegalovirus humanbeta5",
  "preterm", "max", "family Erysipelotrichaceae",
  "preterm", "max", "species Holdemanella biformis",
  "preterm", "max", "species Lactobacillus kitasatonis",
  "preterm", "max", "species Firmicutes bacterium CAG:176",
  "preterm", "min", "species Blautia wexlerae",
  "preterm", "min", "genus Vibrio",
  "preterm", "min", "order Compylobacterales",
  "preterm", "min", "species Pecentumvirus list36",
  "delplace", "max", "species Simplexvirus humanalpha1",
  "delplace", "max", "genus Bordetella",
  "delplace", "max", "species Enterovirus C",
  "delplace", "min", "species Gelderlandvirus s16",
  "delplace", "min", "anno_viruses_bacteriophage mastadenovirusA",
  "delplace", "min", "species Cytomegalovirus humanbeta5",
  "CDrisk", "max", "species; Roseolovirus humanbeta7",
  "CDrisk", "max", "species Staphylococcus aureus",
  "CDrisk", "max", "order Eurotiales",
  "CDrisk", "max", "family Aspergillaceae",
  "lockdown", "max", "genus Faecalibacterium",
  "lockdown", "max", "species Fletchervirus CP30A",
  "lockdown", "max", "genus Mediterraneibacter",
  "lockdown", "max", "species Streptococcus agalactiae",
  "lockdown", "min", "species Streptococcus sp. G148",
  "sex", "max", "species Haemophilus influenzae",
  "sex", "min", "species Actinidia deliciosa"
)

build_custom_feature_plot(custom_feature_spec, suffix = "_custom_ranked")

custom_feature_spec_top3_filtered <- tibble::tribble(
  ~metadata, ~selection_side, ~feature_requested,
  "feeding", "max", "anno_food_item milk",
  "feeding", "max", "anno_food_item pig",
  "feeding", "max", "phylum Cossaviricota",
  "feeding", "min", "species Listeria innocua",
  "feeding", "min", "family Orthoherpesviridae",
  "feeding", "min", "species Lactiplantibacillus plantarum",
  "delmode", "max", "species Human coronavirus NL63",
  "delmode", "max", "genus Triticum",
  "delmode", "max", "genus Aegilops",
  "delmode", "min", "species Bacteroides heparinolyticus",
  "delmode", "min", "species Bacteroides sp. HF-5092",
  "delmode", "min", "species Simplexvirus humanalpha2",
  "siblings", "max", "species Streptococcus pneumoniae",
  "siblings", "max", "species Rhinovirus A",
  "siblings", "max", "species Bacteroides eggerthii",
  "siblings", "min", "species Staphylococcus aureus",
  "siblings", "min", "species Ovis aries",
  "siblings", "min", "order Galliformes",
  "smoking", "max", "species Bifidobacterium breve",
  "smoking", "max", "species Coprococcus comes",
  "smoking", "max", "species Orthopneumovirus hominis",
  "smoking", "min", "species Betainfluenzavirus influenzae",
  "smoking", "min", "species Poales",
  "smoking", "min", "genus Bacillus",
  "preterm", "max", "family Erysipelotrichaceae",
  "preterm", "max", "species Holdemanella biformis",
  "preterm", "max", "species Lactobacillus kitasatonis",
  "preterm", "min", "species Blautia wexlerae",
  "preterm", "min", "genus Vibrio",
  "preterm", "min", "order Compylobacterales",
  "delplace", "max", "species Simplexvirus humanalpha1",
  "delplace", "max", "genus Bordetella",
  "delplace", "max", "species Enterovirus C",
  "delplace", "min", "species Gelderlandvirus s16",
  "delplace", "min", "anno_viruses_bacteriophage mastadenovirusA",
  "delplace", "min", "species Cytomegalovirus humanbeta5",
  "sex", "max", "species Haemophilus influenzae",
  "sex", "min", "species Actinidia deliciosa"
)

build_custom_feature_plot(custom_feature_spec_top3_filtered, suffix = "_custom_ranked_top3_filtered")
