#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Rtsne)
  library(ggalluvial)
  library(plotly)
  library(htmlwidgets)
  library(Cairo)
})

timepoints <- c("T2", "T6", "T8")
z_thr <- 0
root_out_dir <- file.path("results", "sankey_meta")
dir.create(root_out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_specs <- list(
  siblings = c("siblings", "no_siblings"),
  delmode = c("delmode_VG", "delmode_CS"),
  delplace = c("delplace_home", "delplace_hospital"),
  feeding = c("feeding_BF", "feeding_MF"),
  preterm = c("preterm_yes", "preterm_no"),
  CDrisk = c("CDrisk_yes", "CDrisk_no"),
  lockdown = c("lockdown_before", "lockdown_after"),
  smoking = c("smoking_yes", "smoking_no"),
  sex = c("sex_male", "sex_female")
)

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

report_helpers_path <- file.path("results", "report_interesting_metadata", "report_helpers.R")
if (!file.exists(report_helpers_path)) {
  stop("Missing report helper file: ", report_helpers_path)
}

report_env <- new.env(parent = baseenv())
suppressWarnings(source(report_helpers_path, local = report_env))

if (exists("root_results_dir", envir = report_env, inherits = FALSE)) {
  base_dir_default <- get("root_results_dir", envir = report_env, inherits = FALSE)
} else {
  base_dir_default <- file.path("results", "results")
}

report_comparisons <- if (exists("comparisons", envir = report_env, inherits = FALSE)) {
  get("comparisons", envir = report_env, inherits = FALSE)
} else {
  list()
}

message("Default results directory: ", base_dir_default)
message("Root output directory: ", root_out_dir)
message("z_thr: ", z_thr)

resolve_cmp_base_dir <- function(cmp_id) {
  cmp_cfg <- report_comparisons[[cmp_id]]
  cmp_base <- base_dir_default
  if (!is.null(cmp_cfg) && !is.null(cmp_cfg$base_dir) && nzchar(as.character(cmp_cfg$base_dir))) {
    cmp_base <- as.character(cmp_cfg$base_dir)
  }
  cmp_base
}

load_delta_table <- function(cmp_id, tp_label) {
  cmp_base_dir <- resolve_cmp_base_dir(cmp_id)
  delta_path <- file.path(cmp_base_dir, cmp_id, "DELTA_framework", "delta_table.csv")

  if (!file.exists(delta_path)) {
    stop("Missing DELTA table for ", cmp_id, ": ", delta_path)
  }

  d <- readr::read_csv(delta_path, show_col_types = FALSE, name_repair = "unique_quiet")
  d <- d %>% dplyr::select(-dplyr::any_of("...1"))

  required_cols <- c("feature", "p_perm", "T_obs_stand")
  missing_cols <- setdiff(required_cols, names(d))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns for ", cmp_id, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!("rank" %in% names(d))) d$rank <- NA_character_
  if (!("n_peptides_used" %in% names(d))) d$n_peptides_used <- NA_real_

  d %>%
    dplyr::mutate(
      feature = as.character(.data$feature),
      rank = as.character(.data$rank),
      p_perm = as.numeric(.data$p_perm),
      n_peptides_used = as.numeric(.data$n_peptides_used),
      T_obs_stand = as.numeric(.data$T_obs_stand)
    ) %>%
    dplyr::group_by(.data$feature) %>%
    dplyr::slice_max(order_by = abs(.data$T_obs_stand), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(timepoint = tp_label)
}

classify_state <- function(z, thr) {
  dplyr::case_when(
    is.na(z) ~ "0",
    z > thr ~ "+",
    z < -thr ~ "-",
    TRUE ~ "0"
  )
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 1e-4, formatC(x, format = "e", digits = 2), formatC(x, format = "fg", digits = 3))
  )
}

fmt_z <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
}

fmt_n <- function(x) {
  ifelse(is.na(x), "NA", as.character(as.integer(round(x))))
}

run_one_metadata <- function(meta_name, meta_suffixes) {
  lhs <- meta_suffixes[[1]]
  rhs <- meta_suffixes[[2]]
  comparison_ids <- make_comparison_ids(lhs, rhs)
  out_dir <- file.path(root_out_dir, meta_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n=== Metadata: ", meta_name, " ===")
  message("Comparisons:")
  for (tp in timepoints) {
    message(" - ", tp, ": ", comparison_ids[[tp]])
  }

  delta_by_tp <- lapply(timepoints, function(tp) {
    load_delta_table(comparison_ids[[tp]], tp_label = tp)
  })
  names(delta_by_tp) <- timepoints

  sig_features <- sort(unique(unlist(lapply(delta_by_tp, function(d) {
    d[["feature"]][d[["p_perm"]] < 0.05]
  }))))

  if (length(sig_features) == 0L) {
    warning("No significant features for metadata ", meta_name, " (p_perm < 0.05). Skipping.")
    return(invisible(NULL))
  }

  message("Union significant features (p_perm < 0.05): ", length(sig_features))

  z_wide <- tibble::tibble(feature = sig_features)
  for (tp in timepoints) {
    z_wide <- z_wide %>%
      dplyr::left_join(
        delta_by_tp[[tp]] %>%
          dplyr::select(
            "feature",
            !!paste0("rank_", tp) := "rank",
            !!paste0("p_perm_", tp) := "p_perm",
            !!paste0("n_peptides_", tp) := "n_peptides_used",
            !!paste0("z_", tp) := "T_obs_stand"
          ),
        by = "feature"
      )
  }

  z_cols <- paste0("z_", timepoints)
  z_mat <- as.matrix(z_wide[, z_cols, drop = FALSE])
  z_mat[!is.finite(z_mat)] <- NA_real_
  z_mat_for_tsne <- z_mat
  z_mat_for_tsne[is.na(z_mat_for_tsne)] <- 0
  n_features <- nrow(z_mat_for_tsne)

  state_df <- z_wide %>%
    dplyr::mutate(
      state_T2 = classify_state(.data$z_T2, z_thr),
      state_T6 = classify_state(.data$z_T6, z_thr),
      state_T8 = classify_state(.data$z_T8, z_thr),
      trajectory = paste(.data$state_T2, .data$state_T6, .data$state_T8, sep = " -> "),
      class_triplet = paste0(.data$state_T2, .data$state_T6, .data$state_T8),
      primary_rank = dplyr::coalesce(.data$rank_T2, .data$rank_T6, .data$rank_T8)
    )

  traj_counts <- state_df %>%
    dplyr::count(.data$state_T2, .data$state_T6, .data$state_T8, .data$trajectory, name = "n_features") %>%
    dplyr::arrange(dplyr::desc(.data$n_features))

  target_classes <- c("+++", "---")

  summary_triplet_tbl <- state_df %>%
    dplyr::filter(.data$class_triplet %in% target_classes) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(z_cols),
      names_to = "timepoint",
      values_to = "T_obs_stand"
    ) %>%
    dplyr::mutate(timepoint = sub("^z_", "", .data$timepoint)) %>%
    dplyr::group_by(.data$class_triplet) %>%
    dplyr::summarise(
      n_features = dplyr::n_distinct(.data$feature),
      n_values = sum(is.finite(.data$T_obs_stand)),
      median_T_obs_stand = stats::median(.data$T_obs_stand, na.rm = TRUE),
      min_T_obs_stand = min(.data$T_obs_stand, na.rm = TRUE),
      max_T_obs_stand = max(.data$T_obs_stand, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(.data$class_triplet, target_classes))

  summary_triplet_by_tp_tbl <- state_df %>%
    dplyr::filter(.data$class_triplet %in% target_classes) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(z_cols),
      names_to = "timepoint",
      values_to = "T_obs_stand"
    ) %>%
    dplyr::mutate(timepoint = sub("^z_", "", .data$timepoint)) %>%
    dplyr::group_by(.data$class_triplet, .data$timepoint) %>%
    dplyr::summarise(
      n_features = dplyr::n_distinct(.data$feature),
      median_T_obs_stand = stats::median(.data$T_obs_stand, na.rm = TRUE),
      min_T_obs_stand = min(.data$T_obs_stand, na.rm = TRUE),
      max_T_obs_stand = max(.data$T_obs_stand, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(.data$class_triplet, target_classes), .data$timepoint)

  species_triplet_tbl <- state_df %>%
    dplyr::filter(.data$class_triplet %in% target_classes, .data$primary_rank == "species") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      median_z_T2_T8 = stats::median(c(.data$z_T2, .data$z_T6, .data$z_T8), na.rm = TRUE),
      abs_median_z_T2_T8 = abs(.data$median_z_T2_T8)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      "class_triplet",
      "feature",
      "median_z_T2_T8",
      "abs_median_z_T2_T8",
      "rank_T2",
      "rank_T6",
      "rank_T8",
      "p_perm_T2",
      "p_perm_T6",
      "p_perm_T8",
      "n_peptides_T2",
      "n_peptides_T6",
      "n_peptides_T8",
      "z_T2",
      "z_T6",
      "z_T8"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$abs_median_z_T2_T8), .data$class_triplet, .data$feature)

  readr::write_csv(state_df, file.path(out_dir, "feature_states_zthr0.csv"))
  readr::write_csv(traj_counts, file.path(out_dir, "trajectory_counts_zthr0.csv"))
  readr::write_csv(summary_triplet_tbl, file.path(out_dir, "summary_T_obs_stand_triplet_plusplus_minusminusminus.csv"))
  readr::write_csv(summary_triplet_by_tp_tbl, file.path(out_dir, "summary_T_obs_stand_triplet_plusplus_minusminusminus_by_timepoint.csv"))
  readr::write_csv(species_triplet_tbl, file.path(out_dir, "species_triplet_plusplus_minusminusminus.csv"))

  if (n_features >= 2) {
    max_perplexity <- floor((n_features - 1) / 3)
    if (max_perplexity >= 1) {
      perplexity <- min(30, max_perplexity)
      set.seed(632961)
      tsne_fit <- Rtsne::Rtsne(
        z_mat_for_tsne,
        dims = 2,
        perplexity = perplexity,
        check_duplicates = FALSE,
        pca = TRUE,
        verbose = FALSE
      )

      tsne_df <- state_df %>%
        dplyr::mutate(
          tsne1 = tsne_fit$Y[, 1],
          tsne2 = tsne_fit$Y[, 2],
          hover_text = paste0(
            "<b>", .data$feature, "</b>",
            "<br>rank: ", ifelse(is.na(.data$primary_rank), "NA", .data$primary_rank),
            "<br>trajectory: ", .data$trajectory,
            "<br><br><b>T2</b> | p_perm: ", fmt_p(.data$p_perm_T2),
            " | n_peptides: ", fmt_n(.data$n_peptides_T2),
            " | T_obs_stand: ", fmt_z(.data$z_T2),
            "<br><b>T6</b> | p_perm: ", fmt_p(.data$p_perm_T6),
            " | n_peptides: ", fmt_n(.data$n_peptides_T6),
            " | T_obs_stand: ", fmt_z(.data$z_T6),
            "<br><b>T8</b> | p_perm: ", fmt_p(.data$p_perm_T8),
            " | n_peptides: ", fmt_n(.data$n_peptides_T8),
            " | T_obs_stand: ", fmt_z(.data$z_T8)
          )
        )

      readr::write_csv(tsne_df, file.path(out_dir, "tsne_points_zthr0_all.csv"))

      tsne_df <- tsne_df %>%
        dplyr::filter(.data$state_T2 != "0", .data$state_T6 != "0", .data$state_T8 != "0")

      if (nrow(tsne_df) > 0L) {
        readr::write_csv(tsne_df, file.path(out_dir, "tsne_points_zthr0.csv"))

        traj_levels <- tsne_df %>%
          dplyr::count(.data$trajectory, sort = TRUE) %>%
          dplyr::pull(.data$trajectory)

        tsne_df <- tsne_df %>%
          dplyr::mutate(trajectory = factor(.data$trajectory, levels = traj_levels))

        p_tsne <- ggplot(tsne_df, aes(x = .data$tsne1, y = .data$tsne2, color = .data$trajectory)) +
          geom_point(size = 2.2, alpha = 0.85) +
          theme_bw(base_size = 12) +
          labs(
            title = paste0("t-SNE of significant ", meta_name, " features across T2/T6/T8 (no zero-state classes)"),
            subtitle = paste0(
              "Input: T_obs_stand (missing -> 0), union p_perm < 0.05; kept n=",
              nrow(tsne_df), " / ", n_features, ", perplexity=", perplexity
            ),
            x = "t-SNE 1",
            y = "t-SNE 2",
            color = "Trajectory"
          ) +
          theme(
            plot.title = element_text(face = "bold"),
            legend.position = "right"
          )

        Cairo::CairoSVG(
          file.path(out_dir, "tsne_significant_features_T_obs_stand.svg"),
          width = 30, height = 20, unit = "cm", bg = "white", dpi = 300
        )
        print(p_tsne)
        grDevices::dev.off()

        p_tsne_interactive <- plotly::plot_ly(
          data = tsne_df,
          x = ~tsne1,
          y = ~tsne2,
          type = "scatter",
          mode = "markers",
          color = ~trajectory,
          colors = grDevices::hcl.colors(max(3, length(traj_levels)), palette = "Dark 3"),
          text = ~hover_text,
          hoverinfo = "text",
          marker = list(size = 16, opacity = 0.85)
        ) %>%
          plotly::layout(
            title = list(text = paste0("Interactive t-SNE (no zero-state classes): ", meta_name, " (T_obs_stand)")),
            xaxis = list(title = "t-SNE 1"),
            yaxis = list(title = "t-SNE 2"),
            legend = list(title = list(text = "Trajectory"))
          )

        htmlwidgets::saveWidget(
          p_tsne_interactive,
          file = file.path(out_dir, "tsne_significant_features_T_obs_stand_interactive.html"),
          selfcontained = FALSE
        )
      } else {
        warning("No non-zero trajectories left for t-SNE in metadata ", meta_name)
      }
    } else {
      warning("Not enough features for t-SNE perplexity constraints in metadata ", meta_name)
    }
  } else {
    warning("Need at least 2 features for t-SNE in metadata ", meta_name)
  }

  flow_df <- traj_counts %>%
    dplyr::filter(.data$state_T2 != "0", .data$state_T6 != "0", .data$state_T8 != "0") %>%
    dplyr::mutate(
      state_T2 = factor(.data$state_T2, levels = c("+", "-")),
      state_T6 = factor(.data$state_T6, levels = c("+", "-")),
      state_T8 = factor(.data$state_T8, levels = c("+", "-"))
    )

  if (nrow(flow_df) > 0L) {
    lodes_df <- ggalluvial::to_lodes_form(
      data = flow_df,
      axes = c("state_T2", "state_T6", "state_T8"),
      key = "timepoint",
      value = "state",
      id = "trajectory"
    )

    p_alluvial <- ggplot(
      lodes_df,
      aes(
        x = .data$timepoint,
        stratum = .data$state,
        alluvium = .data$trajectory,
        y = .data$n_features
      )
    ) +
      ggalluvial::geom_flow(aes(fill = .data$state), alpha = 0.65) +
      ggalluvial::geom_stratum(
        aes(fill = .data$state),
        width = 1 / 7,
        color = "grey25"
      ) +
      ggplot2::geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 4
      ) +
      scale_x_discrete(
        limits = c("state_T2", "state_T6", "state_T8"),
        labels = c("T2", "T6", "T8"),
        expand = c(0.05, 0.05)
      ) +
      scale_fill_manual(
        values = c("+" = "#2E7D32", "-" = "#C62828"),
        breaks = c("+", "-"),
        drop = FALSE
      ) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0("Alluvial trajectories of feature states (", meta_name, ", no zero-state classes)"),
        subtitle = paste0(
          "States from T_obs_stand with z_thr=", z_thr,
          " using sign bands; kept n=", sum(flow_df$n_features), " / ", n_features, " features"
        ),
        x = "Timepoint",
        y = "Number of features",
        fill = "State"
      ) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )

    Cairo::CairoSVG(
      file.path(out_dir, "alluvial_feature_trajectories_zthr0.svg"),
      width = 30, height = 20, unit = "cm", bg = "white", dpi = 300
    )
    print(p_alluvial)
    grDevices::dev.off()
  } else {
    warning("No non-zero trajectories left for Sankey in metadata ", meta_name)
  }

  message("Saved outputs in: ", out_dir)
  message(" - feature_states_zthr0.csv")
  message(" - trajectory_counts_zthr0.csv")
  message(" - summary_T_obs_stand_triplet_plusplus_minusminusminus.csv")
  message(" - summary_T_obs_stand_triplet_plusplus_minusminusminus_by_timepoint.csv")
  message(" - species_triplet_plusplus_minusminusminus.csv")
  if (file.exists(file.path(out_dir, "tsne_points_zthr0_all.csv"))) {
    message(" - tsne_points_zthr0_all.csv")
  }
  if (file.exists(file.path(out_dir, "tsne_points_zthr0.csv"))) {
    message(" - tsne_points_zthr0.csv")
  }
  if (file.exists(file.path(out_dir, "tsne_significant_features_T_obs_stand.svg"))) {
    message(" - tsne_significant_features_T_obs_stand.svg")
  }
  if (file.exists(file.path(out_dir, "tsne_significant_features_T_obs_stand_interactive.html"))) {
    message(" - tsne_significant_features_T_obs_stand_interactive.html")
  }
  if (file.exists(file.path(out_dir, "alluvial_feature_trajectories_zthr0.svg"))) {
    message(" - alluvial_feature_trajectories_zthr0.svg")
  }

  invisible(NULL)
}

for (meta_name in names(metadata_specs)) {
  meta_suffixes <- metadata_specs[[meta_name]]
  tryCatch(
    run_one_metadata(meta_name, meta_suffixes),
    error = function(e) {
      warning("Metadata ", meta_name, " failed: ", conditionMessage(e))
      NULL
    }
  )
}

message("\nAll metadata runs finished. Output root: ", root_out_dir)
