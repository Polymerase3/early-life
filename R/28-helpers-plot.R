# R/helpers_plot.R

library(plotly)

# comments in English

build_plotly_3d <- function(df, color_var,
                            hover_cols = c("subject_id","big_group","timepoint_factor","dyade"),
                            size = 5, opacity = 0.9, colors = NULL, title = NULL) {

  stopifnot(all(c("tSNE1","tSNE2","tSNE3") %in% names(df)))
  if (!color_var %in% names(df)) stop("color_var not found in df")

  hc <- intersect(hover_cols, names(df))
  # Build a sprintf format like "a: %s<br>b: %s<br>x: %s<br>y: %s<br>z: %s"
  fmt_parts <- c(sprintf("%s: %%s", hc), "x: %s", "y: %s", "z: %s")
  fmt <- paste(fmt_parts, collapse = "<br>")
  args <- c(
    list(fmt),
    lapply(hc, function(v) df[[v]]),
    list(signif(df$tSNE1,4), signif(df$tSNE2,4), signif(df$tSNE3,4))
  )
  hover_text <- do.call(sprintf, args)

  plotly::plot_ly(
    df,
    x = ~tSNE1, y = ~tSNE2, z = ~tSNE3,
    type = "scatter3d", mode = "markers",
    color = ~.data[[color_var]], colors = colors,
    text = hover_text, hoverinfo = "text",
    marker = list(size = size, opacity = opacity)
  ) |>
    plotly::layout(
      title = title %||% paste("t-SNE 3D by", color_var),
      scene = list(
        xaxis = list(title = "t-SNE 1"),
        yaxis = list(title = "t-SNE 2"),
        zaxis = list(title = "t-SNE 3"),
        aspectmode = "data"
      )
    ) |>
    plotly::toWebGL() |>
    plotly::config(displaylogo = FALSE) |>
    plotly::layout(uirevision = "keep-camera")
}


build_plotly_2d <- function(df, color_var, size = 4, opacity = 0.9, colors = NULL, title = NULL) {
  stopifnot(all(c("tSNE1","tSNE2") %in% names(df)))
  if (!color_var %in% names(df)) stop("color_var not found in df")

  plot_ly(
    df, x = ~tSNE1, y = ~tSNE2, type = "scatter", mode = "markers",
    color = ~.data[[color_var]], colors = colors,
    marker = list(size = size, opacity = opacity),
    text = ~paste(
      "x:", signif(tSNE1, 4), "<br>",
      "y:", signif(tSNE2, 4)
    ),
    hoverinfo = "text"
  ) |>
    layout(
      title = title %||% paste("t-SNE 2D by", color_var),
      xaxis = list(title = "t-SNE 1"),
      yaxis = list(title = "t-SNE 2")
    ) |>
    toWebGL() |>
    config(displaylogo = FALSE) |>
    layout(uirevision = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(htmlwidgets)
  library(plotly)
  library(readr)
})

# Build + save a single Plotly heatmap HTML using your precomputed DuckDB table
export_stability_heatmap_html <- function(
    sim_tbl,
    groups,                  # c(g1, g2)
    times,                   # c(t1, t2)
    pairing = c("cross","same"),
    dyads   = c("order","only","no"),
    grid    = c("intersection","union_pad","require_full","observed_only"),
    dedupe  = c("mean","first","median","max","min"),
    z_limits = c(0, 1),      # Kulczynski in [0,1]
    outdir   = "results/heatmaps_html",
    file_stub = NULL,        # optional prefix (e.g., "heat-01-")
    selfcontained = TRUE,    # TRUE = one portable HTML file
    alpha_tbl = NULL         # OPTIONAL: subject_id, group, time, richness, shannon
) {
  pairing <- match.arg(pairing)
  dyads   <- match.arg(dyads)
  grid    <- match.arg(grid)
  dedupe  <- match.arg(dedupe)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # 1) Use your ggplot builder to get axis order + data (no need to render it)
  p_gp <- plot_similarity_heatmap(
    sim_tbl   = sim_tbl,
    groups    = groups,
    times     = times,
    pairing   = pairing,
    dyads     = dyads,
    grid      = grid,
    dedupe    = dedupe,
    limits    = z_limits,
    na_fill   = "#F0F0F0",
    verbose   = FALSE
  )

  # Extract data used for the tiles and the final axis levels/order
  hm_df <- attr(p_gp, "heatmap_df")
  axes  <- attr(p_gp, "axis_levels")
  stopifnot(!is.null(hm_df), !is.null(axes))

  # 2) Wide matrix in the exact on-plot order
  wide <- tidyr::pivot_wider(
    hm_df,
    id_cols   = subject_left,
    names_from = subject_right,
    values_from = similarity
  )

  # Enforce order: rows (top->bottom) = reverse of axes$rows from ggplot,
  # columns left->right = axes$cols
  row_levels_plot <- rev(axes$rows)
  col_levels_plot <- axes$cols

  # keep only columns we have
  keep_cols <- intersect(col_levels_plot, setdiff(names(wide), "subject_left"))

  wide <- wide %>%
    dplyr::mutate(subject_left = as.character(subject_left)) %>%
    dplyr::arrange(factor(subject_left, levels = row_levels_plot)) %>%
    dplyr::select(subject_left, dplyr::all_of(keep_cols))

  # 3) Build Z matrix and IDs
  Z <- as.matrix(wide[, -1, drop = FALSE])  # numeric matrix
  row_ids <- wide$subject_left
  col_ids <- colnames(Z)

  # Labels
  gl <- groups[1]; gr <- groups[2]
  tl <- times[1];  tr <- times[2]

  # 3a) Optional alpha-div lookups (must be already normalized: group/time as in stab_df_full)
  # Expected columns: subject_id, group, time, richness, shannon
  # Prepare default (all NA) and fill where available
  lr <- ls <- rr <- rs <- NULL
  fmt0 <- function(x) ifelse(is.na(x), "NA", sprintf("%.0f", x))
  fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))

  if (!is.null(alpha_tbl)) {
    need_cols <- c("subject_id","group","time","richness","shannon")
    stopifnot(all(need_cols %in% colnames(alpha_tbl)))

    left_info <- alpha_tbl %>%
      dplyr::filter(.data$group == gl, .data$time == tl) %>%
      dplyr::transmute(subject_id, rich = richness, shan = shannon)

    right_info <- alpha_tbl %>%
      dplyr::filter(.data$group == gr, .data$time == tr) %>%
      dplyr::transmute(subject_id, rich = richness, shan = shannon)

    # maps by subject_id
    l_rich_map <- stats::setNames(left_info$rich,  left_info$subject_id)
    l_shan_map <- stats::setNames(left_info$shan,  left_info$subject_id)
    r_rich_map <- stats::setNames(right_info$rich, right_info$subject_id)
    r_shan_map <- stats::setNames(right_info$shan, right_info$subject_id)

    lr <- unname(l_rich_map[row_ids]); ls <- unname(l_shan_map[row_ids])
    rr <- unname(r_rich_map[col_ids]); rs <- unname(r_shan_map[col_ids])
  } else {
    lr <- ls <- rep(NA_real_, length(row_ids))
    rr <- rs <- rep(NA_real_, length(col_ids))
  }

  # 3b) NA color bin: map NA -> -0.001 and add a colorscale stop for gray
  Z_plot <- Z
  Z_plot_na_mask <- is.na(Z_plot)
  Z_plot[Z_plot_na_mask] <- -0.001

  # 3c) Hovertext (with alpha-div metrics)
  # replicate vectors across the grid in the correct orientation
  lr_rep <- rep(fmt0(lr), times = length(col_ids))
  ls_rep <- rep(fmt2(ls), times = length(col_ids))
  rr_rep <- rep(fmt0(rr), each  = length(row_ids))
  rs_rep <- rep(fmt2(rs), each  = length(row_ids))

  hover_mat <- matrix(
    paste0(
      "Left subject: ", rep(row_ids, times = length(col_ids)), "<br>",
      "Right subject: ", rep(col_ids, each = length(row_ids)), "<br>",
      "Left: ", gl, " @ ", tl,
      " | richness=", lr_rep, ", shannon=", ls_rep, "<br>",
      "Right: ", gr, " @ ", tr,
      " | richness=", rr_rep, ", shannon=", rs_rep, "<br>",
      "Kulczynski: ", sprintf("%.3f", as.numeric(Z))
    ),
    nrow = nrow(Z), ncol = ncol(Z), byrow = FALSE
  )
  hover_mat[Z_plot_na_mask] <- matrix(
    paste0(
      "Left subject: ", rep(row_ids, times = length(col_ids)), "<br>",
      "Right subject: ", rep(col_ids, each = length(row_ids)), "<br>",
      "Left: ", gl, " @ ", tl,
      " | richness=", lr_rep, ", shannon=", ls_rep, "<br>",
      "Right: ", gr, " @ ", tr,
      " | richness=", rr_rep, ", shannon=", rs_rep, "<br>",
      "Kulczynski: NA"
    ),
    nrow = nrow(Z), ncol = ncol(Z), byrow = FALSE
  )

  # 4) Tick thinning every 5 labels (work with label strings, not indices)
  pick_ticks_labels <- function(ids) {
    n <- length(ids)
    idx <- unique(c(1, seq(5, n, by = 5)))
    list(vals = ids[idx], text = c("1", 5 * seq_len(length(idx) - 1)))
  }
  tx <- pick_ticks_labels(col_ids)
  ty <- pick_ticks_labels(row_ids)

  # 5) Plotly heatmap — Inferno **reversed**; first stop reserved for NA gray
  inferno_scale <- list(
    list(0.0000, "#F0F0F0"),  # NA bin (mapped to -0.001)
    list(0.0005, "#fcffa4"),
    list(0.1250, "#f9d44a"),
    list(0.2500, "#f98e09"),
    list(0.3750, "#e35933"),
    list(0.5000, "#b83655"),
    list(0.6250, "#88226a"),
    list(0.7500, "#550f6d"),
    list(0.8750, "#1f0c48"),
    list(1.0000, "#000004")
  )

  p <- plot_ly(
    type = "heatmap",
    z = Z_plot,
    x = col_ids,
    y = row_ids,
    zauto = FALSE, zmin = -0.001, zmax = 1.0,
    colorscale = inferno_scale,
    hoverinfo = "text",
    text = hover_mat,
    showscale = TRUE
  ) |>
    colorbar(title = "Kulczynski similarity") |>
    layout(
      xaxis = list(
        title = paste0(gr, " @ ", tr),
        type = "category",
        categoryorder = "array", categoryarray = col_ids,
        tickmode = "array", tickvals = tx$vals, ticktext = tx$text
      ),
      yaxis = list(
        title = paste0(gl, " @ ", tl),
        type = "category",
        categoryorder = "array", categoryarray = row_ids,
        tickmode = "array", tickvals = ty$vals, ticktext = ty$text,
        autorange = "reversed"
      ),
      margin = list(l = 60, r = 20, t = 50, b = 60),
      title = list(text = paste0(gl, " @ ", tl, "  vs  ", gr, " @ ", tr))
    ) |>
    config(displaylogo = FALSE)

  # 6) Filenames
  safe <- function(s) {
    s <- gsub("[^A-Za-z0-9._-]+", "_", s)
    s <- gsub("_+", "_", s)
    tolower(trimws(s, "both", "_"))
  }
  base <- paste0(
    safe(paste0(gl, "-", tl)), "__", safe(paste0(gr, "-", tr)),
    "_kulczynski_plotly"
  )
  if (!is.null(file_stub) && nchar(file_stub)) base <- paste0(file_stub, base)
  file_html <- file.path(outdir, paste0(base, ".html"))
  file_csv  <- file.path(outdir, paste0(base, ".csv"))

  # 7) Save CSV in visual order (row_ids top->bottom, col_ids left->right)
  out_csv <- cbind(subject_left = row_ids, as.data.frame(Z))
  names(out_csv) <- c("subject_left", col_ids)
  readr::write_csv(out_csv, file_csv, na = "")

  # 8) Save HTML (plotly)
  htmlwidgets::saveWidget(p, file_html, selfcontained = selfcontained)

  list(html = file_html, csv = file_csv)
}

