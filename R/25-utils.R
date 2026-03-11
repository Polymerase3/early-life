long_to_wide_exist <- function(
    x,
    row_cols   = c("subject_id","big_group","timepoint_factor"),
    col_col    = "peptide_id",
    value_col  = "exist",
    sep        = "__",
    sparse     = TRUE
) {
  # if it's a dbplyr/duckdb lazy table, pull only the needed columns first
  if (inherits(x, "tbl_sql")) {
    x <- dplyr::select(x, dplyr::all_of(c(row_cols, col_col, value_col))) %>%
      dplyr::collect()
  }

  if (rlang::is_installed("data.table")) {
    DT <- data.table::as.data.table(x)

    # build row_id once; keep only needed cols
    data.table::set(DT, j = "row_id",
                    value = do.call(paste, c(DT[, ..row_cols], list(sep = sep))))
    DT <- DT[, c("row_id", col_col, value_col), with = FALSE]

    # force integer 0/1
    data.table::set(DT, j = value_col,
                    value = as.integer(DT[[value_col]] > 0))

    if (isTRUE(sparse)) {
      # sparse path (recommended for many peptides)
      # duplicates are *summed* by sparseMatrix; threshold back to 0/1
      row_f <- factor(DT$row_id)
      col_f <- factor(DT[[col_col]])
      M <- Matrix::sparseMatrix(
        i = as.integer(row_f),
        j = as.integer(col_f),
        x = DT[[value_col]],
        dims = c(nlevels(row_f), nlevels(col_f)),
        dimnames = list(levels(row_f), levels(col_f))
      )
      return((M > 0) + 0L) # dgCMatrix with 0/1
    } else {
      # dense matrix via dcast (watch memory!)
      form <- reformulate(col_col, response = "row_id")
      wide_dt <- data.table::dcast(
        DT, formula = form,
        value.var = value_col,
        fun.aggregate = max,  # collapse duplicates to 0/1
        fill = 0L
      )
      rn <- wide_dt[["row_id"]]
      mat <- as.matrix(wide_dt[, setdiff(names(wide_dt), "row_id"), with = FALSE])
      rownames(mat) <- rn
      storage.mode(mat) <- "integer"
      return(mat)
    }
  } else {
    # tidyr fallback (will be slower and memory-heavier)
    wide <- x |>
      dplyr::mutate(
        row_id = do.call(paste, c(dplyr::pick(dplyr::all_of(row_cols)), list(sep = sep))),
        !!value_col := as.integer(.data[[value_col]] > 0)
      ) |>
      dplyr::select(row_id, dplyr::all_of(c(col_col, value_col))) |>
      dplyr::group_by(row_id, .data[[col_col]]) |>
      dplyr::summarise(!!value_col := max(.data[[value_col]]), .groups = "drop") |>
      tidyr::pivot_wider(
        id_cols     = row_id,
        names_from  = dplyr::all_of(col_col),
        values_from = dplyr::all_of(value_col),
        values_fill = 0L
      )

    rn <- wide$row_id
    mat <- as.matrix(wide[, setdiff(names(wide), "row_id"), drop = FALSE])
    rownames(mat) <- rn
    storage.mode(mat) <- "integer"
    if (isTRUE(sparse)) {
      return(Matrix::Matrix(mat, sparse = TRUE))
    } else {
      return(mat)
    }
  }
}

kulczynski_dist <- function(M) {
  # ensure sparse binary
  M <- drop0(M)
  M <- as(M != 0, "dgCMatrix")

  # intersections (dense n x n; fine for n = nrow(M))
  A   <- as.matrix(tcrossprod(M))     # a_ij
  r   <- Matrix::rowSums(M)           # r_i
  r_s <- ifelse(r == 0, 1, r)         # safe denom to avoid /0

  # 0.5 * (a/r_i + a/r_j)
  K <- 0.5 * (sweep(A, 1, r_s, "/") + sweep(A, 2, r_s, "/"))

  # distances
  D <- 1 - K
  diag(D) <- 0
  D
}

# ==============================================================================
# Repertoire stability: compute (Jaccard / Sorensen / Kulczynski / Bray-Curtis / Phi) --> tibble
# ==============================================================================

#' @title Compute repertoire stability over time (Jaccard, Sorensen, Kulczynski, Bray-Curtis, Phi)
#'
#' @description For each **subject**, within each **group**, compute the
#' similarity between the set of *present* peptides at each time point and that
#' subject’s **baseline** (their earliest time point). Supported:
#' - **"jaccard"**    : \eqn{|A ∩ B| / |A ∪ B|}
#' - **"sorensen"**   : \eqn{2|A ∩ B| / (|A| + |B|)}  (Dice)
#' - **"kulczynski"** : \eqn{0.5 * [ |A∩B|/(|A∩B|+|A\B|) + |A∩B|/(|A∩B|+|B\A|) ]}
#' - **"bray_curtis"** (binary here): identical to Sørensen/Dice
#' - **"phi"** / **"pearson"** : phi (Pearson corr.) for two binarne wektory; zakres [-1, 1]
#'
#' Baseline self-comparisons (time == baseline) equal 1 by definition for set-based
#' similarities; if `drop_self = TRUE` these rows are removed, a marker is kept.
#'
#' @param x A `<phip_data>` object. **Must be longitudinal**.
#' @param group_col Character; grouping variable in `x$data_long`.
#' @param time_col Character; numeric/Date/POSIX time variable.
#' @param similarity One of: "jaccard","sorensen","kulczynski","bray_curtis","phi","pearson".
#' @param drop_self Logical; drop baseline rows (default `TRUE`).
#' @param auto_expand Logical; if `TRUE`, require n0>=1 & nt>=1.
#'
#' @return tibble: `subject_id`, `group`, `time`, `similarity` (+ attributes).
#' @export
compute_repertoire_stability <- function(
    x,
    group_col,
    time_col,
    similarity  = c("jaccard", "sorensen", "kulczynski", "bray_curtis", "phi", "pearson"),
    drop_self   = TRUE,
    auto_expand = FALSE,  # require n0>=1 & nt>=1 when TRUE
    full_cross  = FALSE   # NEW: all pairwise (group,time) combos per subject
) {
  stopifnot(inherits(x, "phip_data"))
  similarity <- match.arg(similarity)

  group_col_local <- if (is.null(group_col)) ".phip__group" else group_col

  .ph_with_timing(
    headline = "Computing repertoire similarity (<phip_data>)",
    step     = sprintf("similarity: %s; mode: %s",
                       similarity, if (full_cross) "all-pairs" else "baseline"),
    expr = {
      if (!isTRUE(x$meta$longitudinal)) {
        .ph_abort(
          headline = "Longitudinal data required.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      # --- wymogi kolumn ------------------------------------------------------
      need <- c("subject_id", "sample_id", "peptide_id", "exist", time_col)
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(x$data_long))
      .chk_cond(length(miss) > 0,
                sprintf("Missing required columns in data_long: %s",
                        paste(miss, collapse = ", ")))

      .data <- rlang::.data
      tcol  <- rlang::sym(time_col)
      grp   <- rlang::sym(group_col_local)

      # --- presence-only widok (exist==1) + syntetyczna grupa "all" ----------
      present <- if (is.null(group_col)) {
        x$data_long |>
          dplyr::filter(.data$exist == 1L) |>
          dplyr::select(.data$subject_id, .data$sample_id, .data$peptide_id, !!tcol) |>
          dplyr::filter(!is.na(!!tcol)) |>
          dplyr::mutate(!!grp := "all")
      } else {
        x$data_long |>
          dplyr::filter(.data$exist == 1L) |>
          dplyr::select(.data$subject_id, .data$sample_id, .data$peptide_id, !!grp, !!tcol) |>
          dplyr::filter(!is.na(!!tcol))
      }

      # =======================================================================
      # BRANCH A: full_cross = TRUE  --> wszystkie pary (group,time)
      # =======================================================================
      if (isTRUE(full_cross)) {

        # zbiory: unikatowe peptydy per (subject, group, time)
        sets <- present |>
          dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id, .data$sample_id)

        # rozmiary zbiorów
        sizes <- sets |>
          dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n")

        # self-join po peptydzie w obrębie osoby; filtrujemy do par unikalnych (key1 < key2)
        sets_l <- sets |>
          dplyr::rename(group1 = !!grp, time1 = !!tcol,
                        sample_id1 = .data$sample_id, peptide_id1 = .data$peptide_id)
        sets_r <- sets |>
          dplyr::rename(group2 = !!grp, time2 = !!tcol,
                        sample_id2 = .data$sample_id, peptide_id2 = .data$peptide_id)

        inters <- sets_l |>
          dplyr::inner_join(
            sets_r,
            by = c("subject_id", "peptide_id1" = "peptide_id2")
          ) |>
          dplyr::filter(.data$group1 != .data$group2 | .data$time1 != .data$time2) |>
          dplyr::mutate(
            key1 = paste0(.data$group1, "||", .data$time1),
            key2 = paste0(.data$group2, "||", .data$time2)
          ) |>
          dplyr::filter(.data$key1 < .data$key2) |>
          dplyr::count(.data$subject_id, .data$group1, .data$time1,
                       .data$group2, .data$time2, name = "a")

        # dołącz rozmiary
        pairs <- inters |>
          dplyr::left_join(
            sizes |> dplyr::rename(group1 = !!grp, time1 = !!tcol, n1 = .data$n),
            by = c("subject_id","group1","time1")
          ) |>
          dplyr::left_join(
            sizes |> dplyr::rename(group2 = !!grp, time2 = !!tcol, n2 = .data$n),
            by = c("subject_id","group2","time2")
          ) |>
          dplyr::mutate(
            a = dplyr::coalesce(.data$a, 0L),
            n1 = dplyr::coalesce(.data$n1, 0L),
            n2 = dplyr::coalesce(.data$n2, 0L),
            b = pmax(.data$n1 - .data$a, 0L),
            c = pmax(.data$n2 - .data$a, 0L)
          )

        # phi potrzebuje N – policz rozmiary przestrzeni cech per próbka
        if (identical(similarity, "phi") || identical(similarity, "pearson")) {
          timesamples <- x$data_long |>
            dplyr::filter(!is.na(!!tcol)) |>
            dplyr::select(.data$subject_id, !!grp, !!tcol, .data$sample_id) |>
            dplyr::distinct()

          sample_sizes <- x$data_long |>
            dplyr::group_by(.data$sample_id) |>
            dplyr::summarise(N = dplyr::n_distinct(.data$peptide_id), .groups = "drop")

          pairs <- pairs |>
            # próbki odpowiadające parze (group1,time1) / (group2,time2)
            dplyr::left_join(
              timesamples |>
                dplyr::rename(group1 = !!grp, time1 = !!tcol, sample_id1 = .data$sample_id),
              by = c("subject_id","group1","time1")
            ) |>
            dplyr::left_join(
              timesamples |>
                dplyr::rename(group2 = !!grp, time2 = !!tcol, sample_id2 = .data$sample_id),
              by = c("subject_id","group2","time2")
            ) |>
            dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
            dplyr::rename(N1 = .data$N) |>
            dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
            dplyr::rename(N2 = .data$N) |>
            dplyr::mutate(
              N   = pmin(.data$N1, .data$N2, na.rm = TRUE),
              n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
              n1d = .data$a + .data$b,
              nd1 = .data$a + .data$c,
              n0d = .data$N - .data$n1d,
              nd0 = .data$N - .data$nd1,
              den = sqrt(pmax(.data$n1d,0) * pmax(.data$n0d,0) *
                           pmax(.data$nd1,0) * pmax(.data$nd0,0)),
              similarity = dplyr::if_else(.data$den > 0,
                                          (.data$a * .data$n00 - .data$b * .data$c) / .data$den,
                                          NA_real_)
            ) |>
            dplyr::select(-n00,-n1d,-nd1,-n0d,-nd0,-den,-N1,-N2)
          if (!isTRUE(x$meta$full_cross)) {
            .ph_warn("Phi/Pearson requested but meta$full_cross is FALSE; N may be underestimated.")
          }
        } else if (identical(similarity, "jaccard")) {
          pairs <- pairs |>
            dplyr::mutate(
              similarity = dplyr::if_else((.data$n1 + .data$n2 - .data$a) > 0,
                                          .data$a / (.data$n1 + .data$n2 - .data$a),
                                          NA_real_)
            )
        } else if (identical(similarity, "sorensen") || identical(similarity, "bray_curtis")) {
          pairs <- pairs |>
            dplyr::mutate(
              similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                          (2 * .data$a) / (.data$n1 + .data$n2),
                                          NA_real_)
            )
        } else if (identical(similarity, "kulczynski")) {
          pairs <- pairs |>
            dplyr::mutate(
              term1 = dplyr::if_else((.data$a + .data$b) > 0, .data$a / (.data$a + .data$b), NA_real_),
              term2 = dplyr::if_else((.data$a + .data$c) > 0, .data$a / (.data$a + .data$c), NA_real_),
              similarity = dplyr::coalesce(0.5 * (.data$term1 + .data$term2),
                                           dplyr::if_else((.data$a + .data$b + .data$c) == 0, 1.0, 0.0))
            ) |>
            dplyr::select(-term1, -term2)
        }

        out <- pairs |>
          dplyr::select(.data$subject_id, .data$group1, .data$time1, .data$group2, .data$time2, .data$similarity) |>
          dplyr::collect()
        class(out) <- c("phip_repertoire_pairs", class(out))
        return(out)
      }

      # =======================================================================
      # BRANCH B: baseline (oryginalny tryb)
      # =======================================================================

      # ---- baseline (per subject within group) & sets ------------------------
      baseline_tbl <- present |>
        dplyr::group_by(.data$subject_id, !!grp) |>
        dplyr::summarise(baseline_time = min(!!tcol, na.rm = TRUE), .groups = "drop")

      base_set <- present |>
        dplyr::inner_join(baseline_tbl, by = c("subject_id", group_col_local)) |>
        dplyr::filter(!!tcol == .data$baseline_time) |>
        dplyr::distinct(.data$subject_id, !!grp, .data$peptide_id)

      tp_set <- present |>
        dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id)

      bsize <- base_set |>
        dplyr::count(.data$subject_id, !!grp, name = "n0")
      tsize <- tp_set |>
        dplyr::count(.data$subject_id, !!grp, !!tcol, name = "nt")
      ints  <- tp_set |>
        dplyr::inner_join(base_set,
                          by = c("subject_id", group_col_local, "peptide_id")) |>
        dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n_int")

      stab <- tsize |>
        dplyr::left_join(ints,  by = c("subject_id", group_col_local, time_col)) |>
        dplyr::left_join(bsize, by = c("subject_id", group_col_local)) |>
        dplyr::mutate(
          n_int = dplyr::coalesce(.data$n_int, 0L),
          nt    = dplyr::coalesce(.data$nt,    0L),
          n0    = dplyr::coalesce(.data$n0,    0L)
        ) |>
        dplyr::left_join(baseline_tbl, by = c("subject_id", group_col_local))

      if (isTRUE(auto_expand)) {
        n_before <- nrow(stab)
        stab <- stab |> dplyr::filter(.data$n0 >= 1L, .data$nt >= 1L)
        .ph_log_info("Auto-expand filter (require n0 >= 1 & nt >= 1)",
                     bullets = sprintf("removed rows: %d", n_before - nrow(stab)))
      }

      # ---- similarity (baseline) --------------------------------------------
      if (identical(similarity, "jaccard")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax(.data$n_int / pmax(.data$n0 + .data$nt - .data$n_int, 1L), 0), 1
          ))

      } else if (identical(similarity, "sorensen")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax((2 * .data$n_int) / pmax(.data$n0 + .data$nt, 1L), 0), 1
          ))

      } else if (identical(similarity, "bray_curtis")) {
        stab <- stab |>
          dplyr::mutate(similarity = pmin(
            pmax((2 * .data$n_int) / pmax(.data$n0 + .data$nt, 1L), 0), 1
          ))

      } else if (identical(similarity, "kulczynski")) {
        stab <- stab |>
          dplyr::mutate(
            a = .data$n_int,
            b = pmax(.data$n0 - .data$n_int, 0L),
            c = pmax(.data$nt - .data$n_int, 0L),
            term1 = dplyr::if_else((a + b) > 0, a / (a + b), NA_real_),
            term2 = dplyr::if_else((a + c) > 0, a / (a + c), NA_real_),
            similarity = dplyr::coalesce(0.5 * (term1 + term2),
                                         dplyr::if_else((a + b + c) == 0, 1.0, 0.0))
          ) |>
          dplyr::select(-a, -b, -c, -term1, -term2)

      } else { # phi / pearson
        timesamples <- x$data_long |>
          dplyr::filter(!is.na(!!tcol)) |>
          dplyr::select(.data$subject_id, !!grp, !!tcol, .data$sample_id) |>
          dplyr::distinct()

        baseline_sids <- timesamples |>
          dplyr::inner_join(
            baseline_tbl,
            by = c("subject_id", group_col_local, time_col = "baseline_time")
          ) |>
          dplyr::rename(baseline_sample_id = .data$sample_id) |>
          dplyr::select(.data$subject_id, !!grp, .data$baseline_sample_id)

        pair_sids <- timesamples |>
          dplyr::inner_join(baseline_sids, by = c("subject_id", group_col_local)) |>
          dplyr::rename(tp_sample_id = .data$sample_id)

        sample_sizes <- x$data_long |>
          dplyr::group_by(.data$sample_id) |>
          dplyr::summarise(N = dplyr::n_distinct(.data$peptide_id), .groups = "drop")

        N_map <- pair_sids |>
          dplyr::left_join(sample_sizes, by = c("tp_sample_id" = "sample_id")) |>
          dplyr::rename(N_tp = .data$N) |>
          dplyr::left_join(sample_sizes, by = c("baseline_sample_id" = "sample_id")) |>
          dplyr::rename(N_base = .data$N) |>
          dplyr::mutate(N = pmin(.data$N_tp, .data$N_base, na.rm = TRUE)) |>
          dplyr::select(.data$subject_id, !!grp, !!tcol, .data$N)

        stab <- stab |>
          dplyr::left_join(N_map, by = c("subject_id", group_col_local, time_col)) |>
          dplyr::mutate(
            a    = .data$n_int,
            b    = pmax(.data$n0 - .data$n_int, 0L),
            c    = pmax(.data$nt - .data$n_int, 0L),
            n00  = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
            n1d  = .data$a + .data$b,
            nd1  = .data$a + .data$c,
            n0d  = .data$N - .data$n1d,
            nd0  = .data$N - .data$nd1,
            den  = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) *
                          pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
            similarity = dplyr::if_else(.data$den > 0,
                                        (.data$a * .data$n00 - .data$b * .data$c) / .data$den,
                                        NA_real_)
          ) |>
          dplyr::select(-a, -b, -c, -n00, -n1d, -nd1, -n0d, -nd0, -den)

        if (!isTRUE(x$meta$full_cross)) {
          .ph_warn("Phi/Pearson requested but meta$full_cross is FALSE; N may be underestimated. Similarity may be NA for some rows.")
        }
        if (isTRUE(any(is.na(stab$N)))) {
          n_na <- sum(is.na(stab$N))
          .ph_warn(sprintf("Phi/Pearson: missing N for %d row(s); similarity set to NA.", n_na))
        }
      }

      # --- finalize (baseline mode) ------------------------------------------
      out <- stab |>
        dplyr::select(.data$subject_id, !!grp, !!tcol, .data$similarity, .data$baseline_time) |>
        dplyr::collect() |>
        dplyr::rename(group = !!grp, time = !!tcol)

      base_per_group <- out |>
        dplyr::group_by(.data$group) |>
        dplyr::summarise(baseline_time = min(.data$baseline_time, na.rm = TRUE), .groups = "drop")

      if (isTRUE(drop_self)) {
        n_before <- nrow(out)
        out <- dplyr::filter(out, .data$time != .data$baseline_time)
        .ph_log_info("Dropping baseline self-comparisons",
                     bullets = sprintf("removed rows: %d", n_before - nrow(out)))
      }

      # baseline marker heights
      marker_tbl <- NULL
      if (nrow(out)) {
        marker_tbl <- out |>
          dplyr::select(-baseline_time) |>
          dplyr::inner_join(base_per_group, by = "group") |>
          dplyr::group_by(.data$group, .data$baseline_time) |>
          dplyr::summarise(
            next_time = {
              v <- .data$time[.data$time > .data$baseline_time]
              if (length(v)) min(v, na.rm = TRUE) else NA_real_
            },
            .groups = "drop"
          ) |>
          dplyr::left_join(
            out |> dplyr::group_by(.data$group, .data$time) |>
              dplyr::summarise(y_marker = mean(.data$similarity, na.rm = TRUE), .groups = "drop"),
            by = c("group", "next_time" = "time")
          ) |>
          dplyr::mutate(y_marker = dplyr::coalesce(.data$y_marker, 1)) |>
          dplyr::rename(x_marker = baseline_time)
      }

      class(out) <- c("phip_repertoire_stability", class(out))
      attr(out, "baseline_times")  <- base_per_group
      attr(out, "baseline_marker") <- marker_tbl
      attr(out, "drop_self")       <- isTRUE(drop_self)
      out
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# Repertoire stability: plot (takes tibble from compute_*)
# ==============================================================================

#' Plot repertoire stability over time (Jaccard or Sørensen)
#'
#' @description
#' Plot the output of [compute_repertoire_stability()] as group-wise curves
#' over time with confidence bands. If `drop_self = TRUE` was used at compute
#' time, the plot adds **one baseline marker per group**: a larger point at the
#' group baseline time (x), with y positioned at the **mean similarity of the
#' earliest non-baseline time** in that group, filled by group colour and with a
#' black border.
#'
#' @param stab_df A tibble returned by [compute_repertoire_stability()].
#'   Must contain columns `group`, `time`, `similarity`.
#' @param custom_colors Optional named vector of colors for groups.
#'   Defaults to PHIP palette.
#' @param continuous_mode `"gam"` (default), `"binned"`, or `"loess"`.
#' @param gam_k Requested GAM basis size per series (auto-shrinks to safe k).
#' @param ci_method `"model"` (default) or `"bootstrap"`.
#' @param ci_level Confidence level (default `0.95`).
#' @param boot_R Number of bootstrap replicates (default `500`).
#' @param boot_seed Optional integer seed for reproducibility.
#' @param boot_progress Logical; show textual progress bar during bootstrap
#'   (default `TRUE`).
#' @param ci_fill Ribbon fill color (default `"grey70"`).
#' @param ci_alpha Ribbon alpha (default `0.15`).
#' @param point_alpha Alpha for raw points overlay (default `0.25`; set `0` to hide).
#'
#' @return A `ggplot` object.
#' @export
plot_repertoire_stability <- function(
    stab_df,
    custom_colors        = NULL,
    # smoothing / ribbons
    continuous_mode      = c("gam","binned","loess"),
    gam_k                = 7,
    ci_method            = c("model","bootstrap"),
    ci_level             = 0.95,
    boot_R               = 500,
    boot_seed            = NULL,
    boot_progress        = TRUE,
    # visuals
    ci_fill              = "grey70",
    ci_alpha             = 0.15,
    point_alpha          = 0.25,
    # NEW: override fill colors for big baseline markers (named by group)
    marker_fill_colors   = NULL
) {
  stopifnot(is.data.frame(stab_df))
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting repertoire stability (precomputed)",
    step     = sprintf("mode: %s; ci: %s", continuous_mode, ci_method),
    expr = {
      need <- c("group", "time", "similarity")
      miss <- setdiff(need, names(stab_df))
      .chk_cond(length(miss) > 0,
                sprintf("Missing required columns in stability tibble: %s",
                        paste(miss, collapse = ", ")))

      df <- tibble::as_tibble(stab_df)
      marker_tbl <- attr(stab_df, "baseline_marker", exact = TRUE)
      drop_self  <- isTRUE(attr(stab_df, "drop_self",     exact = TRUE))

      # lock group levels
      grp_levels <- if (!is.null(custom_colors) && length(names(custom_colors))) {
        names(custom_colors)
      } else unique(df$group)
      df$group <- factor(df$group, levels = grp_levels)
      if (is.data.frame(marker_tbl) && nrow(marker_tbl)) {
        marker_tbl$group <- factor(marker_tbl$group, levels = grp_levels)
      }

      # smooths per group
      split_list <- split(df, df$group)
      .ph_log_info("Fitting smooths per group",
                   bullets = c(sprintf("k requested: %d", gam_k),
                               sprintf("groups: %d", length(split_list))))
      preds <- purrr::map_dfr(split_list, function(d) {
        d <- as.data.frame(d)
        if (ci_method == "model") {
          out <- .gam_band_one(d, xvar = "time", yvar = "similarity",
                               k_req = gam_k, level = ci_level, nonneg = TRUE)
        } else {
          out <- .bootstrap_gam_one(d, xvar = "time", yvar = "similarity",
                                    k_req = gam_k, R = boot_R, level = ci_level,
                                    seed = boot_seed, progress = boot_progress,
                                    nonneg = TRUE)
        }
        out$group <- unique(d$group)
        out
      })
      preds$group <- factor(preds$group, levels = grp_levels)
      preds$.grp  <- preds$group

      # base plot (phiper style)
      p <- ggplot2::ggplot()
      if (point_alpha > 0) {
        p <- p +
          ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = .data$time, y = .data$similarity, color = .data$group),
            alpha = point_alpha, size = 1.6, show.legend = TRUE
          )
      }
      p <- p +
        ggplot2::geom_ribbon(
          data = preds,
          ggplot2::aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
          fill = ci_fill, alpha = ci_alpha, colour = NA, inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = preds,
          ggplot2::aes(x = .data$.x, y = .data$.y, color = .data$group, group = .data$.grp),
          linewidth = 1
        ) +
        ggplot2::labs(x = "Months since birth",
                      y = "Repertoire stability (Jaccard)",
                      color = "Group") +
        theme_phip()

      # normal palette for lines/points
      p <- .add_phip_scales(p, custom_colors)

      # baseline markers
      if (isTRUE(drop_self) && is.data.frame(marker_tbl) && nrow(marker_tbl)) {
        if (!is.null(marker_fill_colors)) {
          # normalize mapping: allow single color, unnamed vector, or named vector
          if (length(marker_fill_colors) == 1L && is.null(names(marker_fill_colors))) {
            marker_fill_colors <- setNames(rep(marker_fill_colors, length(grp_levels)), grp_levels)
          } else if (is.null(names(marker_fill_colors)) &&
                     length(marker_fill_colors) == length(grp_levels)) {
            names(marker_fill_colors) <- grp_levels
          }
          # new fill scale ONLY for markers
          if (!requireNamespace("ggnewscale", quietly = TRUE)) {
            .ph_abort("ggnewscale is required when marker_fill_colors is used.", step = "markers")
          }
          p <- p +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_point(
              data = marker_tbl,
              ggplot2::aes(x = .data$x_marker, y = .data$y_marker, fill = .data$group),
              shape = 21, size = 2.5, stroke = 0.6, colour = "black",
              inherit.aes = FALSE, show.legend = FALSE
            ) +
            ggplot2::scale_fill_manual(values = marker_fill_colors, guide = "none")
        } else {
          # default: use same group palette; hide fill legend
          p <- p +
            ggplot2::geom_point(
              data = marker_tbl,
              ggplot2::aes(x = .data$x_marker, y = .data$y_marker, fill = .data$group),
              shape = 21, size = 2.0, stroke = 0.6, colour = "black",
              inherit.aes = FALSE, show.legend = FALSE
            ) +
            ggplot2::guides(fill = "none")
        }
      }

      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# Correlate stability trajectories (on <phip_repertoire_stability>) and plot
# ==============================================================================

#' @title Correlation heatmap of repertoire *stability* trajectories (subjects × subjects)
#' @description
#' Given a tibble produced by [compute_repertoire_stability()], compute
#' pairwise correlations between **subjects' stability time series** (0..1)
#' within an optional group and plot a heatmap.
#'
#' Each subject is a vector over time: `similarity(time)` vs their baseline.
#' We correlate those vectors between subjects using Pearson/Spearman with
#' pairwise-complete handling of missing time points.
#'
#' @param stab      A tibble of class `"phip_repertoire_stability"` with cols:
#'                  `subject_id`, `group`, `time`, `similarity`. Baseline rows
#'                  may have been dropped depending on compute settings.
#' @param group     Character scalar or NULL. If provided, filter to this group
#'                  level (matching `stab$group`). If NULL, all groups pooled.
#' @param min_times Integer; require at least this many time points per subject
#'                  to be included (default 2).
#' @param method    "pearson" (default) or "spearman".
#' @param require_shared_times Logical; if TRUE, keep only **timepoints shared by
#'                  at least 2 subjects** (reduces spurious ±1 from singletons).
#' @param sort_by_status Logical; if TRUE, order subjects by status columns.
#' @param subject_meta Optional tibble with `subject_id` plus optional
#'                  `pre_status_col`, `post_status_col` used when sorting.
#' @param pre_status_col,post_status_col Character names in `subject_meta`
#'                  for ordering (defaults: "status_t1","status_t2").
#' @param palette   Diverging palette for fill scale (default "RdYlBu").
#' @param limits    Numeric length-2 for color scale; default c(-1, 1).
#' @return A ggplot object; invisibly a list(list(plot=..., corr_matrix=...,
#'         subjects=..., order_info=...)).
#' @export
ph_plot_stability_corr_heatmap <- function(
    stab,
    group     = NULL,
    min_times = 2,
    method    = c("pearson","spearman"),
    require_shared_times = TRUE,
    sort_by_status = FALSE,
    subject_meta = NULL,
    pre_status_col  = "status_t1",
    post_status_col = "status_t2",
    palette  = "RdYlBu",
    limits   = c(-1, 1)
) {
  method <- match.arg(method)
  .data <- rlang::.data

  if (!inherits(stab, "phip_repertoire_stability")) {
    .ph_abort("Input must be a <phip_repertoire_stability> tibble.",
              step = "class check",
              bullets = "Call compute_repertoire_stability() first.")
  }

  .ph_with_timing(
    headline = "Stability trajectory correlation (subjects × subjects)",
    step     = sprintf("method=%s, min_times=%d", method, min_times),
    expr = {
      # --- 1) Filter to group (optional) -------------------------------------
      stab0 <- stab
      if (!is.null(group)) {
        stab0 <- dplyr::filter(stab0, .data$group == !!group)
        .ph_log_info("Filtering to group", bullets = paste("group =", group))
      }

      # --- 2) Keep subjects with enough time points ---------------------------
      stab1 <- stab0 |>
        dplyr::group_by(.data$subject_id) |>
        dplyr::mutate(.n_t = dplyr::n_distinct(.data$time)) |>
        dplyr::ungroup() |>
        dplyr::filter(.data$.n_t >= !!min_times) |>
        dplyr::select(-.data$.n_t)

      n_keep <- dplyr::n_distinct(stab1$subject_id)
      .chk_cond(n_keep < 2L, "Fewer than 2 subjects remain after filtering; cannot compute correlations.")

      # --- 3) Optionally restrict to shared time points -----------------------
      if (isTRUE(require_shared_times)) {
        shared_times <- stab1 |>
          dplyr::distinct(.data$subject_id, .data$time) |>
          dplyr::count(.data$time, name = "n_subj") |>
          dplyr::filter(.data$n_subj >= 2L) |>
          dplyr::pull(.data$time)
        stab1 <- dplyr::filter(stab1, .data$time %in% shared_times)
        .ph_log_info("Keeping timepoints shared by ≥2 subjects",
                     bullets = sprintf("timepoints kept: %d", length(shared_times)))
      }

      # --- 4) Build subjects × time wide matrix of stability ------------------
      # Rows: subjects; Cols: times; Values: similarity
      wide <- stab1 |>
        dplyr::select(.data$subject_id, .data$time, .data$similarity) |>
        tidyr::pivot_wider(names_from = .data$time,
                           values_from = .data$similarity) |>
        dplyr::arrange(.data$subject_id)

      subj_ids <- wide$subject_id
      mat <- as.matrix(wide[ , setdiff(colnames(wide), "subject_id"), drop = FALSE])

      # --- 5) Correlation across subjects (pairwise complete) -----------------
      # cor between rows -> cor(t(mat))
      if (ncol(mat) == 0L) {
        .ph_abort("No time columns available after filtering.",
                  step = "build matrix")
      }
      C <- try(stats::cor(t(mat), use = "pairwise.complete.obs", method = method), silent = TRUE)
      if (inherits(C, "try-error")) {
        .ph_abort("Correlation failed", step = "stats::cor", bullets = "Check for constant or empty vectors.")
      }
      dimnames(C) <- list(subj_ids, subj_ids)

      # --- 6) Optional ordering by status (metadata) --------------------------
      order_info <- NULL
      if (isTRUE(sort_by_status)) {
        if (is.null(subject_meta)) {
          .ph_warn("sort_by_status=TRUE but subject_meta is NULL; using alphabetical order.")
        } else {
          need_cols <- c("subject_id", pre_status_col, post_status_col)
          miss <- setdiff(need_cols, colnames(subject_meta))
          if (length(miss)) {
            .ph_warn("subject_meta missing columns; falling back to alphabetical.",
                     bullets = paste("missing:", paste(miss, collapse = ", ")))
          } else {
            order_info <- subject_meta |>
              dplyr::filter(.data$subject_id %in% subj_ids) |>
              dplyr::arrange(.data[[pre_status_col]], .data[[post_status_col]], .data$subject_id)
            ord <- rev(order_info$subject_id)  # mimic your colleague’s “slice(rev(row_number()))”
            C <- C[ord, ord, drop = FALSE]
          }
        }
      }

      # --- 7) Melt for plotting ----------------------------------------------
      cor_df <- reshape2::melt(C, varnames = c("Var1","Var2"), value.name = "value")
      cor_df$Var1 <- factor(cor_df$Var1, levels = rownames(C))
      cor_df$Var2 <- factor(cor_df$Var2, levels = colnames(C))

      # --- 8) Plot ------------------------------------------------------------
      p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(colour = "gray60", linewidth = 0.25) +
        ggplot2::scale_fill_distiller(
          palette  = palette,
          direction = -1,
          limits   = limits,
          na.value = "gray90",
          name     = paste0(stringr::str_to_title(method), "\nCorrelation")
        ) +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(size = 6.5, angle = 45, vjust = 1, hjust = 0.5),
          axis.text.y  = ggplot2::element_text(size = 6.5),
          panel.grid   = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          x = "Subjects (stability trajectories)",
          y = "Subjects (stability trajectories)"
        )

      invisible(list(
        plot        = p,
        corr_matrix = C,
        subjects    = subj_ids,
        order_info  = order_info
      ))

    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# phiper logging utilities (ASCII only; based on the chk and cli packages)
# ==============================================================================
# ---- user-tweakable globals (set via options() in .onLoad or zzz.R) ----------
# options(
#   phiper.log.verbose   = TRUE,
#   phiper.log.time_fmt  = "%Y-%m-%d %H:%M:%S",
#   phiper.log.width     = getOption("width", 80)
# )

.ph_opt <- function(key,
                    default = NULL) {
  getOption(paste0("phiper.log.", key), default)
}

.ph_now <- function() {
  format(Sys.time(), .ph_opt("time_fmt", "%H:%M:%S"))
}

.ph_base_prefix <- function(level = "INFO") {
  sprintf("[%s] %-5s ", .ph_now(), toupper(level)[1])
}

# wraps the text nicely, regardless of the console width
.ph_wrap <- function(text,
                     prefix) {
  w <- .ph_opt("width", getOption("width", 80))
  # strwrap: 'initial' for first line, 'prefix' for continuations
  strwrap(text, width = w, initial = prefix, prefix = strrep(
    " ",
    nchar(prefix)
  ))
}

# Compose multi-depth message lines
# currently the maximal supported log depth is 3:
# headline (required), step (optional), bullets (optional chr vec)
.ph_compose_lines <- function(level,
                              headline,
                              step = NULL,
                              bullets = NULL) {
  base <- .ph_base_prefix(level)
  stepP <- paste0(strrep(" ", nchar(base)), "-> ")
  bullP <- paste0(strrep(" ", nchar(base)), "  - ")

  out <- character(0)
  if (!is.null(headline) && nzchar(headline)) {
    out <- c(out, .ph_wrap(headline, base))
  }
  if (!is.null(step) && nzchar(step)) {
    out <- c(out, .ph_wrap(step, stepP))
  }
  if (length(bullets)) {
    for (b in bullets) {
      if (isTRUE(is.na(b)) || !nzchar(b)) next
      out <- c(out, .ph_wrap(b, bullP))
    }
  }
  out
}

# ---- Public logging helpers --------------------------------------------------
## monitor task progress
.ph_log_info <- function(headline,
                         step = NULL,
                         bullets = NULL,
                         verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("INFO", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

## monitor task progression
.ph_log_ok <- function(headline,
                       step = NULL,
                       bullets = NULL,
                       verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("OK", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

# Warnings/errors via chk, but formatted to match the style of the logger
.ph_warn <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("WARN", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::wrn(msg, ...) # single string -> respects \n
  } else {
    warning(msg, call. = FALSE, ...) # fallback if chk not installed
  }
  invisible(lines)
}

.ph_abort <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("ERROR", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::abort_chk(msg, ...) # single string -> respects \n
  } else {
    stop(msg, call. = FALSE, ...) # fallback if chk not installed
  }
}

# ---- Back-compaibility (vecmatch) --------------------------------------------
# drop-in replacement for older .msg() funct
.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) .ph_log_info(headline = paste(...))
  invisible(NULL)
}

# original conditional helper to not break down older code
# upgraded to the unified phiper style
.chk_cond <- function(condition,
                      error_message,
                      error = TRUE,
                      step = NULL,
                      bullets = NULL,
                      ...) {
  # log nopthing
  if (!isTRUE(condition)) {
    return(invisible(FALSE))
  }

  # print error and abort exec
  if (isTRUE(error)) {
    .ph_abort(headline = error_message, step = step, bullets = bullets, ...)
  } else {
    # print warning and go on
    .ph_warn(headline = error_message, step = step, bullets = bullets, ...)
  }
  invisible(TRUE)
}

# ---- timing helper for sections ----------------------------------------------
## many tasks in phiper can be long/take a while; it was important to have the
## infos on timing - this func wraps a task to get a start/end pair in the same
## style
.ph_with_timing <- function(headline,
                            step = NULL,
                            bullets = NULL,
                            expr,
                            verbose = .ph_opt("verbose", TRUE)) {
  t0 <- Sys.time()
  .ph_log_info(headline = headline, step = step, bullets = bullets, verbose = verbose)

  res <- tryCatch(
    {
      force(expr)
    }, # evaluate user's code
    finally = {
      dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
      .ph_log_ok(
        headline = paste0(headline, " - done"),
        step     = sprintf("elapsed: %ss", dt),
        verbose  = verbose
      )
    }
  )

  res
}

# ==============================================================================
# phiper checks + additional helpers (ASCII-only, unified with phiper logger)
# it depends on: .ph_abort(), .ph_warn(), .chk_cond(), word_list(), add_quotes()
# ==============================================================================

# -- check if filename has given extension ------------------------------------
# comes in handy when loading .csv or .parquet; provide filename and vector of
# extensions to check (eg c(".csv", ".parquet"))
.chk_extension <- function(name,
                           x_name,
                           ext_vec) {
  if (is.null(ext_vec) || !length(ext_vec)) {
    return(invisible(TRUE))
  }

  base <- basename(name %||% "") # extracting filename from paths
  parts <- strsplit(base, "\\.", fixed = FALSE)[[1]] # names + complex ext

  # taking last ext after . (eg .tar.gz --> .gz)
  ext <- if (length(parts) > 1L) {
    tolower(paste(parts[-1L], collapse = "."))
  } else {
    ""
  }

  norm <- function(x) sub("^\\.+", "", tolower(x)) # normalize ext
  got <- if (nzchar(ext)) norm(ext) else "<none>"
  ok <- nzchar(ext) && got %in% norm(ext_vec) # final ext check

  if (!ok) {
    .ph_abort(
      headline = sprintf("Invalid file extension for `%s`.", x_name),
      step = sprintf("validating path: %s", name),
      bullets = c(
        sprintf("got: %s", add_quotes(got, 2L)),
        sprintf(
          "allowed: %s",
          word_list(add_quotes(norm(ext_vec), 2L), and_or = "or")
        )
      )
    )
  }
  invisible(TRUE)
}

# -- check if NULL and replace with default when TRUE (warn in unified style) --
.chk_null_default <- function(x,
                              x_name,
                              method,
                              default) {
  if (is.null(x)) {
    # format the default for print
    fmt <- function(v) {
      if (is.character(v) && length(v) == 1L) {
        return(add_quotes(v, 2L))
      }
      if (is.atomic(v) && length(v) == 1L) {
        return(as.character(v))
      }
      sprintf("<%s>", paste(class(v), collapse = "/"))
    }

    # generate warning and the replace
    .ph_warn(
      headline = sprintf("Argument `%s` missing; using default.", x_name),
      step     = sprintf("method: %s", add_quotes(method, 2L)),
      bullets  = sprintf("default: %s", fmt(default))
    )
    x <- default
  }
  x
}

# -- validate path to file -----------------------------------------------------
.chk_path <- function(path,
                      arg_name,
                      extension) {
  ## error when path not a string
  .chk_cond(
    !chk::vld_string(path),
    sprintf("`%s` must be a character scalar.", arg_name),
    step    = "path validation",
    bullets = sprintf("got class: %s", paste(class(path), collapse = "/"))
  )

  ## error when path does not exist
  .chk_cond(
    !chk::vld_file(path),
    sprintf("File for `%s` does not exist.", arg_name),
    step    = "path validation",
    bullets = sprintf("path: %s", path)
  )

  # optionally extension check if provided
  if (!missing(extension) && length(extension)) {
    .chk_extension(
      path,
      arg_name,
      extension
    )
  }

  invisible(TRUE)
}

# -- clean wordlists for message generation ------------------------------------
# for multiple arguments/values
word_list <- function(word_list = NULL,
                      and_or = "and",
                      is_are = FALSE,
                      quotes = FALSE) {
  # Make "a and b" / "a, b, and c"; optionally append "is/are".
  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))
    if (len_wl == 2L) {
      out <- sprintf("%s %s %s", word_list[1L], and_or, word_list[2L])
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or, word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)
  attr(out, "plural") <- TRUE
  out
}

# -- quoting helper (unified error style) --------------------------------------
# define number of quotes you want --> for printing logs/messages/warnings
# or define the quotes itself as a string
add_quotes <- function(x,
                       quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  if (isTRUE(quotes)) quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    .ph_abort(
      headline = "Invalid `quotes` argument.",
      step = "formatting add_quotes()",
      bullets = c(
        "allowed: FALSE, TRUE, 0, 1, 2, or a single-character string",
        sprintf("got class: %s", paste(class(quotes), collapse = "/"))
      )
    )
  }

  if (quotes == 0L) {
    return(x)
  }
  if (quotes == 1L) {
    return(sprintf("'%s'", x))
  }
  sprintf('"%s"', x)
}

# -- not-in operator -----------------------------------------------------------
`%nin%` <- function(x, inx) {
  !(x %in% inx)
}

# -- NULL-coalescing helper ----------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y


# ==============================================================================
# database helpers
# ==============================================================================
# ensure peptide_library is queryable from the SAME connection as data_long
.ensure_peplib_on_main <- function(x, schema_alias = "peplib") {
  main_con <- dbplyr::remote_con(x$data_long)
  pep_con <- if (!is.null(x$meta$peptide_con)) {
    x$meta$peptide_con
  } else {
    dbplyr::remote_con(x$peptide_library)
  }

  # try zero-copy ATTACH when both are DuckDB
  if (inherits(main_con, "duckdb_connection") &&
      inherits(pep_con, "duckdb_connection")) {
    pep_db_path <- try(pep_con@driver@dbdir, silent = TRUE)
    if (!inherits(pep_db_path, "try-error") &&
        is.character(pep_db_path) &&
        nzchar(pep_db_path)) {
      # ATTACH (ignore "already attached")
      try(
        DBI::dbExecute(
          main_con,
          sprintf(
            "ATTACH '%s' AS %s;",
            pep_db_path,
            schema_alias
          )
        ),
        silent = TRUE
      )

      # Detect base table name (fallback: "peptide_meta")
      base_name <- tryCatch(
        {
          nm <- dbplyr::remote_name(x$peptide_library)
          if (is.null(nm) || !nzchar(nm)) {
            "peptide_meta"
          } else {
            sub("^.*\\.", "", nm)
          }
        },
        error = function(e) "peptide_meta"
      )

      # Try both two-part and three-part references
      try_tbl <- function(sql_expr) {
        tryCatch(dplyr::tbl(main_con, dbplyr::sql(sql_expr)),
                 error = function(e) NULL
        )
      }

      candidates <- c(
        sprintf("SELECT * FROM %s.%s", schema_alias, base_name),
        sprintf("SELECT * FROM %s.main.%s", schema_alias, base_name)
      )

      for (q in candidates) {
        out <- try_tbl(q)
        if (!is.null(out)) {
          return(out)
        }
      }
      # If both fail, we’ll fall through to the copy_to() fallback below.
    }
  }

  # Fallback: copy peptidelib into main_con as a TEMP table
  peplib_local <- dplyr::collect(x$peptide_library)
  tmp_name <- paste0("peptide_meta_tmp_", as.integer(Sys.time()))
  dplyr::copy_to(main_con, peplib_local, tmp_name,
                 temporary = TRUE,
                 overwrite = TRUE
  )
}

# ------------------------------------------------------------------------------
# GENERICS
# ------------------------------------------------------------------------------
plot_enrichment_counts <- function(phip_data,
                                   group_cols = NULL,
                                   prevalence_threshold = 0.05,
                                   custom_colors = NULL,
                                   binwidth = 1,
                                   group_interaction = FALSE,
                                   interaction_sep = " * ",
                                   ...) {
  UseMethod("plot_enrichment_counts")
}

plot_alpha_diversity <- function(x, ...) {
  UseMethod("plot_alpha_diversity")
}

# ------------------------------------------------------------------------------
# METHODS (S3)
# ------------------------------------------------------------------------------

# -- plot_enrichment_counts: <phip_data> ---------------------------------------
plot_enrichment_counts.phip_data <- function(phip_data,
                                             group_cols = NULL,
                                             prevalence_threshold = 0.05,
                                             custom_colors = NULL,
                                             binwidth = 1,
                                             group_interaction = FALSE,
                                             interaction_sep = " * ",
                                             ...) {
  x <- phip_data
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data

  .ph_with_timing(
    headline = "Plotting enrichment counts (<phip_data>)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else {
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
                                      collapse = ", "
      ))
    },
    expr = {
      tbl <- x$data_long

      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
        red_txt <- tryCatch(
          {
            ep <- as.numeric(x$meta$exist_prop)
            if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1 / ep) else "<unknown>"
          },
          error = function(e) "<unknown>"
        )
        .ph_log_info(
          "Full-cross detected; pruning non-existing rows before plotting",
          bullets = c(
            "rule: keep exist == 1",
            sprintf("estimated reduction: %s", red_txt)
          )
        )
        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      if (!is.null(group_cols) && length(group_cols)) {
        missing_gcs <- setdiff(group_cols, colnames(tbl))
        if (length(missing_gcs)) {
          .ph_abort(
            headline = "Missing grouping columns in data_long.",
            step = "input validation",
            bullets = sprintf(
              "missing: %s",
              paste(add_quotes(missing_gcs, 1L),
                    collapse = ", "
              )
            )
          )
        }
      }

      if (is.null(group_cols) || !length(group_cols)) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        plots <- lapply(group_cols, function(gc) {
          .plot_enrichment_counts_one(
            tbl,
            group_col = gc,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = gc
          )
        })
        names(plots) <- group_cols

        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep
                                       ))

          plots[[combo_nm]] <- .plot_enrichment_counts_one(
            tbl_inter,
            group_col = inter_col,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = combo_nm
          )
        } else if (isTRUE(group_interaction) && length(group_cols) < 2L) {
          .ph_warn(
            headline = "group_interaction requested but fewer than
            2 group_cols supplied.",
            step     = "interaction plot",
            bullets  = "Interaction plot skipped."
          )
        }

        if (length(plots) == 1L) plots[[1L]] else plots
      }
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_enrichment_counts: data.frame ----------------------------------------
plot_enrichment_counts.data.frame <- function(phip_data,
                                              group_cols = NULL,
                                              prevalence_threshold = 0.05,
                                              custom_colors = NULL,
                                              binwidth = 1,
                                              group_interaction = FALSE,
                                              interaction_sep = " * ",
                                              ...) {
  tbl <- phip_data
  .data <- rlang::.data

  .ph_with_timing(
    headline = "Plotting enrichment counts (data.frame)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else {
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
                                      collapse = ", "
      ))
    },
    expr = {
      need_base <- c("sample_id", "peptide_id", "exist")
      miss_base <- setdiff(need_base, colnames(tbl))
      if (length(miss_base)) {
        .ph_abort(
          headline = "Missing required columns in data.frame.",
          step = "input validation",
          bullets = sprintf(
            "missing: %s",
            paste(add_quotes(miss_base, 1L),
                  collapse = ", "
            )
          )
        )
      }

      if (is.null(group_cols) || !length(group_cols)) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        plots <- lapply(group_cols, function(gc) {
          if (!gc %in% colnames(tbl)) {
            .ph_abort(
              headline = "Grouping column not found in data.frame.",
              step     = "input validation",
              bullets  = sprintf("group_col: %s", add_quotes(gc, 1L))
            )
          }
          .plot_enrichment_counts_one(
            tbl,
            group_col = gc,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = gc
          )
        })
        names(plots) <- group_cols

        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep
                                       ))

          plots[[combo_nm]] <- .plot_enrichment_counts_one(
            tbl_inter,
            group_col = inter_col,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = combo_nm
          )
        }

        if (length(plots) == 1L) plots[[1L]] else plots
      }
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: <phip_data> -----------------------------------------
plot_alpha_diversity.phip_data <- function(
    x,
    metric = c("richness", "shannon_diversity", "simpson_diversity"),
    group_col = "group",
    rank_col = "rank",
    filter_groups = NULL,
    filter_ranks = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam", "binned", "loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model", "bootstrap"),
    ci_level = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...) {
  stopifnot(inherits(x, "phip_data"))
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (<phip_data>)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      if (!is.null(time_col) && isTRUE(!isTRUE(x$meta$longitudinal))) {
        .ph_abort(
          headline = "Longitudinal plotting not allowed for this object.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      plot_alpha_diversity(
        alpha_list,
        metric = metric,
        group_col = group_col,
        rank_col = rank_col,
        filter_groups = filter_groups,
        filter_ranks = filter_ranks,
        custom_colors = custom_colors,
        facet_by_rank = facet_by_rank,
        ncol = ncol,
        facet_scales = facet_scales,
        time_col = time_col,
        continuous_mode = continuous_mode,
        gam_k = gam_k,
        point_alpha = point_alpha,
        ci_method = ci_method,
        ci_level = ci_level,
        boot_R = boot_R,
        boot_seed = boot_seed,
        boot_progress = boot_progress,
        ci_fill = ci_fill,
        ci_alpha = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: data.frame ------------------------------------------
plot_alpha_diversity.data.frame <- function(
    x,
    metric = c("richness", "shannon_diversity", "simpson_diversity"),
    group_col = "group",
    rank_col = "rank",
    filter_groups = NULL,
    filter_ranks = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam", "binned", "loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model", "bootstrap"),
    ci_level = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...) {
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (data.frame)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      plot_alpha_diversity(
        alpha_list,
        metric = metric,
        group_col = group_col,
        rank_col = rank_col,
        filter_groups = filter_groups,
        filter_ranks = filter_ranks,
        custom_colors = custom_colors,
        facet_by_rank = facet_by_rank,
        ncol = ncol,
        facet_scales = facet_scales,
        time_col = time_col,
        continuous_mode = continuous_mode,
        gam_k = gam_k,
        point_alpha = point_alpha,
        ci_method = ci_method,
        ci_level = ci_level,
        boot_R = boot_R,
        boot_seed = boot_seed,
        boot_progress = boot_progress,
        ci_fill = ci_fill,
        ci_alpha = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: precomputed -----------------------------------------
plot_alpha_diversity.phip_alpha_diversity <- function(
    x,
    metric = c("richness", "shannon_diversity", "simpson_diversity"),
    group_col = "group",
    rank_col = "rank",
    filter_groups = NULL,
    filter_ranks = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam", "binned", "loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model", "bootstrap"),
    ci_level = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...) {
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (precomputed)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      alpha_df <- if (is.list(x) && !is.null(group_col) &&
                      group_col %in% names(x)) {
        tibble::as_tibble(x[[group_col]])
      } else if (is.list(x)) {
        inner <- x[vapply(
          x, function(el) inherits(el, "data.frame"),
          logical(1)
        )]
        if (length(inner) == 1L) {
          tibble::as_tibble(inner[[1L]])
        } else if (length(inner) >= 1L) {
          suppressWarnings(dplyr::bind_rows(inner))
        } else {
          if (length(x) && inherits(x[[1L]], "data.frame")) {
            tibble::as_tibble(x[[1L]])
          } else {
            tibble::as_tibble(x)
          }
        }
      } else {
        tibble::as_tibble(x)
      }

      if (!is.null(attr(x, "group_cols", exact = TRUE))) {
        attr(alpha_df, "group_cols") <- attr(x, "group_cols", exact = TRUE)
      }
      if (!is.null(attr(x, "ranks", exact = TRUE))) {
        attr(alpha_df, "ranks") <- attr(x, "ranks", exact = TRUE)
      }

      metrics_norm <- tolower(metric)
      for (m in metrics_norm) {
        if (!m %in% names(alpha_df)) {
          .ph_abort(
            headline = "Metric column not found in alpha data.",
            step     = "input validation",
            bullets  = sprintf("missing: %s", add_quotes(m, 1L))
          )
        }
      }

      if (!is.null(group_col) && !group_col %in% names(alpha_df)) {
        .ph_abort(
          headline = "Grouping column not found in precomputed alpha data.",
          step = "input validation",
          bullets = c(
            sprintf("requested group_col: %s", add_quotes(group_col, 1L)),
            sprintf(
              "available columns: %s",
              paste(add_quotes(colnames(alpha_df), 1L), collapse = ", ")
            )
          )
        )
      }

      if (!is.null(rank_col) && !rank_col %in% names(alpha_df)) {
        .ph_warn(
          headline = "Rank column not found in precomputed alpha data;
          disabling faceting.",
          step     = "input validation",
          bullets  = sprintf("rank_col: %s", add_quotes(rank_col, 1L))
        )
        rank_col <- NULL
      }

      if (!is.null(group_col) && !is.null(filter_groups)) {
        alpha_df <- dplyr::filter(
          alpha_df,
          .data[[group_col]] %in% !!filter_groups
        )
      }
      if (!is.null(rank_col) && !is.null(filter_ranks)) {
        alpha_df <- dplyr::filter(
          alpha_df,
          .data[[rank_col]] %in% !!filter_ranks
        )
      }

      plots <- lapply(metrics_norm, function(m) {
        .build_alpha_plot(
          alpha_df = alpha_df,
          metric = m,
          group_col = group_col,
          rank_col = rank_col,
          filter_groups = filter_groups,
          filter_ranks = filter_ranks,
          custom_colors = custom_colors,
          facet_by_rank = facet_by_rank,
          ncol = ncol,
          facet_scales = facet_scales,
          time_col = if (!is.null(time_col) &&
                         time_col %in% names(alpha_df)) {
            time_col
          } else {
            NULL
          },
          continuous_mode = continuous_mode,
          gam_k = gam_k,
          point_alpha = point_alpha,
          ci_method = ci_method,
          ci_level = ci_level,
          boot_R = boot_R,
          boot_seed = boot_seed,
          boot_progress = boot_progress,
          ci_fill = ci_fill,
          ci_alpha = ci_alpha
        )
      })
      names(plots) <- metrics_norm
      if (length(plots) == 1L) plots[[1L]] else plots
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# HELPERS (internal) – put at the end as requested
# ------------------------------------------------------------------------------

# -- enrichment: shared builder -------------------------------------------------
.plot_enrichment_counts_one <- function(tbl,
                                        group_col = NULL,
                                        prevalence_threshold = 0.05,
                                        custom_colors = NULL,
                                        binwidth = 1,
                                        title_label = if (is.null(group_col)) {
                                          "all samples"
                                        } else {
                                          group_col
                                        }) {
  .data <- rlang::.data
  grp_sym <- if (!is.null(group_col)) rlang::sym(group_col) else NULL

  .ph_with_timing(
    headline = "Building enrichment count plot",
    step = if (is.null(group_col)) "no grouping (aggregate)" else
      sprintf("grouping variable: '%s'", group_col),
    expr = {
      need <- c("sample_id", "peptide_id", "exist")
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(tbl))
      if (length(miss)) {
        .ph_abort(
          headline = "Missing required columns for plotting.",
          step = "input validation",
          bullets = sprintf(
            "missing: %s",
            paste(add_quotes(miss, 1L), collapse = ", ")
          )
        )
      }

      tbl <- if (is.null(group_col)) {
        dplyr::mutate(tbl, Cohort = "All samples")
      } else {
        dplyr::mutate(tbl, Cohort = !!grp_sym)
      }

      tbl_present <- dplyr::filter(tbl, .data$exist > 0)

      cohort_sizes <- tbl |>
        dplyr::distinct(.data$sample_id, .data$Cohort) |>
        dplyr::count(.data$Cohort, name = "n_samples") |>
        dplyr::collect()

      if (nrow(cohort_sizes) == 0L || any(cohort_sizes$n_samples <= 0)) {
        .ph_warn(
          headline = "No samples for at least one cohort (facet may be empty).",
          step     = "cohort sizes"
        )
        cohort_sizes <- dplyr::filter(cohort_sizes, .data$n_samples > 0)
      }

      pep_counts <- tbl_present |>
        dplyr::group_by(.data$Cohort, .data$peptide_id) |>
        dplyr::summarise(
          n_present = dplyr::n_distinct(.data$sample_id),
          .groups = "drop"
        ) |>
        dplyr::filter(.data$n_present > 0) |>
        dplyr::collect()

      thresholds <- cohort_sizes |>
        dplyr::mutate(thresh = ceiling(.data$n_samples * prevalence_threshold))

      n_thresh_tbl <- pep_counts |>
        dplyr::inner_join(dplyr::select(thresholds, .data$Cohort, .data$thresh),
                          by = "Cohort") |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_peptides_thresh = sum(.data$n_present >= .data$thresh),
          .groups = "drop"
        )

      overall_tbl <- pep_counts |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_overall = dplyr::n_distinct(.data$peptide_id),
          .groups = "drop"
        )

      thresholds_df <- thresholds |>
        dplyr::left_join(n_thresh_tbl, by = "Cohort") |>
        dplyr::left_join(overall_tbl, by = "Cohort") |>
        dplyr::mutate(
          n_peptides_thresh = dplyr::coalesce(.data$n_peptides_thresh, 0L),
          n_overall         = dplyr::coalesce(.data$n_overall, 0L),
          y_line            = pmax(.data$n_peptides_thresh, 1L),
          x_mid             = (.data$thresh + .data$n_samples) / 2
        )

      real_order <- thresholds_df$Cohort
      count_df <- pep_counts |>
        dplyr::mutate(Cohort = factor(.data$Cohort, levels = real_order))

      lvl <- levels(count_df$Cohort)
      pal_map <- if (is.null(custom_colors)) {
        .phip_palette_map(lvl)
      } else {
        if (is.null(names(custom_colors))) {
          stats::setNames(rep_len(custom_colors, length(lvl)), lvl)
        } else {
          stats::setNames(unname(custom_colors[lvl]), lvl)
        }
      }

      count_df <- dplyr::mutate(
        count_df,
        fill_col = unname(pal_map[as.character(.data$Cohort)])
      )

      label_size <- 4
      p <- ggplot2::ggplot(count_df, ggplot2::aes(
        x = .data$n_present,
        fill = .data$fill_col
      )) +
        ggplot2::geom_histogram(
          binwidth = binwidth,
          position = "identity",
          alpha = 1,
          colour = NA
        ) +
        ggplot2::scale_y_log10(
          breaks = 10^(0:6),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          expand = ggplot2::expansion(mult = c(0, .15))
        ) +
        ggplot2::annotation_logticks(sides = "l", scaled = TRUE) +
        ggplot2::scale_fill_identity(guide = "none") +
        ggplot2::labs(
          x = "# of observations",
          y = expression("# of significantly bound peptides (" * log[10] * ")")
        ) +
        ggplot2::geom_segment(
          data = thresholds_df,
          ggplot2::aes(x = .data$thresh, xend = .data$n_samples,
                       y = .data$y_line, yend = .data$y_line),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black",
          linewidth = 0.4,
          arrow = ggplot2::arrow(length = grid::unit(0.1, "cm"), ends = "both")
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(
              .data$n_peptides_thresh, " peptides in ≥",
              round(prevalence_threshold * 100), "%"
            )
          ),
          inherit.aes = FALSE, vjust = -0.6, size = label_size
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(.data$n_overall, " peptides overall")
          ),
          inherit.aes = FALSE, vjust = 1.6, size = label_size
        ) +
        theme_phip() +
        ggplot2::ggtitle(sprintf("Enrichment counts by %s", title_label))

      if (!is.null(group_col)) {
        p <- p + ggplot2::facet_wrap(~Cohort, ncol = 2, scales = "free_x")
      }

      .ph_log_ok("Plot built")
      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- alpha-diversity helpers ----------------------------------------------------
.norm_metric <- function(metric) {
  allowed <- c("richness", "shannon_diversity", "simpson_diversity")
  key <- tolower(metric[1])
  if (key %in% allowed) return(key)
  .ph_abort(
    headline = "Unknown metric.",
    step     = "argument validation",
    bullets  = sprintf("allowed: %s", paste(allowed, collapse = ", "))
  )
}

.safe_k <- function(x, k_req) {
  n_ux <- length(unique(stats::na.omit(x)))
  k_eff <- max(3L, min(as.integer(k_req), as.integer(n_ux) - 1L))
  if (!is.finite(k_eff) || k_eff < 3L) 3L else k_eff
}

.gam_band_one <- function(d, xvar, yvar,
                          k_req = 7, level = 0.95, nonneg = FALSE) {
  x <- d[[xvar]]
  rng <- range(x, na.rm = TRUE)
  xgrid <- if (rng[1] == rng[2]) sort(unique(x)) else seq(rng[1], rng[2], length.out = 200)
  k_eff <- .safe_k(x, k_req)
  form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
  fit <- mgcv::gam(form, data = d)
  nd <- data.frame(xgrid); names(nd) <- xvar
  pr <- mgcv::predict.gam(fit, newdata = nd, type = "response", se.fit = TRUE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  mu <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  lwr <- mu - z * se; upr <- mu + z * se
  if (isTRUE(nonneg)) { mu <- pmax(mu, 0); lwr <- pmax(lwr, 0) }
  tibble::tibble(.x = nd[[xvar]], .y = mu, lwr = lwr, upr = upr, .k = k_eff)
}

.bootstrap_gam_one <- function(d, xvar, yvar,
                               k_req = 7, R = 500, level = 0.95,
                               seed = NULL, progress = TRUE, nonneg = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  ctr <- .gam_band_one(d, xvar, yvar, k_req = k_req, level = 0.95, nonneg = nonneg)
  xgrid <- ctr$.x; mu0 <- ctr$.y
  mat <- matrix(NA_real_, nrow = length(xgrid), ncol = R)
  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = R, style = 3) else NULL
  for (r in seq_len(R)) {
    dd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    k_eff <- .safe_k(dd[[xvar]], k_req)
    form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
    fit <- try(mgcv::gam(form, data = dd), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      nd <- data.frame(xgrid); names(nd) <- xvar
      mat[, r] <- as.numeric(mgcv::predict.gam(fit, newdata = nd, type = "response"))
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, r)
  }
  if (!is.null(pb)) close(pb)
  alpha <- 1 - level
  lwr <- apply(mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, type = 6)
  upr <- apply(mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 6)
  if (isTRUE(nonneg)) lwr <- pmax(lwr, 0)
  tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr, .k = unique(ctr$.k))
}

.add_phip_scales <- function(p, custom_colors = NULL) {
  if (is.null(custom_colors)) {
    p + scale_color_phip() + scale_fill_phip()
  } else {
    p + ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::scale_fill_manual(values = custom_colors)
  }
}

.maybe_facet_by_rank <- function(p, df, rank_col, facet_by_rank, ncol, facet_scales) {
  if (!isTRUE(facet_by_rank)) return(p)
  if (is.null(rank_col) || !rank_col %in% names(df)) return(p)
  if (length(unique(df[[rank_col]])) <= 1) return(p)
  p + ggplot2::facet_wrap(stats::as.formula(paste("~", rank_col)),
                          ncol = ncol, scales = facet_scales)
}

.build_alpha_plot <- function(alpha_df,
                              metric,
                              group_col = "group",
                              rank_col = "rank",
                              filter_groups = NULL,
                              filter_ranks = NULL,
                              custom_colors = NULL,
                              facet_by_rank = TRUE,
                              ncol = 2,
                              facet_scales = "fixed",
                              time_col = NULL,
                              continuous_mode = c("gam", "binned", "loess"),
                              gam_k = 7,
                              point_alpha = 0.25,
                              ci_method = c("model", "bootstrap"),
                              ci_level = 0.95,
                              boot_R = 500,
                              boot_seed = NULL,
                              boot_progress = TRUE,
                              ci_fill = "grey70",
                              ci_alpha = 0.15) {
  metric_col <- .norm_metric(metric)
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)
  df <- tibble::as_tibble(alpha_df)

  if (!metric_col %in% names(df)) {
    .ph_abort(
      headline = "Metric column not found in alpha data.",
      step     = "input validation",
      bullets  = sprintf("missing: %s", add_quotes(metric_col, 1L))
    )
  }
  if (!is.null(group_col) && !group_col %in% names(df)) {
    .ph_abort(
      headline = "Grouping column not found in alpha data.",
      step     = "input validation",
      bullets  = sprintf("group_col: %s", add_quotes(group_col, 1L))
    )
  }
  if (!is.null(rank_col) && !rank_col %in% names(df)) {
    .ph_warn(
      headline = "Rank column not found in alpha data; ignoring faceting.",
      step     = "input validation",
      bullets  = sprintf("rank_col: %s", add_quotes(rank_col, 1L))
    )
    rank_col <- NULL
  }

  if (!is.null(filter_groups) && !is.null(group_col)) {
    df <- df[df[[group_col]] %in% filter_groups, , drop = FALSE]
  }
  if (!is.null(filter_ranks) && !is.null(rank_col)) {
    df <- df[df[[rank_col]] %in% filter_ranks, , drop = FALSE]
  }

  keep_cols <- c(metric_col, group_col, rank_col, time_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  if (length(keep_cols)) {
    df <- df[stats::complete.cases(df[, keep_cols, drop = FALSE]), , drop = FALSE]
  }

  ylab <- switch(metric_col,
                 "richness"          = "Richness",
                 "shannon_diversity" = "Shannon diversity",
                 "simpson_diversity" = "Simpson diversity (1 - \u03A3 p^2)")

  # CROSS-SECTIONAL ------------------------------------------------------------
  if (is.null(time_col) || !time_col %in% names(df)) {
    if (is.null(group_col)) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = "", y = .data[[metric_col]])) +
        ggplot2::geom_boxplot(fill = "grey70", colour = "black", outlier.shape = NA) +
        ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 1) +
        ggplot2::labs(x = NULL, y = ylab) +
        theme_phip() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
      p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
      return(p)
    }

    gsym <- rlang::sym(group_col)
    msym <- rlang::sym(metric_col)
    df_counts <- df |>
      dplyr::group_by(!!gsym) |>
      dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
    xlab_map <- stats::setNames(
      paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
      df_counts[[group_col]]
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!gsym, y = !!msym, fill = !!gsym)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::geom_jitter(color = "black", size = 1, width = 0.2,
                           alpha = point_alpha, show.legend = FALSE) +
      ggplot2::scale_x_discrete(labels = xlab_map) +
      ggplot2::labs(x = "Group", y = ylab, fill = group_col) +
      theme_phip()
    p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  # LONGITUDINAL ---------------------------------------------------------------
  if (continuous_mode == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      .ph_abort("mgcv is required for GAM mode.", step = "longitudinal")
    }

    split_keys <- c(group_col, if (!is.null(rank_col) && isTRUE(facet_by_rank)) rank_col else NULL)
    df_split <- if (length(split_keys)) {
      df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
        dplyr::group_split()
    } else list(df)

    .ph_log_info("Fitting GAM smooths (auto-shrinking k per series)",
                 bullets = c(sprintf("k requested: %d", gam_k),
                             sprintf("series: %d", length(df_split))))

    preds <- purrr::map_dfr(df_split, function(dsub) {
      dsub <- as.data.frame(dsub)
      out <- if (ci_method == "model") {
        .gam_band_one(dsub, xvar = time_col, yvar = metric_col,
                      k_req = gam_k, level = ci_level,
                      nonneg = (metric_col == "richness"))
      } else {
        .bootstrap_gam_one(dsub, xvar = time_col, yvar = metric_col,
                           k_req = gam_k, R = boot_R, level = ci_level,
                           seed = boot_seed, progress = boot_progress,
                           nonneg = (metric_col == "richness"))
      }
      if (!is.null(group_col)) out[[group_col]] <- unique(dsub[[group_col]])
      if (!is.null(rank_col) && isTRUE(facet_by_rank)) out[[rank_col]] <- unique(dsub[[rank_col]])
      out
    })

    preds$.grp <- if (!is.null(group_col) && !is.null(rank_col) && isTRUE(facet_by_rank)) {
      interaction(preds[[group_col]], preds[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      preds[[group_col]]
    } else factor("all")

    p <- ggplot2::ggplot()
    if (!is.null(group_col)) {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]], color = .data[[group_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0, colour = "black"
      )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = preds,
        ggplot2::aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
        inherit.aes = FALSE, fill = ci_fill, alpha = ci_alpha, colour = NA
      ) +
      {
        if (!is.null(group_col)) {
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y,
                         color = .data[[group_col]], group = .data$.grp),
            linewidth = 1
          )
        } else {
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.grp),
            linewidth = 1, colour = "black"
          )
        }
      } +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      theme_phip()

    if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  if (continuous_mode == "binned") {
    rng <- range(df[[time_col]], na.rm = TRUE)
    nbins <- 20
    binwidth <- (rng[2] - rng[1]) / nbins
    make_bins <- function(x, width, xmin) floor((x - xmin) / width)
    df$.bin <- make_bins(df[[time_col]], binwidth, rng[1])
    df$.bin_mid <- rng[1] + (df$.bin + 0.5) * binwidth

    group_map <- c(
      if (!is.null(group_col)) group_col else NULL,
      if (!is.null(rank_col) && isTRUE(facet_by_rank)) rank_col else NULL,
      ".bin", ".bin_mid"
    )

    sdat <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_map))) |>
      dplyr::summarise(
        y = mean(.data[[metric_col]], na.rm = TRUE),
        sd = stats::sd(.data[[metric_col]], na.rm = TRUE),
        n = dplyr::n(), .groups = "drop"
      )

    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    sdat <- sdat |>
      dplyr::mutate(
        se = sd / sqrt(pmax(n, 1)),
        lwr = y - z * se, upr = y + z * se
      )

    sdat$.grp <- if (!is.null(group_col) && !is.null(rank_col) && isTRUE(facet_by_rank)) {
      interaction(sdat[[group_col]], sdat[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      sdat[[group_col]]
    } else {
      factor("all")
    }

    p <- ggplot2::ggplot()
    if (!is.null(group_col)) {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]], color = .data[[group_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0, colour = "black"
      )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = sdat,
        ggplot2::aes(x = .data$.bin_mid, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
        fill = ci_fill, alpha = ci_alpha, colour = NA, inherit.aes = FALSE
      ) +
      {
        if (!is.null(group_col)) {
          ggplot2::geom_line(
            data = sdat,
            ggplot2::aes(x = .data$.bin_mid, y = .data$y,
                         color = .data[[group_col]], group = .data$.grp),
            linewidth = 1
          )
        } else {
          ggplot2::geom_line(
            data = sdat,
            ggplot2::aes(x = .data$.bin_mid, y = .data$y, group = .data$.grp),
            linewidth = 1, colour = "black"
          )
        }
      } +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      theme_phip()

    if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  # LOESS fallback
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]])) +
    {
      if (!is.null(group_col)) ggplot2::aes(color = .data[[group_col]])
      else ggplot2::aes()
    } +
    ggplot2::geom_point(alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75, level = ci_level) +
    ggplot2::labs(x = time_col, y = ylab, color = group_col) +
    theme_phip()

  if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
  p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
  p
}

# ==============================================================================
# Repertoire similarity / stability: longitudinal + group-aware (DuckDB-backed)
# ==============================================================================

#' @title Compute repertoire similarity with group/subject/time pairing controls
#'
#' @description
#' Flexible engine for binary repertoire similarity that supports:
#' - Group scope: within-group, between-group, or both (union).
#' - Subject pairing: only self, self-or-dyad/triad, or all subjects.
#' - Time behavior: reference to a baseline (earliest/latest) **or** all pairwise
#'   comparisons across categorical timepoints.
#'
#' By default this reproduces the earlier "stability" behavior (baseline within
#' subject & group), but now also enables subject-by-subject matrices per
#' group×time for categorical designs (dyads/triads supported).
#'
#' **Similarities (binary):**
#' - `"jaccard"`      : |A∩B| / |A∪B|
#' - `"sorensen"`     : 2|A∩B| / (|A|+|B|)
#' - `"bray_curtis"`  : identical to Sørensen in the binary case
#' - `"kulczynski"`   : 0.5 * `a/(a+b) + a/(a+c)` with a=|A∩B|
#' - `"phi"` (≡ Pearson for binary)
#'
#' @param x             A `<phip_data>` object; **must be longitudinal**.
#' @param group_col     Character or NULL; grouping column in `x$data_long`.
#'                      If `NULL`, a synthetic group `"all"` is added.
#' @param time_col      Character; time column in `x$data_long`. May be numeric,
#'                      Date/POSIX (continuous) or factor/character (categorical).
#' @param similarity    One of: "jaccard","sorensen","bray_curtis","kulczynski",
#'                      "phi","pearson". ("pearson" is accepted and mapped to "phi".)
#' @param mode          Group scope: one of "within","between","all".
#' @param dyade_col     Optional character; column that encodes dyad/triad IDs.
#'                      Members of the same dyad/triad share the same ID.
#' @param time_mode     One of: "reference","pairwise".
#'   - "reference": compare each subject/time to a within-subject baseline
#'                 (earliest or latest) **within the same group**.
#'   - "pairwise" : only for **categorical** time; build subject×subject pairs
#'                 for each (time) and the selected group scope.
#' @param ref_time      If `time_mode=="reference"`, baseline selection:
#'                      one of "earliest","latest".
#' @param subject_mode  One of: "self_or_dyad" (default), "self", "all".
#'   - "self"         : only pairs with identical `subject_id`.
#'   - "self_or_dyad" : identical subjects OR same `dyade_col` (if provided).
#'   - "all"          : all subject pairs consistent with `mode` and time alignment.
#' @param drop_self_time Logical; for "reference" mode, drop baseline self rows
#'                      (time == baseline). Default TRUE.
#' @param out_table     Optional character; name for the materialized DuckDB table.
#'                      If NULL, a unique name is generated.
#' @param overwrite     Logical; if TRUE and `out_table` exists, it will be replaced.
#' @param verbose       Logical; phiper logging verbosity (default from `.ph_opt`).
#'
#' @return A lazy `tbl` pointing to the materialized DuckDB table on `x$meta$con`.
#'         Attributes include `table_name`, `time_mode`, `similarity`, and `mode`.
#' @export
compute_repertoire_similarity <- function(
    x,
    group_col,
    time_col,
    similarity = c("jaccard", "sorensen", "bray_curtis", "kulczynski", "phi", "pearson", "beta_sim"),
    mode = c("within", "between", "all"),
    dyade_col = NULL,
    time_mode = c("reference", "pairwise"),
    time_pairing = c("same", "cross"),
    ref_time = c("earliest", "latest"),
    subject_mode = c("self_or_dyad", "self", "all"),
    drop_self_time = TRUE,
    out_table = NULL,
    overwrite = FALSE,
    verbose = .ph_opt("verbose", TRUE)
) {
  .ph_with_timing(
    headline = "Computing repertoire similarity (<phip_data>)",
    step = sprintf(
      "args: similarity=%s; time_mode=%s; mode=%s; subject_mode=%s; ref_time=%s",
      paste(similarity, collapse = "/"), paste(time_mode, collapse = "/"),
      paste(mode, collapse = "/"), paste(subject_mode, collapse = "/"),
      paste(ref_time, collapse = "/")
    ),
    verbose = verbose,
    expr = {
      stopifnot(inherits(x, "phip_data"))

      if (!isTRUE(x$meta$longitudinal)) {
        .ph_abort(
          headline = "Longitudinal data required.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      if (is.null(time_col) || !nzchar(time_col)) {
        .ph_abort(
          headline = "Missing time column.",
          step     = "arguments check",
          bullets  = "time_col must be provided (character, non-empty)."
        )
      }

      similarity <- match.arg(similarity)
      # (no remap of pearson -> phi)

      mode <- match.arg(mode)
      time_mode <- match.arg(time_mode)
      ref_time <- match.arg(ref_time)
      subject_mode <- match.arg(subject_mode)
      time_pairing <- match.arg(time_pairing, several.ok = TRUE)

      con <- x$meta$con %||% NULL
      if (is.null(con)) {
        .ph_abort(
          headline = "No DuckDB connection found.",
          step     = "x$meta$con",
          bullets  = "x$meta$con is NULL; cannot materialize results."
        )
      }

      # -------------------------- columns & inputs -----------------------------
      need <- c("subject_id", "sample_id", "peptide_id", "exist", time_col)
      group_col_local <- if (is.null(group_col)) ".phip__group" else group_col
      if (!is.null(group_col)) need <- c(need, group_col)
      if (!is.null(dyade_col)) need <- c(need, dyade_col)

      miss <- setdiff(need, colnames(x$data_long))
      .chk_cond(
        length(miss) > 0,
        sprintf("Missing required columns in `data_long`: %s", paste(miss, collapse = ", "))
      )

      .data <- rlang::.data
      tcol <- rlang::sym(time_col)
      grp <- rlang::sym(group_col_local)

      # Synthetic group if group_col is NULL
      dl <- if (is.null(group_col)) {
        x$data_long |> dplyr::mutate(!!grp := "all")
      } else {
        x$data_long
      }

      # Only time is OK; but block between/all without a group column
      if (is.null(group_col) && mode %in% c("between", "all")) {
        .ph_abort(
          headline = "Between-group mode requires a group column.",
          step = "arguments check",
          bullets = c(
            "You provided only `time_col` (group_col = NULL).",
            "Use `mode = \"within\"` or supply a valid `group_col`."
          )
        )
      }

      # Detect time scale (categorical vs continuous) using dbplyr-friendly ops
      detect_time_continuous <- function(tbl, tcol) {
        tmp <- tryCatch(
          {
            tbl |>
              dplyr::filter(!is.na(!!tcol)) |>
              dplyr::distinct(!!tcol) |>
              dplyr::slice_min(order_by = !!tcol, n = 100L, with_ties = FALSE) |>
              dplyr::collect()
          },
          error = function(e) {
            .ph_abort(
              headline = "Failed to probe time column for type.",
              step = "detect_time_continuous()",
              bullets = c(
                sprintf("time_col: %s", rlang::as_name(tcol)),
                sprintf("backend: %s", class(dbplyr::remote_con(tbl))[1]),
                sprintf("error: %s", conditionMessage(e)),
                "Hint: ensure the time column is present and queryable."
              )
            )
          }
        )

        if (!nrow(tmp)) {
          return(FALSE)
        }
        v <- tmp[[rlang::as_name(tcol)]]
        is.numeric(v) || inherits(v, "Date") || inherits(v, "POSIXt")
      }
      is_cont_time <- detect_time_continuous(dl, tcol)

      if (identical(time_mode, "pairwise") && is_cont_time) {
        .ph_abort(
          headline = "Categorical time required for pairwise cross.",
          step = "time_mode check",
          bullets = c(
            sprintf("`%s` appears continuous (numeric/Date/POSIX).", time_col),
            "Full subject×subject cross on continuous time does not make sense.",
            "Use `time_mode = \"reference\"` or discretize time."
          )
        )
      }

      # Base presence view
      base_present <- dl |>
        dplyr::filter(.data$exist == 1L, !is.na(!!tcol)) |>
        dplyr::select(
          .data$subject_id, .data$sample_id, .data$peptide_id,
          !!grp, !!tcol,
          dplyr::all_of(dyade_col %||% character())
        )

      # One sample per (subject,group,time) - collapse to the smallest sample_id if dup
      timesamples <- base_present |>
        dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$sample_id)

      dup_count <- timesamples |>
        dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n") |>
        dplyr::filter(.data$n > 1L) |>
        dplyr::tally(name = "n_cells") |>
        dplyr::collect()
      if (nrow(dup_count) && dup_count$n_cells > 0L) {
        .ph_warn(sprintf(
          "Multiple samples per (subject,group,time) in %d cell(s). Using the smallest sample_id per cell.",
          dup_count$n_cells
        ))
        timesamples <- timesamples |>
          dplyr::group_by(.data$subject_id, !!grp, !!tcol) |>
          dplyr::summarise(sample_id = min(.data$sample_id), .groups = "drop")
      }

      # Universe sizes per sample (needed for phi/pearson)
      sample_sizes <- dl |>
        dplyr::group_by(.data$sample_id) |>
        dplyr::summarise(N = dplyr::n_distinct(.data$peptide_id), .groups = "drop")

      # FULL peptide universe size for binary correlations (φ / Pearson)
      libN <- x$data_long |>
        dplyr::distinct(.data$peptide_id) |>
        dplyr::tally(name = "N") |>
        dplyr::collect() |>
        dplyr::pull("N")

      # A tidy dyad map we can use later (ONLY used at the very end)
      has_dyads <- !is.null(dyade_col)
      if (has_dyads) {
        dy_sym <- rlang::sym(dyade_col)
        dy_map <- base_present |>
          dplyr::select(subject_id, !!grp, !!tcol, dyad_id = !!dy_sym) |>
          dplyr::distinct()
      }

      # Output table name & overwrite handling
      unique_name <- function(prefix = "phip_similarity") {
        ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
        rnd <- as.integer(runif(1, 1e5, 1e6))
        paste0(prefix, "_", ts, "_", rnd)
      }
      table_name <- out_table %||% unique_name()

      if (DBI::dbExistsTable(con, table_name)) {
        if (isTRUE(overwrite)) {
          DBI::dbRemoveTable(con, table_name)
        } else {
          .ph_abort(
            headline = "Output table already exists.",
            step     = "out_table / overwrite",
            bullets  = sprintf("Table `%s` exists. Set `overwrite = TRUE` or choose another `out_table`.", table_name)
          )
        }
      }

      # ------------------------------- ENGINE ---------------------------------

      if (identical(time_mode, "reference")) {
        if (mode %in% c("between", "all")) {
          .ph_warn("`time_mode = \"reference\"` compares within subject & group; `mode` ignored.")
        }

        baseline_tbl <- base_present |>
          dplyr::group_by(.data$subject_id, !!grp) |>
          dplyr::summarise(
            baseline_time = if (identical(ref_time, "earliest")) {
              min(!!tcol, na.rm = TRUE)
            } else {
              max(!!tcol, na.rm = TRUE)
            },
            .groups = "drop"
          )

        base_set <- base_present |>
          dplyr::inner_join(baseline_tbl, by = c("subject_id", group_col_local)) |>
          dplyr::filter(!!tcol == .data$baseline_time) |>
          dplyr::distinct(.data$subject_id, !!grp, .data$peptide_id)

        tp_set <- base_present |>
          dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id)

        bsize <- base_set |>
          dplyr::count(.data$subject_id, !!grp, name = "n0")
        tsize <- tp_set |>
          dplyr::count(.data$subject_id, !!grp, !!tcol, name = "nt")
        ints <- tp_set |>
          dplyr::inner_join(base_set, by = c("subject_id", group_col_local, "peptide_id")) |>
          dplyr::count(.data$subject_id, !!grp, !!tcol, name = "a")

        stab <- tsize |>
          dplyr::left_join(ints, by = c("subject_id", group_col_local, time_col)) |>
          dplyr::left_join(bsize, by = c("subject_id", group_col_local)) |>
          dplyr::mutate(
            a  = dplyr::coalesce(.data$a, 0L),
            nt = dplyr::coalesce(.data$nt, 0L),
            n0 = dplyr::coalesce(.data$n0, 0L)
          ) |>
          dplyr::left_join(baseline_tbl, by = c("subject_id", group_col_local))

        # Similarity calc
        stab <- switch(similarity,
                       jaccard = stab |>
                         dplyr::mutate(similarity = dplyr::if_else((.data$n0 + .data$nt - .data$a) > 0,
                                                                   .data$a / (.data$n0 + .data$nt - .data$a), NA_real_
                         )),
                       sorensen = stab |>
                         dplyr::mutate(similarity = dplyr::if_else((.data$n0 + .data$nt) > 0,
                                                                   (2 * .data$a) / (.data$n0 + .data$nt), NA_real_
                         )),
                       bray_curtis = stab |>
                         dplyr::mutate(similarity = dplyr::if_else((.data$n0 + .data$nt) > 0,
                                                                   (2 * .data$a) / (.data$n0 + .data$nt), NA_real_
                         )),
                       kulczynski = stab |>
                         dplyr::mutate(
                           b = pmax(.data$n0 - .data$a, 0L),
                           c = pmax(.data$nt - .data$a, 0L),
                           t1 = dplyr::if_else((.data$a + .data$b) > 0, .data$a / (.data$a + .data$b), NA_real_),
                           t2 = dplyr::if_else((.data$a + .data$c) > 0, .data$a / (.data$a + .data$c), NA_real_),
                           similarity = dplyr::coalesce(
                             0.5 * (t1 + t2),
                             dplyr::if_else((.data$a + .data$b + .data$c) == 0, 1.0, 0.0)
                           )
                         ) |>
                         dplyr::select(-b, -c, -t1, -t2),
                       beta_sim = stab |>
                         dplyr::mutate(
                           similarity = dplyr::if_else(pmin(.data$n0, .data$nt) > 0,
                                                       .data$a / pmin(.data$n0, .data$nt), NA_real_
                           )
                         ),
                       phi = {
                         pair_map <- timesamples |>
                           dplyr::inner_join(
                             baseline_tbl,
                             by = c("subject_id", group_col_local, time_col = "baseline_time")
                           ) |>
                           dplyr::rename(sample_id0 = .data$sample_id) |>
                           dplyr::select(.data$subject_id, !!grp, .data$sample_id0) |>
                           dplyr::inner_join(
                             timesamples |> dplyr::rename(sample_idt = .data$sample_id),
                             by = c("subject_id", group_col_local)
                           )

                         N_map <- pair_map |>
                           dplyr::left_join(sample_sizes, by = c("sample_idt" = "sample_id")) |>
                           dplyr::rename(Nt = .data$N) |>
                           dplyr::left_join(sample_sizes, by = c("sample_id0" = "sample_id")) |>
                           dplyr::rename(N0 = .data$N) |>
                           dplyr::mutate(N = pmin(.data$Nt, .data$N0, na.rm = TRUE)) |>
                           dplyr::select(.data$subject_id, !!grp, !!tcol, .data$N)

                         stab |>
                           dplyr::left_join(N_map, by = c("subject_id", group_col_local, time_col)) |>
                           dplyr::mutate(
                             N = !!libN, # use full universe
                             b = pmax(.data$n0 - .data$a, 0L),
                             c = pmax(.data$nt - .data$a, 0L),
                             n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                             n1d = .data$a + .data$b,
                             nd1 = .data$a + .data$c,
                             n0d = .data$N - .data$n1d,
                             nd0 = .data$N - .data$nd1,
                             den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                             num = (.data$a * .data$n00 - .data$b * .data$c),
                             similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                           ) |>
                           dplyr::select(-b, -c, -n00, -n1d, -nd1, -n0d, -nd0, -den, -num)
                       },
                       pearson = {
                         pair_map <- timesamples |>
                           dplyr::inner_join(
                             baseline_tbl,
                             by = c("subject_id", group_col_local, time_col = "baseline_time")
                           ) |>
                           dplyr::rename(sample_id0 = .data$sample_id) |>
                           dplyr::select(.data$subject_id, !!grp, .data$sample_id0) |>
                           dplyr::inner_join(
                             timesamples |> dplyr::rename(sample_idt = .data$sample_id),
                             by = c("subject_id", group_col_local)
                           )

                         N_map <- pair_map |>
                           dplyr::left_join(sample_sizes, by = c("sample_idt" = "sample_id")) |>
                           dplyr::rename(Nt = .data$N) |>
                           dplyr::left_join(sample_sizes, by = c("sample_id0" = "sample_id")) |>
                           dplyr::rename(N0 = .data$N) |>
                           dplyr::mutate(N = pmin(.data$Nt, .data$N0, na.rm = TRUE)) |>
                           dplyr::select(.data$subject_id, !!grp, !!tcol, .data$N)

                         stab |>
                           dplyr::left_join(N_map, by = c("subject_id", group_col_local, time_col)) |>
                           dplyr::mutate(
                             N = !!libN, # use full universe
                             b = pmax(.data$n0 - .data$a, 0L),
                             c = pmax(.data$nt - .data$a, 0L),
                             n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                             n1d = .data$a + .data$b,
                             nd1 = .data$a + .data$c,
                             n0d = .data$N - .data$n1d,
                             nd0 = .data$N - .data$nd1,
                             den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                             num = (.data$a * .data$n00 - .data$b * .data$c),
                             similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                           ) |>
                           dplyr::select(-b, -c, -n00, -n1d, -nd1, -n0d, -nd0, -den, -num)
                       }
        )

        out <- stab |>
          dplyr::select(.data$subject_id, !!grp, !!tcol, .data$baseline_time, .data$similarity) |>
          dplyr::rename(group = !!grp, time = !!tcol)

        if (isTRUE(drop_self_time)) {
          out <- out |> dplyr::filter(.data$time != .data$baseline_time)
        }

        # Materialize
        tryCatch(
          {
            out |> dplyr::compute(name = table_name, temporary = FALSE)
          },
          error = function(e) {
            .ph_abort(
              headline = "Failed to materialize output table.",
              step = "dplyr::compute()",
              bullets = c(
                sprintf("table: %s", table_name),
                sprintf("error: %s", conditionMessage(e))
              )
            )
          }
        )
      } else {
        # -------------------- PAIRWISE (categorical time matrices) -------------

        sets <- base_present |>
          dplyr::distinct(.data$subject_id, !!grp, !!tcol, .data$peptide_id, .data$sample_id)
        sizes <- sets |>
          dplyr::count(.data$subject_id, !!grp, !!tcol, name = "n")

        L0 <- sets |>
          dplyr::rename(
            subject_id1 = .data$subject_id, group1 = !!grp, time1 = !!tcol,
            peptide_id = .data$peptide_id, sample_id1 = .data$sample_id
          )
        R0 <- sets |>
          dplyr::rename(
            subject_id2 = .data$subject_id, group2 = !!grp, time2 = !!tcol,
            peptide_id = .data$peptide_id, sample_id2 = .data$sample_id
          )

        # helper to compute one arm (no dyad columns here; attach later)
        .compute_arm <- function(kind, L, R) {
          if (kind == "same") {
            pairs_raw <- L |>
              dplyr::rename(time = .data$time1) |>
              dplyr::inner_join(R |> dplyr::rename(time = .data$time2),
                                by = c("time", "peptide_id")
              )

            if (mode == "within") pairs_raw <- pairs_raw |> dplyr::filter(.data$group1 == .data$group2)
            if (mode == "between") pairs_raw <- pairs_raw |> dplyr::filter(.data$group1 != .data$group2)

            if (subject_mode == "self") {
              pairs_scoped <- pairs_raw |> dplyr::filter(.data$subject_id1 == .data$subject_id2)
            } else if (subject_mode == "self_or_dyad") {
              if (has_dyads) {
                pairs_scoped <- pairs_raw
              } else {
                pairs_scoped <- pairs_raw |> dplyr::filter(.data$subject_id1 == .data$subject_id2)
              }
            } else {
              pairs_scoped <- pairs_raw
            }

            # De-duplicate unordered pairs; keep diagonal for SAME when subject_mode=="all"
            keep_diag <- (subject_mode == "all")
            pairs_scoped <- pairs_scoped |>
              dplyr::mutate(
                key_left  = paste0(.data$group1, "|", .data$time, "|", .data$subject_id1),
                key_right = paste0(.data$group2, "|", .data$time, "|", .data$subject_id2)
              ) |>
              dplyr::filter(
                (.data$key_left < .data$key_right) |
                  (keep_diag &
                     .data$subject_id1 == .data$subject_id2 &
                     .data$group1 == .data$group2)
              ) |>
              dplyr::select(-.data$key_left, -.data$key_right)

            inters <- pairs_scoped |>
              dplyr::count(.data$time, .data$group1, .data$subject_id1,
                           .data$group2, .data$subject_id2,
                           name = "a"
              )

            pairs <- inters |>
              dplyr::left_join(
                sizes |> dplyr::rename(
                  group1 = !!grp, time = !!tcol,
                  subject_id1 = .data$subject_id, n1 = .data$n
                ),
                by = c("time", "group1", "subject_id1")
              ) |>
              dplyr::left_join(
                sizes |> dplyr::rename(
                  group2 = !!grp, time = !!tcol,
                  subject_id2 = .data$subject_id, n2 = .data$n
                ),
                by = c("time", "group2", "subject_id2")
              ) |>
              dplyr::mutate(
                a = dplyr::coalesce(.data$a, 0L),
                n1 = dplyr::coalesce(.data$n1, 0L),
                n2 = dplyr::coalesce(.data$n2, 0L),
                b = pmax(.data$n1 - .data$a, 0L),
                c = pmax(.data$n2 - .data$a, 0L)
              )

            pairs <- switch(similarity,
                            jaccard = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2 - .data$a) > 0,
                                                                                         .data$a / (.data$n1 + .data$n2 - .data$a),
                                                                                         NA_real_
                            )),
                            sorensen = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                                                                          (2 * .data$a) / (.data$n1 + .data$n2),
                                                                                          NA_real_
                            )),
                            bray_curtis = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                                                                             (2 * .data$a) / (.data$n1 + .data$n2),
                                                                                             NA_real_
                            )),
                            kulczynski = pairs |>
                              dplyr::mutate(
                                t1 = dplyr::if_else((.data$a + .data$b) > 0, .data$a / (.data$a + .data$b), NA_real_),
                                t2 = dplyr::if_else((.data$a + .data$c) > 0, .data$a / (.data$a + .data$c), NA_real_),
                                similarity = dplyr::coalesce(
                                  0.5 * (t1 + t2),
                                  dplyr::if_else((.data$a + .data$b + .data$c) == 0, 1.0, 0.0)
                                )
                              ) |> dplyr::select(-t1, -t2),
                            beta_sim = pairs |>
                              dplyr::mutate(
                                similarity = dplyr::if_else(pmin(.data$n1, .data$n2) > 0,
                                                            .data$a / pmin(.data$n1, .data$n2), NA_real_
                                )
                              ),
                            phi = {
                              sid_map <- timesamples |> dplyr::rename(group = !!grp, time = !!tcol)
                              pairs |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group1 = .data$group, subject_id1 = .data$subject_id,
                                    sample_id1 = .data$sample_id
                                  ),
                                  by = c("time", "group1", "subject_id1")
                                ) |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group2 = .data$group, subject_id2 = .data$subject_id,
                                    sample_id2 = .data$sample_id
                                  ),
                                  by = c("time", "group2", "subject_id2")
                                ) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
                                dplyr::rename(N1 = .data$N) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
                                dplyr::rename(N2 = .data$N) |>
                                dplyr::mutate(
                                  N = !!libN, # use full universe
                                  n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                                  n1d = .data$a + .data$b,
                                  nd1 = .data$a + .data$c,
                                  n0d = .data$N - .data$n1d,
                                  nd0 = .data$N - .data$nd1,
                                  den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                                  num = (.data$a * .data$n00 - .data$b * .data$c),
                                  similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                                ) |>
                                dplyr::select(-n00, -n1d, -nd1, -n0d, -nd0, -den, -N1, -N2, -num)
                            },
                            pearson = {
                              sid_map <- timesamples |> dplyr::rename(group = !!grp, time = !!tcol)
                              pairs |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group1 = .data$group, subject_id1 = .data$subject_id,
                                    sample_id1 = .data$sample_id
                                  ),
                                  by = c("time", "group1", "subject_id1")
                                ) |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group2 = .data$group, subject_id2 = .data$subject_id,
                                    sample_id2 = .data$sample_id
                                  ),
                                  by = c("time", "group2", "subject_id2")
                                ) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
                                dplyr::rename(N1 = .data$N) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
                                dplyr::rename(N2 = .data$N) |>
                                dplyr::mutate(
                                  N = !!libN, # use full universe
                                  n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                                  n1d = .data$a + .data$b,
                                  nd1 = .data$a + .data$c,
                                  n0d = .data$N - .data$n1d,
                                  nd0 = .data$N - .data$nd1,
                                  den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                                  num = (.data$a * .data$n00 - .data$b * .data$c),
                                  similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                                ) |>
                                dplyr::select(-n00, -n1d, -nd1, -n0d, -nd0, -den, -N1, -N2, -num)
                            }
            )

            pairs |>
              dplyr::transmute(
                time_left      = .data$time,
                group_left     = .data$group1,
                subject_left   = .data$subject_id1,
                time_right     = .data$time,
                group_right    = .data$group2,
                subject_right  = .data$subject_id2,
                similarity     = .data$similarity,
                pairing        = "same"
              )
          } else { # kind == "cross"
            pairs_raw <- dplyr::inner_join(L, R, by = "peptide_id") |>
              dplyr::filter(.data$time1 != .data$time2)

            if (mode == "within") pairs_raw <- pairs_raw |> dplyr::filter(.data$group1 == .data$group2)
            if (mode == "between") pairs_raw <- pairs_raw |> dplyr::filter(.data$group1 != .data$group2)

            if (subject_mode == "self") {
              pairs_scoped <- pairs_raw |> dplyr::filter(.data$subject_id1 == .data$subject_id2)
            } else if (subject_mode == "self_or_dyad") {
              if (has_dyads) {
                pairs_scoped <- pairs_raw
              } else {
                pairs_scoped <- pairs_raw |> dplyr::filter(.data$subject_id1 == .data$subject_id2)
              }
            } else {
              pairs_scoped <- pairs_raw
            }

            pairs_scoped <- pairs_scoped |>
              dplyr::mutate(
                key_left  = paste0(.data$group1, "|", .data$time1, "|", .data$subject_id1),
                key_right = paste0(.data$group2, "|", .data$time2, "|", .data$subject_id2)
              ) |>
              dplyr::filter(.data$key_left < .data$key_right) |>
              dplyr::select(-.data$key_left, -.data$key_right)

            inters <- pairs_scoped |>
              dplyr::count(.data$time1, .data$group1, .data$subject_id1,
                           .data$time2, .data$group2, .data$subject_id2,
                           name = "a"
              )

            pairs <- inters |>
              dplyr::left_join(
                sizes |> dplyr::rename(
                  group1 = !!grp, time1 = !!tcol,
                  subject_id1 = .data$subject_id, n1 = .data$n
                ),
                by = c("time1", "group1", "subject_id1")
              ) |>
              dplyr::left_join(
                sizes |> dplyr::rename(
                  group2 = !!grp, time2 = !!tcol,
                  subject_id2 = .data$subject_id, n2 = .data$n
                ),
                by = c("time2", "group2", "subject_id2")
              ) |>
              dplyr::mutate(
                a = dplyr::coalesce(.data$a, 0L),
                n1 = dplyr::coalesce(.data$n1, 0L),
                n2 = dplyr::coalesce(.data$n2, 0L),
                b = pmax(.data$n1 - .data$a, 0L),
                c = pmax(.data$n2 - .data$a, 0L)
              )

            pairs <- switch(similarity,
                            jaccard = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2 - .data$a) > 0,
                                                                                         .data$a / (.data$n1 + .data$n2 - .data$a),
                                                                                         NA_real_
                            )),
                            sorensen = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                                                                          (2 * .data$a) / (.data$n1 + .data$n2),
                                                                                          NA_real_
                            )),
                            bray_curtis = pairs |> dplyr::mutate(similarity = dplyr::if_else((.data$n1 + .data$n2) > 0,
                                                                                             (2 * .data$a) / (.data$n1 + .data$n2),
                                                                                             NA_real_
                            )),
                            kulczynski = pairs |>
                              dplyr::mutate(
                                t1 = dplyr::if_else((.data$a + .data$b) > 0, .data$a / (.data$a + .data$b), NA_real_),
                                t2 = dplyr::if_else((.data$a + .data$c) > 0, .data$a / (.data$a + .data$c), NA_real_),
                                similarity = dplyr::coalesce(
                                  0.5 * (t1 + t2),
                                  dplyr::if_else((.data$a + .data$b + .data$c) == 0, 1.0, 0.0)
                                )
                              ) |> dplyr::select(-t1, -t2),
                            beta_sim = pairs |>
                              dplyr::mutate(
                                similarity = dplyr::if_else(pmin(.data$n1, .data$n2) > 0,
                                                            .data$a / pmin(.data$n1, .data$n2), NA_real_
                                )
                              ),
                            phi = {
                              sid_map <- timesamples |> dplyr::rename(group = !!grp, time = !!tcol)
                              pairs |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group1 = .data$group, time1 = .data$time,
                                    subject_id1 = .data$subject_id, sample_id1 = .data$sample_id
                                  ),
                                  by = c("time1", "group1", "subject_id1")
                                ) |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group2 = .data$group, time2 = .data$time,
                                    subject_id2 = .data$subject_id, sample_id2 = .data$sample_id
                                  ),
                                  by = c("time2", "group2", "subject_id2")
                                ) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
                                dplyr::rename(N1 = .data$N) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
                                dplyr::rename(N2 = .data$N) |>
                                dplyr::mutate(
                                  N = !!libN, # use full universe
                                  n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                                  n1d = .data$a + .data$b,
                                  nd1 = .data$a + .data$c,
                                  n0d = .data$N - .data$n1d,
                                  nd0 = .data$N - .data$nd1,
                                  den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                                  num = (.data$a * .data$n00 - .data$b * .data$c),
                                  similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                                ) |>
                                dplyr::select(-n00, -n1d, -nd1, -n0d, -nd0, -den, -N1, -N2, -num)
                            },
                            pearson = {
                              sid_map <- timesamples |> dplyr::rename(group = !!grp, time = !!tcol)
                              pairs |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group1 = .data$group, time1 = .data$time,
                                    subject_id1 = .data$subject_id, sample_id1 = .data$sample_id
                                  ),
                                  by = c("time1", "group1", "subject_id1")
                                ) |>
                                dplyr::left_join(
                                  sid_map |> dplyr::rename(
                                    group2 = .data$group, time2 = .data$time,
                                    subject_id2 = .data$subject_id, sample_id2 = .data$sample_id
                                  ),
                                  by = c("time2", "group2", "subject_id2")
                                ) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id1" = "sample_id")) |>
                                dplyr::rename(N1 = .data$N) |>
                                dplyr::left_join(sample_sizes, by = c("sample_id2" = "sample_id")) |>
                                dplyr::rename(N2 = .data$N) |>
                                dplyr::mutate(
                                  N = !!libN, # use full universe
                                  n00 = pmax(.data$N - (.data$a + .data$b + .data$c), 0),
                                  n1d = .data$a + .data$b,
                                  nd1 = .data$a + .data$c,
                                  n0d = .data$N - .data$n1d,
                                  nd0 = .data$N - .data$nd1,
                                  den = sqrt(pmax(.data$n1d, 0) * pmax(.data$n0d, 0) * pmax(.data$nd1, 0) * pmax(.data$nd0, 0)),
                                  num = (.data$a * .data$n00 - .data$b * .data$c),
                                  similarity = dplyr::if_else(.data$den > 0, .data$num / .data$den, NA_real_)
                                ) |>
                                dplyr::select(-n00, -n1d, -nd1, -n0d, -nd0, -den, -N1, -N2, -num)
                            }
            )

            pairs |>
              dplyr::transmute(
                time_left      = .data$time1,
                group_left     = .data$group1,
                subject_left   = .data$subject_id1,
                time_right     = .data$time2,
                group_right    = .data$group2,
                subject_right  = .data$subject_id2,
                similarity     = .data$similarity,
                pairing        = "cross"
              )
          }
        }

        # Build requested arms and UNION them (schema identical)
        arms <- list()
        if ("same" %in% time_pairing) arms <- append(arms, list(.compute_arm("same", L0, R0)))
        if ("cross" %in% time_pairing) arms <- append(arms, list(.compute_arm("cross", L0, R0)))

        .chk_cond(length(arms) == 0L,
                  "`time_pairing` produced no arms.",
                  step = "pairwise assembly",
                  bullets = "Pass at least one of: 'same', 'cross'."
        )

        out <- if (length(arms) == 1L) arms[[1]] else Reduce(dplyr::union_all, arms)

        # subject_mode == "self_or_dyad": if dyads are present, keep self OR same-dyad pairs
        if (has_dyads && subject_mode == "self_or_dyad") {
          out <- out |>
            dplyr::left_join(
              dy_map |> dplyr::rename(dyad_left = .data$dyad_id),
              by = c("subject_left" = "subject_id", "group_left" = group_col_local, "time_left" = time_col)
            ) |>
            dplyr::left_join(
              dy_map |> dplyr::rename(dyad_right = .data$dyad_id),
              by = c("subject_right" = "subject_id", "group_right" = group_col_local, "time_right" = time_col)
            ) |>
            dplyr::filter((.data$subject_left == .data$subject_right) |
                            (!is.na(.data$dyad_left) & .data$dyad_left == .data$dyad_right))
        }

        # Always attach dyad columns when dyads exist (no COALESCE to NULL!)
        if (has_dyads) {
          out <- out |>
            dplyr::left_join(
              dy_map |> dplyr::rename(dyad_left = .data$dyad_id),
              by = c("subject_left" = "subject_id", "group_left" = group_col_local, "time_left" = time_col)
            ) |>
            dplyr::left_join(
              dy_map |> dplyr::rename(dyad_right = .data$dyad_id),
              by = c("subject_right" = "subject_id", "group_right" = group_col_local, "time_right" = time_col)
            )
        } else {
          out <- out |>
            dplyr::mutate(dyad_left = NA_character_, dyad_right = NA_character_)
        }

        # Back-compat: if only SAME requested, add alias column
        if (length(time_pairing) == 1L && identical(time_pairing, "same")) {
          out <- out |> dplyr::mutate(timepoint = .data$time_left, .after = .data$time_left)
        }

        # Strongly type character-ish columns to TEXT for DuckDB union safety
        out <- out |>
          dplyr::mutate(
            dplyr::across(
              tidyselect::any_of(c(
                "time_left", "time_right", "group_left", "group_right",
                "subject_left", "subject_right", "dyad_left", "dyad_right", "pairing"
              )),
              as.character
            )
          )

        # Materialize
        tryCatch(
          {
            out |> dplyr::compute(name = table_name, temporary = FALSE)
          },
          error = function(e) {
            .ph_abort(
              headline = "Failed to materialize output table.",
              step = "dplyr::compute()",
              bullets = c(
                sprintf("table: %s", table_name),
                sprintf("error: %s", conditionMessage(e))
              )
            )
          }
        )
      }

      # Return a lazy handle to the materialized table
      out_lazy <- dplyr::tbl(con, table_name)
      attr(out_lazy, "table_name") <- table_name
      attr(out_lazy, "time_mode") <- time_mode
      attr(out_lazy, "time_pairing") <- time_pairing
      attr(out_lazy, "similarity") <- similarity
      attr(out_lazy, "mode") <- mode
      out_lazy
    }
  )
}

#' @title Plot similarity heatmap for exactly two (group × time) combinations
#'
#' @description
#' Build a subject-by-subject similarity **heatmap** from the DuckDB table
#' created by `compute_repertoire_similarity()` with `time_mode = "pairwise"`
#' and `subject_mode = "all"`.
#'
#' **Hard constraints:**
#' - Only **two** distinct (group × time) combos are allowed.
#' - A **full** matrix is required (every left-subject × right-subject present).
#' - If no group was used, you must provide **two times**.
#'
#' The function validates inputs and aborts early with PHIP-style messages.
#' If dyad columns are present (`dyad_left`, `dyad_right`), you can:
#'   - restrict to dyad-intersection to make the matrix square (`dyads = "only"`),
#'   - just order axes by dyad but keep all pairs (`dyads = "order"`),
#'   - ignore dyads (`dyads = "no"`).
#'
#'
#' @param sim_tbl Lazy DuckDB table from `compute_repertoire_similarity()`.
#'   New schema columns: `time_left`, `group_left`, `subject_left`,
#'   `time_right`, `group_right`, `subject_right`, `similarity`, `pairing`.
#'   Legacy schema columns (auto-adapted): `timepoint`, `group_left`,
#'   `subject_left`, `group_right`, `subject_right`, `similarity`.
#' @param groups Character(2) or NULL. If NULL, there must be **one** group in
#'   the table and you must supply `times = c(t1, t2)`.
#' @param times   Character(2) or NULL. If `groups` is given, `times[1]` pairs
#'   with `groups[1]`, `times[2]` with `groups[2]`. If `groups` is NULL, both
#'   times are used within the single group in the table.
#' @param pairing One of `"cross","same"`. `"same"` = same timepoint on both
#'   sides; `"cross"` = different timepoints.
#' @param cluster One of `"none","rows","cols","both"`. Hierarchical clustering
#'   to reorder axes (performed locally after collecting).
#' @param limits  Length-2 numeric for fill scale limits. Default `c(0,1)`.
#' @param na_fill Colour for NA cells. Default `"#F0F0F0"`.
#' @param verbose Logical; pass-through to PHIP logger.
#'
#' @return A `ggplot2` object (modifiable).
#' @export
plot_similarity_heatmap <- function(sim_tbl,
                                    groups = NULL,
                                    times = NULL,
                                    pairing = c("cross", "same"),
                                    cluster = c("none", "rows", "cols", "both"),
                                    dyads = c("only", "order", "no"),
                                    grid = c("union_pad", "require_full", "intersection", "observed_only"),
                                    dedupe = c("mean", "first", "median", "max", "min"),
                                    limits = c(0, 1),
                                    na_fill = "#F0F0F0",
                                    verbose = .ph_opt("verbose", TRUE)) {
  pairing <- match.arg(pairing)
  cluster <- match.arg(cluster)
  dyads <- match.arg(dyads)
  grid <- match.arg(grid)
  dedupe <- match.arg(dedupe)
  .data <- rlang::.data

  # --- schema detection --------------------------------------------------------
  have <- tryCatch(colnames(sim_tbl), error = function(e) character(0))
  if (!length(have)) have <- tryCatch(dplyr::tbl_vars(sim_tbl), error = function(e) character(0))

  req_new <- c(
    "time_left", "group_left", "subject_left",
    "time_right", "group_right", "subject_right",
    "similarity", "pairing"
  )
  req_legacy <- c(
    "timepoint", "group_left", "subject_left",
    "group_right", "subject_right", "similarity"
  )

  using_new <- all(req_new %in% have)
  using_legacy <- !using_new && all(req_legacy %in% have)

  if (!using_new && !using_legacy) {
    avail <- if (length(have)) paste(have, collapse = ", ") else "<none>"
    .ph_abort(
      headline = "Similarity table has an unexpected schema.",
      step = "column check",
      bullets = c(
        sprintf("required (new):   %s", paste(req_new, collapse = ", ")),
        sprintf("required (legacy):%s", paste(req_legacy, collapse = ", ")),
        sprintf("available:        %s", avail)
      )
    )
  }

  if (using_legacy) {
    .ph_log_info("Detected legacy similarity schema; adapting to new columns.",
                 step = "compat mode", verbose = verbose
    )
    sim_tbl <- sim_tbl |>
      dplyr::mutate(
        time_left = .data$timepoint,
        time_right = .data$timepoint,
        pairing = "same"
      ) |>
      dplyr::select(
        .data$time_left, .data$group_left, .data$subject_left,
        .data$time_right, .data$group_right, .data$subject_right,
        .data$similarity, .data$pairing
      )
    have <- colnames(sim_tbl)
  }

  has_dy <- all(c("dyad_left", "dyad_right") %in% have)

  if (dyads != "no" && !has_dy) {
    .ph_warn(
      headline = "Dyad ordering requested but dyad columns are missing.",
      step = "plot_similarity_heatmap()",
      bullets = c(
        "Expected: dyad_left, dyad_right.",
        "Recompute with `dyade_col` or set `dyads = 'no'`."
      )
    )
    dyads <- "no"
  }
  if (dyads != "no" && cluster != "none") {
    .ph_warn(
      headline = "Clustering disabled to preserve dyad diagonal.",
      step     = "plot_similarity_heatmap()",
      bullets  = sprintf("Requested `cluster = %s`, applying `cluster = 'none'`.", cluster)
    )
    cluster <- "none"
  }

  # --- resolve requested (group,time) combos ----------------------------------
  ddistinct <- function(tbl, ...) dplyr::distinct(tbl, ..., .keep_all = FALSE)

  gt <- ddistinct(sim_tbl, .data$group_left, .data$time_left) |>
    dplyr::rename(group = .data$group_left, time = .data$time_left) |>
    dplyr::union_all(
      ddistinct(sim_tbl, .data$group_right, .data$time_right) |>
        dplyr::rename(group = .data$group_right, time = .data$time_right)
    ) |>
    ddistinct(.data$group, .data$time) |>
    dplyr::collect()

  groups_all <- sort(unique(gt$group))
  times_all <- sort(unique(gt$time))

  if (is.null(groups)) {
    .chk_cond(length(groups_all) != 1L,
              "Argument `groups` is NULL but multiple groups present.",
              step = "input resolution",
              bullets = c(
                sprintf(
                  "unique groups: %s",
                  paste(head(groups_all, 10), collapse = ", ")
                ),
                "Supply `groups` or filter the table."
              )
    )
    .chk_cond(is.null(times) || length(times) != 2L,
              "`times` must be length 2 when `groups` is NULL.",
              step = "input resolution",
              bullets = sprintf(
                "available times: %s",
                paste(head(times_all, 20), collapse = ", ")
              )
    )
    left_group <- groups_all[1]
    right_group <- groups_all[1]
    left_time <- times[1]
    right_time <- times[2]
  } else {
    .chk_cond(length(groups) != 2L,
              "`groups` must be character(2).",
              step    = "input resolution",
              bullets = sprintf("got length: %d", length(groups))
    )
    left_group <- groups[1]
    right_group <- groups[2]
    .chk_cond(is.null(times) || length(times) != 2L,
              "`times` must be character(2).",
              step = "input resolution"
    )
    left_time <- times[1]
    right_time <- times[2]
  }

  if (pairing == "cross" && identical(left_time, right_time)) {
    .ph_abort("`pairing = 'cross'` requires two different timepoints.", step = "input validation")
  }
  if (pairing == "same" && !identical(left_time, right_time)) {
    .ph_abort("`pairing = 'same'` requires the same timepoint on both sides.", step = "input validation")
  }

  # --- filter + normalize orientation -----------------------------------------
  sub <- sim_tbl |>
    dplyr::filter(.data$pairing == pairing) |>
    dplyr::filter(
      (.data$group_left == left_group & .data$time_left == left_time &
         .data$group_right == right_group & .data$time_right == right_time) |
        (.data$group_left == right_group & .data$time_left == right_time &
           .data$group_right == left_group & .data$time_right == left_time)
    ) |>
    dplyr::collect()

  .chk_cond(!nrow(sub), "No rows match the requested combinations.", step = "data filter")

  swap <- !(sub$group_left == left_group & sub$time_left == left_time)
  df <- data.frame(
    time_left = ifelse(swap, sub$time_right, sub$time_left),
    group_left = ifelse(swap, sub$group_right, sub$group_left),
    subject_left = ifelse(swap, sub$subject_right, sub$subject_left),
    dyad_left = if (has_dy) ifelse(swap, sub$dyad_right, sub$dyad_left) else NA_character_,
    time_right = ifelse(swap, sub$time_left, sub$time_right),
    group_right = ifelse(swap, sub$group_left, sub$group_right),
    subject_right = ifelse(swap, sub$subject_left, sub$subject_right),
    dyad_right = if (has_dy) ifelse(swap, sub$dyad_left, sub$dyad_right) else NA_character_,
    similarity = sub$similarity,
    stringsAsFactors = FALSE
  )

  # --- optional dedupe when dyads == "no" -------------------------------------
  if (dyads == "no" && has_dy) {
    agg <- switch(dedupe,
                  mean   = function(z) mean(z, na.rm = TRUE),
                  median = function(z) stats::median(z, na.rm = TRUE),
                  max    = function(z) max(z, na.rm = TRUE),
                  min    = function(z) min(z, na.rm = TRUE),
                  first  = function(z) z[match(TRUE, !is.na(z), nomatch = 1L)]
    )
    df <- df |>
      dplyr::group_by(
        subject_left, subject_right,
        time_left, group_left, time_right, group_right
      ) |>
      dplyr::summarise(similarity = agg(similarity), .groups = "drop")
  }

  axis_title <- sprintf(
    "%s @ %s  vs  %s @ %s",
    df$group_left[1], df$time_left[1],
    df$group_right[1], df$time_right[1]
  )

  # --- helper: duplicate missing reversed pairs (in axis space) ----------------
  symmetrize_pairs <- function(df_axis, same_case) {
    if (!same_case) {
      return(df_axis)
    }
    # work on character to avoid factor level traps
    keys <- df_axis |>
      dplyr::transmute(
        sl = as.character(subject_left),
        sr = as.character(subject_right)
      )
    swap <- df_axis |>
      dplyr::transmute(
        sl = as.character(subject_right),
        sr = as.character(subject_left),
        similarity = similarity
      )
    add <- dplyr::anti_join(swap, keys, by = c("sl" = "sl", "sr" = "sr"))
    if (!nrow(add)) {
      return(df_axis)
    }
    add <- add |>
      dplyr::transmute(
        subject_left  = factor(sl, levels = levels(df_axis$subject_left)),
        subject_right = factor(sr, levels = levels(df_axis$subject_right)),
        similarity    = similarity,
        time_left     = df_axis$time_left[1],
        group_left    = df_axis$group_left[1],
        time_right    = df_axis$time_right[1],
        group_right   = df_axis$group_right[1]
      )
    dplyr::bind_rows(df_axis, dplyr::select(add, colnames(df_axis)))
  }

  same_case <- (pairing == "same" &&
                  identical(left_group, right_group) &&
                  identical(left_time, right_time))

  # ======================= DYAD-AWARE BRANCH ===================================
  if (dyads != "no" && has_dy) {
    left_dy <- unique(df$dyad_left[!is.na(df$dyad_left)])
    right_dy <- unique(df$dyad_right[!is.na(df$dyad_right)])
    common_dy <- sort(intersect(left_dy, right_dy))
    .chk_cond(!length(common_dy), "No common dyad IDs.", step = "dyad intersection")

    if (dyads == "only") {
      df <- df[df$dyad_left %in% common_dy & df$dyad_right %in% common_dy &
                 df$dyad_left == df$dyad_right, , drop = FALSE]
    } else {
      df <- df[df$dyad_left %in% common_dy & df$dyad_right %in% common_dy, , drop = FALSE]
    }

    Llist <- split(
      unique(df[, c("dyad_left", "subject_left")])$subject_left,
      unique(df[, c("dyad_left", "subject_left")])$dyad_left
    )
    Rlist <- split(
      unique(df[, c("dyad_right", "subject_right")])$subject_right,
      unique(df[, c("dyad_right", "subject_right")])$dyad_right
    )

    .build_map <- function(subj_vec, K) {
      subj_vec <- sort(unique(subj_vec))
      if (length(subj_vec) == 0L) {
        return(data.frame(subject = character(0), idx = integer(0), axis = character(0)))
      }
      take <- subj_vec[(seq_len(K) - 1) %% length(subj_vec) + 1L]
      dup_rank <- ave(take, take, FUN = seq_along)
      dup_n <- ave(take, take, FUN = length)
      axis <- ifelse(dup_n > 1L, paste0(take, "\u00B7", dup_rank), take)
      data.frame(subject = take, idx = seq_len(K), axis = axis, stringsAsFactors = FALSE)
    }

    left_map <- list()
    right_map <- list()
    for (d in common_dy) {
      nL <- length(Llist[[d]] %||% character(0))
      nR <- length(Rlist[[d]] %||% character(0))
      K <- max(nL, nR)
      if (K == 0L) next
      Lmap <- .build_map(Llist[[d]], K)
      Lmap$dyad <- d
      Rmap <- .build_map(Rlist[[d]], K)
      Rmap$dyad <- d
      left_map[[d]] <- Lmap
      right_map[[d]] <- Rmap
    }
    left_map <- do.call(rbind, left_map)
    right_map <- do.call(rbind, right_map)

    df <- df |>
      dplyr::inner_join(left_map,
                        by = c("dyad_left" = "dyad", "subject_left" = "subject")
      ) |>
      dplyr::rename(idx_left = .data$idx, axis_left = .data$axis) |>
      dplyr::inner_join(right_map,
                        by = c("dyad_right" = "dyad", "subject_right" = "subject")
      ) |>
      dplyr::rename(idx_right = .data$idx, axis_right = .data$axis)

    y_levels <- left_map |>
      dplyr::arrange(.data$dyad, .data$idx) |>
      dplyr::pull(.data$axis) |>
      unique()
    x_levels <- right_map |>
      dplyr::arrange(.data$dyad, .data$idx) |>
      dplyr::pull(.data$axis) |>
      unique()

    df$subject_left <- factor(df$axis_left, levels = y_levels)
    df$subject_right <- factor(df$axis_right, levels = x_levels)

    # >>> symmetry BEFORE padding (duplicates missing reversed pairs)
    df <- symmetrize_pairs(df, same_case)

    present <- dplyr::distinct(df, .data$subject_left, .data$subject_right)
    want_n <- length(y_levels) * length(x_levels)

    if (grid == "require_full" && nrow(present) < want_n) {
      .ph_abort("Full matrix required and data are incomplete after dyad padding.", step = "matrix validation")
    } else if (grid %in% c("union_pad", "intersection")) {
      base <- df[1, c("time_left", "group_left", "time_right", "group_right")]
      df <- tidyr::complete(df,
                            .data$subject_left, .data$subject_right,
                            fill = list(similarity = NA_real_)
      ) |>
        dplyr::mutate(
          time_left  = base$time_left,  group_left  = base$group_left,
          time_right = base$time_right, group_right = base$group_right
        )
    }
    # ======================= NO-DYAD BRANCH ======================================
  } else {
    left_subs <- sort(unique(df$subject_left))
    right_subs <- sort(unique(df$subject_right))

    if (cluster != "none") {
      mat_df <- tidyr::pivot_wider(
        df,
        id_cols = .data$subject_left,
        names_from = .data$subject_right, values_from = .data$similarity
      )
      rn <- mat_df$subject_left
      mat <- as.matrix(mat_df[, setdiff(names(mat_df), "subject_left"), drop = FALSE])
      rownames(mat) <- rn
      mat[is.na(mat)] <- 0
      if (cluster %in% c("rows", "both")) {
        ord_r <- stats::hclust(stats::dist(mat))$order
        left_subs <- rn[ord_r]
      }
      if (cluster %in% c("cols", "both")) {
        ord_c <- stats::hclust(stats::dist(t(mat)))$order
        right_subs <- colnames(mat)[ord_c]
      }
    }

    if (grid == "intersection") {
      common <- sort(intersect(left_subs, right_subs))
      df <- df[df$subject_left %in% common & df$subject_right %in% common, , drop = FALSE]
      left_subs <- common
      right_subs <- common
    }

    df$subject_left <- factor(df$subject_left, levels = left_subs)
    df$subject_right <- factor(df$subject_right, levels = right_subs)

    # >>> symmetry BEFORE padding (duplicates missing reversed pairs)
    df <- symmetrize_pairs(df, same_case)

    want_n <- length(left_subs) * length(right_subs)
    have_n <- nrow(dplyr::distinct(df, .data$subject_left, .data$subject_right))

    if (grid == "require_full" && have_n < want_n) {
      .ph_abort("Full matrix required and data are incomplete.", step = "matrix validation")
    } else if (grid %in% c("union_pad", "intersection")) {
      base <- df[1, c("time_left", "group_left", "time_right", "group_right")]
      df <- tidyr::complete(df,
                            .data$subject_left, .data$subject_right,
                            fill = list(similarity = NA_real_)
      ) |>
        dplyr::mutate(
          time_left  = base$time_left,  group_left  = base$group_left,
          time_right = base$time_right, group_right = base$group_right
        )
      padded <- want_n - have_n
      if (padded > 0L) {
        .ph_log_info(sprintf("Grid=%s: padded %d tiles with NA.", grid, padded),
                     step = "grid padding", verbose = verbose
        )
      }
    }
  }

  # --- draw --------------------------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$subject_right, y = .data$subject_left, fill = .data$similarity)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradientn(
      colours = c("#F7F7F7", "#A6CEE3", "#1F78B4"),
      limits = limits,
      na.value = na_fill,
      name = "Similarity"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = axis_title,
      x = sprintf("%s @ %s", df$group_right[1], df$time_right[1]),
      y = sprintf("%s @ %s", df$group_left[1], df$time_left[1])
    ) +
    theme_phip() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  .ph_log_ok("Heatmap ready",
             step = sprintf(
               "|L|=%d, |R|=%d",
               nlevels(df$subject_left), nlevels(df$subject_right)
             ),
             verbose = verbose
  )

  # --- attach ordered data used for the heatmap (for CSV export) ---------------
  export_df <- df %>%
    dplyr::select(
      subject_left, subject_right, similarity,
      time_left, group_left, time_right, group_right
    )

  attr(p, "heatmap_df") <- export_df
  attr(p, "axis_levels") <- list(
    rows = levels(df$subject_left),
    cols = levels(df$subject_right)
  )

  p
}
