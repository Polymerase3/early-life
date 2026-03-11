# load required packages
library(phiper)
library(dplyr)
library(tidyr)

set.seed(632961)

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

ps$data_long <- ps$data_long %>%
  dplyr::mutate_if(is.logical, as.integer)

# =============================================================================
# Crosstables: unique subject_id counts
# =============================================================================
count_distinct_subjects <- function(tbl, cols) {
  if (length(cols) == 0) {
    return(tibble::tibble(label = character(), n_subjects = integer()))
  }

  tbl %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(cols),
        ~ dplyr::n_distinct(
          dplyr::if_else(.x == 1L, as.character(subject_id), NA_character_),
          na.rm = TRUE
        )
      )
    ) %>%
    dplyr::collect() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "label", values_to = "n_subjects")
}

tp_levels <- paste0("T", 0:8)
all_cols <- colnames(ps$data_long)

# -----------------------------------------------------------------------------
# 1) sample group crosstables (mom_serum, kid_serum, mom_milk)
# -----------------------------------------------------------------------------
sample_cols <- all_cols[grepl("^(mom_serum|kid_serum|mom_milk)_T\\d+$", all_cols)]
sample_counts_long <- count_distinct_subjects(ps$data_long, sample_cols) %>%
  tidyr::extract(
    "label",
    into = c("group", "timepoint"),
    regex = "^(.*)_(T\\d+)$",
    remove = FALSE
  ) %>%
  dplyr::mutate(timepoint = factor(.data$timepoint, levels = tp_levels))

sample_tables <- split(sample_counts_long, sample_counts_long$group)
sample_tables <- lapply(sample_tables, function(df) {
  df %>%
    dplyr::select(.data$timepoint, .data$n_subjects) %>%
    tidyr::pivot_wider(
      names_from = "timepoint",
      values_from = "n_subjects",
      values_fill = 0
    )
})

# -----------------------------------------------------------------------------
# 2) metadata crosstables (one table per metadata variable)
# -----------------------------------------------------------------------------
meta_cols <- all_cols[grepl("^(kid_serum|mom_kid_serum|mom_milk)_T\\d+_.+", all_cols)]
meta_counts_long <- count_distinct_subjects(ps$data_long, meta_cols) %>%
  tidyr::extract(
    "label",
    into = c("group", "timepoint", "rest"),
    regex = "^(.*)_(T\\d+)_(.*)$",
    remove = FALSE
  )

rest_set <- unique(meta_counts_long$rest)
level_priority <- c("yes", "no", "female", "male", "BF", "MF", "VG", "CS",
  "home", "hospital", "before", "after")

parse_meta_rest <- function(rest, rest_set) {
  # 1) standard suffix levels (e.g., smoking_yes, sex_female)
  if (grepl("_", rest)) {
    last_tok <- sub("^.*_([^_]+)$", "\\1", rest)
    if (last_tok %in% level_priority) {
      return(list(
        meta = sub("_[^_]+$", "", rest),
        level = last_tok
      ))
    }
  }

  # 2) no_* prefix (e.g., no_siblings)
  if (startsWith(rest, "no_")) {
    return(list(
      meta = sub("^no_", "", rest),
      level = "no"
    ))
  }

  # 3) if a no_* counterpart exists, treat this as "yes"
  if (paste0("no_", rest) %in% rest_set) {
    return(list(
      meta = rest,
      level = "yes"
    ))
  }

  # 4) fallback: split last token if present
  if (grepl("_", rest)) {
    return(list(
      meta = sub("_[^_]+$", "", rest),
      level = sub("^.*_", "", rest)
    ))
  }

  # 5) final fallback: single-level
  list(meta = rest, level = "yes")
}

if (nrow(meta_counts_long) > 0) {
  parsed <- lapply(meta_counts_long$rest, parse_meta_rest, rest_set = rest_set)
  meta_counts_long <- meta_counts_long %>%
    dplyr::mutate(
      meta = vapply(parsed, `[[`, character(1), "meta"),
      level = vapply(parsed, `[[`, character(1), "level"),
      timepoint = factor(.data$timepoint, levels = tp_levels)
    )

  # If multiple groups exist, keep them separated in table names.
  group_mult <- length(unique(meta_counts_long$group)) > 1
  meta_counts_long <- meta_counts_long %>%
    dplyr::mutate(
      meta_key = if (group_mult) paste(.data$group, .data$meta, sep = "__") else .data$meta
    )

  meta_tables <- split(meta_counts_long, meta_counts_long$meta_key)
  meta_tables <- lapply(meta_tables, function(df) {
    levels_here <- unique(df$level)
    lvl_order <- c(level_priority, setdiff(levels_here, level_priority))
    df %>%
      dplyr::mutate(level = factor(.data$level, levels = lvl_order)) %>%
      dplyr::select(.data$level, .data$timepoint, .data$n_subjects) %>%
      tidyr::pivot_wider(
        names_from = "timepoint",
        values_from = "n_subjects",
        values_fill = 0
      ) %>%
      dplyr::arrange(.data$level)
  })
} else {
  meta_tables <- list()
}

# Bundle all crosstables into a single list
crosstables <- c(sample_tables, meta_tables)

crosstables
