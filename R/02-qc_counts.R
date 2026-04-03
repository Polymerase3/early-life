# load required packages
library(phiper)
library(dplyr)
library(tidyr)
library(ggplot2)

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

# =============================================================================
# Save tables
# =============================================================================
out_dir  <- "results/sample_sizes"
plot_dir <- file.path(out_dir, "barplots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Timepoint code -> human label
tp_labels       <- c(T0="P12", T1="P28", T2="B",  T3="W2",
                     T4="M1",  T5="M2",  T6="M3",  T7="M6", T8="M12")
tp_label_levels <- unname(tp_labels)

group_display <- c(kid_serum     = "Kid serum",
                   mom_serum     = "Mother serum",
                   mom_milk      = "Mother milk",
                   mom_kid_serum = "Mother-kid serum")

# Sample count tables (one CSV per group)
for (grp in names(sample_tables)) {
  write.csv(sample_tables[[grp]],
            file.path(out_dir, paste0("sample_counts_", grp, ".csv")),
            row.names = FALSE)
}

# Meta count tables (one CSV per meta_key)
for (mk in names(meta_tables)) {
  safe_mk <- gsub("[^A-Za-z0-9_]", "_", mk)
  write.csv(meta_tables[[mk]],
            file.path(out_dir, paste0("meta_counts_", safe_mk, ".csv")),
            row.names = FALSE)
}

# Long-format RDS for the Shiny app
saveRDS(sample_counts_long, file.path(out_dir, "sample_counts_long.rds"))
if (nrow(meta_counts_long) > 0) {
  saveRDS(meta_counts_long, file.path(out_dir, "meta_counts_long.rds"))
}

# =============================================================================
# Barplots
# =============================================================================
phiper::phip_use_montserrat()

# Pretty variable name lookup (covers all likely meta names from parse_meta_rest)
meta_pretty_map <- c(
  risk_CD          = "CDrisk",
  risk_CD_bin      = "CDrisk",
  CDrisk           = "CDrisk",
  cd               = "CDrisk",
  delivery_mode    = "Delivery mode",
  delivery         = "Delivery mode",
  delivery_place   = "Delivery place",
  place            = "Delivery place",
  feedmode_m3_bin  = "Feeding mode at M3",
  feedmode_m3      = "Feeding mode at M3",
  feeding_m3       = "Feeding mode at M3",
  feed             = "Feeding mode at M3",
  lockdown_status  = "Lockdown",
  lockdown         = "Lockdown",
  preterm          = "Preterm",
  sex              = "Sex",
  infant_sex       = "Sex",
  siblings         = "Siblings",
  smoking          = "Smoking"
)

pretty_meta <- function(x) {
  looked_up <- meta_pretty_map[x]
  ifelse(!is.na(looked_up), looked_up, x)
}

# ---- 1) Sample sizes overview -----------------------------------------------
sample_plot_df <- sample_counts_long %>%
  dplyr::mutate(
    timepoint_label = factor(tp_labels[as.character(timepoint)],
                             levels = tp_label_levels),
    group_label = dplyr::recode(group, !!!as.list(group_display))
  ) %>%
  dplyr::filter(!is.na(timepoint_label))

p_overview <- ggplot(sample_plot_df,
                     aes(x = timepoint_label, y = n_subjects, fill = group_label)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = n_subjects),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4.5, fontface = "bold", color = "black") +
  labs(title = "Sample sizes by group and timepoint",
       x = "Timepoint", y = "N subjects", fill = "Group") +
  phiper::scale_fill_phip() +
  phiper::theme_phip() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "sample_sizes_overview.svg"),
       p_overview, width = 8, height = 8, bg = "white")

# ---- 2) Metadata distribution over time (one plot per meta_key) -------------
if (nrow(meta_counts_long) > 0) {

  for (mk in sort(unique(meta_counts_long$meta_key))) {

    df_mk <- meta_counts_long %>%
      dplyr::filter(meta_key == mk) %>%
      dplyr::mutate(
        timepoint_label = factor(tp_labels[as.character(timepoint)],
                                 levels = tp_label_levels),
        # preserve level_priority ordering
        level = factor(level, levels = c(
          level_priority[level_priority %in% unique(level)],
          setdiff(unique(level), level_priority)
        ))
      ) %>%
      dplyr::filter(!is.na(timepoint_label))

    if (nrow(df_mk) == 0) next

    meta_name  <- unique(df_mk$meta)[1]
    meta_label <- pretty_meta(meta_name)
    grp_label  <- dplyr::recode(unique(df_mk$group)[1], !!!as.list(group_display))
    plot_title <- paste0(meta_label, ": ", grp_label)

    # Compute within-timepoint percentage for labels
    df_mk <- df_mk %>%
      dplyr::group_by(timepoint_label) %>%
      dplyr::mutate(pct = 100 * n_subjects / sum(n_subjects, na.rm = TRUE)) %>%
      dplyr::ungroup()

    p_meta <- ggplot(df_mk, aes(x = timepoint_label, y = n_subjects, fill = level)) +
      geom_col(position = "dodge", color = "black", linewidth = 0.3, width = 0.75) +
      geom_text(
        aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
        position = position_dodge(width = 0.75),
        vjust = -0.5, size = 4.5, fontface = "bold", color = "black"
      ) +
      labs(title = plot_title,
           x = "Timepoint", y = "N subjects", fill = meta_label) +
      phiper::scale_fill_phip() +
      phiper::theme_phip() +
      theme(text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1))

    safe_mk    <- gsub("[^A-Za-z0-9_]", "_", mk)
    safe_label <- gsub("[^A-Za-z0-9_]", "_", meta_label)
    ggsave(file.path(plot_dir, paste0(safe_label, "_", safe_mk, ".svg")),
           p_meta, width = 8, height = 8, bg = "white")
  }
}

message("Done. Results written to: ", out_dir)

invisible(crosstables)
