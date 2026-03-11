#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

data_dir <- file.path("data", "alex_flagellins")

matrix_paths <- list.files(
  data_dir,
  pattern = "^(phipseq4analysis|microbe4analysis)_.+\\.txt$",
  full.names = TRUE
)

manifest <- tibble(path = matrix_paths) %>%
  mutate(
    file = basename(.data$path),
    matrix_name = str_remove(.data$file, "\\.txt$"),
    type = case_when(
      str_detect(.data$file, "^phipseq4analysis_") ~ "phipseq",
      str_detect(.data$file, "^microbe4analysis_") ~ "microbiome",
      TRUE ~ NA_character_
    ),
    panel = str_match(.data$file, "analysis_(.+)\\.txt$")[, 2]
  ) %>%
  arrange(.data$type, .data$panel)

read_wide_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)

  # Files are stored with an unnamed first column that contains sample IDs.
  lines[[1]] <- paste0("sample_id\t", lines[[1]])

  read_tsv(
    I(lines),
    col_types = cols(
      sample_id = col_character(),
      .default = col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
}

# 1) Import all wide matrices into a list.
wide_matrices <- manifest %>%
  mutate(data = map(.data$path, read_wide_matrix)) %>%
  select("matrix_name", "data") %>%
  deframe()

phipseq_wide <- wide_matrices[str_detect(names(wide_matrices), "^phipseq4analysis_")]
microbiome_wide <- wide_matrices[str_detect(names(wide_matrices), "^microbe4analysis_")]

phipseq_panels <- str_remove(names(phipseq_wide), "^phipseq4analysis_")
microbiome_panels <- str_remove(names(microbiome_wide), "^microbe4analysis_")
common_panels <- intersect(phipseq_panels, microbiome_panels)

# 2) Save first-column IDs from one phipseq + one microbiome file, then build remapping.
reference_panel <- if ("combined" %in% common_panels) "combined" else common_panels[[1]]

reference_phipseq_ids <- phipseq_wide[[paste0("phipseq4analysis_", reference_panel)]] %>%
  transmute(row_index = row_number(), phipseq_id = .data$sample_id)

reference_microbiome_ids <- microbiome_wide[[paste0("microbe4analysis_", reference_panel)]] %>%
  transmute(row_index = row_number(), microbiome_id = .data$sample_id)

id_remap_reference <- reference_microbiome_ids %>%
  left_join(reference_phipseq_ids, by = "row_index")

build_panel_remap <- function(panel) {
  phip_name <- paste0("phipseq4analysis_", panel)
  micro_name <- paste0("microbe4analysis_", panel)

  phip_ids <- phipseq_wide[[phip_name]] %>%
    transmute(row_index = row_number(), phipseq_id = .data$sample_id)

  micro_ids <- microbiome_wide[[micro_name]] %>%
    transmute(row_index = row_number(), microbiome_id = .data$sample_id)

  if (nrow(phip_ids) != nrow(micro_ids)) {
    stop("Row count mismatch for panel: ", panel)
  }

  remap <- micro_ids %>%
    left_join(phip_ids, by = "row_index")

  if (anyNA(remap$phipseq_id)) {
    stop("Remapping failed for panel: ", panel)
  }

  remap
}

panel_remaps <- set_names(map(common_panels, build_panel_remap), common_panels)

# 3) Recode microbiome sample IDs to phipseq IDs; first column is already named sample_id.
microbiome_wide_recoded <- map2(
  microbiome_wide,
  names(microbiome_wide),
  function(df, matrix_name) {
    panel <- str_remove(matrix_name, "^microbe4analysis_")
    remap <- panel_remaps[[panel]]

    if (is.null(remap)) {
      stop("No remap found for microbiome panel: ", panel)
    }
    if (nrow(df) != nrow(remap)) {
      stop("Row count mismatch while recoding microbiome panel: ", panel)
    }

    df %>%
      mutate(sample_id = remap$phipseq_id)
  }
)
names(microbiome_wide_recoded) <- str_remove(names(microbiome_wide_recoded), "^microbe4analysis_")

# 4) Add timepoint column to each phipseq dataset.
phipseq_wide_with_timepoint <- map2(
  phipseq_wide,
  names(phipseq_wide),
  function(df, matrix_name) {
    timepoint <- str_remove(matrix_name, "^phipseq4analysis_")
    df %>%
      mutate(timepoint = timepoint, .after = sample_id)
  }
)
names(phipseq_wide_with_timepoint) <- str_remove(names(phipseq_wide_with_timepoint), "^phipseq4analysis_")

# 5) Pivot phipseq wide -> long with columns: sample_id, timepoint, peptide_id, exist.
phipseq_long_by_timepoint <- map(
  phipseq_wide_with_timepoint,
  function(df) {
    peptide_cols <- names(df)[str_detect(names(df), "^agilent_")]
    if (length(peptide_cols) == 0L) {
      stop("No agilent_* columns found in phipseq matrix.")
    }

    df %>%
      pivot_longer(
        cols = all_of(peptide_cols),
        names_to = "peptide_id",
        values_to = "exist"
      ) %>%
      select("sample_id", "timepoint", "peptide_id", "exist")
  }
)

phipseq_long <- bind_rows(phipseq_long_by_timepoint)

# 6) Left-join microbiome to phipseq per timepoint/panel on sample_id.
phipseq_with_microbiome_by_timepoint <- map2(
  phipseq_long_by_timepoint,
  names(phipseq_long_by_timepoint),
  function(phip_long, panel) {
    micro_df <- microbiome_wide_recoded[[panel]]
    if (is.null(micro_df)) {
      stop("No microbiome matrix for panel: ", panel)
    }

    left_join(phip_long, micro_df, by = "sample_id")
  }
)

# 7) Row-bind all joined datasets, keeping all unique microbiome columns.
# Missing microbiome variables in a given panel are filled with NA.
phipseq_with_microbiome_all <- bind_rows(phipseq_with_microbiome_by_timepoint)

# 8) Sanitize sample_id and extract subject_id (substring before first '-').
phipseq_with_microbiome_all <- phipseq_with_microbiome_all %>%
  mutate(subject_id = str_extract(.data$sample_id, "^[^-]+"), .before = "sample_id")

# 9) Save the final merged table to the same source directory.
final_rds_path <- file.path(data_dir, "phipseq_with_microbiome_all.rds")
saveRDS(phipseq_with_microbiome_all, final_rds_path)
