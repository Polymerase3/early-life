## =============================================================================
## Early Life – Lifelines NEXT: minimal reproducible demo
##
## Purpose: demonstrate that phiper is correctly installed and that the
##          demo dataset can be loaded and analysed.
##
## Runtime:  < 5 minutes on a standard laptop (tested on Linux x86_64, R 4.4)
## Output:   demo/demo_richness.png + a richness summary printed to the console
##
## Run from the repository root:
##   source("demo/demo.R")
## =============================================================================

## ---- 0. dependencies --------------------------------------------------------
required_pkgs <- c("phiper", "dplyr", "ggplot2", "withr")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                       quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0L) {
  stop(
    "The following packages are required but not installed:\n  ",
    paste(missing_pkgs, collapse = ", "),
    "\nSee the Installation section in README.md."
  )
}

suppressPackageStartupMessages({
  library(phiper)
  library(dplyr)
  library(ggplot2)
})

## ---- 1. load demo dataset ---------------------------------------------------
demo_path <- file.path("demo", "early_life_demo.parquet")
if (!file.exists(demo_path)) {
  stop("Demo parquet not found at '", demo_path, "'.\n",
       "Make sure you are running this script from the repository root.")
}

message("Loading demo dataset …")
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path  = demo_path,
    sample_id       = "sample_id",
    peptide_id      = "peptide_id",
    subject_id      = "subject_id",
    exist           = "exist",
    fold_change     = "fold_change",
    peptide_library = FALSE,
    materialise_table = TRUE,
    auto_expand     = FALSE,
    n_cores         = 1L
  )
})
message("phip_data object created successfully.")

## ---- 2. compute per-sample richness -----------------------------------------
message("Computing per-sample PhIP-Seq richness …")

richness_df <- ps |>
  dplyr::filter(exist == 1L) |>
  dplyr::group_by(sample_id, group_char) |>
  dplyr::summarise(richness = dplyr::n(), .groups = "drop") |>
  dplyr::collect()

## ---- 3. print summary -------------------------------------------------------
message("\n--- Richness summary (hits per sample) ---")
summary_tbl <- richness_df |>
  dplyr::group_by(group_char) |>
  dplyr::summarise(
    n_samples   = dplyr::n(),
    mean_rich   = round(mean(richness), 1),
    median_rich = round(median(richness), 1),
    sd_rich     = round(sd(richness), 1),
    .groups     = "drop"
  )
print(summary_tbl, n = Inf)

## ---- 4. plot ----------------------------------------------------------------
message("\nGenerating richness plot …")

p <- ggplot(richness_df,
            aes(x = reorder(group_char, richness, median),
                y = richness,
                colour = group_char)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  coord_flip() +
  labs(
    title    = "Per-sample PhIP-Seq richness - demo subset",
    subtitle = "Number of peptides with exist == 1 per sample",
    x        = NULL,
    y        = "Richness (# reactive peptides)",
    colour   = "Group"
  ) +
  theme_phip() +
  theme(legend.position = "none")

out_png <- file.path("demo", "demo_richness.png")
ggsave(out_png, plot = p, width = 9, height = 5, dpi = 150, bg = "white")
message("Richness plot saved to: ", out_png)
message("\nDemo complete.")
