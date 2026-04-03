# Early Life - Lifelines NEXT Analysis

[![phiper](https://img.shields.io/badge/phiper-Polymerase3%2Fphiper-blue)](https://github.com/Polymerase3/phiper)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: inactive](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

This repository contains a frozen analysis pipeline for PhIP-Seq early-life antibody repertoire data from the Lifelines NEXT cohort.
The goal is reproducibility, not ongoing development.

This repository contains the analysis code written by **Mateusz Franciszek Kolek (MFK)** for a scientific article based on the **Lifelines NEXT** cohort. The study examines the PhIP-Seq antibody repertoires of mothers and their infants using serum and breast milk samples across multiple timepoints.

---

## System Requirements

**Operating system:** Linux (tested on Ubuntu 22.04 / kernel 6.17, x86\_64). The pipeline should run on any modern Linux distribution; macOS is untested but likely compatible. Windows is not supported.

**RAM:** ≥ 64 GB recommended for the full Lifelines NEXT dataset. The demo script runs comfortably in < 8 GB.

**R version:** ≥ 4.2.0 (tested on R 4.4.x).

**Required CRAN / Bioconductor packages:**

| Package | Purpose |
|---------|---------|
| `stringr`, `stringi` | String manipulation |
| `magrittr` | Pipe operator |
| `data.table` | Fast tabular data |
| `MASS` | Statistical methods |
| `ggplot2`, `ggpubr`, `ggExtra`, `ggsignif`, `ggnewscale`, `ggbreak` | Visualisation |
| `dplyr`, `tidyr`, `purrr`, `forcats`, `tibble` | Tidyverse data wrangling |
| `lme4`, `lmerTest`, `multilevelTools`, `JWileymisc` | Mixed-effects modelling |
| `extraoperators` | Extended operators |
| `showtext` | Custom fonts |
| `DBI`, `duckdb`, `dbplyr` | In-process SQL / parquet I/O |
| `knitr`, `kableExtra`, `htmltools`, `webshot2` | Table rendering |
| `magick` | Image processing |
| `mclust` | Gaussian mixture models |
| `readxl`, `readr`, `writexl` | Spreadsheet I/O |
| `glue`, `fs` | String interpolation / filesystem |
| `DescTools`, `corrplot` | Descriptive statistics |
| `wesanderson`, `RColorBrewer`, `randomcoloR` | Colour palettes |
| `tidytext` | Text mining |
| `Rtsne` | t-SNE |
| `Matrix` | Sparse matrices |
| `plotly`, `htmlwidgets`, `rgl` | Interactive graphics |
| `scales` | Scale helpers |
| `withr` | Scoped side effects |
| `pak` | Fast package installation |

**In-house package:** `phiper` 0.2.6 (source tarball included; see Installation below).

---

## Installation

Estimated time: **< 5 minutes** on a standard internet connection.

### 1. Install phiper from the included tarball

```r
install.packages("phiper_0.2.6.tar.gz", repos = NULL, type = "source")
```

Alternatively, install the latest development version from GitHub (although it is not recommended, as breaking changes occured since v. 0.2.6):

```r
# install.packages("remotes")
remotes::install_github("Polymerase3/phiper")
```

### 2. Install all CRAN dependencies

```r
# install.packages("pak")   # run once if pak is not available
pak::pkg_install(c(
  "stringr", "magrittr", "data.table", "MASS", "ggplot2",
  "multilevelTools", "dplyr", "ggpubr", "tidyr", "lmerTest",
  "lme4", "extraoperators", "JWileymisc", "showtext",
  "DBI", "duckdb", "dbplyr", "knitr", "kableExtra", "htmltools",
  "webshot2", "magick", "ggExtra", "mclust", "readxl",
  "stringi", "glue", "purrr", "forcats", "DescTools",
  "corrplot", "tibble", "wesanderson", "tidytext", "Rtsne",
  "Matrix", "plotly", "randomcoloR", "scales", "readr",
  "ggsignif", "ggnewscale", "ggbreak", "RColorBrewer",
  "rgl", "htmlwidgets", "writexl", "fs", "withr", "pak"
))
```

---

## Demo

A small simulated PhIP-Seq dataset is provided in `demo/early_life_demo.parquet`. It contains a subset of peptides and samples representative of the full Lifelines NEXT data structure and runs in **< 5 minutes** on a standard laptop.

### Run the demo

```r
source("demo/demo.R")
```

This will:
1. Load `demo/early_life_demo.parquet` into a `phip_data` object using `phip_convert()`.
2. Compute per-sample PhIP-Seq richness.
3. Print a summary table and save a richness plot to `demo/demo_richness.png`.

Expected output: a PNG file `demo/demo_richness.png` and a printed richness summary table in the R console.

---

## Instructions for Use

To run the full analysis on your own PhIP-Seq data:

1. Place your long-format parquet file in `data/`. The expected columns (names are configurable in `phip_convert()`) are: `sampleID` (or `sample_id`), `peptideID` (or `peptide_id`), `exist`, `fold_change`, `input`, `count`.
2. Update the path in `R/01-data-prepare.R` (line 92):
   ```r
   data_long_path = "data/your_file.parquet"
   ```
3. Adjust `subject_id`, `peptide_id`, and other column names in `phip_convert()` to match your dataset.
4. Run the numbered scripts in order (`01-data-prepare.R` → `02-qc_counts.R` → …). Each script sources `R/25-utils.R` and `R/26-zzz.R` for shared utilities and comparison definitions.
5. Outputs (figures, tables) are written to `results/` and `main_figures/`.

---

## phiper

The analyses rely on the in-house R package **phiper**. The version used in this project is **0.2.6**; the source tarball is included in the root of this repository:

- `phiper_0.2.6.tar.gz`

The latest version of phiper is maintained at: https://github.com/Polymerase3/phiper

`R/05-richness-figure1.R` also loads **phiperio**, a companion I/O helper package. It is not required for the demo or any other script in this repository.

---

## R/ - analysis scripts

| File | Description |
|------|-------------|
| `01-data-prepare.R` | Loads packages, imports PhIP-Seq count data from parquet, processes ages and metadata, and exports prepared objects for downstream scripts. |
| `02-qc_counts.R` | Builds QC crosstables counting unique subjects and samples stratified by timepoint, group, sequencing run, and plate. |
| `03-analysis.R` | Main analysis script: implements POP (prevalence) and DELTA (differential abundance) frameworks with permutation testing, forest plots, and scatter plots. |
| `04-tsne-multiple-groups.R` | Computes 2D and 3D t-SNE ordinations on Kulczynski distance matrices; accepts command-line arguments (N_CORES, CMP_INDEX, CMP_KEY) for HPC array job execution. |
| `05-richness-figure1.R` | Generates richness and repertoire-stability plots with time-gradient colour coding per group (mom serum, kid serum, mom milk) using GAM smoothing and bootstrap confidence intervals. |
| `06-run-plate.R` | Counts samples per run/plate/group/timepoint and computes alpha-diversity (richness) stratified by sequencing run and plate with `rstatix` statistical comparisons. |
| `07-heatmaps.R` | Computes full pairwise Kulczynski repertoire similarity across all group × timepoint combinations and exports static ggplot2 similarity heatmaps with dyad-swapped ordering. |
| `07-outlier-detection.R` | Flags per-group richness outliers using Mahalanobis distance; produces histogram, QQ, and Mahalanobis distribution plots (SVG) and outlier CSV tables in `results/outliers/`. |
| `08-custom-scatters.R` | Produces custom scatter plots comparing prevalence percentages between groups. |
| `09-quick-look.R` | Preliminary visualization of DELTA results (forest plots, enrichment summaries) for a single comparison. |
| `10-fdr.R` | Applies Benjamini-Hochberg FDR correction and empirical-null FDR via locfdr on standardized test statistics across all comparisons. |
| `11-pcoa.R` | Performs PCoA (and optional NMDS) on Kulczynski distance matrices with ellipses and centroid plotting. |
| `12-descriptives.R` | Summary statistics for maternal age, gestational age, and categorical metadata variables. |
| `13-breast-milk-meta.R` | Significance counts and locfdr visualization grid for breast milk metadata comparisons. |
| `14-normal-meta.R` | FDR summary for selected kid-serum metadata comparisons (siblings, delivery mode, feeding, preterm status, etc.). |
| `15-replot-snippet.R` | Minimal script to reproduce a specific scatter plot (mom P12 vs. birth) from existing POP output. |
| `16-sankey-meta.R` | Alluvial/Sankey plots and t-SNE visualizations tracking feature significance trajectories across three timepoints. |
| `17-sankey-summary.R` | Heatmap summarizing metadata variable effects with preferential feature ordering and custom colour schemes. |
| `18-alex-flagellins.R` | POP and DELTA analysis of flagellin-specific peptides comparing presence/absence across sample groups. |
| `19-preterm-t2-bh-summary.R` | Focused Benjamini-Hochberg FDR analysis on preterm status at T2 with detailed significance counting. |
| `20-alex-flagellins-analysis.R` | Comprehensive flagellin analysis merged with microbiome data: richness, distance, POP, and DELTA comparisons. |
| `21-alex-fdr.R` | locfdr correction across all flagellin comparisons with summary tables and feature visualizations. |
| `22-alex-bubble.R` | Bubble plot of microbiome feature associations across timepoints (point size = effect magnitude, shape = significance). |
| `23-flags-summary.R` | Summary table of nominal and FDR-corrected significance counts for flagellin analyses by timepoint. |
| `24-sankey-fig2.R` | Sankey/alluvial plots, heatmaps, scatter plots, and t-SNE visualizations comparing feature significance across two specific comparisons. |
| `25-utils.R` | Shared utility functions used across scripts. |
| `26-zzz.R` | Centralised definition of all comparison pairs (group/metadata combinations) and their longitudinal status flags. |
| `27-richness.R` | Exports alpha-diversity tables (peptide to phylum ranks) and enrichment plots for `big_group × timepoint_recoded`; also computes kid_serum subgroup alpha-diversity and enrichment plots for 9 categorical metadata variables; outputs to `results/richness/` for the Shiny app. |
| `28-helpers-plot.R` | Shared helper module: `build_plotly_3d`/`build_plotly_2d` for interactive t-SNE scatter plots and `export_stability_heatmap_html` for exporting Kulczynski similarity matrices as self-contained Plotly HTML heatmaps. |
| `29-heatmaps-html.R` | Iterates over all group × timepoint state pairs, exports each as an interactive Plotly Kulczynski heatmap (with per-subject alpha-diversity tooltips) and a CSV to `results/heatmaps_html/`; sources `28-helpers-plot.R`. |
| `30-stability-meta.R` | Tests associations between kid_serum repertoire stability (Kulczynski, M0/M3/M12 contrasts) and categorical metadata (siblings, delivery mode/place, feeding, preterm, smoking, sex) using bootstrap mean-difference tests; exports an Excel table and violin+boxplots per variable. |

---

## slurm/

SLURM job scripts for running array jobs on the HPC cluster.
