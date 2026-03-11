# Early Life — Lifelines NEXT Analysis

This repository contains the analysis code written by **Mateusz Franciszek Kolek (MFK)** for a scientific article based on the **Lifelines NEXT** cohort. The study examines the PhIP-Seq antibody repertoires of mothers and their infants using serum and breast milk samples across multiple timepoints.

---

## phiper

The analyses rely on the in-house R package **phiper**. The version used in this project is **0.2.6**; the source tarball is included in the root of this repository:

- `phiper_0.2.6.tar.gz`

The latest version of phiper is maintained at: https://github.com/Polymerase3/phiper

---

## R/ - analysis scripts

| File | Description |
|------|-------------|
| `01-data-prepare.R` | Loads packages, imports PhIP-Seq count data from parquet, processes ages and metadata, and exports prepared objects for downstream scripts. |
| `02-qc_counts.R` | Generates crosstables of sample and metadata counts by timepoint and group for quality control. |
| `03-analysis.R` | Main analysis script: implements POP (prevalence) and DELTA (differential abundance) frameworks with permutation testing, forest plots, and scatter plots. |
| `04-tsne-multiple-groups.R` | Computes t-SNE ordinations on Kulczynski distance matrices for mom serum, kid serum, and mom milk across timepoint comparisons. |
| `05-richness-figure1.R` | Generates richness and repertoire-stability plots with GAM smoothing and bootstrap confidence intervals across groups. |
| `06-run-plate.R` | Computes alpha/beta diversity diagnostics stratified by sequencing run and plate, including statistical comparisons. |
| `07-heatmaps.R` | Creates heatmaps from DELTA results showing taxa enrichment patterns across comparisons. |
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

---

## slurm/

SLURM job scripts for running array jobs on the HPC cluster.
