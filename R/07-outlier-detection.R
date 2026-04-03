# =============================================================================
# setup
# =============================================================================
library(phiper)
library(dplyr)
library(ggplot2)
library(tibble)
library(withr)

# =============================================================================
# 1) data import  (same call as 03-analysis.R)
# =============================================================================
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

phiper::phip_use_montserrat()

# =============================================================================
# 2) alpha diversity per sample, grouped by big_group
# =============================================================================
out_alpha <- phiper::compute_alpha_diversity(
  ps,
  group_cols        = "big_group",
  group_interaction = FALSE,
  interaction_only  = FALSE,
  carry_cols        = c("timepoint_recoded", "sample_id")
)

# filter to peptide-level rank (same pattern as 06-run-plate.R)
out_tbl <- out_alpha[[1]] %>%
  dplyr::filter(.data$rank == "peptide_id")

# =============================================================================
# settings
# =============================================================================
vars    <- c("richness")          # variables for Mahalanobis distance
small_eps <- 1e-8                 # covariance regularization
bins    <- 30                     # histogram bins
alpha   <- 0.01                   # significance level for Mahalanobis cutoff
p_vars  <- length(vars)

out_dir <- "results/outliers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# helpers
# =============================================================================

# compute Mahalanobis distance per group
compute_group_maha <- function(df_group, vars, small_eps = 1e-8) {
  df     <- df_group
  cc_idx <- complete.cases(df[, vars])
  X_cc   <- as.matrix(df[cc_idx, vars])
  p      <- ncol(X_cc)

  maha2 <- rep(NA_real_, nrow(df))
  maha  <- rep(NA_real_, nrow(df))
  pval  <- rep(NA_real_, nrow(df))

  if (nrow(X_cc) == 0) {
    df$maha2     <- maha2
    df$maha      <- maha
    df$maha_pval <- pval
    return(df)
  }

  mu <- colMeans(X_cc)
  S  <- tryCatch(cov(X_cc), error = function(e) NA)

  if (any(is.na(S)) || det(S) <= .Machine$double.eps) {
    S <- S + diag(small_eps, p)
  }

  d2 <- mahalanobis(X_cc, center = mu, cov = S)

  maha2[cc_idx] <- d2
  maha[cc_idx]  <- sqrt(d2)
  pval[cc_idx]  <- pchisq(d2, df = p, lower.tail = FALSE)

  df$maha2     <- maha2
  df$maha      <- maha
  df$maha_pval <- pval
  df
}

# safe file-name helper
clean_name <- function(x) {
  x2 <- gsub("\\s+", "_", as.character(x))
  x2 <- gsub("[^A-Za-z0-9_\\-]", "_", x2)
  substr(x2, 1, 80)
}

# make and save all plots + summary for one group
make_and_save_for_group <- function(df_group, group_name, suffix = "") {
  df_rich <- df_group %>% dplyr::filter(!is.na(richness))
  safe_g  <- clean_name(group_name)

  # --- histogram + density ---
  p_hist <- ggplot(df_rich, aes(x = richness)) +
    geom_histogram(aes(y = ..density..),
                   bins = bins, alpha = 0.5, color = "black", fill = "grey70") +
    geom_density(linewidth = 1) +
    labs(x = "Richness", y = "Density") +
    phiper::theme_phip() +
    theme(text = element_text(size = 25))
  ggsave(file.path(out_dir, paste0("hist_richness_", safe_g, suffix, ".svg")),
         p_hist, width = 12, height = 12, dpi = 300, bg = "white")

  # --- QQ plot ---
  p_qq <- ggplot(df_rich, aes(sample = richness)) +
    stat_qq(size = 0.8) +
    stat_qq_line(color = "blue") +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") +
    phiper::theme_phip() +
    theme(text = element_text(size = 25))
  ggsave(file.path(out_dir, paste0("qq_richness_", safe_g, suffix, ".svg")),
         p_qq, width = 12, height = 12, dpi = 300, bg = "white")

  # --- Mahalanobis ---
  df_maha      <- compute_group_maha(df_group, vars = vars, small_eps = small_eps)
  df_maha_plot <- df_maha %>% dplyr::filter(!is.na(maha))

  cut_d2   <- qchisq(1 - alpha, df = p_vars)
  cut_sqrt <- sqrt(cut_d2)
  n_flagged <- sum(df_maha$maha2 > cut_d2, na.rm = TRUE)

  if (nrow(df_maha_plot) > 0) {
    p_maha <- ggplot(df_maha_plot, aes(x = maha)) +
      geom_histogram(aes(y = ..density..),
                     bins = bins, alpha = 0.5, color = "black", fill = "grey80") +
      geom_density(linewidth = 1) +
      geom_vline(xintercept = cut_sqrt, color = "red", linewidth = 1) +
      annotate("text",
               x = cut_sqrt, y = Inf,
               label = paste0("alpha=", alpha, "\ncut=", round(cut_sqrt, 3)),
               vjust = 2.7, hjust = 1.15, color = "red", size = 7.5, lineheight = 0.32) +
      labs(x = "Mahalanobis distance (sqrt)", y = "Density") +
      phiper::theme_phip() +
      theme(text = element_text(size = 25))
    ggsave(file.path(out_dir, paste0("maha_", safe_g, suffix, ".svg")),
           p_maha, width = 12, height = 12, dpi = 300, bg = "white")
  }

  summary <- tibble::tibble(
    big_group            = group_name,
    n_total              = nrow(df_group),
    n_rich_nonNA         = sum(!is.na(df_group$richness)),
    n_complete_for_maha  = sum(complete.cases(df_group[, vars])),
    n_flagged_alpha      = n_flagged,
    cut_d2               = cut_d2,
    cut_sqrt             = cut_sqrt
  )

  list(summary = summary, df_maha = df_maha)
}

# =============================================================================
# 3) run per-group
# =============================================================================
groups <- unique(as.character(out_tbl$big_group))

all_summaries <- list()
all_maha      <- list()

for (g in groups) {
  df_g   <- out_tbl %>% dplyr::filter(big_group == g)
  if (nrow(df_g) == 0) next
  safe_g <- clean_name(g)

  res_full <- make_and_save_for_group(df_g, g, suffix = "")
  all_summaries[[paste0(safe_g, "_full")]] <- res_full$summary
  all_maha[[paste0(safe_g, "_full")]]      <- res_full$df_maha

  # mom_milk variant: richness < 4000
  if (g == "mom_milk") {
    df_g_f <- df_g %>% dplyr::filter(!is.na(richness) & richness < 4000)
    if (nrow(df_g_f) > 0) {
      res_lt <- make_and_save_for_group(df_g_f, g, suffix = "_lt4000")
      all_summaries[[paste0(safe_g, "_lt4000")]] <- res_lt$summary
      all_maha[[paste0(safe_g, "_lt4000")]]      <- res_lt$df_maha
    }
  }
}

# =============================================================================
# 4) flag outliers
# =============================================================================
th_lt   <- 2.576   # mom_milk_lt4000
th_full <- 10      # mom_milk_full

# collect sample metadata from ps for mapping
sample_meta <- ps$data_long %>%
  dplyr::select("sample_id", "subject_id", "timepoint_recoded", "big_group") %>%
  dplyr::distinct() %>%
  dplyr::collect()

# mom_milk_lt4000 outliers
hits_lt <- all_maha[["mom_milk_lt4000"]] %>%
  dplyr::filter(maha > th_lt) %>%
  dplyr::select("sample_id", "maha", "maha2", dplyr::everything())

map_lt <- sample_meta %>%
  dplyr::filter(sample_id %in% hits_lt$sample_id)

# mom_milk_full outliers
hits_full <- all_maha[["mom_milk_full"]] %>%
  dplyr::filter(maha > th_full) %>%
  dplyr::select("sample_id", "maha", "maha2", dplyr::everything())

map_full <- sample_meta %>%
  dplyr::filter(sample_id %in% hits_full$sample_id)

# =============================================================================
# 5) save summaries and full Mahalanobis table
# =============================================================================
summary_tbl <- dplyr::bind_rows(all_summaries)
write.csv(summary_tbl,
          file = file.path(out_dir, "summary_by_group.csv"),
          row.names = FALSE)

maha_all <- dplyr::bind_rows(all_maha, .id = "source_group_variant")
write.csv(maha_all,
          file = file.path(out_dir, "maha_all_groups.csv"),
          row.names = FALSE)

write.csv(hits_lt,
          file = file.path(out_dir, "outliers_mom_milk_lt4000.csv"),
          row.names = FALSE)

write.csv(hits_full,
          file = file.path(out_dir, "outliers_mom_milk_full.csv"),
          row.names = FALSE)

message("Done. Results written to: ", out_dir)
