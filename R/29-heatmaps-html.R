# load required packages
library(phiper)
library(rlang)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(withr)
source("R/28-helpers-plot.R")
source("R/25-utils.R")

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

################################################################################
########################### COMPUTING STABILITIES ##############################
################################################################################
con <- ps$meta$con
ps$meta$longitudinal <- TRUE

# 1) compute stability using dyade_recoded and keep only kid_serum vs kid_serum
stab_recoded_tbl <- compute_repertoire_similarity(
  ps,
  group_col    = "big_group",
  time_col     = "timepoint_recoded",
  similarity   = "kulczynski",
  mode         = "all",
  dyade_col    = "dyade_recoded",
  time_mode    = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "all",
  drop_self_time = FALSE,
  overwrite    = TRUE,
  verbose      = FALSE
)

# materialise kid vs kid subset
if (DBI::dbExistsTable(con, "stab_kid_kid")) {
  DBI::dbExecute(con, "DROP TABLE IF EXISTS stab_kid_kid")
}
stab_recoded_tbl %>%
  filter(group_left == "kid_serum", group_right == "kid_serum") %>%
  compute(name = "stab_kid_kid", temporary = FALSE)

# 2) compute full stability using original dyade
stab_orig_tbl <- compute_repertoire_similarity(
  ps,
  group_col    = "big_group",
  time_col     = "timepoint_recoded",
  similarity   = "kulczynski",
  mode         = "all",
  dyade_col    = "dyade",
  time_mode    = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "all",
  drop_self_time = FALSE,
  overwrite    = TRUE,
  verbose      = FALSE
)

# 3) union: drop kid-kid from dyade-based result, append kid-kid from recoded
if (DBI::dbExistsTable(con, "stab_df_full")) {
  DBI::dbExecute(con, "DROP TABLE IF EXISTS stab_df_full")
}
stab_orig_tbl %>%
  filter(!(group_left == "kid_serum" & group_right == "kid_serum")) %>%
  union_all(tbl(con, "stab_kid_kid")) %>%
  compute(name = "stab_df_full", temporary = FALSE)

stab_df_full <- dplyr::tbl(con, "stab_df_full")

# 4) recode timepoint columns T0-T8 -> P12/P28/M0/M1/M3/M6/M12
con      <- dbplyr::remote_con(stab_df_full)
tbl_name <- dbplyr::remote_name(stab_df_full)

map_t2p <- tibble::tibble(
  from = c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"),
  to   = c("P12", "P28", "M0", "M1", "M1", "M3", "M3", "M6", "M12")
)
copy_to(con, map_t2p, name = "map_t2p_tmp", temporary = TRUE, overwrite = TRUE)

DBI::dbExecute(con, sprintf(
  'UPDATE "%1$s" AS s
      SET time_left = COALESCE(m.to, s.time_left)
    FROM map_t2p_tmp AS m
   WHERE s.time_left = m.from;', tbl_name
))
DBI::dbExecute(con, sprintf(
  'UPDATE "%1$s" AS s
      SET time_right = COALESCE(m.to, s.time_right)
    FROM map_t2p_tmp AS m
   WHERE s.time_right = m.from;', tbl_name
))
DBI::dbExecute(con, sprintf('ANALYZE "%s";', tbl_name))

# 5) swap high dyad IDs towards the middle for better plot ordering
max_dyad <- stab_df_full %>%
  transmute(
    dyad_left  = as.integer(dyad_left),
    dyad_right = as.integer(dyad_right)
  ) %>%
  summarise(
    m1 = max(dyad_left,  na.rm = TRUE),
    m2 = max(dyad_right, na.rm = TRUE)
  ) %>%
  collect() %>%
  { max(c(.$m1, .$m2), na.rm = TRUE) }

hi_range  <- seq(max_dyad - 39L, max_dyad)
mid_range <- 40L:79L
map_df <- tibble::tibble(
  old = as.character(c(hi_range, mid_range)),
  new = as.character(c(mid_range, hi_range))
)

stab_df_full <- stab_df_full %>%
  left_join(map_df, by = c("dyad_left"  = "old"), copy = TRUE) %>%
  mutate(dyad_left  = dplyr::coalesce(new, dyad_left))  %>%
  select(-new) %>%
  left_join(map_df, by = c("dyad_right" = "old"), copy = TRUE) %>%
  mutate(dyad_right = dplyr::coalesce(new, dyad_right)) %>%
  select(-new)

################################################################################
########################### ALPHA DIVERSITY ####################################
################################################################################
alpha_div <- compute_alpha_diversity(
  ps,
  group_cols        = c("big_group", "timepoint_recoded"),
  ranks             = "peptide_id",
  group_interaction = TRUE,
  carry_cols        = "subject_id"
)

alpha_tbl <- alpha_div$`big_group * timepoint_recoded` %>%
  filter(grepl("^big_group \\* timepoint_recoded$", `phip_interaction`, perl = TRUE) |
           !is.na(`phip_interaction`)) %>%
  rename(interaction = `phip_interaction`) %>%
  filter(grepl(" \\* ", interaction)) %>%
  tidyr::separate(interaction,
                  into   = c("big_group", "timepoint_recoded"),
                  sep    = " \\* ",
                  remove = FALSE) %>%
  mutate(time_mapped = dplyr::recode(timepoint_recoded,
                                     !!!setNames(map_t2p$to, map_t2p$from))) %>%
  group_by(subject_id, group = big_group, time = time_mapped) %>%
  summarise(
    richness = mean(richness,          na.rm = TRUE),
    shannon  = mean(shannon_diversity, na.rm = TRUE),
    .groups  = "drop"
  )

################################################################################
########################## HTML HEATMAP EXPORT #################################
################################################################################
dir.create("results/heatmaps_html", recursive = TRUE, showWarnings = FALSE)

# cartesian of available states
left_states  <- stab_df_full %>% distinct(g1 = group_left,  t1 = time_left)  %>% collect()
right_states <- stab_df_full %>% distinct(g2 = group_right, t2 = time_right) %>% collect()
pairs_df <- tidyr::crossing(left_states, right_states)

# deterministic ordering
time_order <- c("P12", "P28", "M0", "M1", "M2", "M3", "M6", "M12")
ord_time <- function(x) factor(x, levels = unique(c(time_order, setdiff(x, time_order))), ordered = TRUE)

pairs_df <- pairs_df %>%
  mutate(
    g1o = as.character(g1), g2o = as.character(g2),
    t1o = ord_time(t1),     t2o = ord_time(t2)
  ) %>%
  arrange(g1o, t1o, g2o, t2o) %>%
  select(-g1o, -g2o, -t1o, -t2o)

# main loop
n   <- nrow(pairs_df)
pad <- max(2, nchar(as.character(n)))
log_list <- vector("list", n)

cat(sprintf("Exporting %d Plotly heatmap(s) ...\n", n))

for (i in seq_len(n)) {
  g1 <- pairs_df$g1[[i]]; t1 <- pairs_df$t1[[i]]
  g2 <- pairs_df$g2[[i]]; t2 <- pairs_df$t2[[i]]

  pairing <- if (identical(t1, t2)) "same" else "cross"
  dyads   <- "order"
  grid    <- "intersection"

  res <- tryCatch(
    export_stability_heatmap_html(
      sim_tbl       = stab_df_full,
      groups        = c(g1, g2),
      times         = c(t1, t2),
      pairing       = pairing,
      dyads         = dyads,
      grid          = grid,
      dedupe        = "mean",
      z_limits      = c(0, 1),
      outdir        = "results/heatmaps_html",
      file_stub     = sprintf(paste0("heat-%0", pad, "d-"), i),
      selfcontained = TRUE,
      alpha_tbl     = alpha_tbl
    ),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    cat(sprintf("[%-*d] FAIL: %s\n", pad, i, conditionMessage(res)))
    log_list[[i]] <- data.frame(
      idx = i, g1 = g1, t1 = t1, g2 = g2, t2 = t2,
      status = "fail", file_html = NA, file_csv = NA,
      error = conditionMessage(res)
    )
  } else {
    cat(sprintf("[%-*d] OK  -> %s\n", pad, i, basename(res$html)))
    log_list[[i]] <- data.frame(
      idx = i, g1 = g1, t1 = t1, g2 = g2, t2 = t2,
      status = "ok", file_html = res$html, file_csv = res$csv, error = NA
    )
  }
}

log_df <- dplyr::bind_rows(log_list)
readr::write_csv(log_df, "results/heatmaps_html/plotly_heatmap_export_log.csv")
cat("Done.\n")
