#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

metas <- readRDS("data/metas_clean.rds")

# Descriptives
age_desc <- metas %>%
  summarise(
    n = sum(!is.na(mother_birthcard_age_at_delivery)),
    mean = mean(mother_birthcard_age_at_delivery, na.rm = TRUE),
    sd = sd(mother_birthcard_age_at_delivery, na.rm = TRUE),
    min = min(mother_birthcard_age_at_delivery, na.rm = TRUE),
    max = max(mother_birthcard_age_at_delivery, na.rm = TRUE)
  )

ga_desc <- metas %>%
  summarise(
    n = sum(!is.na(birth_birthcard_reg_gestational_age_weeks)),
    mean = mean(birth_birthcard_reg_gestational_age_weeks, na.rm = TRUE),
    sd = sd(birth_birthcard_reg_gestational_age_weeks, na.rm = TRUE),
    min = min(birth_birthcard_reg_gestational_age_weeks, na.rm = TRUE),
    max = max(birth_birthcard_reg_gestational_age_weeks, na.rm = TRUE)
  )

cat("\nMother birthcard age at delivery (years):\n")
print(age_desc)

cat("\nGestational age at birth (weeks):\n")
print(ga_desc)

# ------------------------------------------------------------------
# Marginal counts per subject (any timepoint), using cleaned metas
# ------------------------------------------------------------------

count_unique <- function(df, var, label = deparse(substitute(var))) {
  df %>%
    distinct(id_infant, {{ var }}) %>%
    filter(!is.na({{ var }})) %>%
    group_by(level = {{ var }}) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(variable = label, .before = 1)
}

sex_counts <- count_unique(metas, infant_sex, "sex")
delmode_counts <- count_unique(metas, delivery_mode, "delivery_mode")
delplace_counts <- count_unique(metas, delivery_place, "delivery_place")
preterm_counts <- count_unique(metas, preterm, "preterm")
feeding_counts <- count_unique(metas, feedmode_m3_bin, "feeding_m3")
siblings_counts <- count_unique(metas, siblings, "siblings")
smoking_counts <- count_unique(metas, smoking, "smoking")
cdrisk_counts <- count_unique(metas, risk_CD_bin, "CDrisk")
lockdown_counts <- count_unique(metas, lockdown_status, "lockdown")

all_counts <- bind_rows(
  sex_counts,
  delmode_counts,
  delplace_counts,
  preterm_counts,
  feeding_counts,
  siblings_counts,
  smoking_counts,
  cdrisk_counts,
  lockdown_counts
) %>%
  group_by(variable) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

cat("\nMarginal counts by subject (all timepoints collapsed):\n")
print(all_counts)

table(metas$infant_misc_sex)
