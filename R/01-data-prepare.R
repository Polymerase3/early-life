## ----------------------------READING PACKAGES---------------------------------
# setting random generator seed
set.seed(16748991)
seed_before <- .Random.seed

# creating vector of necessary packages
packages <- c(
  "stringr",
  "magrittr",
  "data.table",
  "MASS",
  "ggplot2",
  "multilevelTools",
  "dplyr",
  "ggpubr",
  "tidyr",
  "lmerTest",
  "lme4",
  "extraoperators",
  "JWileymisc",
  "showtext",
  "DBI",
  "duckdb",
  "knitr",
  "kableExtra",
  "htmltools",
  "webshot2",
  "magick",
  "ggExtra",
  "mclust",
  "dplyr",
  "DBI",
  "readxl",
  "stringi",
  "glue",
  "purrr",
  "forcats",
  "DescTools",
  "corrplot",
  "tibble",
  "wesanderson",
  "tidytext",
  "Rtsne",
  "Matrix",
  "plotly",
  "randomcoloR",
  "stringr",
  "scales",
  "readr",
  "ggsignif",
  "ggnewscale",
  "scales",
  "ggbreak",
  "RColorBrewer",
  "rgl",
  "htmlwidgets",
  "forcats",
  "writexl",
  "fs",
  "glue",
  "plotly"
)

# load phiper (installed from the tarball)
library(phiper)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  pak::pkg_install(packages[!installed_packages])
}

# packages loading
invisible(lapply(packages, library, character.only = TRUE))

# add font + phiper use font in all plots
font_add_google("Montserrat", "monte")
phip_use_montserrat()
showtext_auto()


# removing unnecessary variables
rm(list = c("installed_packages", "packages"))

## ------------------ PHIP DATA + ages -----------------------------------------
## create phip_data
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path = "data/babies_counts.parquet",
    peptide_library = TRUE,
    subject_id = "subjectID",
    peptide_id = "peptideID",
    materialise_table = TRUE,
    auto_expand = FALSE,
    n_cores = 10
  )
})

## load the ages data.frame
ages <- read_excel("data/time_collection_blood_LLNEXT.xlsx",
                   col_types = c(
                     "text", "text", "text",
                     "text", "text", "numeric", "text",
                     "text", "text", "text", "text", "numeric",
                     "text", "numeric", "numeric", "numeric"
                   )
)

ages$next_id <- gsub("LLNEXT", "", ages$next_id, fixed = TRUE)

## Converting TIMEPOINTS to factors
ages$timepoint_factor <- dplyr::recode_factor(ages$timepoint,
                                              "P12" = "T0",
                                              "P28" = "T1",
                                              "B" = "T2",
                                              "W2" = "T3",
                                              "M1" = "T4",
                                              "M2" = "T5",
                                              "M3" = "T6",
                                              "M6" = "T7",
                                              "M12" = "T8"
)

## NOTE: after i performed the standardization, there were some observations
## missing for the following subjects:
# subject_id timepoint big_group
# <chr>          <dbl> <chr>
# 1 010577            38 kid_serum
# 2 008101            86 kid_serum
# 3 008525            38 mom_serum
# 4 206131            38 mom_serum
#
# The first two are obvious missing data in the ages data.frame, they have NA's
# in the place of exact age; for the subject 010577 there was no exact_age at
# birth, but we can approximate it with 0; for the subject 008101 exact age at
# 12 months was missing, so i approximated it with adding ~ 9 months to the M3
# value: 122.3229167 + 9 * 30
# For the other two mums, the exact_ages at birth were also missing, i resolved
# it manually using the last given value.

ages$exact_age[ages$next_id == "010577" & ages$timepoint == "B"] <- 0
ages$exact_age[ages$next_id == "008101" &
                 ages$timepoint == "M12"] <- 122.3229167 + 9 * 30
ages$exact_age[ages$next_id == "008525" &
                 ages$timepoint == "B"] <- 11148.52 + 22 * 7
ages$exact_age[ages$next_id == "206131" &
                 ages$timepoint == "B"] <- 13839.48 + 22 * 7

## ----- standardizing the ages to days_since_birth on a common scale ----------
## calculating the mom age at birth
mom_birth_age <- ages %>%
  filter(type == "M") %>%
  group_by(next_id) %>%
  arrange(
    !(timepoint == "B" & !is.na(exact_age)),
    exact_age
  ) %>%
  slice(1) %>%
  transmute(next_id,
            mom_minimal_age = exact_age,
            timepoint
  ) %>%
  ungroup()

## for some moms there is no data at birth --> i have to extrapolate from the
## timepoints
mom_birth_age <- mom_birth_age %>%
  mutate(
    mom_birth_age = case_when(
      timepoint == "B" ~ mom_minimal_age, # nothin changes
      timepoint == "P12" ~ mom_minimal_age + 26 * 7, # 38 weeks = 266 days
      timepoint == "P28" ~ mom_minimal_age + 10 * 7, # 22 weeks = 154 days
      TRUE ~ NA_real_ # fallback
    )
  ) %>%
  dplyr::select(next_id, mom_birth_age)

## merge the data with the ages
withr::with_preserve_seed({
  ps_merged <- left_join(ps,
                         ages[, c("next_id", "timepoint_factor", "exact_age")],
                         by = dplyr::join_by(
                           subject_id == next_id,
                           timepoint_factor == timepoint_factor
                         ),
                         copy = TRUE
  )
})

## selecting only necessary variables
ps_merged %<>% dplyr::select(
  sample_id, subject_id, big_group,
  timepoint, timepoint_factor, peptide_id, relative,
  dyade, exact_age, fold_change
)

## adding mom_birth_age
withr::with_seed(1653, {
  ps_merged <- left_join(ps_merged, mom_birth_age,
                         by = dplyr::join_by(subject_id == next_id),
                         copy = TRUE
  )
})

# converting the exact_age to days_since_birth variable
ps_merged %<>%
  mutate(
    months_since_birth = case_when(
      big_group == "kid_serum" ~ exact_age / 30,
      big_group == "mom_serum" ~ (exact_age - mom_birth_age) / 30,
      big_group == "mom_milk" & timepoint == 40 ~ 0.5,
      big_group == "mom_milk" & timepoint == 42 ~ 1,
      big_group == "mom_milk" & timepoint == 46 ~ 2,
      big_group == "mom_milk" & timepoint == 50 ~ 3,
      big_group == "mom_milk" & timepoint == 62 ~ 6,
      big_group == "mom_milk" & timepoint == 86 ~ 12,
      # dont have measurments --> approx
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(-mom_birth_age, -timepoint, -exact_age)

## sanity check --> any missing values in the timepoint?
ps_merged %>%
  summarise(
    n_missing = sum(is.na(months_since_birth))
  )

## remove unnecessary
rm(list = c("ps", "ages", "mom_birth_age"))
gc()

## ------------------- APPENDING METADATA --------------------------------------
## we will need them anyway for the GLMMs and tSNEs, so better to append them rn
## the metadata are actually only important for the kids, so i will append them
## only to the kids --> we dont do any grouping analyses for moms/milks -->
## wouldn't make much sense actually
metas <- read_tsv(
  file = "data/Metadata_LLNEXT_including_feedingmode.txt",
  na = c("NA", "", "n/a", "N/A"),
  trim_ws = TRUE,
  locale = locale(encoding = "UTF-8"),
  col_types = cols(
    NG_ID = col_character(),
    next_id_infant = col_character(),
    next_id_mother = col_character(),
    FAMILY.x = col_character(),
    infant_relations = col_character(),
    sibling_number = col_character(),
    twin_pair = col_character(),
    birth_deliverybirthcard_mode_binary = col_character(),
    infant_birthcard_feeding_mode_after_birth = col_character(),
    infant_birthcard_med_antibiotics_T195_J01CA_B = col_character(),
    mother_birthcard_age_at_delivery = col_double(),
    infant_misc_sex = col_character(),
    infant_ffq_ever_never_breastfed = col_character(),
    birth_birthcard_reg_gestational_age_days = col_double(),
    birth_birthcard_reg_gestational_age_weeks = col_double(),
    birth_birthcard_reg_gestational_age_week_string = col_character(),
    birth_birthcard_reg_source = col_character(),
    feedmode_m3 = col_character(),
    feedmode_m3_txt = col_character(),
    baby_birthcard_delivery_mode = col_character(),
    Feeding_mode_Birth = col_character(),
    mother_health_smoked_one_whole_year_p18 = col_character(),
    mother_birthcard_parity = col_double(),
    birth_deliverybirthcard_place_delivery_simple = col_character(),
    FAMILY.y = col_character(),
    Timepoint_categorical = col_character(),
    infant_ffq_feeding_mode_simple = col_character(),
    infant_ffq_feeding_mode_complex = col_character(),
    infant_ffq_breastfeeding_freq_weighted = col_double(),
    infant_ffq_formula_follow_on_freq_weighted = col_double(),
    infant_ffq_stopped_breastfeeding = col_character()
  )
) %>%
  # trim any leftover whitespace in character columns
  mutate(across(where(is.character), ~ trimws(.))) %>%
  # tidy duplicated family columns: prefer FAMILY.x, if missing use FAMILY.y
  mutate(
    FAMILY = coalesce(FAMILY.x, FAMILY.y)
  ) %>%
  select(-FAMILY.x, -FAMILY.y) # drop originals

# select the vars of interest: other are either duplicates or not-necessary or
# doesn't add any value (antibiotics)
metas %<>% select(
  NG_ID,
  next_id_infant, # each infant is unique
  Timepoint_categorical,
  infant_ffq_feeding_mode_simple,
  infant_ffq_feeding_mode_complex
)

## the problem in the new metadata is, that there are multiple rows per subject;
## it is because somebody appended the feeding mode at the end of the old meta
## table with multiple timepoints per subject. We now have to pivot them to wide
id_col <- "next_id_infant"
time_col <- "Timepoint_categorical"
pivot_cols <- c("infant_ffq_feeding_mode_simple",
                "infant_ffq_feeding_mode_complex")

# pivot timepoint feeding cols wide
time_wide <- metas %>%
  select(all_of(c(id_col, time_col, pivot_cols))) %>%
  distinct() %>%
  pivot_wider(
    id_cols = all_of(id_col),
    names_from = all_of(time_col),
    values_from = all_of(pivot_cols),
    names_sep = "_"
  )

## the thing is: we want to fill M1 with W2 and M3 with M2; the merging rule is:
## when M1 or M3 present, then leave it, when not present, then try to fill with
## W2 or M2 respectively; at the end remove the redundant cols
# --- fill M1 from W2 when M1 is missing ---
if (all(c("infant_ffq_feeding_mode_simple_W2",
          "infant_ffq_feeding_mode_simple_M1") %in% names(time_wide))) {
  time_wide <- time_wide %>%
    mutate(infant_ffq_feeding_mode_simple_M1 =
             coalesce(infant_ffq_feeding_mode_simple_M1,
                      infant_ffq_feeding_mode_simple_W2))
}
if (all(c("infant_ffq_feeding_mode_complex_W2",
          "infant_ffq_feeding_mode_complex_M1") %in% names(time_wide))) {
  time_wide <- time_wide %>%
    mutate(infant_ffq_feeding_mode_complex_M1 =
             coalesce(infant_ffq_feeding_mode_complex_M1,
                      infant_ffq_feeding_mode_complex_W2))
}

# --- fill M3 from M2 when M3 is missing ---
if (all(c("infant_ffq_feeding_mode_simple_M2",
          "infant_ffq_feeding_mode_simple_M3") %in% names(time_wide))) {
  time_wide <- time_wide %>%
    mutate(infant_ffq_feeding_mode_simple_M3 =
             coalesce(infant_ffq_feeding_mode_simple_M3,
                      infant_ffq_feeding_mode_simple_M2))
}
if (all(c("infant_ffq_feeding_mode_complex_M2",
          "infant_ffq_feeding_mode_complex_M3") %in% names(time_wide))) {
  time_wide <- time_wide %>%
    mutate(infant_ffq_feeding_mode_complex_M3 =
             coalesce(infant_ffq_feeding_mode_complex_M3,
                      infant_ffq_feeding_mode_complex_M2))
}

# --- drop the source columns (W2 and M2), if present ---
time_wide <- time_wide %>%
  select(-any_of(c(
    "infant_ffq_feeding_mode_simple_W2",
    "infant_ffq_feeding_mode_complex_W2",
    "infant_ffq_feeding_mode_simple_M2",
    "infant_ffq_feeding_mode_complex_M2"
  )))

# --- now read the other file with metadata ---
# it was really a pain for me, cause the two files have different number of
# subjects... firstly i used only the file with different timepoints for feeding
# mode and run all the tests, but then i realised, that i have a lot less
# samples in the analyses, cause somebody didn't match the files at all...
# so right now i import the second file and merge it with the time_wide to get
# the metadata for all samples and dont lose any samples
# and CAVE: in the metadata file thats on slack, one infantID is literally "s";
# i replaced it with the right one
other <- read_excel("data/Metadata_LLNEXT_early_life_paper.xlsx")

other %>%
  left_join(time_wide, by = "next_id_infant")

# join -> one row per subject
metas <- other %>%
  left_join(time_wide, by = id_col) %>%
  rename(
    id_infant = next_id_infant,
    delivery_mode = birth_deliverybirthcard_mode_binary,
    feedmode_birth = infant_birthcard_feeding_mode_after_birth,
    mother_delivery_age = mother_birthcard_age_at_delivery,
    infant_sex = infant_misc_sex,
    ever_breastfed = infant_ffq_ever_never_breastfed,
    gestage_birth_days = birth_birthcard_reg_gestational_age_days,
    gestage_birth_weeks = birth_birthcard_reg_gestational_age_weeks,
    feedmode_m3 = feedmode_m3,
    smoking = mother_health_smoked_one_whole_year_p18,
    parity = mother_birthcard_parity,
    delivery_place = birth_deliverybirthcard_place_delivery_simple
  ) %>%
  dplyr::mutate(
    parity_num = as.numeric(dplyr::na_if(trimws(parity), "NA")),
    ga_weeks   = as.numeric(dplyr::na_if(trimws(gestage_birth_weeks), "NA")),
    siblings   = dplyr::case_when(
      !is.na(parity_num) & parity_num > 0 ~ "yes",   # per your rule: >1
      !is.na(parity_num)                  ~ "no",
      TRUE                                ~ NA_character_
    ),
    preterm    = dplyr::case_when(
      !is.na(ga_weeks) & ga_weeks < 37 ~ "yes",
      !is.na(ga_weeks)                 ~ "no",
      TRUE                             ~ NA_character_
    )
  ) %>%
  dplyr::select(-parity_num, -ga_weeks)

# shorten the long ffq column prefixes:
names(metas) <- names(metas) %>%
  gsub("^infant_ffq_feeding_mode_simple_", "ffq_simple_", .) %>%
  gsub("^infant_ffq_feeding_mode_complex_", "ffq_complex_", .)

# delete the LLNEXT prefix from the ID
metas$id_infant <- str_sub(metas$id_infant, 7)

# convert the character NAs to real NAs for all variables
metas <- metas %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA")))

# cd risk, getting rid of the LLNEXT prefix and merging RiskCD
cd <- read_excel("data/CD_risk_share.xlsx")
cd$LLNEXTsample.id <- str_sub(cd$LLNEXTsample.id , 7)
colnames(cd)[1] <- "id_infant"

metas <- metas %>%
  left_join(
    cd %>%
      select(id_infant, RiskCD) %>%
      rename(risk_CD = RiskCD),
    by = "id_infant"
  ) %>%
  mutate(
    risk_CD_bin = case_when(
      is.na(risk_CD) ~ NA_character_,
      tolower(risk_CD) %in% c("high", "medium") ~ "yes",
      TRUE ~ "no"
    )
  )

# lockdown_status --> read only the required columns
lk <- read_tsv(
  "data/Metadata_LLNEXT_including_feedingmode_Covid.txt",
  col_types = cols(.default = col_guess())
) %>%
  select(next_id_infant, Timepoint_categorical, lockdown_status) %>%
  # normalize formatting (optional but robust)
  mutate(
    lockdown_status = na_if(str_squish(str_to_lower(lockdown_status)), "")
  )

# per-subject consistency check
check_lockdown <- lk %>%
  distinct(next_id_infant, Timepoint_categorical, lockdown_status) %>%
  group_by(next_id_infant) %>%
  summarise(
    n_timepoints = n_distinct(Timepoint_categorical),
    n_status_non_na = n_distinct(lockdown_status[!is.na(lockdown_status)]),
    all_na = all(is.na(lockdown_status)),
    consistent = n_status_non_na <= 1,
    status = case_when(
      all_na ~ NA_character_,
      n_status_non_na == 1 ~ first(lockdown_status[!is.na(lockdown_status)]),
      TRUE ~ "mixed"
    ),
    .groups = "drop"
  )

# collapsed subject × lockdown_status data frame
subject_lockdown <- check_lockdown %>%
  transmute(next_id_infant, lockdown_status = status)
subject_lockdown$next_id_infant <- str_sub(subject_lockdown$next_id_infant, 7)
colnames(subject_lockdown)[1] <- "id_infant"

# merge with metas
metas <- metas %>%
  left_join(subject_lockdown, by = "id_infant")

# recoding the breastfeeding m3, according to Sasha's comment
metas <- metas %>%
  mutate(
    feedmode_m3_clean = feedmode_m3 %>%
      str_squish() %>% str_to_upper() %>% na_if(""),

    feedmode_m3_bin = case_when(
      is.na(feedmode_m3_clean) ~ NA_character_,
      feedmode_m3_clean %in% c("BF", "OTHERWISE_NAMELY:") ~ "BF",
      feedmode_m3_clean %in% c("FF", "FF/BF") ~ "MF+FF",
      TRUE ~ "MF+FF"
    )) %>%
  select(-feedmode_m3_clean)  # drop helper

## ------------------- BINARIZING METADATA -------------------------------------
## converting the metadata to binary dummies
metas_out <- metas %>%
  transmute(
    id_infant,
    siblings_yes      = as.logical(siblings == "yes"),
    siblings_no       = as.logical(siblings == "no"),
    delmode_VG        = as.logical(delivery_mode == "VG"),
    delmode_CS        = as.logical(delivery_mode == "CS"),
    delplace_home     = as.logical(delivery_place == "home"),
    delplace_hospital = as.logical(delivery_place == "hospital"),
    feeding_BF        = as.logical(feedmode_m3_bin == "BF"),
    feeding_MF        = as.logical(feedmode_m3_bin == "MF+FF"),
    preterm_yes       = as.logical(preterm == "yes"),
    preterm_no        = as.logical(preterm == "no"),
    CDrisk_yes        = as.logical(risk_CD_bin == "yes"),
    CDrisk_no         = as.logical(risk_CD_bin == "no"),
    lockdown_before   = as.logical(lockdown_status == "voor_lockdown"),
    lockdown_after    = as.logical(lockdown_status == "na_lockdown"),
    smoking_yes       = as.logical(smoking == "yes"),
    smoking_no        = as.logical(smoking == "no"),
    sex_male          = as.logical(infant_sex == "male"),
    sex_female        = as.logical(infant_sex == "female")
  )

## ------------------- JOINING METADATA -------------------------------------
# duckdb-optimised merge of the metadata to the main ps_merged object;
# only for the infants --> all other data will have NAs in the place of the
# metadata!!!

# same DuckDB connection as the big phiper object
con <- dbplyr::remote_con(ps_merged$data_long)

# 1) copy metas into duckdb (TEMP) so no auto-copy happens
tmp <- paste0("metas_tmp_", as.integer(Sys.time()))
copy_to(con, metas_out, name = tmp, temporary = TRUE, overwrite = TRUE)
metas_db <- tbl(con, tmp)

# 2) LEFT JOIN directly on the tbl_lazy then MATERIALIZE that tbl
joined_tbl <- left_join(
  ps_merged$data_long,
  metas_db,
  by = dplyr::join_by(subject_id == id_infant)
)

joined_tbl_mat <- compute(
  joined_tbl,
  name = paste0("ps_merged_long_", as.integer(Sys.time())),
  temporary = TRUE
)

# 3) put the materialized table back into phip_data
ps_merged$data_long <- joined_tbl_mat
ps_merged_meta <- ps_merged

# 4) now its safe to drop the temp metas table
DBI::dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', tmp))

# 5) clean other mess
rm(
  "ps_merged", "metas", "con", "joined_tbl", "joined_tbl_mat", "metas_db",
  "tmp", "cd", "subject_lockdown", "check_lockdown", "metas_out"
)
gc()

## ------------------- BINARIZING THE DATA -------------------------------------
# right now the data are usable only for fold_change-specific analyses - it
# consists only of the enriched peptides and it is not possible to directly
# compute prevalences and do the repertoire analyses; we have to expand the data
# to the binary form, which causes the nrow() to explode, as we are storing the
# data in the long format; hence the decision to use duckdb and phiper - the
# long format is necessary for plotting/GLMMs and many other analyses
# expand and tag existence
ps_merged_meta_bin <- expand_phip_data(ps_merged_meta,
                                       add_exist = TRUE,
                                       exist_col = "exist"
)

# then some cosmetics: add new column for the timepoints of breast milk, we
# wanted to recode W2 --> M1 and M2 --> M3; then also recode the dyades,
# originally i hardcoded the same dyades for siblings/twins and mother, so all
# 3 subjects had the same dyade ID; it was a little bit problematic in the
# heatmaps the, as the siblings were plotted against each other, and not the
# same subjects, so i recoded the dyade column to a new one, so we have another
# option to try in the heatmaps

# --- SETTINGS: (each element length-2) ---
pairs_list <- list(
  c("120677", "150788"),
  c("202903", "202915"),
  c("204122", "205217")
)
subjects <- unlist(pairs_list)
new_tbl <- "ps_merged_meta_bin_materialized"
match_col <- "subject_id"

# --- connection + detect source table ---
con <- ps_merged_meta_bin$meta$con
orig_tbl <- if (is.character(ps_merged_meta_bin$data_long)) {
  ps_merged_meta_bin$data_long
} else {
  dbplyr::remote_name(ps_merged_meta_bin$data_long)
}
ps_tbl <- tbl(con, orig_tbl)

# --- get current subject -> dyade mapping (small table) ---
df_map <- ps_tbl %>%
  select(!!sym(match_col), dyade) %>%
  distinct() %>%
  collect()

df_map <- df_map %>%
  mutate(dyade_num = suppressWarnings(as.integer(trimws(as.character(dyade)))))

# --- compute start = current max numeric dyade (0 if none) ---
start <- if (all(is.na(df_map$dyade_num))) {
  0L
} else {
  max(df_map$dyade_num, na.rm = TRUE)
}

start <- ifelse(is.infinite(start) | is.na(start), 0L, as.integer(start))

# --- allocate new dyade numbers by pair, appended at the end ---
# for pair i -> numbers: start + (2*(i-1)) + 1 and +2
pair_allocs <- lapply(seq_along(pairs_list), function(i) {
  base <- start + 2 * (i - 1)
  tibble(
    subject_id = pairs_list[[i]],
    new_dyade  = as.character(base + seq_len(length(pairs_list[[i]])))
  )
})
mapping <- bind_rows(pair_allocs)

# --- write mapping to DB temporary table ---
DBI::dbWriteTable(con, "tmp_pairs_dyade_map", mapping, temporary = TRUE,
                  overwrite = TRUE)

# --- build lazy recoded table and materialize it ---
# create a new table with timepoint_recoded (same rules) and dyade_recoded
# (coalesce)
ps_tbl %>%
  mutate(
    timepoint_recoded = case_when(
      big_group == "mom_milk" & timepoint_factor == "T3" ~ "T4",
      big_group == "mom_milk" & timepoint_factor == "T5" ~ "T6",
      TRUE ~ timepoint_factor
    )
  ) %>%
  left_join(tbl(con, "tmp_pairs_dyade_map"), by = match_col) %>%
  mutate(dyade_recoded = coalesce(new_dyade, dyade)) %>%
  select(everything(), -new_dyade) %>%
  compute(name = new_tbl, temporary = FALSE)

# --- safe swap: rename original -> backup, new -> original ---
backup_name <- paste0(orig_tbl, "_backup_", format(Sys.time(), "%Y%m%d%H%M%S"))
DBI::dbExecute(con, paste0(
  "ALTER TABLE ", DBI::dbQuoteIdentifier(con, orig_tbl),
  " RENAME TO ", DBI::dbQuoteIdentifier(con, backup_name)
))
DBI::dbExecute(con, paste0(
  "ALTER TABLE ", DBI::dbQuoteIdentifier(con, new_tbl),
  " RENAME TO ", DBI::dbQuoteIdentifier(con, orig_tbl)
))

# --- drop other temp tables (keep only the new data_long) ---
all_tables <- DBI::dbListTables(con)
to_keep <- orig_tbl
to_drop <- setdiff(all_tables, to_keep)
if (length(to_drop) > 0) {
  for (t in to_drop) {
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ",
                               DBI::dbQuoteIdentifier(con, t)))
  }
}

# refresh R reference
ps_merged_meta_bin$data_long <- dplyr::tbl(con, orig_tbl)

# -------------------___CREATEING THE DUMMY COMP COLUMNS-----------------------------
con <- ps_merged_meta_bin$meta$con
orig_tbl <- dbplyr::remote_name(ps_merged_meta_bin$data_long)

core_cols <- c(
  "sample_id", "peptide_id", "subject_id",
  "big_group", "timepoint_recoded", "dyade_recoded", "exist"
)

# define list of comparisons (pairs of group labels) to analyze
comparisons <- list(
  c("mom_serum_T0", "mom_serum_T1"),
  c("mom_serum_T0", "mom_serum_T2"),
  #c("mom_serum_T1", "mom_serum_T2"),

  c("mom_serum_T2", "kid_serum_T2"),
  c("kid_serum_T2", "kid_serum_T6"),
  c("kid_serum_T2", "kid_serum_T8"),
  c("kid_serum_T6", "kid_serum_T8"),
  c("kid_serum_T2", "kid_serum_T8"),

  c("mom_milk_T4", "mom_milk_T6"),
  c("mom_milk_T4", "mom_milk_T7"),
  c("mom_milk_T4", "mom_milk_T8"),

  c("kid_serum_T6", "mom_milk_T6"),
  c("kid_serum_T8", "mom_milk_T8"),

  c("kid_serum_T2_siblings", "kid_serum_T2_no_siblings"),
  c("kid_serum_T6_siblings", "kid_serum_T6_no_siblings"),
  c("kid_serum_T8_siblings", "kid_serum_T8_no_siblings"),

  c("kid_serum_T2_delmode_VG", "kid_serum_T2_delmode_CS"),
  c("kid_serum_T6_delmode_VG", "kid_serum_T6_delmode_CS"),
  c("kid_serum_T8_delmode_VG", "kid_serum_T8_delmode_CS"),

  c("kid_serum_T2_delplace_home", "kid_serum_T2_delplace_hospital"),
  c("kid_serum_T6_delplace_home", "kid_serum_T6_delplace_hospital"),
  c("kid_serum_T8_delplace_home", "kid_serum_T8_delplace_hospital"),

  c("kid_serum_T2_feeding_BF", "kid_serum_T2_feeding_MF"),
  c("kid_serum_T6_feeding_BF", "kid_serum_T6_feeding_MF"),
  c("kid_serum_T8_feeding_BF", "kid_serum_T8_feeding_MF"),

  c("kid_serum_T2_preterm_yes", "kid_serum_T2_preterm_no"),
  c("kid_serum_T6_preterm_yes", "kid_serum_T6_preterm_no"),
  c("kid_serum_T8_preterm_yes", "kid_serum_T8_preterm_no"),

  c("kid_serum_T2_CDrisk_yes", "kid_serum_T2_CDrisk_no"),
  c("kid_serum_T6_CDrisk_yes", "kid_serum_T6_CDrisk_no"),
  c("kid_serum_T8_CDrisk_yes", "kid_serum_T8_CDrisk_no"),

  c("kid_serum_T2_lockdown_before", "kid_serum_T2_lockdown_after"),
  c("kid_serum_T6_lockdown_before", "kid_serum_T6_lockdown_after"),
  c("kid_serum_T8_lockdown_before", "kid_serum_T8_lockdown_after"),

  c("kid_serum_T2_smoking_yes", "kid_serum_T2_smoking_no"),
  c("kid_serum_T6_smoking_yes", "kid_serum_T6_smoking_no"),
  c("kid_serum_T8_smoking_yes", "kid_serum_T8_smoking_no"),

  c("kid_serum_T2_sex_male", "kid_serum_T2_sex_female"),
  c("kid_serum_T6_sex_male", "kid_serum_T6_sex_female"),
  c("kid_serum_T8_sex_male", "kid_serum_T8_sex_female")
)

# comparisons must exist already in your env
labels <- unique(unlist(comparisons))

# --- parse labels like "kid_serum_T2_smoking_yes" ---
parse_label <- function(label) {
  m <- regexec("^(.*)_(T\\d+)(?:_(.*))?$", label)
  x <- regmatches(label, m)[[1]]
  if (length(x) == 0) stop("Bad label format: ", label)

  list(
    big_group = x[2],
    tp        = x[3],
    cov       = if (length(x) >= 5 && nzchar(x[4])) x[4] else NA_character_
  )
}

# --- build SQL expression for a label (BOOLEAN column) ---
make_label_expr_sql <- function(label, con) {
  p <- parse_label(label)

  base <- sprintf(
    "(big_group = %s AND timepoint_recoded = %s)",
    DBI::dbQuoteString(con, p$big_group),
    DBI::dbQuoteString(con, p$tp)
  )

  if (!is.na(p$cov)) {
    # cov is an existing 0/1 dummy; treat NA as 0
    cov_sql <- sprintf(
      "(COALESCE(%s, 0) = 1)",
      DBI::dbQuoteIdentifier(con, p$cov)
    )
    base <- sprintf("(%s AND %s)", base, cov_sql)
  }

  sprintf("%s AS %s", base, DBI::dbQuoteIdentifier(con, label))
}

# Optional: fail fast if any covariate columns are missing
desc <- DBI::dbGetQuery(con, sprintf("DESCRIBE %s;", DBI::dbQuoteIdentifier(con, orig_tbl)))
colname_field <- if ("column_name" %in% names(desc)) "column_name" else names(desc)[1]
existing_cols <- desc[[colname_field]]

covs_needed <- unique(na.omit(vapply(labels, function(x) parse_label(x)$cov, character(1))))
missing_covs <- setdiff(covs_needed, existing_cols)
if (length(missing_covs) > 0) {
  stop("Missing covariate columns needed for labels: ", paste(missing_covs, collapse = ", "))
}

# Build SELECT list
label_exprs_sql <- vapply(labels, make_label_expr_sql, character(1), con = con)

select_sql <- paste(
  c(
    paste(DBI::dbQuoteIdentifier(con, core_cols), collapse = ",\n  "),
    paste(label_exprs_sql, collapse = ",\n  ")
  ),
  collapse = ",\n  "
)

new_tbl <- paste0(orig_tbl, "__core_plus_labels")

sql_create <- sprintf(
  "CREATE TABLE %s AS
   SELECT
     %s
   FROM %s;",
  DBI::dbQuoteIdentifier(con, new_tbl),
  select_sql,
  DBI::dbQuoteIdentifier(con, orig_tbl)
)

DBI::dbExecute(con, sql_create)

# --- swap tables (backup old) ---
backup_tbl <- paste0(orig_tbl, "_backup_", format(Sys.time(), "%Y%m%d%H%M%S"))

DBI::dbExecute(con, sprintf(
  "ALTER TABLE %s RENAME TO %s;",
  DBI::dbQuoteIdentifier(con, orig_tbl),
  DBI::dbQuoteIdentifier(con, backup_tbl)
))

DBI::dbExecute(con, sprintf(
  "ALTER TABLE %s RENAME TO %s;",
  DBI::dbQuoteIdentifier(con, new_tbl),
  DBI::dbQuoteIdentifier(con, orig_tbl)
))

# refresh lazy ref
ps_merged_meta_bin$data_long <- dplyr::tbl(con, orig_tbl)

## ---------------------------------- DATA EXPORT ------------------------------
# export an existing DuckDB table to .parquet
export_path <- "data/phiper_validated_babies_full.parquet"

con <- ps_merged_meta_bin$meta$con
tblname <- dbplyr::remote_name(ps_merged_meta_bin$data_long)

DBI::dbExecute(
  con,
  sprintf(
    "COPY %s TO '%s' (FORMAT PARQUET);",
    DBI::dbQuoteIdentifier(con, tblname),
    export_path
  )
)

#
rm(list = c(
  "con", "mapping", "meta_kids", "other", "ps_merged_meta", "ps_tbl",
  "time_wide", "all_tables", "backup_name", "cache_path", "id_col",
  "match_col", "new_tbl", "orig_tbl", "pivot_cols", "t",
  "targets", "time_col", "to_drop", "to_keep", "df_map"
))

gc()
gc()
