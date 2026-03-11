# Centralized comparison definitions shared across scripts.
comparisons <- list(
  c("mom_serum_T0", "mom_serum_T1"),
  c("mom_serum_T0", "mom_serum_T2"),
  c("mom_serum_T1", "mom_serum_T2"),

  c("mom_serum_T2", "kid_serum_T2"),
  c("kid_serum_T2", "kid_serum_T6"),
  c("kid_serum_T2", "kid_serum_T8"),
  c("kid_serum_T6", "kid_serum_T8"),
  c("mom_milk_T4", "mom_milk_T6"),
  c("mom_milk_T4", "mom_milk_T7"),
  c("mom_milk_T4", "mom_milk_T8"),
  c("kid_serum_T6", "mom_milk_T6"),
  c("kid_serum_T8", "mom_milk_T8"),
  c("mom_serum_T2", "mom_milk_T4"),
  c("mom_serum_T2", "mom_milk_T6"),
  c("mom_serum_T2", "mom_milk_T7"),
  c("kid_serum_T6", "mom_milk_T4"),
  c("kid_serum_T6", "mom_milk_T7"),
  c("kid_serum_T8", "mom_milk_T4"),
  c("kid_serum_T8", "mom_milk_T7"),

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
  c("kid_serum_T8_sex_male", "kid_serum_T8_sex_female"),

  c("mom_milk_T4_siblings", "mom_milk_T4_no_siblings"),
  c("mom_milk_T6_siblings", "mom_milk_T6_no_siblings"),
  c("mom_milk_T7_siblings", "mom_milk_T7_no_siblings"),
  c("mom_milk_T8_siblings", "mom_milk_T8_no_siblings"),
  c("mom_milk_T4_delmode_VG", "mom_milk_T4_delmode_CS"),
  c("mom_milk_T6_delmode_VG", "mom_milk_T6_delmode_CS"),
  c("mom_milk_T7_delmode_VG", "mom_milk_T7_delmode_CS"),
  c("mom_milk_T8_delmode_VG", "mom_milk_T8_delmode_CS"),
  c("mom_milk_T4_delplace_home", "mom_milk_T4_delplace_hospital"),
  c("mom_milk_T6_delplace_home", "mom_milk_T6_delplace_hospital"),
  c("mom_milk_T7_delplace_home", "mom_milk_T7_delplace_hospital"),
  c("mom_milk_T8_delplace_home", "mom_milk_T8_delplace_hospital"),
  c("mom_milk_T4_feeding_BF", "mom_milk_T4_feeding_MF"),
  c("mom_milk_T6_feeding_BF", "mom_milk_T6_feeding_MF"),
  c("mom_milk_T7_feeding_BF", "mom_milk_T7_feeding_MF"),
  c("mom_milk_T8_feeding_BF", "mom_milk_T8_feeding_MF"),
  c("mom_milk_T4_preterm_yes", "mom_milk_T4_preterm_no"),
  c("mom_milk_T6_preterm_yes", "mom_milk_T6_preterm_no"),
  c("mom_milk_T7_preterm_yes", "mom_milk_T7_preterm_no"),
  c("mom_milk_T8_preterm_yes", "mom_milk_T8_preterm_no"),
  c("mom_milk_T4_CDrisk_yes", "mom_milk_T4_CDrisk_no"),
  c("mom_milk_T6_CDrisk_yes", "mom_milk_T6_CDrisk_no"),
  c("mom_milk_T7_CDrisk_yes", "mom_milk_T7_CDrisk_no"),
  c("mom_milk_T8_CDrisk_yes", "mom_milk_T8_CDrisk_no"),
  c("mom_milk_T4_lockdown_before", "mom_milk_T4_lockdown_after"),
  c("mom_milk_T6_lockdown_before", "mom_milk_T6_lockdown_after"),
  c("mom_milk_T7_lockdown_before", "mom_milk_T7_lockdown_after"),
  c("mom_milk_T8_lockdown_before", "mom_milk_T8_lockdown_after"),
  c("mom_milk_T4_smoking_yes", "mom_milk_T4_smoking_no"),
  c("mom_milk_T6_smoking_yes", "mom_milk_T6_smoking_no"),
  c("mom_milk_T7_smoking_yes", "mom_milk_T7_smoking_no"),
  c("mom_milk_T8_smoking_yes", "mom_milk_T8_smoking_no"),
  c("mom_milk_T4_sex_male", "mom_milk_T4_sex_female"),
  c("mom_milk_T6_sex_male", "mom_milk_T6_sex_female"),
  c("mom_milk_T7_sex_male", "mom_milk_T7_sex_female"),
  c("mom_milk_T8_sex_male", "mom_milk_T8_sex_female")
)

longitudinal <- c(
  TRUE,  # mom_serum_T0 vs mom_serum_T1
  TRUE,  # mom_serum_T0 vs mom_serum_T2
  TRUE,  # mom_serum_T1 vs mom_serum_T2
  TRUE,  # mom_serum_T2 vs kid_serum_T2
  TRUE,  # kid_serum_T2 vs kid_serum_T6
  TRUE,  # kid_serum_T2 vs kid_serum_T8
  TRUE,  # kid_serum_T6 vs kid_serum_T8
  TRUE,  # mom_milk_T4 vs mom_milk_T6
  TRUE,  # mom_milk_T4 vs mom_milk_T7
  TRUE,  # mom_milk_T4 vs mom_milk_T8
  TRUE,  # kid_serum_T6 vs mom_milk_T6
  TRUE,  # kid_serum_T8 vs mom_milk_T8
  TRUE,  # mom_serum_T2 vs mom_milk_T4
  TRUE,  # mom_serum_T2 vs mom_milk_T6
  TRUE,  # mom_serum_T2 vs mom_milk_T7
  TRUE,  # kid_serum_T6 vs mom_milk_T4
  TRUE,  # kid_serum_T6 vs mom_milk_T7
  TRUE,  # kid_serum_T8 vs mom_milk_T4
  TRUE,  # kid_serum_T8 vs mom_milk_T7

  FALSE, # kid_serum_T2_siblings vs kid_serum_T2_no_siblings
  FALSE, # kid_serum_T6_siblings vs kid_serum_T6_no_siblings
  FALSE, # kid_serum_T8_siblings vs kid_serum_T8_no_siblings
  FALSE, # kid_serum_T2_delmode_VG vs kid_serum_T2_delmode_CS
  FALSE, # kid_serum_T6_delmode_VG vs kid_serum_T6_delmode_CS
  FALSE, # kid_serum_T8_delmode_VG vs kid_serum_T8_delmode_CS
  FALSE, # kid_serum_T2_delplace_home vs kid_serum_T2_delplace_hospital
  FALSE, # kid_serum_T6_delplace_home vs kid_serum_T6_delplace_hospital
  FALSE, # kid_serum_T8_delplace_home vs kid_serum_T8_delplace_hospital
  FALSE, # kid_serum_T2_feeding_BF vs kid_serum_T2_feeding_MF
  FALSE, # kid_serum_T6_feeding_BF vs kid_serum_T6_feeding_MF
  FALSE, # kid_serum_T8_feeding_BF vs kid_serum_T8_feeding_MF
  FALSE, # kid_serum_T2_preterm_yes vs kid_serum_T2_preterm_no
  FALSE, # kid_serum_T6_preterm_yes vs kid_serum_T6_preterm_no
  FALSE, # kid_serum_T8_preterm_yes vs kid_serum_T8_preterm_no
  FALSE, # kid_serum_T2_CDrisk_yes vs kid_serum_T2_CDrisk_no
  FALSE, # kid_serum_T6_CDrisk_yes vs kid_serum_T6_CDrisk_no
  FALSE, # kid_serum_T8_CDrisk_yes vs kid_serum_T8_CDrisk_no
  FALSE, # kid_serum_T2_lockdown_before vs kid_serum_T2_lockdown_after
  FALSE, # kid_serum_T6_lockdown_before vs kid_serum_T6_lockdown_after
  FALSE, # kid_serum_T8_lockdown_before vs kid_serum_T8_lockdown_after
  FALSE, # kid_serum_T2_smoking_yes vs kid_serum_T2_smoking_no
  FALSE, # kid_serum_T6_smoking_yes vs kid_serum_T6_smoking_no
  FALSE, # kid_serum_T8_smoking_yes vs kid_serum_T8_smoking_no
  FALSE, # kid_serum_T2_sex_male vs kid_serum_T2_sex_female
  FALSE, # kid_serum_T6_sex_male vs kid_serum_T6_sex_female
  FALSE, # kid_serum_T8_sex_male vs kid_serum_T8_sex_female

  FALSE, # mom_milk_T4_siblings vs mom_milk_T4_no_siblings
  FALSE, # mom_milk_T6_siblings vs mom_milk_T6_no_siblings
  FALSE, # mom_milk_T7_siblings vs mom_milk_T7_no_siblings
  FALSE, # mom_milk_T8_siblings vs mom_milk_T8_no_siblings
  FALSE, # mom_milk_T4_delmode_VG vs mom_milk_T4_delmode_CS
  FALSE, # mom_milk_T6_delmode_VG vs mom_milk_T6_delmode_CS
  FALSE, # mom_milk_T7_delmode_VG vs mom_milk_T7_delmode_CS
  FALSE, # mom_milk_T8_delmode_VG vs mom_milk_T8_delmode_CS
  FALSE, # mom_milk_T4_delplace_home vs mom_milk_T4_delplace_hospital
  FALSE, # mom_milk_T6_delplace_home vs mom_milk_T6_delplace_hospital
  FALSE, # mom_milk_T7_delplace_home vs mom_milk_T7_delplace_hospital
  FALSE, # mom_milk_T8_delplace_home vs mom_milk_T8_delplace_hospital
  FALSE, # mom_milk_T4_feeding_BF vs mom_milk_T4_feeding_MF
  FALSE, # mom_milk_T6_feeding_BF vs mom_milk_T6_feeding_MF
  FALSE, # mom_milk_T7_feeding_BF vs mom_milk_T7_feeding_MF
  FALSE, # mom_milk_T8_feeding_BF vs mom_milk_T8_feeding_MF
  FALSE, # mom_milk_T4_preterm_yes vs mom_milk_T4_preterm_no
  FALSE, # mom_milk_T6_preterm_yes vs mom_milk_T6_preterm_no
  FALSE, # mom_milk_T7_preterm_yes vs mom_milk_T7_preterm_no
  FALSE, # mom_milk_T8_preterm_yes vs mom_milk_T8_preterm_no
  FALSE, # mom_milk_T4_CDrisk_yes vs mom_milk_T4_CDrisk_no
  FALSE, # mom_milk_T6_CDrisk_yes vs mom_milk_T6_CDrisk_no
  FALSE, # mom_milk_T7_CDrisk_yes vs mom_milk_T7_CDrisk_no
  FALSE, # mom_milk_T8_CDrisk_yes vs mom_milk_T8_CDrisk_no
  FALSE, # mom_milk_T4_lockdown_before vs mom_milk_T4_lockdown_after
  FALSE, # mom_milk_T6_lockdown_before vs mom_milk_T6_lockdown_after
  FALSE, # mom_milk_T7_lockdown_before vs mom_milk_T7_lockdown_after
  FALSE, # mom_milk_T8_lockdown_before vs mom_milk_T8_lockdown_after
  FALSE, # mom_milk_T4_smoking_yes vs mom_milk_T4_smoking_no
  FALSE, # mom_milk_T6_smoking_yes vs mom_milk_T6_smoking_no
  FALSE, # mom_milk_T7_smoking_yes vs mom_milk_T7_smoking_no
  FALSE, # mom_milk_T8_smoking_yes vs mom_milk_T8_smoking_no
  FALSE, # mom_milk_T4_sex_male vs mom_milk_T4_sex_female
  FALSE, # mom_milk_T6_sex_male vs mom_milk_T6_sex_female
  FALSE, # mom_milk_T7_sex_male vs mom_milk_T7_sex_female
  FALSE  # mom_milk_T8_sex_male vs mom_milk_T8_sex_female
)

# Extended comparisons for analyses: add mom_kid_serum_T2_* fallbacks alongside kid_serum_T2_*.
kid_t2_mask <- vapply(comparisons, function(x) any(grepl("^kid_serum_T2", x)), logical(1))
comparisons2 <- comparisons
if (any(kid_t2_mask)) {
  comparisons2_extra <- lapply(
    comparisons[kid_t2_mask],
    function(x) sub("^kid_serum_T2", "mom_kid_serum_T2", x)
  )
  drop_mask <- vapply(
    comparisons2_extra,
    function(x) identical(x, c("mom_serum_T2", "mom_kid_serum_T2")),
    logical(1)
  )
  comparisons2_extra <- comparisons2_extra[!drop_mask]
  comparisons2 <- c(comparisons, comparisons2_extra)
}
if (any(kid_t2_mask)) {
  longitudinal2_extra <- longitudinal[kid_t2_mask]
  longitudinal2_extra <- longitudinal2_extra[!drop_mask]
  longitudinal2 <- c(longitudinal, longitudinal2_extra)
} else {
  longitudinal2 <- longitudinal
}
