#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
table_dir <- file.path(paths$results, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))

# Keep displayed p values consistent across all generated tables and the manuscript.
format_p_value <- function(p) {
  vapply(p, function(value) {
    if (is.na(value)) return(NA_character_)
    if (value < 0.001) format(value, scientific = TRUE, digits = 2, trim = TRUE)
    else sprintf("%.3f", value)
  }, character(1))
}

pfas_design <- make_domain_design(
  data, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
overlap_design <- make_domain_design(
  data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
)
metals_design <- make_domain_design(
  data, "domain_metals_source", "cohort_lead_full", "weight_metals_all"
)

mean_row <- function(design, variable, label, cohort, digits = 1) {
  estimate <- survey::svymean(
    stats::as.formula(paste0("~", variable)), design, na.rm = TRUE
  )
  value <- unname(stats::coef(estimate)[[1]])
  se <- sqrt(unname(stats::vcov(estimate)[1, 1]))
  data.frame(
    cohort = cohort, characteristic = label, level = "Mean (SE)",
    estimate = value, robust_se = se,
    display = sprintf(paste0("%.", digits, "f (%.", digits, "f)"), value, se),
    stringsAsFactors = FALSE
  )
}

proportion_rows <- function(design, variable, label, cohort) {
  levels <- levels(design$variables[[variable]])
  do.call(rbind, lapply(levels, function(level) {
    indicator <- as.numeric(design$variables[[variable]] == level)
    estimate <- survey::svymean(~indicator, design, na.rm = TRUE)
    value <- 100 * unname(stats::coef(estimate)[[1]])
    se <- 100 * sqrt(unname(stats::vcov(estimate)[1, 1]))
    data.frame(
      cohort = cohort, characteristic = label, level = level,
      estimate = value, robust_se = se,
      display = sprintf("%.1f (%.1f)", value, se), stringsAsFactors = FALSE
    )
  }))
}

count_rows <- function(design, cohort) {
  data.frame(
    cohort = cohort,
    characteristic = c("Participants", "All-cause deaths"),
    level = "Unweighted n",
    estimate = c(nrow(design$variables), sum(design$variables$death_all)),
    robust_se = NA_real_,
    display = format(c(nrow(design$variables), sum(design$variables$death_all)),
                     big.mark = ",", scientific = FALSE),
    stringsAsFactors = FALSE
  )
}

baseline_for_design <- function(design, cohort) {
  rbind(
    count_rows(design, cohort),
    mean_row(design, "age_years", "Age, years", cohort),
    proportion_rows(design, "sex", "Sex, %", cohort),
    proportion_rows(design, "race_ethnicity", "Race and ethnicity, %", cohort),
    proportion_rows(design, "education", "Education, %", cohort),
    mean_row(design, "pir", "Poverty-income ratio", cohort, 2),
    proportion_rows(design, "smoking", "Smoking status, %", cohort),
    mean_row(design, "bmi", "Body mass index, kg/m²", cohort),
    mean_row(design, "albumin_g_dl", "Serum albumin, g/dL", cohort, 2),
    mean_row(design, "egfr", "eGFR, mL/min/1.73 m²", cohort),
    mean_row(design, "log2_sum4_pfas", "log2 Sigma4PFAS, nmol/L", cohort, 2)
  )
}

baseline <- rbind(
  baseline_for_design(pfas_design, "Primary PFAS cohort"),
  baseline_for_design(overlap_design, "PFAS-metals overlap cohort")
)
metals_baseline <- rbind(
  count_rows(metals_design, "Full metals cohort"),
  mean_row(metals_design, "age_years", "Age, years", "Full metals cohort"),
  proportion_rows(metals_design, "sex", "Sex, %", "Full metals cohort"),
  proportion_rows(metals_design, "race_ethnicity", "Race and ethnicity, %", "Full metals cohort"),
  proportion_rows(metals_design, "education", "Education, %", "Full metals cohort"),
  mean_row(metals_design, "pir", "Poverty-income ratio", "Full metals cohort", 2),
  proportion_rows(metals_design, "smoking", "Smoking status, %", "Full metals cohort"),
  mean_row(metals_design, "bmi", "Body mass index, kg/m²", "Full metals cohort"),
  mean_row(metals_design, "albumin_g_dl", "Serum albumin, g/dL", "Full metals cohort", 3),
  mean_row(metals_design, "egfr", "eGFR, mL/min/1.73 m²", "Full metals cohort")
)

row_key <- unique(baseline[, c("characteristic", "level")])
row_key$row_order <- seq_len(nrow(row_key))
table1 <- merge(row_key, baseline[, c("cohort", "characteristic", "level", "display")],
                by = c("characteristic", "level"), all.x = TRUE, sort = FALSE)
table1 <- reshape(
  table1, idvar = c("row_order", "characteristic", "level"),
  timevar = "cohort", direction = "wide"
)
table1 <- table1[order(table1$row_order), ]
names(table1) <- sub("display\\.", "", names(table1))

cox <- utils::read.csv(file.path(paths$results, "primary_cox_models.csv"))
cox$model <- factor(cox$model, levels = c("M0", "M1", "M2", "M3"))
table2 <- cox[, c(
  "cohort", "analyte", "model", "n", "events", "beta", "robust_se",
  "hr", "ci_low", "ci_high", "p_value"
)]
table2$hr_95_ci <- sprintf("%.3f (%.3f-%.3f)", table2$hr, table2$ci_low, table2$ci_high)
table2$p_display <- format_p_value(table2$p_value)

crude <- utils::read.csv(file.path(paths$results, "crude_cox_models.csv"))
crude$hr_95_ci <- sprintf("%.3f (%.3f-%.3f)", crude$hr, crude$ci_low, crude$ci_high)
crude$p_display <- format_p_value(crude$p_value)

contrasts <- utils::read.csv(file.path(paths$results, "coefficient_contrasts.csv"))
between <- utils::read.csv(file.path(paths$results, "between_analyte_tests.csv"))
structure <- utils::read.csv(file.path(paths$results, "albumin_exposure_structure.csv"))
contrasts$estimate_95_ci <- sprintf(
  "%.4f (%.4f to %.4f)", contrasts$estimate, contrasts$ci_low, contrasts$ci_high
)
between$estimate_95_ci <- sprintf(
  "%.4f (%.4f to %.4f)", between$estimate, between$ci_low, between$ci_high
)
table3 <- merge(
  contrasts, structure[, c("analyte", "beta_albumin", "deviance_incremental_r2")],
  by = "analyte", all.x = TRUE, sort = FALSE
)

linear_models <- utils::read.csv(file.path(paths$results, "linear_overlap_models.csv"))
linear_models$hr_95_ci <- sprintf(
  "%.3f (%.3f-%.3f)",
  linear_models$hr, linear_models$ci_low, linear_models$ci_high
)
linear_models$p_display <- format_p_value(linear_models$p_value)

linear_contrasts <- utils::read.csv(file.path(
  paths$results, "linear_overlap_contrasts.csv"
))
linear_between <- utils::read.csv(file.path(
  paths$results, "linear_overlap_between_analyte_tests.csv"
))
linear_contrasts$comparison_type <- "Within-analyte coefficient change"
linear_contrasts$p_holm <- NA_real_
linear_between$comparison_type <- "PFAS-minus-metal difference"
linear_table_contrasts <- rbind(
  linear_contrasts[, c(
    "comparison_type", "analyte", "contrast", "estimate", "robust_se",
    "ci_low", "ci_high", "p_value", "p_holm"
  )],
  linear_between[, c(
    "comparison_type", "analyte", "contrast", "estimate", "robust_se",
    "ci_low", "ci_high", "p_value", "p_holm"
  )]
)
linear_table_contrasts$estimate_95_ci <- sprintf(
  "%.4f (%.4f to %.4f)",
  linear_table_contrasts$estimate,
  linear_table_contrasts$ci_low,
  linear_table_contrasts$ci_high
)

congener_models <- utils::read.csv(file.path(paths$results, "congener_models.csv"))
congener_models$hr_95_ci <- sprintf(
  "%.3f (%.3f-%.3f)",
  congener_models$hr, congener_models$ci_low, congener_models$ci_high
)
congener_models$p_display <- format_p_value(congener_models$p_value)

congener_contrasts <- utils::read.csv(file.path(
  paths$results, "congener_contrasts.csv"
))
congener_contrasts$estimate_95_ci <- sprintf(
  "%.4f (%.4f to %.4f)",
  congener_contrasts$estimate,
  congener_contrasts$ci_low,
  congener_contrasts$ci_high
)
congener_structure <- utils::read.csv(file.path(
  paths$results, "congener_albumin_structure.csv"
))
congener_structure$beta_95_ci <- sprintf(
  "%.3f (%.3f to %.3f)",
  congener_structure$beta_albumin,
  congener_structure$ci_low,
  congener_structure$ci_high
)

structural <- utils::read.csv(file.path(
  paths$results, "structural_contrast_sensitivity.csv"
))
scenario_labels <- c(
  primary_per_doubling = "Primary, per doubling",
  weighted_iqr_log2 = "Per weighted IQR of log2 exposure",
  reduced_demographic_base = "Reduced demographic base",
  pre_2015_overlap_cycles = "2003–2012 overlap cycles"
)
scenario_domains <- c(
  primary_per_doubling = "Reference",
  weighted_iqr_log2 = "Exposure scaling",
  reduced_demographic_base = "Adjustment specification",
  pre_2015_overlap_cycles = "Cycle restriction"
)
format_contrast <- function(row) {
  sprintf("%.4f (%.4f to %.4f)", row$estimate, row$ci_low, row$ci_high)
}
structural_table_rows <- lapply(names(scenario_labels), function(scenario) {
  rows <- structural[structural$scenario == scenario, ]
  get_row <- function(contrast) {
    result <- rows[rows$contrast == contrast, ]
    if (nrow(result) != 1L) stop("Missing structural contrast: ", scenario, "/", contrast)
    result
  }
  pfas <- get_row("delta_pfas")
  lead <- get_row("delta_lead")
  cadmium <- get_row("delta_cadmium")
  pfas_lead <- get_row("pfas_minus_lead")
  pfas_cadmium <- get_row("pfas_minus_cadmium")
  data.frame(
    domain = unname(scenario_domains[[scenario]]),
    specification = unname(scenario_labels[[scenario]]),
    n = pfas$n, events = pfas$events,
    pfas_delta = format_contrast(pfas),
    lead_delta = format_contrast(lead),
    cadmium_delta = format_contrast(cadmium),
    pfas_minus_lead = format_contrast(pfas_lead),
    pfas_minus_lead_holm_p = pfas_lead$p_holm_within_scenario,
    pfas_minus_cadmium = format_contrast(pfas_cadmium),
    pfas_minus_cadmium_holm_p = pfas_cadmium$p_holm_within_scenario,
    stringsAsFactors = FALSE
  )
})
table_s12 <- do.call(rbind, structural_table_rows)

# Freeze the supplementary numeric sources that are assembled from diagnostics.
leave_one_cycle <- utils::read.csv(file.path(paths$results, "leave_one_cycle_out.csv"))
leave_one_cycle$hr_95_ci <- sprintf(
  "%.3f (%.3f, %.3f)",
  leave_one_cycle$hr, leave_one_cycle$ci_low, leave_one_cycle$ci_high
)
leave_one_cycle$p_display <- format_p_value(leave_one_cycle$p_value)

mi_missingness <- utils::read.csv(file.path(paths$results, "mi_missingness.csv"))

covariance_jkn <- utils::read.csv(
  file.path(paths$results, "core_covariance_jkn_check.csv")
)
covariance_jkn$jkn_p_display <- format_p_value(covariance_jkn$jkn_p_value)

targeted_sensitivities <- utils::read.csv(
  file.path(paths$results, "targeted_contrast_sensitivity.csv")
)
targeted_sensitivities$p_display <- format_p_value(targeted_sensitivities$p_value)
targeted_sensitivities$p_holm_display <- format_p_value(
  targeted_sensitivities$p_holm_within_scenario
)

secondary_sensitivity <- utils::read.csv(
  file.path(paths$results, "secondary_sensitivity_models.csv")
)
secondary_sensitivity$hr_95_ci <- sprintf(
  "%.3f (%.3f, %.3f)",
  secondary_sensitivity$hr,
  secondary_sensitivity$ci_low,
  secondary_sensitivity$ci_high
)
secondary_sensitivity$p_display <- format_p_value(secondary_sensitivity$p_value)

utils::write.csv(baseline, file.path(paths$results, "baseline_characteristics_long.csv"), row.names = FALSE)
utils::write.csv(table1, file.path(table_dir, "table1_baseline.csv"), row.names = FALSE)
utils::write.csv(table2, file.path(table_dir, "table2_primary_models.csv"), row.names = FALSE)
utils::write.csv(crude, file.path(table_dir, "table_s13_unadjusted_models.csv"), row.names = FALSE)
utils::write.csv(table3, file.path(table_dir, "table3_coefficient_contrasts.csv"), row.names = FALSE)
utils::write.csv(between, file.path(table_dir, "table3_between_analyte_tests.csv"), row.names = FALSE)
utils::write.csv(metals_baseline, file.path(table_dir, "table_s2_full_metals_baseline.csv"), row.names = FALSE)
utils::write.csv(
  linear_models,
  file.path(table_dir, "table_s10a_linear_overlap_models.csv"),
  row.names = FALSE
)
utils::write.csv(
  linear_table_contrasts,
  file.path(table_dir, "table_s10b_linear_overlap_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  congener_models,
  file.path(table_dir, "table_s11a_congener_models.csv"),
  row.names = FALSE
)
utils::write.csv(
  congener_contrasts,
  file.path(table_dir, "table_s11b_congener_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  congener_structure,
  file.path(table_dir, "table_s11c_congener_structure.csv"),
  row.names = FALSE
)
utils::write.csv(
  table_s12,
  file.path(table_dir, "table_s12_structural_contrast_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  leave_one_cycle,
  file.path(table_dir, "table_s7b_leave_one_cycle_out.csv"),
  row.names = FALSE
)
utils::write.csv(
  mi_missingness,
  file.path(table_dir, "table_s8a_mi_missingness.csv"),
  row.names = FALSE
)
utils::write.csv(
  covariance_jkn,
  file.path(table_dir, "table_s8b_covariance_validation.csv"),
  row.names = FALSE
)
utils::write.csv(
  targeted_sensitivities,
  file.path(table_dir, "table_s8b_targeted_sensitivities.csv"),
  row.names = FALSE
)
utils::write.csv(
  secondary_sensitivity,
  file.path(table_dir, "table_s3_secondary_sensitivity_models.csv"),
  row.names = FALSE
)

message("Main and supplementary table sources complete")
