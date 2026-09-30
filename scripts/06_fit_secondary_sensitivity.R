#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))

collect_models <- function(fits, exposure, analyte, cohort, outcome, analysis) {
  output <- do.call(rbind, lapply(names(fits), function(model) {
    extract_cox_result(fits[[model]], exposure, analyte, cohort, outcome, model)
  }))
  output$analysis <- analysis
  output
}

pfas_design <- make_domain_design(
  data, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
overlap_design <- make_domain_design(
  data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
)
lead_design <- make_domain_design(
  data, "domain_metals_source", "cohort_lead_full", "weight_metals_all"
)
cadmium_design <- make_domain_design(
  data, "domain_metals_source", "cohort_cadmium_full", "weight_metals_all"
)
pfas_2013_design <- make_domain_design(
  data, "domain_pfas_2013_source", "cohort_pfas_2013", "weight_pfas_2013"
)

fit_unadjusted_cox <- function(design, exposure) {
  survey::svycoxph(
    stats::as.formula(paste0(
      "survival::Surv(followup_years, death_all) ~ ", exposure
    )),
    design = design, model = TRUE, x = TRUE
  )
}

crude_models <- list(
  fit_unadjusted_cox(pfas_design, "log2_sum4_pfas"),
  fit_unadjusted_cox(overlap_design, "log2_sum4_pfas"),
  fit_unadjusted_cox(overlap_design, "log2_lead"),
  fit_unadjusted_cox(overlap_design, "log2_cadmium")
)
crude <- do.call(rbind, list(
  extract_cox_result(crude_models[[1]], "log2_sum4_pfas", "pfas", "pfas_primary", "all_cause", "Unadjusted"),
  extract_cox_result(crude_models[[2]], "log2_sum4_pfas", "pfas", "overlap_primary", "all_cause", "Unadjusted"),
  extract_cox_result(crude_models[[3]], "log2_lead", "lead", "overlap_primary", "all_cause", "Unadjusted"),
  extract_cox_result(crude_models[[4]], "log2_cadmium", "cadmium", "overlap_primary", "all_cause", "Unadjusted")
))

results <- list()
results[["lead_full"]] <- collect_models(
  fit_symmetric_cox(lead_design, "log2_lead"), "log2_lead", "lead",
  "lead_full", "all_cause", "full_metals"
)
results[["cadmium_full"]] <- collect_models(
  fit_symmetric_cox(cadmium_design, "log2_cadmium"), "log2_cadmium", "cadmium",
  "cadmium_full", "all_cause", "full_metals"
)
for (outcome in c("death_cvd", "death_cancer")) {
  outcome_label <- sub("death_", "", outcome)
  results[[outcome]] <- collect_models(
    fit_symmetric_cox(
      pfas_design, "log2_sum4_pfas", outcome = outcome, bmi_shape = "linear"
    ),
    "log2_sum4_pfas", "pfas", "pfas_primary", outcome_label, "secondary_outcome"
  )
}

results[["linear_physiology"]] <- collect_models(
  fit_symmetric_cox(
    pfas_design, "log2_sum4_pfas", physiology = "linear"
  ), "log2_sum4_pfas", "pfas", "pfas_primary", "all_cause",
  "linear_albumin_egfr"
)

healthy_design <- subset(
  pfas_design,
  !is.na(restriction_no_baseline_disease) & restriction_no_baseline_disease
)
results[["healthy_restriction"]] <- collect_models(
  fit_symmetric_cox(healthy_design, "log2_sum4_pfas"),
  "log2_sum4_pfas", "pfas", "pfas_primary_healthy", "all_cause",
  "exclude_baseline_diabetes_cvd_cancer"
)

not_pregnant_design <- subset(pfas_design, restriction_not_pregnant)
results[["pregnancy_restriction"]] <- collect_models(
  fit_symmetric_cox(not_pregnant_design, "log2_sum4_pfas"),
  "log2_sum4_pfas", "pfas", "pfas_primary_not_pregnant", "all_cause",
  "exclude_known_pregnancy"
)

secondary <- do.call(rbind, results)
rownames(secondary) <- NULL

loo <- list()
cycles <- levels(droplevels(pfas_design$variables$cycle_factor))
for (cycle in cycles) {
  keep <- as.character(pfas_design$variables$cycle_factor) != cycle
  design_loo <- pfas_design[keep, ]
  design_loo$variables$cycle_factor <- droplevels(design_loo$variables$cycle_factor)
  fits <- fit_symmetric_cox(design_loo, "log2_sum4_pfas")
  item <- do.call(rbind, lapply(c("M0", "M3"), function(model) {
    extract_cox_result(
      fits[[model]], "log2_sum4_pfas", "pfas", "pfas_leave_one_cycle_out",
      "all_cause", model
    )
  }))
  item$excluded_cycle <- cycle
  loo[[cycle]] <- item
}
loo <- do.call(rbind, loo)
rownames(loo) <- NULL

ci_null <- secondary$ci_low <= 1 & secondary$ci_high >= 1
if (any(ci_null != (secondary$p_value >= 0.05))) {
  stop("Secondary CI and p-value consistency check failed")
}

utils::write.csv(secondary, file.path(paths$results, "secondary_sensitivity_models.csv"), row.names = FALSE)
utils::write.csv(crude, file.path(paths$results, "crude_cox_models.csv"), row.names = FALSE)
utils::write.csv(loo, file.path(paths$results, "leave_one_cycle_out.csv"), row.names = FALSE)

mean_2013 <- survey::svymean(~log2_sum4_pfas, pfas_2013_design, na.rm = TRUE)
quantile_2013 <- survey::svyquantile(
  ~log2_sum4_pfas, pfas_2013_design, quantiles = c(0.25, 0.5, 0.75),
  ci = FALSE, na.rm = TRUE
)
descriptive_2013 <- data.frame(
  cohort = "pfas_2013_surplus_specimen", n = nrow(pfas_2013_design$variables),
  deaths_all = sum(pfas_2013_design$variables$death_all),
  weighted_mean_log2_sum4 = unname(stats::coef(mean_2013)),
  robust_se_mean = sqrt(unname(stats::vcov(mean_2013)[1, 1])),
  q25_log2_sum4 = as.numeric(quantile_2013[[1]])[[1]],
  median_log2_sum4 = as.numeric(quantile_2013[[1]])[[2]],
  q75_log2_sum4 = as.numeric(quantile_2013[[1]])[[3]]
)
utils::write.csv(
  descriptive_2013, file.path(paths$results, "pfas_2013_descriptive.csv"), row.names = FALSE
)
message("Secondary and sensitivity models complete")
