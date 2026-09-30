#!/usr/bin/env Rscript

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()

steps <- c(
  "scripts/00_environment_specification.R",
  "scripts/01_audit_source_fields.R",
  "scripts/02_extract_source_data.R",
  "scripts/03_harmonize_variables.R",
  "scripts/04_build_cohorts_designs.R",
  "scripts/05_fit_primary_models.R",
  "scripts/05b_validate_core_covariance.R",
  "scripts/06_fit_secondary_sensitivity.R",
  "scripts/06b_targeted_sensitivity.R",
  "scripts/06c_linear_overlap_contrasts.R",
  "scripts/06d_validate_extended_covariance.R",
  "scripts/06e_congener_analysis.R",
  "scripts/07_multiple_imputation.R",
  "scripts/08_diagnostics.R",
  "scripts/09_build_tables.R",
  "scripts/10_build_figures.R"
)

for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", step))
  if (!identical(status, 0L)) stop("Pipeline stopped at: ", step)
}

message("Completed available rebuild stages")
