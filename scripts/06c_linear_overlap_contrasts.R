#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))
design <- make_domain_design(
  data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
)
analytes <- c(
  pfas = "log2_sum4_pfas", lead = "log2_lead", cadmium = "log2_cadmium"
)

fits <- lapply(analytes, function(exposure) {
  fit_symmetric_cox(design, exposure, physiology = "linear")
})
model_results <- do.call(rbind, unlist(lapply(names(fits), function(analyte) {
  lapply(names(fits[[analyte]]), function(model) {
    extract_cox_result(
      fits[[analyte]][[model]], analytes[[analyte]], analyte,
      "overlap_primary", "all_cause", model
    )
  })
}), recursive = FALSE))
model_results$physiology_specification <- "linear"

fit_vector <- list()
exposure_vector <- labels <- character()
for (analyte in names(fits)) {
  for (model in names(fits[[analyte]])) {
    fit_vector[[length(fit_vector) + 1L]] <- fits[[analyte]][[model]]
    exposure_vector <- c(exposure_vector, analytes[[analyte]])
    labels <- c(labels, paste(analyte, model, sep = "_"))
  }
}
covariance <- joint_target_covariance(fit_vector, exposure_vector, labels)
estimates <- stats::setNames(vapply(seq_along(fit_vector), function(i) {
  unname(stats::coef(fit_vector[[i]])[[exposure_vector[[i]]]])
}, numeric(1)), labels)

contrasts <- coefficient_contrasts(estimates, covariance, names(analytes))
contrasts$cohort <- "overlap_primary"
contrasts$physiology_specification <- "linear"

between <- do.call(rbind, lapply(c("lead", "cadmium"), function(metal) {
  weights <- stats::setNames(rep(0, length(estimates)), names(estimates))
  weights[c("pfas_M1", "pfas_M3")] <- c(-1, 1)
  weights[paste0(metal, c("_M1", "_M3"))] <- c(1, -1)
  result <- linear_contrast(
    estimates, covariance, weights,
    paste0("albumin_conditional_pfas_minus_", metal)
  )
  result$analyte <- paste("pfas", metal, sep = "_vs_")
  result$cohort <- "overlap_primary"
  result$physiology_specification <- "linear"
  result
}))
between$p_holm <- stats::p.adjust(between$p_value, method = "holm")

if (nrow(design$variables) != 8705L || sum(design$variables$death_all) != 1266L) {
  stop("Linear overlap analysis did not retain the frozen cohort and events")
}
if (any((model_results$ci_low <= 1 & model_results$ci_high >= 1) !=
        (model_results$p_value >= 0.05))) {
  stop("Linear overlap CI and p-value consistency check failed")
}

utils::write.csv(
  model_results, file.path(paths$results, "linear_overlap_models.csv"),
  row.names = FALSE
)
utils::write.csv(
  contrasts, file.path(paths$results, "linear_overlap_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  between, file.path(paths$results, "linear_overlap_between_analyte_tests.csv"),
  row.names = FALSE
)
utils::write.csv(
  as.data.frame(covariance),
  file.path(paths$results, "linear_overlap_joint_covariance.csv")
)
message("Linear-physiology overlap analysis complete")
