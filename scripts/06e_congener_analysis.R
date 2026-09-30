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
congeners <- c(
  pfoa = "log2_pfoa", pfos = "log2_pfos",
  pfhxs = "log2_pfhxs", pfna = "log2_pfna"
)

fits <- lapply(congeners, function(exposure) fit_symmetric_cox(design, exposure))
model_results <- do.call(rbind, unlist(lapply(names(fits), function(congener) {
  lapply(names(fits[[congener]]), function(model) {
    extract_cox_result(
      fits[[congener]][[model]], congeners[[congener]], congener,
      "overlap_primary", "all_cause", model
    )
  })
}), recursive = FALSE))

fit_vector <- list()
exposure_vector <- labels <- character()
for (congener in names(fits)) {
  for (model in names(fits[[congener]])) {
    fit_vector[[length(fit_vector) + 1L]] <- fits[[congener]][[model]]
    exposure_vector <- c(exposure_vector, congeners[[congener]])
    labels <- c(labels, paste(congener, model, sep = "_"))
  }
}
covariance <- joint_target_covariance(fit_vector, exposure_vector, labels)
estimates <- stats::setNames(vapply(seq_along(fit_vector), function(i) {
  unname(stats::coef(fit_vector[[i]])[[exposure_vector[[i]]]])
}, numeric(1)), labels)
contrasts <- coefficient_contrasts(estimates, covariance, names(congeners))
contrasts$cohort <- "overlap_primary"
contrasts$p_holm <- NA_real_
albumin_rows <- contrasts$contrast == "albumin_conditional_egfr"
contrasts$p_holm[albumin_rows] <- stats::p.adjust(
  contrasts$p_value[albumin_rows], method = "holm"
)

structure_results <- do.call(rbind, lapply(names(congeners), function(congener) {
  fit_exposure_albumin_structure(design, congeners[[congener]], congener)
}))

if (nrow(design$variables) != 8705L || sum(design$variables$death_all) != 1266L) {
  stop("Congener analysis did not retain the frozen overlap cohort and events")
}
if (any((model_results$ci_low <= 1 & model_results$ci_high >= 1) !=
        (model_results$p_value >= 0.05))) {
  stop("Congener CI and p-value consistency check failed")
}
if (!all(model_results$n == 8705L) || !all(model_results$events == 1266L)) {
  stop("Congener models do not share the frozen analysis sample")
}

utils::write.csv(
  model_results, file.path(paths$results, "congener_models.csv"), row.names = FALSE
)
utils::write.csv(
  contrasts, file.path(paths$results, "congener_contrasts.csv"), row.names = FALSE
)
utils::write.csv(
  structure_results, file.path(paths$results, "congener_albumin_structure.csv"),
  row.names = FALSE
)
utils::write.csv(
  as.data.frame(covariance),
  file.path(paths$results, "congener_joint_covariance.csv")
)
message("PFAS component analysis complete.")
