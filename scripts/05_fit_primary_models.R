#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
source_file <- file.path(paths$processed, "analysis_cohorts.rds")
if (!file.exists(source_file)) stop("Run 04_build_cohorts_designs.R first")
data <- readRDS(source_file)

pfas_design <- make_domain_design(
  data, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
overlap_design <- make_domain_design(
  data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
)

pfas_primary <- fit_symmetric_cox(pfas_design, "log2_sum4_pfas")
analytes <- c(pfas = "log2_sum4_pfas", lead = "log2_lead", cadmium = "log2_cadmium")
overlap_fits <- lapply(analytes, function(exposure) {
  fit_symmetric_cox(overlap_design, exposure)
})

results <- do.call(rbind, c(
  lapply(names(pfas_primary), function(model) {
    extract_cox_result(
      pfas_primary[[model]], "log2_sum4_pfas", "pfas", "pfas_primary",
      "all_cause", model
    )
  }),
  unlist(lapply(names(overlap_fits), function(analyte) {
    lapply(names(overlap_fits[[analyte]]), function(model) {
      extract_cox_result(
        overlap_fits[[analyte]][[model]], analytes[[analyte]], analyte,
        "overlap_primary", "all_cause", model
      )
    })
  }), recursive = FALSE)
))

fit_vector <- list()
exposure_vector <- character()
labels <- character()
for (analyte in names(overlap_fits)) {
  for (model in names(overlap_fits[[analyte]])) {
    fit_vector[[length(fit_vector) + 1L]] <- overlap_fits[[analyte]][[model]]
    exposure_vector <- c(exposure_vector, analytes[[analyte]])
    labels <- c(labels, paste(analyte, model, sep = "_"))
  }
}
joint_covariance <- joint_target_covariance(fit_vector, exposure_vector, labels)

estimates <- stats::setNames(vapply(seq_along(fit_vector), function(i) {
  unname(stats::coef(fit_vector[[i]])[[exposure_vector[[i]]]])
}, numeric(1)), labels)

contrasts <- coefficient_contrasts(estimates, joint_covariance, names(overlap_fits))
contrasts$cohort <- "overlap_primary"

between <- list()
for (metal in c("lead", "cadmium")) {
  weights <- stats::setNames(rep(0, length(estimates)), names(estimates))
  weights[c("pfas_M1", "pfas_M3")] <- c(-1, 1)
  weights[paste0(metal, c("_M1", "_M3"))] <- c(1, -1)
  item <- linear_contrast(
    estimates, joint_covariance, weights,
    paste0("albumin_conditional_pfas_minus_", metal)
  )
  item$analyte <- paste("pfas", metal, sep = "_vs_")
  item$cohort <- "overlap_primary"
  between[[metal]] <- item
}
between <- do.call(rbind, between)
between$p_holm <- stats::p.adjust(between$p_value, method = "holm")

structure_results <- do.call(rbind, lapply(names(analytes), function(analyte) {
  fit_exposure_albumin_structure(overlap_design, analytes[[analyte]], analyte)
}))

ci_null <- results$ci_low <= 1 & results$ci_high >= 1
p_nonsignificant <- results$p_value >= 0.05
if (any(ci_null != p_nonsignificant)) {
  stop("CI and p-value consistency check failed")
}

utils::write.csv(results, file.path(paths$results, "primary_cox_models.csv"), row.names = FALSE)
utils::write.csv(contrasts, file.path(paths$results, "coefficient_contrasts.csv"), row.names = FALSE)
utils::write.csv(between, file.path(paths$results, "between_analyte_tests.csv"), row.names = FALSE)
utils::write.csv(structure_results, file.path(paths$results, "albumin_exposure_structure.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(joint_covariance),
                 file.path(paths$results, "joint_exposure_coefficient_covariance.csv"))

message("Primary models complete: ", nrow(results), " model estimates")
