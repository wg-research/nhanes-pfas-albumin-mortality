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
models <- c("M1", "M3")

# First verify that the lightweight fitter used within multiple imputation is
# numerically identical to svycoxph on the same complete data.
agreement <- list()
for (analyte in names(analytes)) {
  reference <- fit_selected_cox(design, analytes[[analyte]], models)
  for (model in models) {
    light <- fit_target_survey_cox(design, analytes[[analyte]], model)
    agreement[[paste(analyte, model, sep = "_")]] <- data.frame(
      analyte = analyte,
      model = model,
      beta_svycoxph = unname(stats::coef(reference[[model]])[[analytes[[analyte]]]]),
      beta_mi_fitter = light$beta,
      variance_svycoxph = unname(stats::vcov(reference[[model]])[
        analytes[[analyte]], analytes[[analyte]]
      ]),
      variance_mi_fitter = light$variance
    )
  }
}
agreement <- do.call(rbind, agreement)
rownames(agreement) <- NULL
agreement$absolute_beta_difference <- abs(
  agreement$beta_svycoxph - agreement$beta_mi_fitter
)
agreement$relative_variance_difference <- abs(
  agreement$variance_svycoxph / agreement$variance_mi_fitter - 1
)
if (max(agreement$absolute_beta_difference) > 1e-12 ||
    max(agreement$relative_variance_difference) > 1e-10) {
  stop("The lightweight survey-Cox fitter does not reproduce svycoxph")
}

# Independently refit all six models under a stratified delete-one-PSU
# jackknife.  This validates the off-diagonal covariance that determines the
# PFAS-minus-metal tests rather than only rechecking marginal model variances.
replicate_design <- survey::as.svrepdesign(design, type = "JKn", mse = TRUE)
replicate_call <- 0L
warning_records <- list()
replicate_statistic <- function(weights, variables) {
  replicate_call <<- replicate_call + 1L
  variables$.replicate_weight <- as.numeric(weights)
  estimates <- numeric()
  for (analyte in names(analytes)) {
    for (model in models) {
      formula <- stats::as.formula(paste0(
        "survival::Surv(followup_years, death_all) ~ ",
        cox_rhs(analytes[[analyte]], model, include_cycle = TRUE)
      ))
      audited_fit <- fit_audited_replicate_cox(
        formula, variables, target_term = analytes[[analyte]]
      )
      if (nrow(audited_fit$warnings)) {
        warning_records[[length(warning_records) + 1L]] <<- cbind(
          data.frame(
            replicate_call = replicate_call, analyte = analyte, model = model,
            stringsAsFactors = FALSE
          ),
          audited_fit$warnings
        )
      }
      estimates[[paste(analyte, model, sep = "_")]] <- unname(
        audited_fit$coefficients[[analytes[[analyte]]]]
      )
    }
  }
  deltas <- c(
    delta_pfas = estimates[["pfas_M3"]] - estimates[["pfas_M1"]],
    delta_lead = estimates[["lead_M3"]] - estimates[["lead_M1"]],
    delta_cadmium = estimates[["cadmium_M3"]] - estimates[["cadmium_M1"]]
  )
  c(
    estimates, deltas,
    pfas_minus_lead = deltas[["delta_pfas"]] - deltas[["delta_lead"]],
    pfas_minus_cadmium = deltas[["delta_pfas"]] - deltas[["delta_cadmium"]]
  )
}

jkn <- survey::withReplicates(replicate_design, replicate_statistic)

warning_audit <- if (length(warning_records)) {
  do.call(rbind, warning_records)
} else {
  data.frame(
    replicate_call = integer(), analyte = character(), model = character(),
    warning = character(), coefficient_index = integer(),
    affected_term = character(), target_estimate = numeric(),
    stringsAsFactors = FALSE
  )
}
jkn_estimate <- stats::coef(jkn)
jkn_se <- sqrt(diag(stats::vcov(jkn)))
jkn_output <- data.frame(
  estimand = names(jkn_estimate),
  estimate = unname(jkn_estimate),
  jkn_se = unname(jkn_se),
  jkn_ci_low = unname(jkn_estimate - stats::qnorm(0.975) * jkn_se),
  jkn_ci_high = unname(jkn_estimate + stats::qnorm(0.975) * jkn_se),
  jkn_p_value = 2 * stats::pnorm(abs(jkn_estimate / jkn_se), lower.tail = FALSE),
  design_degrees_freedom = survey::degf(design),
  n_design_strata = length(unique(as.character(design$strata[, 1]))),
  n_design_psu = length(unique(interaction(
    design$strata[, 1], design$cluster[, 1], drop = TRUE
  )))
)

primary <- utils::read.csv(file.path(paths$results, "primary_cox_models.csv"))
primary <- primary[
  primary$cohort == "overlap_primary" & primary$model %in% models,
]
primary$estimand <- paste(primary$analyte, primary$model, sep = "_")
linearized <- primary[, c("estimand", "beta", "robust_se", "p_value")]
names(linearized)[2:4] <- c(
  "linearized_estimate", "linearized_se", "linearized_p_value"
)

contrasts <- utils::read.csv(file.path(paths$results, "coefficient_contrasts.csv"))
contrasts <- contrasts[contrasts$contrast == "albumin_conditional_egfr",]
contrasts$estimand <- paste0("delta_", contrasts$analyte)
contrast_linearized <- contrasts[, c("estimand", "estimate", "robust_se", "p_value")]
names(contrast_linearized)[2:4] <- c(
  "linearized_estimate", "linearized_se", "linearized_p_value"
)

between <- utils::read.csv(file.path(paths$results, "between_analyte_tests.csv"))
between$estimand <- sub("pfas_vs_", "pfas_minus_", between$analyte)
between_linearized <- between[, c("estimand", "estimate", "robust_se", "p_value")]
names(between_linearized)[2:4] <- c(
  "linearized_estimate", "linearized_se", "linearized_p_value"
)

linearized <- rbind(linearized, contrast_linearized, between_linearized)
jkn_output <- merge(jkn_output, linearized, by = "estimand", all.x = TRUE, sort = FALSE)
jkn_output$relative_se_difference <-
  jkn_output$jkn_se / jkn_output$linearized_se - 1
jkn_output$jkn_p_holm <- NA_real_
between_rows <- jkn_output$estimand %in% c(
  "pfas_minus_lead", "pfas_minus_cadmium"
)
jkn_output$jkn_p_holm[between_rows] <- stats::p.adjust(
  jkn_output$jkn_p_value[between_rows], method = "holm"
)

utils::write.csv(
  agreement, file.path(paths$results, "lightweight_fitter_agreement.csv"), row.names = FALSE
)
utils::write.csv(
  jkn_output, file.path(paths$results, "core_covariance_jkn_check.csv"), row.names = FALSE
)
utils::write.csv(
  warning_audit,
  file.path(paths$audit, "core_covariance_jkn_warning_audit.csv"),
  row.names = FALSE
)
message(
  "Core covariance validation complete; largest absolute relative SE difference = ",
  signif(max(abs(jkn_output$relative_se_difference), na.rm = TRUE), 4)
)
