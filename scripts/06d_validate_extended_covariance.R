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
        cox_rhs(
          analytes[[analyte]], model, physiology = "linear",
          include_cycle = TRUE
        )
      ))
      audited_fit <- fit_audited_replicate_cox(
        formula, variables, target_term = analytes[[analyte]],
        allowed_separation_terms = c(
          "cycle_factor2009-2010", "cycle_factor2017-2018",
          "splines::ns(bmi, df = 3)3"
        )
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
  delta_pfas <- estimates[["pfas_M3"]] - estimates[["pfas_M1"]]
  delta_lead <- estimates[["lead_M3"]] - estimates[["lead_M1"]]
  delta_cadmium <- estimates[["cadmium_M3"]] - estimates[["cadmium_M1"]]
  c(
    pfas_minus_lead = delta_pfas - delta_lead,
    pfas_minus_cadmium = delta_pfas - delta_cadmium
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
output <- data.frame(
  estimand = names(jkn_estimate),
  estimate = unname(jkn_estimate),
  jkn_se = unname(jkn_se),
  jkn_ci_low = unname(jkn_estimate - stats::qnorm(0.975) * jkn_se),
  jkn_ci_high = unname(jkn_estimate + stats::qnorm(0.975) * jkn_se),
  jkn_p_value = 2 * stats::pnorm(abs(jkn_estimate / jkn_se), lower.tail = FALSE),
  stringsAsFactors = FALSE
)

linearized <- utils::read.csv(file.path(
  paths$results, "linear_overlap_between_analyte_tests.csv"
))
linearized$estimand <- sub(
  "albumin_conditional_", "", linearized$contrast, fixed = TRUE
)
linearized <- linearized[, c("estimand", "estimate", "robust_se", "p_value")]
names(linearized)[2:4] <- c(
  "linearized_estimate", "linearized_se", "linearized_p_value"
)
output <- merge(output, linearized, by = "estimand", all.x = TRUE, sort = FALSE)
output$relative_se_difference <- output$jkn_se / output$linearized_se - 1
output$jkn_p_holm <- stats::p.adjust(output$jkn_p_value, method = "holm")

if (any(sign(output$estimate) != sign(output$linearized_estimate)) ||
    any(abs(output$relative_se_difference) >= 0.10) ||
    any(output$jkn_p_holm >= 0.05)) {
  stop("Linear-physiology jackknife validation did not reproduce the core conclusions")
}

utils::write.csv(
  output, file.path(paths$results, "linear_overlap_jkn_check.csv"),
  row.names = FALSE
)
utils::write.csv(
  warning_audit,
  file.path(paths$audit, "linear_overlap_jkn_warning_audit.csv"),
  row.names = FALSE
)
message(
  "Linear-physiology covariance validation complete; largest absolute relative ",
  "SE difference = ", signif(max(abs(output$relative_se_difference)), 4)
)
