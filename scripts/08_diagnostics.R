#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))
library(survival)

paths <- project_paths()
ensure_output_directories(paths)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))

pfas_design <- make_domain_design(
  data, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
overlap_design <- make_domain_design(
  data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
)

weighted_quantiles <- function(design, variable, probabilities) {
  formula <- stats::as.formula(paste0("~", variable))
  result <- survey::svyquantile(
    formula, design, quantiles = probabilities, ci = FALSE, na.rm = TRUE
  )
  as.numeric(result[[1]])
}

fit_rcs_diagnostic <- function(design, exposure, analyte) {
  knot_probabilities <- c(0.05, 0.275, 0.5, 0.725, 0.95)
  knots <- weighted_quantiles(design, exposure, knot_probabilities)
  basis <- Hmisc::rcspline.eval(
    design$variables[[exposure]], knots = knots, inclx = TRUE
  )
  nonlinear_names <- paste0("rcs_nonlinear_", seq_len(ncol(basis) - 1L))
  for (j in seq_along(nonlinear_names)) {
    design$variables[[nonlinear_names[[j]]]] <- basis[, j + 1L]
  }
  rhs <- paste(
    c(cox_rhs(exposure, "M3"), nonlinear_names), collapse = " + "
  )
  formula <- stats::as.formula(paste0(
    "survival::Surv(followup_years, death_all) ~ ", rhs
  ))
  fit <- survey::svycoxph(formula, design = design, model = TRUE, x = TRUE)
  beta_nonlinear <- stats::coef(fit)[nonlinear_names]
  covariance_nonlinear <- stats::vcov(fit)[nonlinear_names, nonlinear_names, drop = FALSE]
  wald <- as.numeric(t(beta_nonlinear) %*%
                       solve(covariance_nonlinear, beta_nonlinear))
  p_nonlinearity <- stats::pchisq(
    wald, df = length(nonlinear_names), lower.tail = FALSE
  )

  limits <- weighted_quantiles(design, exposure, c(0.01, 0.99))
  reference <- weighted_quantiles(design, exposure, 0.5)
  grid <- seq(limits[[1]], limits[[2]], length.out = 150)
  grid_basis <- Hmisc::rcspline.eval(grid, knots = knots, inclx = TRUE)
  reference_basis <- Hmisc::rcspline.eval(reference, knots = knots, inclx = TRUE)
  target_names <- c(exposure, nonlinear_names)
  target_beta <- stats::coef(fit)[target_names]
  target_covariance <- stats::vcov(fit)[target_names, target_names, drop = FALSE]
  differences <- sweep(grid_basis, 2, as.numeric(reference_basis), "-")
  log_hr <- as.numeric(differences %*% target_beta)
  variance <- rowSums((differences %*% target_covariance) * differences)
  curve <- data.frame(
    analyte = analyte, exposure_value = grid, reference = reference,
    hr = exp(log_hr), ci_low = exp(log_hr - stats::qnorm(0.975) * sqrt(variance)),
    ci_high = exp(log_hr + stats::qnorm(0.975) * sqrt(variance))
  )
  test <- data.frame(
    analyte = analyte, n = fit$n, events = fit$nevent,
    df_nonlinear = length(nonlinear_names), wald_chisq = wald,
    p_nonlinearity = p_nonlinearity,
    knot_1 = knots[[1]], knot_2 = knots[[2]], knot_3 = knots[[3]],
    knot_4 = knots[[4]], knot_5 = knots[[5]], reference = reference
  )
  list(test = test, curve = curve)
}

rcs <- list(
  pfas_primary = fit_rcs_diagnostic(pfas_design, "log2_sum4_pfas", "pfas_primary"),
  pfas_overlap = fit_rcs_diagnostic(overlap_design, "log2_sum4_pfas", "pfas_overlap"),
  lead = fit_rcs_diagnostic(overlap_design, "log2_lead", "lead"),
  cadmium = fit_rcs_diagnostic(overlap_design, "log2_cadmium", "cadmium")
)
rcs_tests <- do.call(rbind, lapply(rcs, `[[`, "test"))
rcs_curves <- do.call(rbind, lapply(rcs, `[[`, "curve"))

fit_ph_diagnostic <- function(design, exposure, analyte, weight_field, cohort) {
  rows <- droplevels(design$variables)
  split_data <- survival::survSplit(
    Surv(followup_years, death_all) ~ ., data = rows,
    cut = c(5, 10), start = "interval_start", end = "interval_end",
    event = "interval_event", episode = "followup_interval"
  )
  split_data$followup_interval <- factor(
    split_data$followup_interval, levels = 1:3,
    labels = c("0 to <5 years", "5 to <10 years", "10 or more years")
  )
  split_design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = stats::as.formula(paste0("~", weight_field)),
    nest = TRUE, data = split_data
  )
  ph_rhs <- paste0(
    cox_rhs(exposure, "M3"),
    " + survival::strata(followup_interval) + ", exposure, ":followup_interval"
  )
  ph_fit <- survey::svycoxph(
    stats::as.formula(paste0(
      "survival::Surv(interval_start, interval_end, interval_event) ~ ", ph_rhs
    )),
    design = split_design, model = TRUE, x = TRUE
  )
  interaction_names <- grep(
    paste0("^", exposure, ":followup_interval"),
    names(stats::coef(ph_fit)), value = TRUE
  )
  interaction_beta <- stats::coef(ph_fit)[interaction_names]
  interaction_covariance <- stats::vcov(ph_fit)[
    interaction_names, interaction_names, drop = FALSE
  ]
  ph_wald <- as.numeric(t(interaction_beta) %*%
                          solve(interaction_covariance, interaction_beta))
  test <- data.frame(
    cohort = cohort, analyte = analyte, model = "M3",
    df = length(interaction_names), wald_chisq = ph_wald,
    p_time_interaction = stats::pchisq(
      ph_wald, df = length(interaction_names), lower.tail = FALSE
    )
  )
  period_results <- lapply(levels(split_data$followup_interval), function(period) {
    weights <- stats::setNames(rep(0, length(stats::coef(ph_fit))),
                               names(stats::coef(ph_fit)))
    weights[[exposure]] <- 1
    if (period != levels(split_data$followup_interval)[[1]]) {
      interaction <- paste0(exposure, ":followup_interval", period)
      weights[[interaction]] <- 1
    }
    beta <- sum(weights * stats::coef(ph_fit))
    se <- sqrt(as.numeric(t(weights) %*% stats::vcov(ph_fit) %*% weights))
    data.frame(
      cohort = cohort, analyte = analyte, interval = period,
      events = sum(split_data$interval_event[
        split_data$followup_interval == period
      ]), beta = beta, robust_se = se, hr = exp(beta),
      ci_low = exp(beta - stats::qnorm(0.975) * se),
      ci_high = exp(beta + stats::qnorm(0.975) * se),
      p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
    )
  })
  list(test = test, intervals = do.call(rbind, period_results))
}

ph_results <- list(
  pfas_primary = fit_ph_diagnostic(
    pfas_design, "log2_sum4_pfas", "pfas", "weight_pfas_primary", "pfas_primary"
  ),
  pfas_overlap = fit_ph_diagnostic(
    overlap_design, "log2_sum4_pfas", "pfas", "weight_overlap", "overlap_primary"
  ),
  lead = fit_ph_diagnostic(
    overlap_design, "log2_lead", "lead", "weight_overlap", "overlap_primary"
  ),
  cadmium = fit_ph_diagnostic(
    overlap_design, "log2_cadmium", "cadmium", "weight_overlap", "overlap_primary"
  )
)
ph_tests <- do.call(rbind, lapply(ph_results, `[[`, "test"))
ph_intervals <- do.call(rbind, lapply(ph_results, `[[`, "intervals"))
primary_intervals <- ph_intervals[ph_intervals$cohort == "pfas_primary", ]

utils::write.csv(rcs_tests, file.path(paths$results, "rcs_nonlinearity_tests.csv"), row.names = FALSE)
utils::write.csv(rcs_curves, file.path(paths$results, "rcs_curves.csv"), row.names = FALSE)
utils::write.csv(ph_tests, file.path(paths$results, "ph_time_interaction_tests.csv"), row.names = FALSE)
utils::write.csv(ph_tests[ph_tests$cohort == "pfas_primary", ],
                 file.path(paths$results, "ph_time_interaction_test.csv"), row.names = FALSE)
utils::write.csv(primary_intervals, file.path(paths$results, "ph_interval_estimates.csv"), row.names = FALSE)
utils::write.csv(ph_intervals, file.path(paths$results, "ph_interval_estimates_all.csv"), row.names = FALSE)
message("Nonlinearity and proportional-hazards diagnostics complete.")
