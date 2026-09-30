make_domain_design <- function(data, source_flag, cohort_flag, weight_field) {
  source <- droplevels(data[data[[source_flag]], , drop = FALSE])
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = stats::as.formula(paste0("~", weight_field)),
    nest = TRUE, data = source
  )
  subset(design, design$variables[[cohort_flag]])
}

# Keep the survey package's internal variance call in one audited location.
# The locked package version and replicate-weight checks guard this interface.
survey_recvar <- function(influence, design) {
  namespace <- asNamespace("survey")
  if (!exists("svyrecvar", envir = namespace, inherits = FALSE)) {
    stop("The locked survey package no longer provides svyrecvar")
  }
  get("svyrecvar", envir = namespace)(
    influence, design$cluster, design$strata, design$fpc,
    postStrata = design$postStrata
  )
}

cox_rhs <- function(
    exposure, model, physiology = c("spline", "linear"), include_cycle = TRUE,
    bmi_shape = c("spline", "linear"),
    adjustment = c("main", "reduced_demographic")) {
  physiology <- match.arg(physiology)
  bmi_shape <- match.arg(bmi_shape)
  adjustment <- match.arg(adjustment)
  base <- c(exposure, "splines::ns(age_years, df = 4)", "sex", "race_ethnicity")
  if (adjustment == "main") {
    base <- c(
      base, "education", "pir", "smoking",
      if (bmi_shape == "spline") "splines::ns(bmi, df = 3)" else "bmi"
    )
  }
  if (include_cycle) base <- c(base, "cycle_factor")
  if (model %in% c("M1", "M3")) {
    base <- c(base, if (physiology == "spline") {
      "splines::ns(egfr, df = 3)"
    } else "egfr")
  }
  if (model %in% c("M2", "M3")) {
    base <- c(base, if (physiology == "spline") {
      "splines::ns(albumin_g_dl, df = 3)"
    } else "albumin_g_dl")
  }
  paste(base, collapse = " + ")
}

fit_selected_cox <- function(
    design, exposure, models, outcome = "death_all",
    physiology = c("spline", "linear"), bmi_shape = c("spline", "linear"),
    adjustment = c("main", "reduced_demographic")) {
  physiology <- match.arg(physiology)
  bmi_shape <- match.arg(bmi_shape)
  adjustment <- match.arg(adjustment)
  include_cycle <- length(unique(as.character(
    design$variables$cycle_factor
  ))) > 1L
  stats::setNames(lapply(models, function(model) {
    formulaString <- paste0(
      "survival::Surv(followup_years, ", outcome, ") ~ ",
      cox_rhs(
        exposure, model, physiology = physiology, include_cycle = include_cycle,
        bmi_shape = bmi_shape, adjustment = adjustment
      )
    )
    formula <- stats::as.formula(formulaString)
    survey::svycoxph(formula, design = design, model = TRUE, x = TRUE)
  }), models)
}

fit_symmetric_cox <- function(
    design, exposure, outcome = "death_all", physiology = c("spline", "linear"),
    bmi_shape = c("spline", "linear"),
    adjustment = c("main", "reduced_demographic")) {
  fit_selected_cox(
    design = design, exposure = exposure, models = c("M0", "M1", "M2", "M3"),
    outcome = outcome, physiology = physiology, bmi_shape = bmi_shape,
    adjustment = adjustment
  )
}

fit_target_survey_cox <- function(
    design, exposure, model, outcome = "death_all",
    physiology = c("spline", "linear"), bmi_shape = c("spline", "linear"),
    adjustment = c("main", "reduced_demographic")) {
  physiology <- match.arg(physiology)
  bmi_shape <- match.arg(bmi_shape)
  adjustment <- match.arg(adjustment)
  include_cycle <- length(unique(as.character(
    design$variables$cycle_factor
  ))) > 1L
  formula <- stats::as.formula(paste0(
    "survival::Surv(followup_years, ", outcome, ") ~ ",
    cox_rhs(
      exposure, model, physiology = physiology, include_cycle = include_cycle,
      bmi_shape = bmi_shape, adjustment = adjustment
    )
  ))
  model_data <- design$variables
  model_data$.survey_probability_weight <-
    (1 / design$prob) / mean(1 / design$prob)
  fit <- survival::coxph(
    formula, data = model_data, weights = .survey_probability_weight,
    subset = .survey_probability_weight > 0, model = TRUE, x = TRUE
  )
  influence <- as.matrix(stats::residuals(fit, type = "dfbeta", weighted = TRUE))
  if (is.null(colnames(influence))) colnames(influence) <- names(stats::coef(fit))
  target_influence <- influence[, exposure, drop = FALSE]
  target_variance <- survey_recvar(target_influence, design)
  list(
    beta = unname(stats::coef(fit)[[exposure]]),
    variance = unname(target_variance[1, 1]), influence = target_influence,
    design = design, n = fit$n, events = fit$nevent,
    degf_resid = survey::degf(design) - length(stats::coef(fit)) + 1
  )
}

joint_light_covariance <- function(fits, labels) {
  if (length(fits) != length(labels)) stop("fits and labels must have equal length")
  joint <- do.call(cbind, lapply(fits, `[[`, "influence"))
  colnames(joint) <- labels
  design <- fits[[1]]$design
  covariance <- survey_recvar(joint, design)
  dimnames(covariance) <- list(labels, labels)
  expected <- vapply(fits, `[[`, numeric(1), "variance")
  if (!isTRUE(all.equal(unname(diag(covariance)), unname(expected), tolerance = 1e-7))) {
    stop("Lightweight joint covariance does not reproduce target variances")
  }
  covariance
}

extract_cox_result <- function(fit, exposure, analyte, cohort, outcome, model) {
  beta <- unname(stats::coef(fit)[[exposure]])
  robust_se <- sqrt(stats::vcov(fit)[exposure, exposure])
  z <- beta / robust_se
  data.frame(
    cohort = cohort, analyte = analyte, outcome = outcome, model = model,
    n = fit$n, events = fit$nevent, beta = beta, robust_se = robust_se,
    hr = exp(beta), ci_low = exp(beta - stats::qnorm(0.975) * robust_se),
    ci_high = exp(beta + stats::qnorm(0.975) * robust_se),
    z = z, p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

target_dfbeta <- function(fit, exposure) {
  influence <- as.matrix(stats::residuals(fit, type = "dfbeta", weighted = TRUE))
  if (is.null(colnames(influence))) colnames(influence) <- names(stats::coef(fit))
  if (!(exposure %in% colnames(influence))) stop("Exposure absent from dfbeta: ", exposure)
  influence[, exposure, drop = FALSE]
}

joint_target_covariance <- function(fits, exposures, labels) {
  if (length(fits) != length(exposures) || length(fits) != length(labels)) {
    stop("fits, exposures, and labels must have equal length")
  }
  row_counts <- vapply(fits, function(fit) nrow(fit$survey.design$variables), integer(1))
  if (length(unique(row_counts)) != 1L) stop("Joint models do not share a sample")
  influences <- Map(target_dfbeta, fits, exposures)
  joint <- do.call(cbind, influences)
  colnames(joint) <- labels
  design <- fits[[1]]$survey.design
  covariance <- survey_recvar(joint, design)
  dimnames(covariance) <- list(labels, labels)
  expected <- vapply(seq_along(fits), function(i) {
    stats::vcov(fits[[i]])[exposures[[i]], exposures[[i]]]
  }, numeric(1))
  if (!isTRUE(all.equal(unname(diag(covariance)), unname(expected),
                        tolerance = 1e-7))) {
    stop("Joint influence covariance does not reproduce model variances")
  }
  covariance
}

linear_contrast <- function(estimates, covariance, weights, label) {
  weights <- weights[names(estimates)]
  estimate <- sum(weights * estimates)
  variance <- as.numeric(t(weights) %*% covariance %*% weights)
  se <- sqrt(max(variance, 0))
  z <- estimate / se
  data.frame(
    contrast = label, estimate = estimate, robust_se = se,
    ci_low = estimate - stats::qnorm(0.975) * se,
    ci_high = estimate + stats::qnorm(0.975) * se,
    z = z, p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

fit_audited_replicate_cox <- function(
    formula, data, target_term,
    allowed_separation_terms = "cycle_factor2017-2018") {
  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      formula, data = data, weights = .replicate_weight,
      subset = .replicate_weight > 0, model = FALSE, x = FALSE
    ),
    warning = function(condition) {
      fit_warnings <<- c(fit_warnings, trimws(conditionMessage(condition)))
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- stats::coef(fit)
  if (!(target_term %in% names(coefficients)) || !is.finite(coefficients[[target_term]])) {
    stop("Replicate Cox fit returned a non-finite target coefficient: ", target_term)
  }
  
  warning_audit <- data.frame(
    warning = character(), coefficient_index = integer(),
    affected_term = character(), target_estimate = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(fit_warnings)) {
    expected_pattern <- paste0(
      "^Loglik converged before variable\\s+[0-9]+\\s*; ",
      "coefficient may be infinite\\.$"
    )
    if (any(!grepl(expected_pattern, fit_warnings, perl = TRUE))) {
      stop("Unexpected warning from replicate Cox fit: ",
           paste(unique(fit_warnings), collapse = " | "))
    }
    coefficient_index <- as.integer(sub(
      ".*variable\\s+([0-9]+)\\s*;.*", "\\1", fit_warnings, perl = TRUE
    ))
    affected_term <- names(coefficients)[coefficient_index]
    if (anyNA(affected_term) || any(!(affected_term %in% allowed_separation_terms))) {
      stop(
        "Replicate Cox separation warning affected an unexpected term: ",
        paste(unique(affected_term), collapse = ", ")
      )
    }
    warning_audit <- data.frame(
      warning = fit_warnings,
      coefficient_index = coefficient_index,
      affected_term = affected_term,
      target_estimate = unname(coefficients[[target_term]]),
      stringsAsFactors = FALSE
    )
  }
  list(coefficients = coefficients, warnings = warning_audit)
}

coefficient_contrasts <- function(estimates, covariance, analytes) {
  output <- list()
  index <- 0L
  contrast_defs <- list(
    albumin_conditional_egfr = c(M1 = -1, M3 = 1),
    egfr_conditional_albumin = c(M2 = -1, M3 = 1),
    joint_albumin_egfr = c(M0 = -1, M3 = 1)
  )
  for (analyte in analytes) {
    for (contrast_name in names(contrast_defs)) {
      weights <- stats::setNames(rep(0, length(estimates)), names(estimates))
      model_weights <- contrast_defs[[contrast_name]]
      weights[paste(analyte, names(model_weights), sep = "_")] <- model_weights
      index <- index + 1L
      item <- linear_contrast(estimates, covariance, weights, contrast_name)
      item$analyte <- analyte
      output[[index]] <- item
    }
  }
  do.call(rbind, output)
}

fit_exposure_albumin_structure <- function(design, exposure, analyte) {
  reduced <- stats::as.formula(paste(
    exposure, "~ splines::ns(age_years, df = 4) + sex + race_ethnicity +",
    "education + pir + smoking + splines::ns(bmi, df = 3) + cycle_factor +",
    "splines::ns(egfr, df = 3)"
  ))
  full <- stats::update(reduced, . ~ . + albumin_g_dl)
  fit_reduced <- survey::svyglm(reduced, design = design, family = stats::gaussian())
  fit_full <- survey::svyglm(full, design = design, family = stats::gaussian())
  beta <- stats::coef(fit_full)[["albumin_g_dl"]]
  se <- sqrt(stats::vcov(fit_full)["albumin_g_dl", "albumin_g_dl"])
  incremental_r2 <- 1 - stats::deviance(fit_full) / stats::deviance(fit_reduced)
  data.frame(
    analyte = analyte, beta_albumin = beta, robust_se = se,
    ci_low = beta - stats::qnorm(0.975) * se,
    ci_high = beta + stats::qnorm(0.975) * se,
    p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
    deviance_incremental_r2 = incremental_r2,
    n = nrow(fit_full$model), stringsAsFactors = FALSE
  )
}
