prepare_imputation_data <- function(data, cohort = c("pfas", "overlap")) {
  cohort <- match.arg(cohort)
  valid_outcome <- data$age_years >= 20 & data$mortality_eligible &
    !is.na(data$followup_years) & data$followup_years > 0 & !is.na(data$death_all)
  if (cohort == "pfas") {
    eligible <- data$domain_pfas_source & valid_outcome &
      stats::complete.cases(data[, c(
        "log2_sum4_pfas", "albumin_g_dl", "egfr", "SDMVSTRA", "SDMVPSU",
        "weight_pfas_primary"
      )])
    weight <- "weight_pfas_primary"
    exposures <- "log2_sum4_pfas"
  } else {
    eligible <- data$domain_overlap_source & valid_outcome &
      stats::complete.cases(data[, c(
        "log2_sum4_pfas", "log2_lead", "log2_cadmium", "albumin_g_dl", "egfr",
        "SDMVSTRA", "SDMVPSU", "weight_overlap"
      )])
    weight <- "weight_overlap"
    exposures <- c("log2_sum4_pfas", "log2_lead", "log2_cadmium")
  }
  fields <- c(
    "SEQN", "SDMVSTRA", "SDMVPSU", weight, "followup_years", "death_all",
    exposures, "age_years", "sex", "race_ethnicity", "education", "pir",
    "smoking", "bmi", "cycle_factor", "albumin_g_dl", "egfr"
  )
  output <- droplevels(data[eligible, fields, drop = FALSE])
  output$nelson_aalen <- mice::nelsonaalen(output, followup_years, death_all)
  output$log_survey_weight <- log(output[[weight]])
  attr(output, "weight_field") <- weight
  output
}

run_mice <- function(data, m = 20L, maxit = 10L, seed = 20260131L) {
  initial <- mice::mice(data, maxit = 0, printFlag = FALSE)
  method <- initial$method
  method[] <- ""
  method[c("education", "smoking")] <- "polyreg"
  method[c("pir", "bmi")] <- "pmm"

  predictor <- initial$predictorMatrix
  predictor[,] <- 0
  targets <- names(method)[method != ""]
  auxiliaries <- setdiff(
    names(data), c("SEQN", "followup_years")
  )
  predictor[targets, auxiliaries] <- 1
  diag(predictor) <- 0

  mice::mice(
    data, m = m, maxit = maxit, method = method,
    predictorMatrix = predictor, seed = seed, printFlag = FALSE
  )
}

split_rhat <- function(draws) {
  # Auxiliary rank-normalized split/folded diagnostic applied to MICE chain
  # means; this is not a MICE-native R-hat statistic.
  draws <- as.matrix(draws)
  if (nrow(draws) < 4L || ncol(draws) < 2L || anyNA(draws)) return(NA_real_)
  half <- floor(nrow(draws) / 2L)
  split_draws <- cbind(
    draws[seq_len(half), , drop = FALSE],
    draws[seq.int(nrow(draws) - half + 1L, nrow(draws)), , drop = FALSE]
  )
  basic_rhat <- function(values) {
    within <- mean(apply(values, 2L, stats::var))
    between <- nrow(values) * stats::var(colMeans(values))
    if (!is.finite(within) || within <= 0 || !is.finite(between)) return(NA_real_)
    sqrt(
      (((nrow(values) - 1) / nrow(values)) * within +
         between / nrow(values)) / within
    )
  }
  rank_normalize <- function(values) {
    ranks <- rank(as.vector(values), ties.method = "average")
    matrix(
      stats::qnorm((ranks - 0.5) / length(ranks)),
      nrow = nrow(values), ncol = ncol(values)
    )
  }
  bulk_rhat <- basic_rhat(rank_normalize(split_draws))
  folded_rhat <- basic_rhat(rank_normalize(
    abs(split_draws - stats::median(split_draws))
  ))
  diagnostics <- c(bulk_rhat, folded_rhat)
  diagnostics <- diagnostics[is.finite(diagnostics)]
  if (!length(diagnostics)) NA_real_ else max(diagnostics)
}

mids_convergence <- function(data) {
  if (!inherits(data, "mids")) stop("data must be a mice mids object")
  if (data$m < 2L || data$iteration < 3L || is.null(data$chainMean)) {
    stop("mids object does not contain sufficient chains for convergence diagnostics")
  }
  variables <- names(data$data)
  iterations <- seq_len(data$iteration)
  output <- do.call(rbind, lapply(variables, function(variable) {
    chain_mean <- data$chainMean[variable, , , drop = FALSE][1, , ]
    adjacent <- vapply(iterations[-1L], function(iteration) {
      value <- suppressWarnings(stats::cor(
        chain_mean[iteration - 1L, ], chain_mean[iteration, ],
        use = "pairwise.complete.obs"
      ))
      if (is.na(value)) 0 else value
    }, numeric(1))
    autocorrelation <- c(NA_real_, cumsum(adjacent) / seq_along(adjacent)) +
      0 * chain_mean[, 1L]
    data.frame(
      .it = iterations,
      vrb = variable,
      ac = autocorrelation,
      psrf = vapply(iterations, function(iteration) {
        split_rhat(chain_mean[seq_len(iteration), , drop = FALSE])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output
}

completed_survey_design <- function(imputation, index, weight_field) {
  completed <- droplevels(mice::complete(imputation, index))
  survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = stats::as.formula(paste0("~", weight_field)),
    nest = TRUE, data = completed
  )
}

pool_scalar_survey <- function(q, u, df_complete, label) {
  pooled <- mice::pool.scalar(q, u, n = df_complete + 1, k = 1)
  se <- sqrt(pooled$t)
  critical <- stats::qt(0.975, df = pooled$df)
  data.frame(
    term = label, estimate = pooled$qbar, robust_se = se,
    ci_low = pooled$qbar - critical * se,
    ci_high = pooled$qbar + critical * se,
    df = pooled$df,
    p_value = 2 * stats::pt(abs(pooled$qbar / se), df = pooled$df,
                            lower.tail = FALSE),
    fraction_missing_information = pooled$fmi,
    between_imputation_variance = pooled$b,
    stringsAsFactors = FALSE
  )
}
