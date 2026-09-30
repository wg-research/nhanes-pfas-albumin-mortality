#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))
source(file.path("scripts", "utils", "mi_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))
m <- 20L
random_seed <- 20260131L

missingness_table <- function(data, cohort) {
  data.frame(
    cohort = cohort, field = names(data), n = nrow(data),
    n_missing = vapply(data, function(x) sum(is.na(x)), integer(1)),
    percent_missing = 100 * vapply(data, function(x) mean(is.na(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
}

pfas_data <- prepare_imputation_data(data, "pfas")
overlap_data <- prepare_imputation_data(data, "overlap")

missingness <- rbind(
  missingness_table(pfas_data, "pfas_pre_imputation"),
  missingness_table(overlap_data, "overlap_pre_imputation")
)

load_or_impute <- function(cache_file, data, seed, label) {
  input_sha256 <- digest::digest(data, algo = "sha256")
  cache_signature <- digest::digest(
    list(
      input_sha256 = input_sha256,
      m = m,
      seed = seed,
      prepare_imputation_data = deparse(body(prepare_imputation_data)),
      run_mice = deparse(body(run_mice))
    ),
    algo = "sha256"
  )
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    same_signature <- identical(
      attr(cached, "cache_signature", exact = TRUE), cache_signature
    )
    if (same_signature) {
      message("Using cached ", label, " imputations (n=", nrow(data), ", m=", m, ")")
      return(cached)
    }
    message("Ignoring stale ", label, " imputation cache: signature mismatch")
  }
  message("Imputing ", label, " cohort (n=", nrow(data), ")")
  imputed <- run_mice(data, m = m, seed = seed)
  attr(imputed, "input_sha256") <- input_sha256
  attr(imputed, "cache_signature") <- cache_signature
  imputed
}

pfas_cache <- file.path(
  paths$processed, sprintf("mi_pfas_m20_seed_%d.rds", random_seed)
)
overlap_cache <- file.path(
  paths$processed, sprintf("mi_overlap_m20_seed_%d.rds", random_seed)
)
pfas_mids <- load_or_impute(pfas_cache, pfas_data, random_seed, "PFAS")
overlap_mids <- load_or_impute(overlap_cache, overlap_data, random_seed, "overlap")
saveRDS(pfas_mids, pfas_cache, compress = "xz")
saveRDS(overlap_mids, overlap_cache, compress = "xz")

mi_pfas_models <- c("M0", "M1", "M3")
pfas_q <- stats::setNames(lapply(mi_pfas_models, function(x) numeric(m)), mi_pfas_models)
pfas_u <- pfas_q
delta_q <- numeric(m)
delta_u <- numeric(m)
df_pfas <- numeric(m)

for (i in seq_len(m)) {
  design <- completed_survey_design(pfas_mids, i, "weight_pfas_primary")
  fits <- stats::setNames(lapply(mi_pfas_models, function(model) {
    fit_target_survey_cox(design, "log2_sum4_pfas", model)
  }), mi_pfas_models)
  for (model in names(fits)) {
    pfas_q[[model]][[i]] <- fits[[model]]$beta
    pfas_u[[model]][[i]] <- fits[[model]]$variance
  }
  covariance <- joint_light_covariance(fits[c("M1", "M3")], c("M1", "M3"))
  weights <- c(M1 = -1, M3 = 1)
  delta_q[[i]] <- sum(weights * c(
    M1 = pfas_q$M1[[i]], M3 = pfas_q$M3[[i]]
  ))
  delta_u[[i]] <- as.numeric(t(weights) %*% covariance %*% weights)
  df_pfas[[i]] <- fits$M3$degf_resid
}

pfas_pooled <- do.call(rbind, lapply(names(pfas_q), function(model) {
  output <- pool_scalar_survey(
    pfas_q[[model]], pfas_u[[model]], min(df_pfas), paste0("pfas_", model)
  )
  output$model <- model
  output$hr <- exp(output$estimate)
  output$hr_ci_low <- exp(output$ci_low)
  output$hr_ci_high <- exp(output$ci_high)
  output
}))
pfas_delta_pooled <- pool_scalar_survey(
  delta_q, delta_u, min(df_pfas), "albumin_conditional_egfr"
)

analytes <- c(pfas = "log2_sum4_pfas", lead = "log2_lead", cadmium = "log2_cadmium")
mi_overlap_models <- c("M1", "M3")
labels <- as.vector(outer(names(analytes), mi_overlap_models, paste, sep = "_"))
overlap_q <- stats::setNames(lapply(labels, function(x) numeric(m)), labels)
overlap_u <- overlap_q
contrast_names <- c(
  "delta_pfas", "delta_lead", "delta_cadmium", "pfas_minus_lead", "pfas_minus_cadmium"
)
contrast_q <- stats::setNames(lapply(contrast_names, function(x) numeric(m)), contrast_names)
contrast_u <- contrast_q
df_overlap <- numeric(m)

for (i in seq_len(m)) {
  design <- completed_survey_design(overlap_mids, i, "weight_overlap")
  fit_vector <- list()
  exposure_vector <- character()
  fit_labels <- character()
  for (analyte in names(analytes)) {
    fits <- stats::setNames(lapply(mi_overlap_models, function(model) {
      fit_target_survey_cox(design, analytes[[analyte]], model)
    }), mi_overlap_models)
    for (model in names(fits)) {
      label <- paste(analyte, model, sep = "_")
      fit_vector[[length(fit_vector) + 1L]] <- fits[[model]]
      exposure_vector <- c(exposure_vector, analytes[[analyte]])
      fit_labels <- c(fit_labels, label)
      overlap_q[[label]][[i]] <- fits[[model]]$beta
      overlap_u[[label]][[i]] <- fits[[model]]$variance
    }
  }
  covariance <- joint_light_covariance(fit_vector, fit_labels)
  estimates <- stats::setNames(vapply(fit_vector, `[[`, numeric(1), "beta"), fit_labels)

  contrast_weights <- list()
  for (analyte in names(analytes)) {
    w <- stats::setNames(rep(0, length(estimates)), names(estimates))
    w[paste0(analyte, c("_M1", "_M3"))] <- c(-1, 1)
    contrast_weights[[paste0("delta_", analyte)]] <- w
  }
  contrast_weights$pfas_minus_lead <- contrast_weights$delta_pfas - contrast_weights$delta_lead
  contrast_weights$pfas_minus_cadmium <- contrast_weights$delta_pfas - contrast_weights$delta_cadmium
  for (name in names(contrast_weights)) {
    w <- contrast_weights[[name]]
    contrast_q[[name]][[i]] <- sum(w * estimates)
    contrast_u[[name]][[i]] <- as.numeric(t(w) %*% covariance %*% w)
  }
  df_overlap[[i]] <- min(vapply(fit_vector, `[[`, numeric(1), "degf_resid"))
}

overlap_pooled <- do.call(rbind, lapply(labels, function(label) {
  output <- pool_scalar_survey(
    overlap_q[[label]], overlap_u[[label]], min(df_overlap), label
  )
  parts <- strsplit(label, "_", fixed = TRUE)[[1]]
  output$analyte <- parts[[1]]
  output$model <- parts[[2]]
  output$hr <- exp(output$estimate)
  output$hr_ci_low <- exp(output$ci_low)
  output$hr_ci_high <- exp(output$ci_high)
  output
}))
contrast_pooled <- do.call(rbind, lapply(contrast_names, function(name) {
  pool_scalar_survey(
    contrast_q[[name]], contrast_u[[name]], min(df_overlap), name
  )
}))
between_rows <- contrast_pooled$term %in% c("pfas_minus_lead", "pfas_minus_cadmium")
contrast_pooled$p_holm <- NA_real_
contrast_pooled$p_holm[between_rows] <- stats::p.adjust(
  contrast_pooled$p_value[between_rows], method = "holm"
)

logged <- rbind(
  if (!is.null(pfas_mids$loggedEvents)) transform(pfas_mids$loggedEvents, cohort = "pfas") else NULL,
  if (!is.null(overlap_mids$loggedEvents)) transform(overlap_mids$loggedEvents, cohort = "overlap") else NULL
)
convergence <- rbind(
  transform(mids_convergence(pfas_mids), cohort = "pfas"),
  transform(mids_convergence(overlap_mids), cohort = "overlap")
)

utils::write.csv(missingness, file.path(paths$results, "mi_missingness.csv"), row.names = FALSE)
utils::write.csv(pfas_pooled, file.path(paths$results, "mi_pfas_models.csv"), row.names = FALSE)
utils::write.csv(pfas_delta_pooled, file.path(paths$results, "mi_pfas_contrast.csv"), row.names = FALSE)
utils::write.csv(overlap_pooled, file.path(paths$results, "mi_overlap_models.csv"), row.names = FALSE)
utils::write.csv(contrast_pooled, file.path(paths$results, "mi_overlap_contrasts.csv"), row.names = FALSE)
utils::write.csv(convergence, file.path(paths$audit, "mi_convergence.csv"), row.names = FALSE)
mi_run_metadata <- data.frame(
  cohort = c("pfas", "overlap"),
  random_seed = rep(random_seed, 2),
  imputations = rep(m, 2),
  iterations = rep(10L, 2),
  cache_file = basename(c(pfas_cache, overlap_cache)),
  input_sha256 = c(
    attr(pfas_mids, "input_sha256", exact = TRUE),
    attr(overlap_mids, "input_sha256", exact = TRUE)
  ),
  cache_signature = c(
    attr(pfas_mids, "cache_signature", exact = TRUE),
    attr(overlap_mids, "cache_signature", exact = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  mi_run_metadata, file.path(paths$audit, "mi_run_metadata.csv"), row.names = FALSE
)
if (!is.null(logged)) {
  utils::write.csv(logged, file.path(paths$audit, "mi_logged_events.csv"), row.names = FALSE)
}
message("Multiple-imputation analysis complete (seed=", random_seed, ").")
