#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))
analytes <- c(
  pfas = "log2_sum4_pfas", lead = "log2_lead", cadmium = "log2_cadmium"
)

fit_overlap_system <- function(
    data, selected = names(analytes), adjustment = "main") {
  design <- make_domain_design(
    data, "domain_overlap_source", "cohort_overlap_primary", "weight_overlap"
  )
  fields <- analytes[selected]
  fits <- lapply(fields, function(exposure) {
    fit_selected_cox(
      design, exposure, c("M1", "M3"), adjustment = adjustment
    )
  })
  fit_vector <- list()
  exposure_vector <- labels <- character()
  for (analyte in names(fits)) {
    for (model in c("M1", "M3")) {
      fit_vector[[length(fit_vector) + 1L]] <- fits[[analyte]][[model]]
      exposure_vector <- c(exposure_vector, fields[[analyte]])
      labels <- c(labels, paste(analyte, model, sep = "_"))
    }
  }
  covariance <- joint_target_covariance(fit_vector, exposure_vector, labels)
  estimates <- stats::setNames(vapply(seq_along(fit_vector), function(i) {
    unname(stats::coef(fit_vector[[i]])[[exposure_vector[[i]]]])
  }, numeric(1)), labels)

  list(
    design = design, fits = fits, estimates = estimates,
    covariance = covariance, selected = names(fits)
  )
}

albumin_contrasts <- function(system, scenario, scale = NULL) {
  estimates <- system$estimates
  covariance <- system$covariance
  selected <- system$selected
  if (is.null(scale)) scale <- stats::setNames(rep(1, length(selected)), selected)
  if (!setequal(names(scale), selected) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Contrast scales must be finite, positive, and named for each analyte")
  }

  output <- list()
  for (analyte in selected) {
    weights <- stats::setNames(rep(0, length(estimates)), names(estimates))
    weights[paste0(analyte, c("_M1", "_M3"))] <-
      scale[[analyte]] * c(-1, 1)
    item <- linear_contrast(
      estimates, covariance, weights, paste0("delta_", analyte)
    )
    item$scenario <- scenario
    item$n <- nrow(system$design$variables)
    item$events <- sum(system$design$variables$death_all)
    item$scale <- scale[[analyte]]
    output[[length(output) + 1L]] <- item
  }
  if ("pfas" %in% selected) {
    for (metal in intersect(c("lead", "cadmium"), selected)) {
      weights <- stats::setNames(rep(0, length(estimates)), names(estimates))
      weights[paste0("pfas", c("_M1", "_M3"))] <-
        scale[["pfas"]] * c(-1, 1)
      weights[paste0(metal, c("_M1", "_M3"))] <-
        scale[[metal]] * c(1, -1)
      item <- linear_contrast(
        estimates, covariance, weights, paste0("pfas_minus_", metal)
      )
      item$scenario <- scenario
      item$n <- nrow(system$design$variables)
      item$events <- sum(system$design$variables$death_all)
      item$scale <- NA_real_
      output[[length(output) + 1L]] <- item
    }
  }
  do.call(rbind, output)
}

fit_overlap_scenario <- function(
    data, scenario, selected = names(analytes), adjustment = "main") {
  system <- fit_overlap_system(data, selected, adjustment)
  albumin_contrasts(system, scenario)
}

# The documented common top code is the primary analysis.  Recreate the former
# raw-age/eGFR specification as a sensitivity analysis on exactly the same
# participant flags.
raw_age <- data
raw_age$age_years <- raw_age$age_years_raw
raw_age$egfr <- raw_age$egfr_raw_age
raw_overlap <- fit_overlap_scenario(raw_age, "raw_age_and_raw_age_egfr")
raw_primary_design <- make_domain_design(
  raw_age, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
raw_primary_fits <- fit_selected_cox(
  raw_primary_design, "log2_sum4_pfas", c("M0", "M3")
)
age_primary <- do.call(rbind, lapply(names(raw_primary_fits), function(model) {
  result <- extract_cox_result(
    raw_primary_fits[[model]], "log2_sum4_pfas", "pfas", "pfas_primary",
    "all_cause", model
  )
  result$scenario <- "raw_age_and_raw_age_egfr"
  result
}))

# Public NHANES values below the LLOD equal LLOD/sqrt(2).  Replacing only
# flagged values with LLOD/2 lowers those observations by 0.5 on the log2 scale.
llod_half <- data
cadmium_flagged <- !is.na(llod_half$LBDBCDLC) & llod_half$LBDBCDLC == 1 &
  !is.na(llod_half$log2_cadmium)
llod_half$log2_cadmium[cadmium_flagged] <-
  llod_half$log2_cadmium[cadmium_flagged] - 0.5
llod_half_result <- fit_overlap_scenario(
  llod_half, "cadmium_llod_over_2", c("pfas", "cadmium")
)

# The highest applicable cadmium LLOD proportion occurred in 2009-2010.
# Excluding that cycle leaves a five-cycle national estimate, so the remaining
# two-year weights are rescaled from 1/6 to 1/5.
exclude_high_llod <- data
keep_cycle <- exclude_high_llod$cycle != "2009-2010"
exclude_high_llod$domain_overlap_source <-
  exclude_high_llod$domain_overlap_source & keep_cycle
exclude_high_llod$cohort_overlap_primary <-
  exclude_high_llod$cohort_overlap_primary & keep_cycle
exclude_high_llod$weight_overlap[keep_cycle] <-
  exclude_high_llod$weight_overlap[keep_cycle] * 6 / 5
exclude_high_llod$weight_overlap[!keep_cycle] <- NA_real_
exclude_cycle_result <- fit_overlap_scenario(
  exclude_high_llod, "exclude_2009_2010", c("pfas", "cadmium")
)

# Structural sensitivity analyses retain the primary same-participant cohort
# unless the scenario explicitly restricts survey cycles.
primary_system <- fit_overlap_system(data)
weighted_quartiles <- do.call(rbind, lapply(names(analytes), function(analyte) {
  exposure <- analytes[[analyte]]
  quartiles <- survey::svyquantile(
    stats::as.formula(paste0("~", exposure)), primary_system$design,
    quantiles = c(0.25, 0.75), ci = FALSE, na.rm = TRUE
  )
  values <- unname(stats::coef(quartiles))
  if (length(values) != 2L) stop("Weighted quartile extraction failed: ", exposure)
  data.frame(
    analyte = analyte, q25_log2 = values[[1]], q75_log2 = values[[2]],
    log2_weighted_iqr = values[[2]] - values[[1]],
    stringsAsFactors = FALSE
  )
}))
weighted_iqr <- stats::setNames(
  weighted_quartiles$log2_weighted_iqr, weighted_quartiles$analyte
)

structural_primary <- albumin_contrasts(
  primary_system, "primary_per_doubling"
)
structural_iqr <- albumin_contrasts(
  primary_system, "weighted_iqr_log2", weighted_iqr
)
reduced_system <- fit_overlap_system(data, adjustment = "reduced_demographic")
if (nrow(reduced_system$design$variables) != nrow(primary_system$design$variables) ||
    sum(reduced_system$design$variables$death_all) !=
      sum(primary_system$design$variables$death_all)) {
  stop("Reduced-base analysis did not retain the primary overlap cohort")
}
structural_reduced <- albumin_contrasts(
  reduced_system, "reduced_demographic_base"
)

pre_2015 <- data
pre_2015_cycles <- c(
  "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012"
)
keep_pre_2015 <- pre_2015$cycle %in% pre_2015_cycles
pre_2015$domain_overlap_source <-
  pre_2015$domain_overlap_source & keep_pre_2015
pre_2015$cohort_overlap_primary <-
  pre_2015$cohort_overlap_primary & keep_pre_2015
pre_2015$weight_overlap[keep_pre_2015] <-
  pre_2015$weight_overlap[keep_pre_2015] * 6 / 5
pre_2015$weight_overlap[!keep_pre_2015] <- NA_real_
pre_2015_system <- fit_overlap_system(pre_2015)
observed_pre_2015_cycles <- sort(unique(as.character(
  pre_2015_system$design$variables$cycle
)))
if (!identical(observed_pre_2015_cycles, sort(pre_2015_cycles))) {
  stop("2003–2012 overlap design contains unexpected cycles")
}
expected_pre_2015 <- data$cohort_overlap_primary & data$cycle %in% pre_2015_cycles
if (nrow(pre_2015_system$design$variables) != sum(expected_pre_2015) ||
    sum(pre_2015_system$design$variables$death_all) !=
      sum(data$death_all[expected_pre_2015])) {
  stop("2003–2012 model sample does not match the scenario cohort")
}
structural_pre_2015 <- albumin_contrasts(
  pre_2015_system, "pre_2015_overlap_cycles"
)

structural <- rbind(
  structural_primary, structural_iqr, structural_reduced,
  structural_pre_2015
)
structural$scaling <- ifelse(
  structural$scenario == "weighted_iqr_log2", "weighted_iqr_log2",
  "per_doubling"
)
structural$p_holm_within_scenario <- NA_real_
for (scenario in unique(structural$scenario)) {
  rows <- structural$scenario == scenario & grepl("^pfas_minus_", structural$contrast)
  structural$p_holm_within_scenario[rows] <- stats::p.adjust(
    structural$p_value[rows], method = "holm"
  )
}

iqr_output <- data.frame(
  weighted_quartiles,
  cohort = "overlap_primary", n = nrow(primary_system$design$variables),
  stringsAsFactors = FALSE
)

scenario_audit <- data.frame(
  scenario = c(
    "primary_per_doubling", "weighted_iqr_log2",
    "reduced_demographic_base", "pre_2015_overlap_cycles"
  ),
  n = c(
    nrow(primary_system$design$variables), nrow(primary_system$design$variables),
    nrow(reduced_system$design$variables), nrow(pre_2015_system$design$variables)
  ),
  events = c(
    sum(primary_system$design$variables$death_all),
    sum(primary_system$design$variables$death_all),
    sum(reduced_system$design$variables$death_all),
    sum(pre_2015_system$design$variables$death_all)
  ),
  cycles = c(
    paste(sort(unique(as.character(primary_system$design$variables$cycle))), collapse = ";"),
    paste(sort(unique(as.character(primary_system$design$variables$cycle))), collapse = ";"),
    paste(sort(unique(as.character(reduced_system$design$variables$cycle))), collapse = ";"),
    paste(observed_pre_2015_cycles, collapse = ";")
  ),
  weight_rescale = c(1, 1, 1, 6 / 5),
  stringsAsFactors = FALSE
)

targeted <- rbind(raw_overlap, llod_half_result, exclude_cycle_result)
targeted$p_holm_within_scenario <- NA_real_
for (scenario in unique(targeted$scenario)) {
  rows <- targeted$scenario == scenario & grepl("^pfas_minus_", targeted$contrast)
  if (any(rows)) {
    targeted$p_holm_within_scenario[rows] <- stats::p.adjust(
      targeted$p_value[rows], method = "holm"
    )
  }
}

utils::write.csv(
  age_primary, file.path(paths$results, "age_definition_primary_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  targeted, file.path(paths$results, "targeted_contrast_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  structural, file.path(paths$results, "structural_contrast_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  iqr_output, file.path(paths$results, "weighted_iqr_log2.csv"),
  row.names = FALSE
)
utils::write.csv(
  scenario_audit,
  file.path(paths$audit, "structural_sensitivity_scenario_audit.csv"),
  row.names = FALSE
)
message("Targeted and structural contrast sensitivity analyses complete")
