#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "source_map.R"))
source(file.path("scripts", "utils", "model_spec.R"))
source(file.path("scripts", "utils", "data_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
source_file <- file.path(paths$processed, "harmonized_analysis_data.rds")
if (!file.exists(source_file)) stop("Run 03_harmonize_variables.R first")
data <- readRDS(source_file)

primary_cycles <- primary_pfas_cycles()
overlap_cycles <- primary_overlap_cycles()
base_fields <- base_model_fields()
physiology <- physiology_fields()
design <- design_fields()

data$weight_pfas_primary <- ifelse(
  data$cycle %in% primary_cycles, data$pfas_weight_2yr / length(primary_cycles), NA_real_
)
data$weight_overlap <- ifelse(
  data$cycle %in% overlap_cycles, data$pfas_weight_2yr / length(overlap_cycles), NA_real_
)
data$weight_metals_all <- data$metal_weight_2yr / length(unique(data$cycle))
data$weight_pfas_2013 <- ifelse(data$cycle == "2013-2014", data$pfas_weight_2yr, NA_real_)

valid_outcome <- valid_mortality_followup(data)
complete_base <- complete_fields(data, c(base_fields, physiology, design))
positive_pfas_primary <- !is.na(data$weight_pfas_primary) & data$weight_pfas_primary > 0
positive_overlap <- !is.na(data$weight_overlap) & data$weight_overlap > 0
positive_metals <- !is.na(data$weight_metals_all) & data$weight_metals_all > 0
positive_2013 <- !is.na(data$weight_pfas_2013) & data$weight_pfas_2013 > 0

data$domain_pfas_source <- data$cycle %in% primary_cycles & positive_pfas_primary
data$domain_overlap_source <- data$cycle %in% overlap_cycles & positive_overlap
data$domain_metals_source <- positive_metals
data$domain_pfas_2013_source <- data$cycle == "2013-2014" & positive_2013

data$cohort_pfas_primary <- data$domain_pfas_source & valid_outcome &
  !is.na(data$log2_sum4_pfas) & complete_base
data$cohort_overlap_primary <- data$domain_overlap_source & valid_outcome &
  complete_fields(data, c("log2_sum4_pfas", "log2_lead", "log2_cadmium")) &
  complete_base
data$cohort_lead_full <- data$domain_metals_source & valid_outcome &
  !is.na(data$log2_lead) & complete_base
data$cohort_cadmium_full <- data$domain_metals_source & valid_outcome &
  !is.na(data$log2_cadmium) & complete_base
data$cohort_pfas_2013 <- data$domain_pfas_2013_source & valid_outcome &
  !is.na(data$log2_sum4_pfas) & complete_base

data$restriction_no_baseline_disease <-
  data$diabetes_nonfasting == 0 & data$cvd_history == 0 & data$cancer_history == 0
data$restriction_not_pregnant <- is.na(data$pregnant) | data$pregnant == 0

flow_branch <- function(data, branch, cycle_filter, weight, exposures) {
  in_cycle <- data$cycle %in% cycle_filter
  steps <- list(
    "NHANES participants in included cycles" = in_cycle,
    "Adults aged 20 years or older" = in_cycle & !is.na(data$age_years) & data$age_years >= 20,
    "Eligible for mortality linkage with follow-up" = in_cycle & valid_outcome,
    "Positive analyte-specific survey weight" = in_cycle & valid_outcome &
      !is.na(weight) & weight > 0,
    "Complete exposure measurement" = in_cycle & valid_outcome &
      !is.na(weight) & weight > 0 & complete_fields(data, exposures),
    "Complete albumin and eGFR" = in_cycle & valid_outcome &
      !is.na(weight) & weight > 0 & complete_fields(data, c(exposures, physiology)),
    "Complete primary adjustment set" = in_cycle & valid_outcome &
      !is.na(weight) & weight > 0 &
      complete_fields(data, c(exposures, physiology, base_fields, design))
  )
  counts <- vapply(steps, sum, integer(1), na.rm = TRUE)
  data.frame(
    branch = branch, step_order = seq_along(steps), step = names(steps),
    n_remaining = unname(counts),
    n_excluded_since_previous = c(NA_integer_, head(counts, -1) - tail(counts, -1)),
    stringsAsFactors = FALSE
  )
}

flow <- rbind(
  flow_branch(data, "Primary PFAS cohort", primary_cycles,
              data$weight_pfas_primary, "log2_sum4_pfas"),
  flow_branch(data, "Primary PFAS-metals overlap cohort", overlap_cycles,
              data$weight_overlap, c("log2_sum4_pfas", "log2_lead", "log2_cadmium")),
  flow_branch(data, "Full metals cohort", unique(data$cycle),
              data$weight_metals_all, c("log2_lead", "log2_cadmium")),
  flow_branch(data, "2013-2014 surplus-specimen PFAS cohort", "2013-2014",
              data$weight_pfas_2013, "log2_sum4_pfas")
)

cohort_flags <- c(
  pfas_primary = "cohort_pfas_primary",
  overlap_primary = "cohort_overlap_primary",
  lead_full = "cohort_lead_full",
  cadmium_full = "cohort_cadmium_full",
  pfas_2013 = "cohort_pfas_2013"
)
cohort_summary <- do.call(rbind, lapply(names(cohort_flags), function(cohort) {
  summarize_cohort_followup(data, data[[cohort_flags[[cohort]]]], cohort)
}))
rownames(cohort_summary) <- NULL

if (!identical(sum(data$cohort_pfas_primary), tail(
  flow$n_remaining[flow$branch == "Primary PFAS cohort"], 1
))) stop("PFAS flow count does not match cohort flag")
if (!identical(sum(data$cohort_overlap_primary), tail(
  flow$n_remaining[flow$branch == "Primary PFAS-metals overlap cohort"], 1
))) stop("Overlap flow count does not match cohort flag")
congener_fields <- c("log2_pfoa", "log2_pfos", "log2_pfhxs", "log2_pfna")
if (any(!complete_fields(data[data$cohort_overlap_primary, ], congener_fields))) {
  stop("Fixed overlap cohort is incomplete for at least one PFAS component")
}

utils::write.csv(flow, file.path(paths$results, "participant_flow.csv"), row.names = FALSE)
utils::write.csv(cohort_summary, file.path(paths$results, "cohort_summary.csv"), row.names = FALSE)
saveRDS(data, file.path(paths$processed, "analysis_cohorts.rds"), compress = "xz")

message("Cohorts built: primary PFAS n=", sum(data$cohort_pfas_primary),
        "; overlap n=", sum(data$cohort_overlap_primary))
