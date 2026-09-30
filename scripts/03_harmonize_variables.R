#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "source_map.R"))
source(file.path("scripts", "utils", "data_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
source_file <- file.path(paths$processed, "extracted_source_data.rds")
if (!file.exists(source_file)) stop("Run 02_extract_source_data.R first")
data <- readRDS(source_file)

for (field in names(data)) {
  if (inherits(data[[field]], "haven_labelled")) data[[field]] <- as.numeric(data[[field]])
}

# NHANES public result fields already contain LLOD/sqrt(2) below detection.
# Preserve those values and construct cycle-specific PFOA and PFOS totals.
early <- data$cycle %in% c(
  "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012"
)
isomer <- data$cycle %in% c("2015-2016", "2017-2018")
surplus <- data$cycle == "2013-2014"

data$pfoa_ng_ml <- NA_real_
data$pfos_ng_ml <- NA_real_
data$pfoa_ng_ml[early] <- data$LBXPFOA[early]
data$pfos_ng_ml[early] <- data$LBXPFOS[early]
data$pfoa_ng_ml[isomer] <- complete_sum(data$LBXNFOA, data$LBXBFOA)[isomer]
data$pfos_ng_ml[isomer] <- complete_sum(data$LBXNFOS, data$LBXMFOS)[isomer]
data$pfoa_ng_ml[surplus] <- complete_sum(data$SSNPFOA, data$SSBPFOA)[surplus]
data$pfos_ng_ml[surplus] <- complete_sum(data$SSNPFOS, data$SSMPFOS)[surplus]
data$pfhxs_ng_ml <- data$LBXPFHS
data$pfna_ng_ml <- data$LBXPFNA

molecular_weight <- c(pfoa = 414.07, pfos = 500.13, pfhxs = 400.12, pfna = 464.08)
for (compound in names(molecular_weight)) {
  concentration <- paste0(compound, "_ng_ml")
  molar <- paste0(compound, "_nmol_l")
  log2_molar <- paste0("log2_", compound)
  data[[molar]] <- data[[concentration]] * 1000 / molecular_weight[[compound]]
  data[[log2_molar]] <- ifelse(data[[molar]] > 0, log2(data[[molar]]), NA_real_)
}
data$sum4_pfas_nmol_l <- rowSums(data[, paste0(
  names(molecular_weight), "_nmol_l"
)])
sum4_fields <- c("pfoa_ng_ml", "pfos_ng_ml", "pfhxs_ng_ml", "pfna_ng_ml")
data$sum4_pfas_nmol_l[rowSums(is.na(data[, sum4_fields])) > 0] <- NA_real_
data$log2_sum4_pfas <- ifelse(
  data$sum4_pfas_nmol_l > 0, log2(data$sum4_pfas_nmol_l), NA_real_
)
data$log2_lead <- ifelse(data$LBXBPB > 0, log2(data$LBXBPB), NA_real_)
data$log2_cadmium <- ifelse(data$LBXBCD > 0, log2(data$LBXBCD), NA_real_)

data$creatinine_uncalibrated <- data$LBXSCR
data$creatinine_mg_dl <- data$LBXSCR
calibration_cycle <- data$cycle == "2005-2006" & !is.na(data$LBXSCR)
data$creatinine_mg_dl[calibration_cycle] <-
  -0.016 + 0.978 * data$LBXSCR[calibration_cycle]

# NHANES used an 85-year top code in 2003-2006 and an 80-year top code
# thereafter.  The pooled analysis uses the documented common 80-year code;
# retain the raw-age version for the corresponding sensitivity audit.
data$age_years_raw <- data$RIDAGEYR
data$age_years <- pmin(data$RIDAGEYR, 80)
data$egfr <- ckd_epi_2021(
  data$creatinine_mg_dl, data$age_years, female = data$RIAGENDR == 2
)
data$egfr_raw_age <- ckd_epi_2021(
  data$creatinine_mg_dl, data$age_years_raw, female = data$RIAGENDR == 2
)

data$sex <- factor(
  ifelse(data$RIAGENDR == 1, "Male", ifelse(data$RIAGENDR == 2, "Female", NA)),
  levels = c("Male", "Female")
)
race <- c(
  `1` = "Mexican American", `2` = "Other Hispanic",
  `3` = "Non-Hispanic White", `4` = "Non-Hispanic Black",
  `5` = "Other or multiracial"
)
data$race_ethnicity <- factor(
  unname(race[as.character(data$RIDRETH1)]), levels = unname(race)
)
education <- rep(NA_character_, nrow(data))
education[data$DMDEDUC2 %in% c(1, 2)] <- "Less than high school"
education[data$DMDEDUC2 == 3] <- "High school or GED"
education[data$DMDEDUC2 == 4] <- "Some college or associate degree"
education[data$DMDEDUC2 == 5] <- "College graduate or above"
data$education <- factor(education, levels = c(
  "Less than high school", "High school or GED",
  "Some college or associate degree", "College graduate or above"
))
data$smoking <- derive_smoking(data$SMQ020, data$SMQ040)
data$diabetes_nonfasting <- derive_diabetes_nonfasting(data$DIQ010, data$LBXGH)
data$cvd_history <- derive_cvd_history(data)
data$cancer_history <- derive_binary_history(data$MCQ220)
data$pregnant <- ifelse(
  data$RIDEXPRG == 1, 1L,
  ifelse(data$RIAGENDR == 1 | data$RIDEXPRG %in% c(2, 3), 0L, NA_integer_)
)

data$mortality_eligible <- data$ELIGSTAT == 1
data$followup_years <- data$PERMTH_EXM / 12
# Zero denotes death in the examination month; half a month retains the event.
data$followup_years[data$mortality_eligible & data$PERMTH_EXM == 0] <- 0.5 / 12
data$death_all <- ifelse(data$mortality_eligible, data$MORTSTAT, NA_integer_)
data$death_cvd <- ifelse(
  data$mortality_eligible,
  as.integer(data$MORTSTAT == 1 & data$UCOD_LEADING %in% c(1, 5)), NA_integer_
)
data$death_cancer <- ifelse(
  data$mortality_eligible,
  as.integer(data$MORTSTAT == 1 & data$UCOD_LEADING == 2), NA_integer_
)

# Public-use cause categories are reported in full through 2013-2014 but are
# collapsed to heart disease, cancer, and all other causes in 2015-2018.
cause_labels <- c(
  "Heart disease", "Malignant neoplasms", "Chronic lower respiratory diseases",
  "Accidents", "Cerebrovascular diseases", "Alzheimer disease", "Diabetes mellitus",
  "Influenza and pneumonia", "Nephritis and related kidney diseases", "All other causes"
)
cause_code_availability_audit <- expand.grid(
  cycle = nhanes_source_map()$cycle,
  ucod_leading = seq_along(cause_labels),
  stringsAsFactors = FALSE
)
cause_code_availability_audit$cause_group <- cause_labels[
  cause_code_availability_audit$ucod_leading
]
cause_code_availability_audit$deaths <- vapply(
  seq_len(nrow(cause_code_availability_audit)),
  function(i) sum(
    data$MORTSTAT == 1 &
      data$cycle == cause_code_availability_audit$cycle[[i]] &
      data$UCOD_LEADING == cause_code_availability_audit$ucod_leading[[i]],
    na.rm = TRUE
  ),
  integer(1)
)
late_cycles <- c("2015-2016", "2017-2018")
late_death_codes <- unique(data$UCOD_LEADING[
  data$MORTSTAT == 1 & data$cycle %in% late_cycles & !is.na(data$UCOD_LEADING)
])
if (!setequal(late_death_codes, c(1L, 2L, 10L))) {
  stop("Unexpected public-use cause categories in 2015-2018 mortality data")
}

data$pfas_weight_2yr <- NA_real_
for (field in unique(stats::na.omit(data$pfas_weight_source))) {
  rows <- data$pfas_weight_source == field
  data$pfas_weight_2yr[rows] <- data[[field]][rows]
}
data$metal_weight_2yr <- ifelse(
  data$metal_weight_source == "WTSH2YR", data$WTSH2YR, data$WTMEC2YR
)

data$bmi <- data$BMXBMI
data$pir <- data$INDFMPIR
data$albumin_g_dl <- data$LBXSAL
data$cycle_factor <- factor(data$cycle, levels = nhanes_source_map()$cycle)

audit_missingness <- function(x, fields) {
  do.call(rbind, lapply(fields, function(field) data.frame(
    field = field,
    n = nrow(x),
    n_observed = sum(!is.na(x[[field]])),
    n_missing = sum(is.na(x[[field]])),
    percent_missing = 100 * mean(is.na(x[[field]])),
    minimum = if (all(is.na(x[[field]]))) NA_real_ else min(x[[field]], na.rm = TRUE),
    maximum = if (all(is.na(x[[field]]))) NA_real_ else max(x[[field]], na.rm = TRUE)
  )))
}

key_fields <- c(
  "age_years", "bmi", "pir", "albumin_g_dl", "creatinine_mg_dl", "egfr",
  paste0(names(molecular_weight), "_nmol_l"),
  paste0("log2_", names(molecular_weight)),
  "sum4_pfas_nmol_l", "log2_sum4_pfas", "LBXBPB", "LBXBCD",
  "log2_lead", "log2_cadmium", "followup_years", "death_all"
)
missingness <- do.call(rbind, lapply(split(data, data$cycle), function(x) {
  output <- audit_missingness(x, key_fields)
  output$cycle <- unique(x$cycle)
  output[, c("cycle", setdiff(names(output), "cycle"))]
}))

comment_fields <- list(
  pfoa = c("LBDPFOAL", "LBDNFOAL", "LBDBFOAL", "SDNPFOAL", "SDBPFOAL"),
  pfos = c("LBDPFOSL", "LBDNFOSL", "LBDMFOSL", "SDNPFOSL", "SDMPFOSL"),
  pfhxs = "LBDPFHSL", pfna = "LBDPFNAL",
  lead = "LBDBPBLC", cadmium = "LBDBCDLC"
)
lod_audit <- do.call(rbind, lapply(names(comment_fields), function(analyte) {
  fields <- intersect(comment_fields[[analyte]], names(data))
  do.call(rbind, lapply(split(data, data$cycle), function(x) {
    indicators <- as.matrix(x[, fields, drop = FALSE])
    any_comment <- rowSums(!is.na(indicators)) > 0
    below <- rowSums(indicators == 1, na.rm = TRUE) > 0
    data.frame(
      cycle = unique(x$cycle), analyte = analyte,
      n_with_comment = sum(any_comment), n_below_lod = sum(below & any_comment),
      percent_below_lod = if (sum(any_comment)) {
        100 * sum(below & any_comment) / sum(any_comment)
      } else NA_real_
    )
  }))
}))

isomer_audit <- do.call(rbind, lapply(
  c("2013-2014", "2015-2016", "2017-2018"), function(cycle) {
    x <- data[data$cycle == cycle, ]
    fields <- if (cycle == "2013-2014") {
      c("SSNPFOA", "SSBPFOA", "SSNPFOS", "SSMPFOS")
    } else c("LBXNFOA", "LBXBFOA", "LBXNFOS", "LBXMFOS")
    observed <- rowSums(!is.na(x[, fields]))
    data.frame(
      cycle = cycle, n = nrow(x), n_any_isomer = sum(observed > 0),
      n_complete_isomers = sum(observed == length(fields)),
      n_partial_isomers = sum(observed > 0 & observed < length(fields)),
      n_complete_sum4 = sum(!is.na(x$sum4_pfas_nmol_l))
    )
  }
))

creatinine_audit <- do.call(rbind, lapply(split(data, data$cycle), function(x) {
  data.frame(
    cycle = unique(x$cycle), n = sum(!is.na(x$creatinine_uncalibrated)),
    mean_uncalibrated = mean(x$creatinine_uncalibrated, na.rm = TRUE),
    mean_calibrated = mean(x$creatinine_mg_dl, na.rm = TRUE),
    n_changed = sum(x$creatinine_uncalibrated != x$creatinine_mg_dl, na.rm = TRUE)
  )
}))

weight_audit <- do.call(rbind, lapply(split(data, data$cycle), function(x) {
  data.frame(
    cycle = unique(x$cycle), pfas_weight_source = unique(x$pfas_weight_source),
    metal_weight_source = unique(x$metal_weight_source),
    n_positive_pfas_weight = sum(x$pfas_weight_2yr > 0, na.rm = TRUE),
    n_positive_metal_weight = sum(x$metal_weight_2yr > 0, na.rm = TRUE)
  )
}))

age_topcode_audit <- do.call(rbind, lapply(split(data, data$cycle), function(x) {
  data.frame(
    cycle = unique(x$cycle),
    n = nrow(x),
    n_age_above_80 = sum(x$age_years_raw > 80, na.rm = TRUE),
    percent_age_above_80 = 100 * mean(x$age_years_raw > 80, na.rm = TRUE),
    max_raw_age = max(x$age_years_raw, na.rm = TRUE),
    max_analysis_age = max(x$age_years, na.rm = TRUE)
  )
}))

utils::write.csv(missingness, file.path(paths$audit, "harmonized_missingness_ranges.csv"), row.names = FALSE)
utils::write.csv(lod_audit, file.path(paths$audit, "lod_fill_audit.csv"), row.names = FALSE)
utils::write.csv(isomer_audit, file.path(paths$audit, "pfas_isomer_completeness.csv"), row.names = FALSE)
utils::write.csv(creatinine_audit, file.path(paths$audit, "creatinine_calibration_audit.csv"), row.names = FALSE)
utils::write.csv(weight_audit, file.path(paths$audit, "cycle_weight_audit.csv"), row.names = FALSE)
utils::write.csv(age_topcode_audit, file.path(paths$audit, "age_topcode_audit.csv"), row.names = FALSE)
utils::write.csv(
  cause_code_availability_audit,
  file.path(paths$audit, "cause_code_availability_audit.csv"),
  row.names = FALSE
)

saveRDS(data, file.path(paths$processed, "harmonized_analysis_data.rds"), compress = "xz")
message("Harmonization complete: ", nrow(data), " records")
