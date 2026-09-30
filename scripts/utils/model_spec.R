base_model_fields <- function() {
  c(
    "age_years", "sex", "race_ethnicity", "education", "pir", "smoking",
    "bmi", "cycle_factor"
  )
}

physiology_fields <- function() c("albumin_g_dl", "egfr")

design_fields <- function() c("SDMVSTRA", "SDMVPSU")

complete_fields <- function(data, fields) {
  stats::complete.cases(data[, fields, drop = FALSE])
}

valid_mortality_followup <- function(data) {
  data$age_years >= 20 & data$mortality_eligible &
    !is.na(data$followup_years) & data$followup_years > 0 &
    !is.na(data$death_all)
}
