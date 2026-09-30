select_available <- function(data, variables) {
  data[, intersect(variables, names(data)), drop = FALSE]
}

assert_unique_key <- function(data, key = "SEQN", label = deparse(substitute(data))) {
  if (!(key %in% names(data))) stop(label, " does not contain ", key)
  if (anyNA(data[[key]])) stop(label, " contains missing ", key)
  if (anyDuplicated(data[[key]])) stop(label, " contains duplicated ", key)
  invisible(TRUE)
}

left_join_one_to_one <- function(x, y, key = "SEQN", label = "source") {
  assert_unique_key(x, key, "merge anchor")
  assert_unique_key(y, key, label)
  result <- merge(x, y, by = key, all.x = TRUE, sort = FALSE)
  if (nrow(result) != nrow(x)) stop("Row count changed after merging ", label)
  result
}

read_xpt_fields <- function(path, fields, label = basename(path)) {
  if (!file.exists(path)) stop("Missing XPT source: ", path)
  data <- haven::read_xpt(path)
  missing <- setdiff(fields, names(data))
  if (length(missing)) {
    stop(label, " is missing required fields: ", paste(missing, collapse = ", "))
  }
  data <- as.data.frame(data[, fields, drop = FALSE])
  assert_unique_key(data, label = label)
  data
}

read_mortality_file <- function(path) {
  if (!file.exists(path)) stop("Missing mortality file: ", path)
  widths <- c(6L, 8L, 1L, 1L, 3L, 1L, 1L, 21L, 3L, 3L)
  names <- c(
    "SEQN", "skip_1", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING",
    "DIABETES_MCOD", "HYPERTEN_MCOD", "skip_2", "PERMTH_INT", "PERMTH_EXM"
  )
  data <- utils::read.fwf(
    path, widths = widths, col.names = names,
    colClasses = c("integer", "character", rep("integer", 5), "character", "integer", "integer"),
    strip.white = TRUE, blank.lines.skip = TRUE, na.strings = c("", ".")
  )
  data <- data[, !grepl("^skip_", names(data)), drop = FALSE]
  assert_unique_key(data, label = basename(path))
  if (!all(stats::na.omit(data$ELIGSTAT) %in% 1:3)) {
    stop("Invalid ELIGSTAT values in ", basename(path))
  }
  if (!all(stats::na.omit(data$MORTSTAT) %in% 0:1)) {
    stop("Invalid MORTSTAT values in ", basename(path))
  }
  if (!all(stats::na.omit(data$PERMTH_EXM) >= 0)) {
    stop("Invalid PERMTH_EXM values in ", basename(path))
  }
  data
}

derive_smoking <- function(ever_100, current) {
  result <- rep(NA_character_, length(ever_100))
  result[ever_100 == 2] <- "Never"
  result[ever_100 == 1 & current == 3] <- "Former"
  result[ever_100 == 1 & current %in% c(1, 2)] <- "Current"
  factor(result, levels = c("Never", "Former", "Current"))
}

derive_diabetes_nonfasting <- function(self_report, hba1c) {
  result <- rep(NA_integer_, length(self_report))
  result[self_report == 1 | (!is.na(hba1c) & hba1c >= 6.5)] <- 1L
  result[self_report %in% c(2, 3) & !is.na(hba1c) & hba1c < 6.5] <- 0L
  result
}

derive_cvd_history <- function(data) {
  fields <- c("MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F")
  values <- as.matrix(data[, fields, drop = FALSE])
  result <- rep(NA_integer_, nrow(values))
  result[rowSums(values == 1, na.rm = TRUE) > 0] <- 1L
  valid_response <- values == 1 | values == 2
  fully_observed <- rowSums(valid_response, na.rm = TRUE) == ncol(values)
  result[fully_observed & rowSums(values == 1, na.rm = TRUE) == 0] <- 0L
  result
}

derive_binary_history <- function(response) {
  result <- rep(NA_integer_, length(response))
  result[response == 1] <- 1L
  result[response == 2] <- 0L
  result
}

ckd_epi_2021 <- function(creatinine_mg_dl, age_years, female) {
  kappa <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  ratio <- creatinine_mg_dl / kappa
  142 * pmin(ratio, 1)^alpha * pmax(ratio, 1)^(-1.200) *
    0.9938^age_years * ifelse(female, 1.012, 1)
}

complete_sum <- function(...) {
  values <- cbind(...)
  result <- rowSums(values)
  result[rowSums(is.na(values)) > 0] <- NA_real_
  result
}

summarize_cohort_followup <- function(data, cohort_flag, cohort_name) {
  required <- c("followup_years", "death_all")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Follow-up summary is missing required fields: ",
         paste(missing, collapse = ", "))
  }
  if (length(cohort_flag) != nrow(data) || !is.logical(cohort_flag)) {
    stop("cohort_flag must be a logical vector with one value per row")
  }
  if (anyNA(cohort_flag)) stop("cohort_flag contains missing values")
  if (length(cohort_name) != 1L || is.na(cohort_name) || !nzchar(cohort_name)) {
    stop("cohort_name must be one non-empty string")
  }

  rows <- which(cohort_flag)
  if (!length(rows)) stop("Cohort contains no participants: ", cohort_name)
  followup <- data$followup_years[rows]
  deaths <- data$death_all[rows]
  if (any(!is.finite(followup)) || any(followup <= 0)) {
    stop("Follow-up must be finite and positive in cohort: ", cohort_name)
  }
  if (anyNA(deaths) || !all(deaths %in% 0:1)) {
    stop("death_all must be complete and binary in cohort: ", cohort_name)
  }

  quartiles <- stats::quantile(
    followup, probs = c(0.25, 0.50, 0.75), names = FALSE, type = 7
  )
  data.frame(
    cohort = cohort_name,
    n = length(rows),
    deaths_all = sum(deaths),
    observed_person_years = sum(followup),
    followup_q1_years = quartiles[[1]],
    followup_median_years = quartiles[[2]],
    followup_q3_years = quartiles[[3]],
    stringsAsFactors = FALSE
  )
}
