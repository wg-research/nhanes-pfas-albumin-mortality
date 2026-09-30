#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
paths <- project_paths()
ensure_output_directories(paths)

analysis_packages <- c(
  "digest", "ggplot2", "haven", "Hmisc", "jsonlite", "mice",
  "patchwork", "ragg", "splines", "survey", "survival"
)

package_status <- data.frame(
  package = analysis_packages,
  available = vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    analysis_packages,
    function(pkg) if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    },
    character(1)
  )
)

if (!all(package_status$available)) {
  stop("Missing required R packages: ",
       paste(package_status$package[!package_status$available], collapse = ", "))
}

expected_inputs <- c(
  paths$raw,
  paths$mortality
)
if (!all(file.exists(expected_inputs))) {
  stop("Missing required project inputs: ",
       paste(expected_inputs[!file.exists(expected_inputs)], collapse = ", "))
}

environment_record <- list(
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  project_root = paths$root,
  r_version = R.version.string,
  platform = R.version$platform,
  packages = package_status,
  survey_lonely_psu = "adjust",
  primary_pfas_cycles = c(
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2015-2016", "2017-2018"
  ),
  secondary_surplus_cycle = "2013-2014",
  random_seed = 20260131L
)

jsonlite::write_json(
  environment_record,
  file.path(paths$audit, "environment.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)
utils::write.csv(
  package_status,
  file.path(paths$audit, "package_versions.csv"),
  row.names = FALSE
)
writeLines(capture.output(utils::sessionInfo()), file.path(paths$results, "session_info.txt"))

message("Environment specification passed: ", paths$root)
