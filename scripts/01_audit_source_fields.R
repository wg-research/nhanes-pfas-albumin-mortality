#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
paths <- project_paths()
ensure_output_directories(paths)

cycles <- data.frame(
  cycle = c(
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2013-2014", "2015-2016", "2017-2018"
  ),
  suffix = LETTERS[3:10],
  stringsAsFactors = FALSE
)

source_map <- data.frame(
  cycle = cycles$cycle,
  demo = paste0("DEMO_", cycles$suffix),
  biochemistry = c("L40_C", paste0("BIOPRO_", LETTERS[4:10])),
  metals = c("L06BMT_C", paste0("PBCD_", LETTERS[4:10])),
  pfas_standard = c("L24PFC_C", "PFC_D", "PFC_E", "PFC_F", "PFC_G",
                    "PFAS_H", "PFAS_I", "PFAS_J"),
  pfas_surplus = c(rep(NA_character_, 5), "SSPFAS_H", NA_character_, NA_character_),
  stringsAsFactors = FALSE
)

component_dir <- c(
  demo = "Demographics",
  biochemistry = "Laboratory",
  metals = "Laboratory",
  pfas_standard = "Laboratory",
  pfas_surplus = "Laboratory"
)

read_metadata <- function(cycle, component, stem) {
  if (is.na(stem)) return(NULL)
  folder <- component_dir[[component]]
  xpt <- file.path(paths$raw, cycle, folder, paste0(stem, ".xpt"))
  if (!file.exists(xpt)) {
    stop("Missing source XPT: ", xpt)
  }
  data <- haven::read_xpt(xpt, n_max = 0)
  labels <- vapply(data, function(x) {
    value <- attr(x, "label", exact = TRUE)
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1))
  data.frame(
    cycle = cycle,
    component = component,
    source_stem = stem,
    xpt_path = xpt,
    variable = names(data),
    label = unname(labels),
    stringsAsFactors = FALSE
  )
}

metadata <- list()
for (i in seq_len(nrow(source_map))) {
  for (component in names(component_dir)) {
    item <- read_metadata(
      source_map$cycle[[i]], component, source_map[[component]][[i]]
    )
    if (!is.null(item)) metadata[[length(metadata) + 1L]] <- item
  }
}
metadata <- do.call(rbind, metadata)

target_patterns <- paste(
  c(
    "^SEQN$", "^SDMVSTRA$", "^SDMVPSU$", "^WT.*2YR$",
    "^RIDAGEYR$", "^RIAGENDR$", "^RIDRETH[13]$", "^RIDEXPRG$",
    "^DMDEDUC2$", "^INDFMPIR$", "^LBXSAL$", "^LBXSCR$",
    "^LBXBPB$", "^LBXBCD$", "^LBDBPBLC$", "^LBDBCDLC$",
    "PFOA", "PFOS", "PFHS", "PFNA", "PFHxS"
  ),
  collapse = "|"
)
target_metadata <- metadata[grepl(target_patterns, metadata$variable, ignore.case = TRUE), ]

utils::write.csv(source_map, file.path(paths$audit, "source_file_map.csv"), row.names = FALSE)
utils::write.csv(metadata, file.path(paths$audit, "source_variable_inventory.csv"), row.names = FALSE)
utils::write.csv(target_metadata, file.path(paths$audit, "target_field_inventory.csv"), row.names = FALSE)

required_by_component <- list(
  demo = c("SEQN", "RIDAGEYR", "RIAGENDR", "SDMVSTRA", "SDMVPSU", "WTMEC2YR"),
  biochemistry = c("SEQN", "LBXSAL", "LBXSCR"),
  metals = c("SEQN", "LBXBPB", "LBXBCD"),
  pfas_standard = "SEQN",
  pfas_surplus = c("SEQN", "WTSSBH2Y")
)

checks <- list()
for (i in seq_len(nrow(source_map))) {
  for (component in names(required_by_component)) {
    stem <- source_map[[component]][[i]]
    if (is.na(stem)) next
    available <- metadata$variable[
      metadata$cycle == source_map$cycle[[i]] & metadata$component == component
    ]
    for (field in required_by_component[[component]]) {
      checks[[length(checks) + 1L]] <- data.frame(
        cycle = source_map$cycle[[i]],
        component = component,
        source_stem = stem,
        field = field,
        present = field %in% available,
        stringsAsFactors = FALSE
      )
    }
  }
}
checks <- do.call(rbind, checks)
utils::write.csv(checks, file.path(paths$audit, "required_field_checks.csv"), row.names = FALSE)

if (any(!checks$present)) {
  missing <- checks[!checks$present, ]
  stop("Required source fields are missing. See required_field_checks.csv: ",
       paste(paste(missing$cycle, missing$component, missing$field, sep = "/"),
             collapse = ", "))
}

message("Source field audit passed: ", nrow(metadata), " variables verified and inventoried")
