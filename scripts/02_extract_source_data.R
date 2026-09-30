#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "source_map.R"))
source(file.path("scripts", "utils", "data_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
map <- nhanes_source_map()

file_path <- function(cycle, component, stem) {
  folder <- switch(
    component,
    demo = "Demographics",
    bmi = "Examination",
    smoking = "Questionnaire",
    diabetes = "Questionnaire",
    medical = "Questionnaire",
    "Laboratory"
  )
  file.path(paths$raw, cycle, folder, paste0(stem, ".xpt"))
}

extract_cycle <- function(spec) {
  cycle <- spec$cycle
  message("Extracting ", cycle)

  demo_fields <- c(
    "SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "RIDEXPRG",
    "DMDEDUC2", "INDFMPIR", "SDMVSTRA", "SDMVPSU", "WTMEC2YR"
  )
  data <- read_xpt_fields(file_path(cycle, "demo", spec$demo), demo_fields)

  biochemistry <- read_xpt_fields(
    file_path(cycle, "biochemistry", spec$biochemistry),
    c("SEQN", "LBXSAL", "LBXSCR")
  )
  data <- left_join_one_to_one(data, biochemistry, label = spec$biochemistry)

  bmi <- read_xpt_fields(file_path(cycle, "bmi", spec$bmi), c("SEQN", "BMXBMI"))
  data <- left_join_one_to_one(data, bmi, label = spec$bmi)

  smoking <- read_xpt_fields(
    file_path(cycle, "smoking", spec$smoking), c("SEQN", "SMQ020", "SMQ040")
  )
  data <- left_join_one_to_one(data, smoking, label = spec$smoking)

  diabetes <- read_xpt_fields(
    file_path(cycle, "diabetes", spec$diabetes), c("SEQN", "DIQ010")
  )
  data <- left_join_one_to_one(data, diabetes, label = spec$diabetes)

  medical <- read_xpt_fields(
    file_path(cycle, "medical", spec$medical),
    c("SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ220")
  )
  data <- left_join_one_to_one(data, medical, label = spec$medical)

  hba1c <- read_xpt_fields(
    file_path(cycle, "hba1c", spec$hba1c), c("SEQN", "LBXGH")
  )
  data <- left_join_one_to_one(data, hba1c, label = spec$hba1c)

  metals_path <- file_path(cycle, "metals", spec$metals)
  metal_fields <- c("SEQN", "LBXBPB", "LBXBCD")
  metal_raw <- haven::read_xpt(metals_path, n_max = 0)
  if (spec$metals_weight == "WTSH2YR") metal_fields <- c(metal_fields, "WTSH2YR")
  for (comment in c("LBDBPBLC", "LBDBCDLC")) {
    if (comment %in% names(metal_raw)) metal_fields <- c(metal_fields, comment)
  }
  metals <- read_xpt_fields(metals_path, metal_fields)
  data <- left_join_one_to_one(data, metals, label = spec$metals)

  pfas_path <- file_path(cycle, "pfas", spec$pfas)
  pfas_raw <- haven::read_xpt(pfas_path)
  if (cycle %in% c("2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012")) {
    pfas_fields <- c(
      "SEQN", spec$pfas_weight, "LBXPFOA", "LBXPFOS", "LBXPFHS", "LBXPFNA",
      "LBDPFOAL", "LBDPFOSL", "LBDPFHSL", "LBDPFNAL"
    )
    pfas <- as.data.frame(pfas_raw[, pfas_fields, drop = FALSE])
  } else if (cycle %in% c("2015-2016", "2017-2018")) {
    pfas_fields <- c(
      "SEQN", spec$pfas_weight, "LBXNFOA", "LBXBFOA", "LBXNFOS", "LBXMFOS",
      "LBXPFHS", "LBXPFNA", "LBDNFOAL", "LBDBFOAL", "LBDNFOSL",
      "LBDMFOSL", "LBDPFHSL", "LBDPFNAL"
    )
    pfas <- as.data.frame(pfas_raw[, pfas_fields, drop = FALSE])
  } else {
    pfas <- as.data.frame(pfas_raw[, c("SEQN", "WTSB2YR", "LBXPFHS", "LBXPFNA",
                                      "LBDPFHSL", "LBDPFNAL"), drop = FALSE])
    surplus <- read_xpt_fields(
      file_path(cycle, "pfas_surplus", spec$pfas_surplus),
      c("SEQN", "WTSSBH2Y", "SSNPFOA", "SSBPFOA", "SSNPFOS", "SSMPFOS",
        "SDNPFOAL", "SDBPFOAL", "SDNPFOSL", "SDMPFOSL")
    )
    pfas <- merge(surplus, pfas, by = "SEQN", all.x = TRUE, sort = FALSE)
  }
  assert_unique_key(pfas, label = paste(cycle, "PFAS"))
  data <- left_join_one_to_one(data, pfas, label = paste(cycle, "PFAS"))

  mortality <- read_mortality_file(file.path(paths$mortality, spec$mortality))
  data <- left_join_one_to_one(data, mortality, label = spec$mortality)

  data$cycle <- cycle
  data$pfas_weight_source <- spec$pfas_weight
  data$metal_weight_source <- spec$metals_weight
  data
}

cycles <- lapply(seq_len(nrow(map)), function(i) extract_cycle(as.list(map[i, ])))
combined <- do.call(rbind, lapply(cycles, function(data) {
  missing <- setdiff(unique(unlist(lapply(cycles, names))), names(data))
  for (field in missing) data[[field]] <- NA
  data[, unique(unlist(lapply(cycles, names))), drop = FALSE]
}))

assert_unique_key(combined, label = "combined extraction")
saveRDS(combined, file.path(paths$processed, "extracted_source_data.rds"), compress = "xz")

cycle_counts <- aggregate(SEQN ~ cycle, combined, length)
names(cycle_counts)[2] <- "n_demo"
utils::write.csv(cycle_counts, file.path(paths$audit, "extracted_cycle_counts.csv"), row.names = FALSE)

message("Extraction complete: ", nrow(combined), " unique participants")

