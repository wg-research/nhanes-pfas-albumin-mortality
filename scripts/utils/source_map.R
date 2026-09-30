nhanes_source_map <- function() {
  data.frame(
    cycle = c(
      "2003-2004", "2005-2006", "2007-2008", "2009-2010",
      "2011-2012", "2013-2014", "2015-2016", "2017-2018"
    ),
    suffix = LETTERS[3:10],
    demo = paste0("DEMO_", LETTERS[3:10]),
    biochemistry = c("L40_C", paste0("BIOPRO_", LETTERS[4:10])),
    pfas = c("L24PFC_C", "PFC_D", "PFC_E", "PFC_F", "PFC_G",
             "PFAS_H", "PFAS_I", "PFAS_J"),
    pfas_surplus = c(rep(NA_character_, 5), "SSPFAS_H", NA_character_, NA_character_),
    pfas_weight = c("WTSA2YR", "WTSA2YR", "WTSC2YR", "WTSC2YR",
                    "WTSA2YR", "WTSSBH2Y", "WTSB2YR", "WTSB2YR"),
    metals = c("L06BMT_C", paste0("PBCD_", LETTERS[4:10])),
    metals_weight = c(rep("WTMEC2YR", 5), "WTSH2YR", "WTSH2YR", "WTMEC2YR"),
    bmi = paste0("BMX_", LETTERS[3:10]),
    smoking = paste0("SMQ_", LETTERS[3:10]),
    diabetes = paste0("DIQ_", LETTERS[3:10]),
    medical = paste0("MCQ_", LETTERS[3:10]),
    hba1c = c("L10_C", paste0("GHB_", LETTERS[4:10])),
    mortality = paste0(
      "NHANES_", gsub("-", "_", c(
        "2003-2004", "2005-2006", "2007-2008", "2009-2010",
        "2011-2012", "2013-2014", "2015-2016", "2017-2018"
      )), "_MORT_2019_PUBLIC.dat"
    ),
    stringsAsFactors = FALSE
  )
}

primary_pfas_cycles <- function() {
  c(
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2015-2016", "2017-2018"
  )
}

primary_overlap_cycles <- function() {
  c(
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2017-2018"
  )
}

