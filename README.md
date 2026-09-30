# Albumin conditioning differentially changes PFAS and blood-metal mortality estimates: a same-participant analysis of NHANES 2003–2018

The scripts in `scripts/` (run via `scripts/run_all.R`) build the analytic cohorts, fit the models, and write the CSV and figure sources behind the tables and figures. Inputs are public NHANES files linked to the 2019 public-use mortality files; generated outputs go to `data_rebuild/` and `results/`.

## Data

Download the NHANES examination, laboratory, questionnaire, and demographic files, and their documentation, from <https://www.cdc.gov/nchs/nhanes/> and the linked mortality files from <https://www.cdc.gov/nchs/data-linkage/mortality.htm>. Arrange the files under `Database/` like this:

```text
Database/
`-- nhanes/
    |-- 2003-2004/
    |   |-- Demographics/
    |   |-- Examination/
    |   |-- Laboratory/
    |   `-- Questionnaire/
    |-- 2005-2006/
    |-- 2007-2008/
    |-- 2009-2010/
    |-- 2011-2012/
    |-- 2013-2014/
    |-- 2015-2016/
    |-- 2017-2018/
    `-- linked_mortality/
```

Only these eight cycles (2003–2018) are read. One file per cycle and component:

| Cycle | Demographics | Biochemistry | PFAS | Blood metals | Body measures | Smoking | Diabetes | Medical conditions | HbA1c |
|---|---|---|---|---|---|---|---|---|---|
| 2003-2004 | `DEMO_C` | `L40_C` | `L24PFC_C` | `L06BMT_C` | `BMX_C` | `SMQ_C` | `DIQ_C` | `MCQ_C` | `L10_C` |
| 2005-2006 | `DEMO_D` | `BIOPRO_D` | `PFC_D` | `PBCD_D` | `BMX_D` | `SMQ_D` | `DIQ_D` | `MCQ_D` | `GHB_D` |
| 2007-2008 | `DEMO_E` | `BIOPRO_E` | `PFC_E` | `PBCD_E` | `BMX_E` | `SMQ_E` | `DIQ_E` | `MCQ_E` | `GHB_E` |
| 2009-2010 | `DEMO_F` | `BIOPRO_F` | `PFC_F` | `PBCD_F` | `BMX_F` | `SMQ_F` | `DIQ_F` | `MCQ_F` | `GHB_F` |
| 2011-2012 | `DEMO_G` | `BIOPRO_G` | `PFC_G` | `PBCD_G` | `BMX_G` | `SMQ_G` | `DIQ_G` | `MCQ_G` | `GHB_G` |
| 2013-2014 | `DEMO_H` | `BIOPRO_H` | `PFAS_H` and `SSPFAS_H` | `PBCD_H` | `BMX_H` | `SMQ_H` | `DIQ_H` | `MCQ_H` | `GHB_H` |
| 2015-2016 | `DEMO_I` | `BIOPRO_I` | `PFAS_I` | `PBCD_I` | `BMX_I` | `SMQ_I` | `DIQ_I` | `MCQ_I` | `GHB_I` |
| 2017-2018 | `DEMO_J` | `BIOPRO_J` | `PFAS_J` | `PBCD_J` | `BMX_J` | `SMQ_J` | `DIQ_J` | `MCQ_J` | `GHB_J` |

`Demographics` holds `DEMO_*`, `Examination` holds `BMX_*`, `Laboratory` holds the biochemistry, PFAS, metals, and HbA1c files, and `Questionnaire` holds `SMQ_*`, `DIQ_*`, and `MCQ_*`. For example, the 2003–2004 demographics file goes at `Database/nhanes/2003-2004/Demographics/DEMO_C.xpt`.

`linked_mortality/` needs the eight files `NHANES_2003_2004_MORT_2019_PUBLIC.dat` through `NHANES_2017_2018_MORT_2019_PUBLIC.dat`. The full file mapping is also stored in `scripts/utils/source_map.R`.

## Requirements

R 4.5.3, plus `digest`, `ggplot2`, `haven`, `Hmisc`, `jsonlite`, `mice`, `patchwork`, `ragg`, `splines`, `survey`, `survival`. The exact versions used for the analysis are pinned in `renv.lock`; install them with:

```r
renv::restore()
```

or, without renv:

```r
install.packages(c(
  "digest", "ggplot2", "haven", "Hmisc", "jsonlite", "mice",
  "patchwork", "ragg", "survey", "survival"
))
```

The covariance implementation calls the non-exported `survey:::svyrecvar` function. Restore the locked environment before rebuilding so that the recorded `survey` version is used. In the sensitivity outputs, `pre_2015_overlap_cycles` denotes the five-cycle 2003–2012 overlap analysis.

PDF figures need a working Cairo graphics device.

## Running

From the repository root:

```sh
Rscript --vanilla scripts/run_all.R
```

On Windows you may need the full path, e.g. `& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" --vanilla scripts/run_all.R`.

The stages run in order and stop at the first failed check; a successful run prints `Completed available rebuild stages`. You can run `Rscript --vanilla scripts/00_environment_specification.R` first to confirm the packages and data are in place. Stage 07 (multiple imputation) is the slowest, roughly 1–2 hours. Outputs go to `results/latest_release/` and intermediates to `data_rebuild/`.

The full proportional-hazards interaction output is `ph_time_interaction_tests.csv`. The singular `ph_time_interaction_test.csv` is retained as a one-row PFAS-primary compatibility file for older downstream checks.

As a sanity check, the three main cohorts should come out to 10,169, 8,705, and 31,469 participants (1,328, 1,266, and 4,110 deaths).

## License

MIT; see `LICENSE`. The NHANES source files themselves remain subject to NCHS terms.

## Citation

Please cite the associated article and the archived reproducibility package:

https://doi.org/10.5281/zenodo.23050962
