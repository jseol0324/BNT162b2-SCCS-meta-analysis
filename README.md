# BNT162b2-SCCS-meta-analysis

Reproducible code and data for a systematic review and meta-analysis of stroke risk
following BNT162b2 vaccination using self-controlled case series (SCCS) studies.

## Reproducibility

All analyses were conducted in R (version ≥ 4.2.0) using the metafor package.
Study-level data and analysis scripts are publicly available in this repository.

Study-level data are provided in the `data/` directory:
- `data/dat.csv`: study-level effect estimates
- `data/meta_info.csv`: moderator and subgroup coding

## How to run

```r
source("R/01_meta_analysis.R")
source("R/02_meta_regression.R")
source("R/03_subgroup_analysis.R")
