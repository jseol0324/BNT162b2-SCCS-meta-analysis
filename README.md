# BNT162b2-SCCS-meta-analysis

Reproducibility repository for:

> Jang S, Pyo S, Kim EJ, Lee W, Kwon S-S.  
> **Stroke risk following BNT162b2 vaccination: A systematic review and meta-analysis of self-controlled case series studies.**  
> *Scientific Reports* (submitted 2026).

---

## Overview

This repository contains the extracted dataset, data dictionary, and R analysis script used to reproduce all meta-analytic results reported in the manuscript. The analysis pools eight self-controlled case series (SCCS) studies examining the association between BNT162b2 (Pfizer–BioNTech) vaccination and stroke risk using a random-effects model with restricted maximum likelihood (REML) estimation and the Knapp–Hartung adjustment.

---

## Repository structure

```
BNT162b2-SCCS-meta-analysis/
├── extracted_data.csv      # Study-level pooled estimates used in meta-analysis
├── data_dictionary.csv     # Variable definitions and coding notes
├── meta_analysis.R         # Full R analysis script (Sections A–K)
└── README.md               # This file
```

---

## Data

### `extracted_data.csv`

Each row corresponds to one included study. For studies that reported only stratum-specific IRRs (Jabagi 2022, Botton 2022, Makkar 2025, Salmaggi 2025), a single within-study estimate was derived by inverse-variance weighting on the log scale prior to meta-analysis (see Section A of `meta_analysis.R`).

Key variables:

| Variable | Description |
|----------|-------------|
| `study_id` | Unique study identifier |
| `yi` | Log-transformed pooled within-study IRR |
| `sei` | Standard error of `yi` |
| `IRR` | Incidence rate ratio (original scale) |
| `LCL_95` / `UCL_95` | 95% confidence interval bounds |
| `risk_window_cat` | `<=21d` or `>21d` |
| `dose_cat` | `primary` (doses 1–2 only) or `booster` (includes dose ≥3) |
| `sccs_binary` | `Standard` or `Non-standard` |
| `outcome_type` | `Both` (IS+HS), `IS` (ischemic only), `HS` (hemorrhagic only) |
| `baseline_type` | `Post-vax`, `Mixed`, or `Pre-vax` |

See `data_dictionary.csv` for full variable definitions.

---

## Code

### `meta_analysis.R`

The script is organized into sections:

| Section | Content |
|---------|---------|
| A | Within-study inverse-variance pooling |
| B | Dataset construction and factor coding |
| C | Primary meta-analysis (REML + Knapp–Hartung) |
| D | Forest plot |
| E | Leave-one-out analysis |
| F | Subgroup analyses (continent, risk window, dose, SCCS variant, age group, outcome subtype, baseline type) |
| G | Univariable meta-regression |
| H | Sensitivity analyses (Sens i–vi) |
| I | Publication bias (Egger's test, funnel plot) |
| J | Influence diagnostics (Cook's distance, rstudent, DFBETAs, Baujat plot) |
| K | Output tables saved to `output/` |

---

## Requirements

```r
R >= 4.4.1
metafor >= 4.8.0
ggplot2
ggrepel
```

Install missing packages:

```r
install.packages(c("metafor", "ggplot2", "ggrepel"))
```

---

## Reproducing the analysis

```r
# Clone the repository and set working directory
setwd("BNT162b2-SCCS-meta-analysis")

# Run the full analysis
source("meta_analysis.R")

# Output tables are saved to output/
```

Results in `output/`:

- `Table_Main.csv` — Study-level and pooled IRR estimates with weights  
- `Table_Subgroups.csv` — All subgroup analysis results  
- `Table_LOO.csv` — Leave-one-out analysis results  

---

## Included studies

| ID | Author | Year | Country | k events |
|----|--------|------|---------|----------|
| HC2021 | Hippisley-Cox et al. | 2021 | UK | 6,439 |
| Chui2022 | Chui et al. | 2022 | Hong Kong | 117 |
| Jabagi2022 | Jabagi et al. | 2022 | France | 17,014 |
| Botton2022 | Botton et al. | 2022 | France | 28,717 |
| AbRahman2024 | Ab Rahman et al. | 2024 | Malaysia | 74,700 |
| Xu2024 | Xu et al. | 2024 | USA | 1,057 |
| Makkar2025 | Makkar et al. | 2025 | USA | 119,275 |
| Salmaggi2025 | Salmaggi et al. | 2025 | Italy | NA |

---

## Pooled result

| Model | IRR | 95% CI | I² | τ² |
|-------|-----|--------|-----|-----|
| Random-effects (REML + KH), k=8 | 0.967 | 0.892–1.049 | 69.2% | 0.0026 |

---

## License

The code in this repository is released under the MIT License.  
The extracted data are derived from previously published studies; please cite the original sources accordingly.

---

## Contact

Soon-Sun Kwon, PhD  
Department of Mathematics, Ajou University  
qrio1010@ajou.ac.kr

