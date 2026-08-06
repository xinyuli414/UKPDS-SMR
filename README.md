# Healthy Participant Effect: Standardised Mortality Ratios in the UK Prospective Diabetes Study (UKPDS)

Analysis code and meta-analysis data for the Brief Report:

> **Healthy Participant Effect: Standardised Mortality Ratios Over 44 Years in the UK Prospective Diabetes Study.**
> Xinyu Li, Jose Leal, Ruth L. Coleman, James Altunkaya, Helen Dakin, Rury R. Holman, Philip Clarke.
> University of Oxford. *(Manuscript in preparation — targeted at* Diabetes Care.)

## Overview

This study asks whether participants in a randomised trial have lower mortality than people with type 2 diabetes (T2D) in routine care — a **healthy-participant effect (HPE)**. We estimated **standardised mortality ratios (SMRs)** — observed deaths divided by the deaths expected if participants had died at UK general-population rates — across the 44-year follow-up of UKPDS (1977–2021), and benchmarked them against a meta-analysis of published SMRs from routine-care T2D cohorts.

Over a median 17.8 years of follow-up (99,571 person-years; 3,380 deaths among 5,102 participants), the overall UKPDS SMR was **1.63 (95% CI 1.58–1.69)**. Mortality was *below* the general population in the first three years after enrolment (SMR **0.84, 95% CI 0.68–1.03**) and rose progressively thereafter. The pooled SMR from routine-care cohorts in high-income countries with follow-up overlapping the UKPDS era was **1.64 (95% CI 1.39–1.94)**, which UKPDS did not reach until around year 11 — consistent with a healthy-participant effect early in the trial.

Expected deaths are derived from age-, sex-, and calendar-year–matched UK life tables from the [Human Mortality Database (HMD)](https://www.mortality.org). All SMR confidence intervals use the exact Poisson method; the meta-analysis uses a multilevel random-effects model with cluster-robust confidence intervals at the cohort level.

## An invitation to collaborate

One trial cannot tell us how general the healthy-participant effect is.  **We would like to repeat this analysis across as many trials as possible, and we are actively looking for collaborators.**

The analysis requires only the following, one row per participant:

| Field | Notes |
| --- | --- |
| Sex | Selects the matching life table |
| Date of birth, or age at entry | Year and month are sufficient |
| Date of entry | Randomisation or enrolment |
| Date of last follow-up | Death, or date last known alive |
| Vital status | All-cause death |

No biomarkers, treatment allocation, cause of death, or record linkage are needed; expected deaths are taken from national life tables that are already public (the [HMD](https://www.mortality.org) covers around 40 countries). Baseline age, smoking status, HbA1c, BMI, ethnicity and treatment arm are optional, and allow the stratified analyses reported in the paper to be reproduced.

**Individual-level data need not leave your institution.** The scripts here run locally against your own data and emit only aggregate SMRs by follow-up band, which are what the cross-trial comparison actually needs. We are happy to adapt the code to your data structure.

If this is of interest — whether you hold a trial dataset, or simply think a particular trial ought to be in scope — please get in touch (contact details at the bottom of this page). 


## Repository contents

| File | What it is |
| --- | --- |
| `UKPDS SMR 44y Main analysis 5102_20260723.r` | Base-case SMR analysis of all **5,102** enrolled participants |
| `UKPDS SMR 44y sensitivity analysis 4209_20260723.r` | Sensitivity analyses in the **4,209** randomised participants, including treatment arms |
| `Meta-analysis.R` | Meta-analysis of published routine-care T2D SMRs |
| `SMR_T2DM_input_final.xlsx` | Extracted literature data used by `Meta-analysis.R` |
| `README.md` | This file |

### `UKPDS SMR 44y Main analysis 5102_20260723.r`

The primary (base-case) analysis. It builds the analytic cohort from the UKPDS longitudinal dataset, splits each participant's follow-up into calendar-year person-time slices (time origin = **enrolment**), attaches age-, sex-, and calendar-year–matched HMD mortality rates, and computes the overall SMR and SMRs stratified by baseline sex, age, smoking, HbA1c, BMI, and ethnicity. It produces the descriptive baseline table, the "SMR by years since enrolment" and "SMR by calendar year" figures, and the six-panel subgroup figure — corresponding to the manuscript **Table 1** and **Figure 1** and to **ESM Table 5.1** and **ESM Figures 5.1 and 5.8**. Outputs are written to `Output/UKPDS_SMR_tables_5102.xlsx` (four sheets) and `Output/Figure/`.

### `UKPDS SMR 44y sensitivity analysis 4209_20260723.r`

Restricts the analysis to the **4,209 participants randomised after the three-month dietary run-in** and re-estimates SMRs on **two time origins** — time since enrolment and time since randomisation. It also reproduces the original UKPDS trial design to compare SMRs by treatment allocation: the glycaemic study (UKPDS 33 — conventional vs intensive sulfonylurea/insulin) and the overweight metformin substudy (UKPDS 34 — conventional vs metformin). It maps to **ESM Tables 5.2 and 5.3** and **ESM Figures 5.2–5.7 and 5.9**. Outputs are written to `Output/UKPDS_SMR_tables_4209.xlsx` (two tabs) and `Output/Figure/` (filenames suffixed `_4209`).

### `Meta-analysis.R`

Pools published all-cause SMRs for routine-care T2D populations to provide the external benchmark. Each SMR is analysed on the log scale, with within-study variance taken from the reported 95% CI (or approximated as 1/deaths when no CI is available). Estimates are pooled with a **multilevel random-effects meta-analysis** (`metafor::rma.mv`) with random intercepts for cohort / study / estimate, and headline confidence intervals are made **cluster-robust by cohort** (`clubSandwich`). It also runs the pre-specified subgroups (WHO region, World Bank income group, OECD membership, and overlap with the UKPDS era) and produces the forest plots behind manuscript **Figure 3** and the **ESM Appendix 4** figures. Outputs are a results workbook (`SMR_T2DM_results.xlsx`), combined forest-plot PDFs, and per-figure PNGs.

### `SMR_T2DM_input_final.xlsx`

The extracted literature-review dataset that `Meta-analysis.R` reads. It has two sheets (data begins on row 4):

- **`1_Paper_overall`** — one row per study, with study/cohort metadata (author, year, country, WHO region, income/OECD flags, study period) and overall, male, and female SMRs with confidence limits, death counts, and cohort sizes.
- **`2_Sub_period`** — additional SMR estimates reported for calendar sub-periods within a study.

The included studies (23 in the meta-analysis; 25 in the review) span 348,228 individuals with T2D and 92,811 deaths, and are listed in **ESM Table 4.1** of the paper.

## Data requirements

The UKPDS scripts depend on data files that are **not included in this repository**:

- **UKPDS individual-patient data** —  These contain confidential participant data and are not redistributed here (see *Data availability* below).
- **HMD UK life tables (1×1)** — `UK fltper_1x1.txt` (female) and `UK mltper_1x1.txt` (male), freely available from the [Human Mortality Database](https://www.mortality.org) after free registration.

The meta-analysis input (`SMR_T2DM_input_final.xlsx`) is included, so `Meta-analysis.R` can be run from this repository alone.

## Requirements and running

The scripts were developed in **R 4.5.1**. Required packages:

- **UKPDS scripts:** `dplyr`, `tidyr`, `lubridate`, `ggplot2`, `readr`, `haven`, `stringr`, `patchwork`, `grid`, `openxlsx`, `knitr`
- **Meta-analysis:** `metafor`, `clubSandwich`, `dplyr`, `openxlsx` (plus `haven`, `tidyr`, `lubridate`, `readr`, `stringr`, `ggplot2`, `patchwork`)

```r
install.packages(c(
  "dplyr", "tidyr", "lubridate", "ggplot2", "readr", "haven", "stringr",
  "patchwork", "openxlsx", "knitr", "metafor", "clubSandwich"
))
```

Each script sets its input/output locations near the top with absolute paths (e.g. `setwd("J:/3rd project_SMR")` and a `ROOT` path in `Meta-analysis.R`). **Edit these paths** to point at your own copy of the data before running, then run a script end to end, for example:

```bash
Rscript "UKPDS SMR 44y Main analysis 5102_20260723.r"
```

Outputs are written to an `Output/` folder (tables and `Output/Figure/`) for the UKPDS scripts and to the configured `outputs/` folder for the meta-analysis.

## Data availability

UKPDS individual-participant data are not publicly available owing to participant confidentiality and data-governance restrictions; they may be requested through the UKPDS/Diabetes Trials Unit, University of Oxford, subject to approval. UK general-population life tables are freely available from the Human Mortality Database. The data underlying the meta-analysis were extracted from published studies and are provided in `SMR_T2DM_input_final.xlsx`.

## Citation

If you use this code, please cite the accompanying paper (details to be added on publication) and this repository. A `CITATION.cff` file can be added once the paper is published.

## Contact

Xinyu Li — <xinyu.li@ndph.ox.ac.uk>
Economics of Population Health Research Centre, Nuffield Department of Population Health, University of Oxford.
