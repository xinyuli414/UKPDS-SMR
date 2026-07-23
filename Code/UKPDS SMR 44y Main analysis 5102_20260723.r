# UKPDS Standardised Mortality Ratio (SMR) analysis (44 years of follow-up) ====
# Author: Xinyu Li (xinyu.li@ndph.ox.ac.uk)
# Date:   2026/June/20
#
# Purpose
#   For the UKPDS cohort, compare the number of deaths observed during
#   follow-up against the number that would be expected if patients had died
#   at UK general-population rates. The ratio of the two is the Standardised
#   Mortality Ratio (SMR = observed deaths / expected deaths). Expected deaths
#   are calculated from the Human Mortality Database (HMD) 1x1 UK life tables.
#   This version works with the 44-year UKPDS dataset.
#
# Time scale
#   Time since enrolment, where the enrolment date is gr1date.
#
# Outputs
#   Output/UKPDS_SMR_tables.xlsx (one file, four sheets):
#       - "Baseline characteristics" : Table 1 for the full cohort
#       - "Cohort summary"           : follow-up, person-years, deaths, SMR
#       - "SMR by subgroup (3-year)" : SMR overall and by 3-year band
#       - "SMR by subgroup (1-year)" : SMR overall and by single-year band
#   Output/Figure/SMR_years_since_enrolment.png
#   Output/Figure/SMR_calendar_year.png
#   Output/Figure/SMR_subgroups_combined_3yr.png


# 1. SETUP =====================================================================
# Aim:    Set the working directory, create the output folders, load packages,
#         and list the source file paths.
# Input:  None.
# Output: Output/ and Output/Figure/ folders (created if they do not exist).

## 1.1 Working directory and output folders ----
setwd("K:/HERC_UKPDS/Xinyu Li Risk equations/Xinyu UKPDS SMR project")


# Output folders (created if they do not already exist)
out_dir <- "Output"
dir.create(out_dir,                      showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "Figure"), showWarnings = FALSE, recursive = TRUE)

## 1.2 Packages ----
suppressPackageStartupMessages({
  library(dplyr)      # data manipulation
  library(tidyr)      # reshaping between long and wide
  library(lubridate)  # working with dates
  library(ggplot2)    # plotting
  library(readr)      # reading text files (life tables)
  library(haven)      # reading Stata .dta files
  library(stringr)    # string handling
  library(patchwork)  # combining plots into panel layouts
  library(grid)       # low-level graphics for shared axis labels
  library(openxlsx)   # writing a multi-sheet Excel workbook
})


## 1.3 File paths ----
path_extended_stata <- "K:/QNAP/UKPDS/DTU data/XL Death Data.dta"             # extra death/censoring dates
path_stata          <- "K:/HERC_UKPDS/UKPDS-OM2 model/Data/locf_anal1Aug11.dta"           # main longitudinal dataset
path_lt_female      <- "Data/Life tables/UK fltper_1x1.txt" # HMD female life table
path_lt_male        <- "Data/Life tables/UK mltper_1x1.txt" # HMD male life table


# --- Software versions: UKPDS SMR main + sensitivity analysis ---
pkgs_main <- c("dplyr", "tidyr", "lubridate", "ggplot2", "readr",
               "haven", "stringr", "patchwork", "openxlsx", "knitr")

pkg_versions <- sapply(pkgs_main, function(p) as.character(packageVersion(p)))

# Echo to console
cat(R.version.string, "\n")
print(pkg_versions)

# Write R version + package versions + full sessionInfo() to one file
writeLines(
  c(R.version.string,
    "",
    "Package versions:",
    sprintf("  %-10s %s", names(pkg_versions), pkg_versions),
    "",
    "sessionInfo():",
    capture.output(sessionInfo())),
  file.path(out_dir, "session_info.txt")
)


# 2. DATA LOADING ==============================================================
# Aim:    Read the two source datasets and attach the extra death and
#         censoring dates to the main longitudinal dataset.
# Input:  path_stata (main dataset), path_extended_stata (extra dates).
# Output: locf, the main dataset with the extra date columns merged in by
#         patient number.

locf           <- read_dta(path_stata)
extended_UKPDS <- read_dta(path_extended_stata)

locf <- merge(
  locf,
  extended_UKPDS[, c("ukpdsno", "anydeath_event", "anydeath_dt",
                     "cens_dt2", "GL_date")],
  by = "ukpdsno"
)


# 3. BASELINE FRAMES ===========================================================
# Aim:    Build the three tidy data frames the rest of the script relies on.
# Input:  locf.
# Output: UKPDS          : one row per patient for the survival/SMR analysis.
#         baseline_t1    : baseline values for the descriptive Table 1.
#         ukpds_baseline : baseline categories used for the subgroup SMRs.

## 3.1 Analytic cohort (one row per patient) ----
# Aim:    For each patient, define the dates that bound their follow-up and
#         whether they died.
# Output: id, sex, date of birth, enrolment date, end-of-follow-up date, death
#         indicator and randomisation date. Patients with missing key dates,
#         or whose follow-up ends before enrolment, are dropped.
UKPDS <- locf %>%
  select(ukpdsno, year, female, dob_fu, gr1date, anydeath_event, anydeath_dt,
         cens_dt2, GL_date) %>%
  arrange(ukpdsno, year) %>%
  select(-year) %>%
  distinct() %>%
  transmute(
    id            = ukpdsno,
    sex           = if_else(female == 1, "Female", "Male"),
    dob           = as.Date(dob_fu),
    enrolment     = as.Date(gr1date),
    end           = as.Date(if_else(is.na(anydeath_dt), cens_dt2, anydeath_dt)),
    event         = anydeath_event,
    randomization = as.Date(GL_date)
  ) %>%
  filter(!is.na(dob), !is.na(enrolment), !is.na(end), end >= enrolment)

## 3.2 Descriptive Table 1 frame ----
# Aim:    Take each patient's first (baseline) record and keep the continuous
#         variables in their original units plus three categorical variables
#         (sex, smoking, ethnicity) with explicit factor levels.
# Output: baseline_t1, one row per patient.
baseline_t1 <- locf %>%
  arrange(ukpdsno, year) %>%
  group_by(ukpdsno) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    id    = ukpdsno,
    sex   = factor(if_else(female == 1, "Female", "Male"),
                   levels = c("Male", "Female")),
    age   = as.numeric(curr_age),
    bmi   = as.numeric(bmi),
    sbp   = as.numeric(sbp),
    dbp   = as.numeric(dbp),
    hba1c = as.numeric(hba1c),
    chol  = as.numeric(chol),
    ldl   = as.numeric(ldl),
    hdl   = as.numeric(hdl),
    
    # 2-level smoking (smoker: 0/1). Switch to smoker012 for never/ex/current.
    smoking = factor(smoker, levels = c(0, 1),
                     labels = c("Non-current smoker", "Current smoker")),
    
    # 4-level ethnicity from race / indian / afro flags.
    ethnicity = factor(case_when(
      race   == 1 ~ "White",
      indian == 1 ~ "Indian Asian",
      afro   == 1 ~ "Afro-Caribbean",
      TRUE        ~ "Other"),
      levels = c("White", "Indian Asian", "Afro-Caribbean", "Other"))
  )

## 3.3 Subgroup categories frame ----
# Aim:    Take each patient's baseline record and bin it into the subgroup
#         categories used in Figure 3 and the SMR table. Categories match the
#         manuscript's subgroup definitions.
# Output: ukpds_baseline, one row per patient.
ukpds_baseline <- locf %>%
  arrange(ukpdsno, year) %>%
  group_by(ukpdsno) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    id          = ukpdsno,
    age_cat     = cut(curr_age, breaks = c(-Inf, 50, 60, Inf),
                      labels = c("<50 years", "50-60 years", ">=60 years"),
                      right = FALSE),
    smoking_cat = factor(smoker, levels = c(0, 1),
                         labels = c("Non-current smoker", "Current smoker")),
    hba1c_cat   = cut(hba1c, breaks = c(-Inf, 7, Inf),
                      labels = c("<7% (<53 mmol/mol)", ">=7% (>=53 mmol/mol)"),
                      right = FALSE),
    bmi_cat     = cut(bmi, breaks = c(-Inf, 25, 30, Inf),
                      labels = c("<25 kg/m2", "25-30 kg/m2", ">=30 kg/m2"),
                      right = FALSE),
    ethnicity_cat = factor(
      ifelse(race == 1, "White", "Afro-Caribbean/Asian Indian/Other"),
      levels = c("White", "Afro-Caribbean/Asian Indian/Other"))
  )


# 4. DESCRIPTIVE TABLE 1: BASELINE CHARACTERISTICS =============================
# Aim:    Summarise baseline characteristics of the cohort: mean (SD) for
#         continuous variables and n (%) for categorical variables, each with
#         a missing-data row.
# Input:  baseline_t1.
# Output: table1_descriptive, a two-column (Variable, Value) data frame that
#         becomes the "Baseline characteristics" sheet in Section 12.

## 4.1 Formatting helpers ----
# fmt_n_pct():    format a count as "n (xx.x%)".
# fmt_mean_sd():  format a numeric vector as "mean (SD)", ignoring NAs.
# row_continuous(): build the mean(SD) row + a missing row for one variable.
# row_categorical(): build a header row + one row per level + a missing row.
N_total <- n_distinct(baseline_t1$id)

fmt_n_pct <- function(n, d = N_total) {
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / d)
}

fmt_mean_sd <- function(x) {
  ok <- !is.na(x)
  sprintf("%.1f (%.1f)", mean(x[ok]), sd(x[ok]))
}

row_continuous <- function(x, label) {
  bind_rows(
    data.frame(Variable = paste0(label, ", mean (SD)"),
               Value    = fmt_mean_sd(x),
               stringsAsFactors = FALSE),
    data.frame(Variable = "  Missing, n (%)",
               Value    = fmt_n_pct(sum(is.na(x))),
               stringsAsFactors = FALSE)
  )
}

row_categorical <- function(x, label) {
  if (!is.factor(x)) x <- factor(x)
  tab <- table(x, useNA = "no")
  cat_rows <- data.frame(
    Variable = paste0("  ", names(tab)),
    Value    = vapply(as.integer(tab), fmt_n_pct, character(1)),
    stringsAsFactors = FALSE
  )
  bind_rows(
    data.frame(Variable = paste0(label, ", n (%)"),
               Value    = "",
               stringsAsFactors = FALSE),
    cat_rows,
    data.frame(Variable = "  Missing, n (%)",
               Value    = fmt_n_pct(sum(is.na(x))),
               stringsAsFactors = FALSE)
  )
}

## 4.2 Assemble Table 1 ----
table1_descriptive <- bind_rows(
  data.frame(Variable = "N (cohort)",
             Value    = format(N_total, big.mark = ","),
             stringsAsFactors = FALSE),
  
  row_categorical(baseline_t1$sex,       "Sex"),
  row_continuous (baseline_t1$age,       "Age, years"),
  row_categorical(baseline_t1$ethnicity, "Ethnicity"),
  row_categorical(baseline_t1$smoking,   "Smoking status"),
  row_continuous (baseline_t1$bmi,       "BMI, kg/m2"),
  row_continuous (baseline_t1$sbp,       "Systolic BP, mmHg"),
  row_continuous (baseline_t1$dbp,       "Diastolic BP, mmHg"),
  row_continuous (baseline_t1$hba1c,     "HbA1c, %"),
  row_continuous (baseline_t1$chol,      "Total cholesterol, mmol/L"),
  row_continuous (baseline_t1$ldl,       "LDL cholesterol, mmol/L"),
  row_continuous (baseline_t1$hdl,       "HDL cholesterol, mmol/L")
)

print(knitr::kable(table1_descriptive, format = "pipe", align = "lr"))


# 5. PERSON-TIME SPLIT =========================================================
# Aim:    Split each patient's follow-up into one row per calendar year so
#         that every row can be matched to the correct age- and year-specific
#         general-population mortality rate. Person-years in a year =
#         days observed in that year / total days in that year.
# Input:  UKPDS.
# Output: ukpds_long, one row per person-year, carrying calendar year,
#         attained age, person-years (py) and a death flag.

## 5.1 split_person_time() helper ----
# Aim:    Expand one patient row into many person-year rows.
# Input:  a one-row-per-person frame with entry date, end date, date of birth
#         and event indicator; entry_col names the column to start the clock.
# Output: a long frame with year, age, days_observed, py and death per row.
split_person_time <- function(data, entry_col = "enrolment") {
  data %>%
    mutate(
      entry_date = .data[[entry_col]],
      entry_year = year(entry_date),
      end_year   = year(end)
    ) %>%
    rowwise() %>%
    mutate(cal_year = list(seq(entry_year, end_year))) %>%
    ungroup() %>%
    unnest(cal_year) %>%
    mutate(
      year_start    = as.Date(sprintf("%d-01-01", cal_year)),
      year_end      = as.Date(sprintf("%d-12-31", cal_year)),
      t0            = pmax(entry_date, year_start),
      t1            = pmin(end, year_end),
      days_observed = as.numeric(difftime(t1, t0, units = "days")) + 1,
      days_in_year  = as.numeric(difftime(year_end, year_start, units = "days")) + 1,
      py            = days_observed / days_in_year,
      age           = floor(as.numeric(difftime(t0, dob, units = "days")) / 365.25)
    ) %>%
    group_by(id) %>%
    mutate(death = if_else(event == 1 & t1 == max(t1), 1L, 0L)) %>%
    ungroup() %>%
    rename(year = cal_year) %>%
    select(-year_start, -year_end, -entry_date, -entry_year, -end_year)
}

## 5.2 Apply to cohort ----
ukpds_long <- split_person_time(UKPDS, entry_col = "enrolment") %>%
  mutate(enrolment_year = year(enrolment))


# 6. LIFE TABLES (HMD 1x1) =====================================================
# Aim:    Read the UK general-population mortality rates and keep only the
#         years and ages present in the cohort.
# Input:  path_lt_female, path_lt_male.
# Output: life_1x1, with sex / year / single age / mx (central death rate).

## 6.1 read_hmd_1x1() helper ----
# Aim:    Read one HMD 1x1 life table, find the data header automatically and
#         return only the columns needed.
# Input:  path to the .txt file and a sex label ("Female"/"Male").
# Output: a frame with sex, Year, Age and mx.
read_hmd_1x1 <- function(path, sex_label) {
  header_row <- which(str_detect(read_lines(path), "^\\s*Year\\s+Age\\s+mx\\s+"))
  read_table(
    file = path, skip = max(header_row - 1, 0),
    col_types = cols(
      Year = col_integer(), Age = col_character(), mx = col_double(),
      qx = col_double(), ax = col_double(), lx = col_double(),
      dx = col_double(), Lx = col_double(), Tx = col_double(), ex = col_double()
    )
  ) %>%
    mutate(Age = as.integer(str_replace(Age, "110\\+", "110")),
           sex = sex_label) %>%
    select(sex, Year, Age, mx)
}

## 6.2 Load and trim to cohort years/ages ----
life_1x1 <- bind_rows(
  read_hmd_1x1(path_lt_female, "Female"),
  read_hmd_1x1(path_lt_male,   "Male")
) %>%
  filter(between(Year, min(ukpds_long$year, na.rm = TRUE),
                 max(ukpds_long$year, na.rm = TRUE)),
         between(Age,  min(ukpds_long$age,  na.rm = TRUE),
                 max(ukpds_long$age,  na.rm = TRUE))) %>%
  arrange(sex, Year, Age)


# 7. JOIN EXPECTED DEATHS ======================================================
# Aim:    Attach the matching mortality rate to each person-year and compute
#         the expected number of deaths for that cell (expected = mx * py).
# Input:  ukpds_long, life_1x1.
# Output: ukpds_joined, person-year rows with an expected-deaths column.

ukpds_joined <- ukpds_long %>%
  left_join(life_1x1, by = c("sex", "year" = "Year", "age" = "Age")) %>%
  mutate(expected = mx * py)


# 8. OVERALL SMR AND COHORT SUMMARY ============================================
# Aim:    Compute the overall SMR with an exact Poisson 95% CI, and assemble a
#         small summary table (sample size, deaths, person-years, follow-up,
#         overall SMR).
# Input:  ukpds_joined, UKPDS.
# Output: smr_overall (list) and cohort_summary (data frame, used in Sect. 12).

## 8.1 Overall SMR ----
# Exact Poisson 95% CI (chi-square form). Used because expected counts can be
# small in some strata, where a normal-approximation CI is unreliable.
O_all <- sum(ukpds_joined$death,    na.rm = TRUE)
E_all <- sum(ukpds_joined$expected, na.rm = TRUE)

smr_overall <- list(
  smr   = O_all / E_all,
  lower = ifelse(O_all == 0, 0, 0.5 * qchisq(0.025, df = 2 * O_all))     / E_all,
  upper =                       0.5 * qchisq(0.975, df = 2 * (O_all + 1)) / E_all
)

cat(sprintf("Overall SMR = %.2f (95%% CI %.2f-%.2f); O = %d, E = %.1f\n",
            smr_overall$smr, smr_overall$lower, smr_overall$upper, O_all, E_all))

## 8.2 Cohort summary table ----
# Follow-up time per person (administrative, on the time-since-enrolment scale),
# using one row per person from the UKPDS frame.
fu_years <- as.numeric(UKPDS$end - UKPDS$enrolment) / 365.25

cohort_summary <- data.frame(
  Metric = c("N (cohort)",
             "Deaths (observed)",
             "Total person-years",
             "Median follow-up, years (IQR)",
             "Median follow-up, years (range)",
             "Overall SMR (95% CI)"),
  Value  = c(
    format(nrow(UKPDS), big.mark = ","),
    format(O_all, big.mark = ","),
    formatC(sum(ukpds_joined$py, na.rm = TRUE), format = "f",
            big.mark = ",", digits = 1),
    sprintf("%.1f (%.1f-%.1f)",
            median(fu_years),
            quantile(fu_years, 0.25),
            quantile(fu_years, 0.75)),
    sprintf("%.1f (%.1f-%.1f)",
            median(fu_years), min(fu_years), max(fu_years)),
    sprintf("%.2f (%.2f-%.2f)",
            smr_overall$smr, smr_overall$lower, smr_overall$upper)
  ),
  stringsAsFactors = FALSE
)

## 8.3 smr_fmt() helper for SMR table cells ----
# Aim:    Format an SMR as "X.XX (X.XX-X.XX)" with an exact Poisson CI; this is
#         reused for every cell of the subgroup table in Section 11.
# Input:  observed (O) and expected (E) counts.
# Output: a formatted string, or a placeholder when E is missing or zero.
smr_fmt <- function(O, E) {
  if (is.na(E) || E == 0) return("—")
  smr   <- O / E
  lower <- if (O == 0) 0 else 0.5 * qchisq(0.025, df = 2 * O)
  upper <-                  0.5 * qchisq(0.975, df = 2 * (O + 1))
  sprintf("%.2f (%.2f-%.2f)", smr, lower / E, upper / E)
}


# 9. BUILD ANALYSIS FRAME ======================================================
# Aim:    Add the baseline subgroup categories to the person-year data and a
#         "years since enrolment" column. This frame drives Figure 3 and the
#         SMR table.
# Input:  ukpds_joined, ukpds_baseline.
# Output: ukpds_analysis.

ukpds_analysis <- ukpds_joined %>%
  left_join(ukpds_baseline, by = "id") %>%
  mutate(years_since_enrolment = year - enrolment_year)


# 10. FIGURES ==================================================================
# Aim:    Produce the three figures kept for the manuscript.
# Input:  ukpds_joined, ukpds_analysis, smr_overall.
# Output: three .png files in Output/Figure/.

## 10.1 Figure 1: SMR by years since enrolment ----
# Aim:    Show how the SMR changes with time since enrolment, alongside the
#         reference line (SMR = 1), the overall UKPDS SMR and a published
#         meta-analysis estimate. A strip table below shows the underlying
#         population, deaths, sex mix and mean age.
# Output: Output/Figure/SMR_years_since_enrolment.png
smr_by_ysd <- ukpds_joined %>%
  mutate(years_since_enrolment = year - enrolment_year) %>%
  group_by(years_since_enrolment) %>%
  summarise(
    n_pop    = n_distinct(id),
    observed = sum(death,    na.rm = TRUE),
    expected = sum(expected, na.rm = TRUE),
    mean_age = mean(age,     na.rm = TRUE),
    P_male   = round(mean(sex == "Male", na.rm = TRUE) * 100, 1),
    .groups  = "drop"
  ) %>%
  mutate(
    smr   = if_else(expected > 0, observed / expected, NA_real_),
    lower = if_else(observed == 0, 0,
                    0.5 * qchisq(0.025, df = 2 * observed)) / expected,
    upper =          0.5 * qchisq(0.975, df = 2 * (observed + 1)) / expected
  ) %>%
  arrange(years_since_enrolment)

x_vals_ysd   <- sort(unique(smr_by_ysd$years_since_enrolment))
x_limits_ysd <- c(min(x_vals_ysd) - 0.5, 39.5)

p_main_ysd <- ggplot(smr_by_ysd, aes(x = years_since_enrolment, y = smr)) +
  # Meta-analysis pooled estimate ribbon
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.39, ymax = 1.94,
           fill = "#1b9e77", alpha = 0.15) +
  geom_hline(yintercept = 1.64, linetype = "dashed",
             colour = "#1b9e77", linewidth = 0.9) +
  annotate("text", x = min(x_vals_ysd) - 0.2, y = 1.8,
           label = "Meta-analysis SMR: 1.64 (95% CI 1.39–1.94)",
           hjust = 0, vjust = 0, size = 3.5,
           colour = "#1b9e77", fontface = "italic") +
  # SMR series
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_text(aes(label = sprintf("%.2f", smr)),
            vjust = -1.5, size = 3.5, fontface = "bold") +
  # Reference SMR = 1
  geom_hline(yintercept = 1, linewidth = 0.7, colour = "grey30") +
  annotate("text", x = 39, y = 0.7,
           label = "Reference (SMR = 1)",
           hjust = 1, vjust = -0.6, size = 3, colour = "grey30") +
  # Overall UKPDS SMR
  geom_hline(yintercept = smr_overall$smr, linetype = "dotted",
             linewidth = 0.9, colour = "grey50") +
  annotate("text", x = min(x_vals_ysd) + 20, y = smr_overall$smr - 0.25,
           label = sprintf("Overall UKPDS SMR %.2f (95%% CI %.2f–%.2f)",
                           smr_overall$smr, smr_overall$lower, smr_overall$upper),
           hjust = 1, vjust = -0.6, size = 3, colour = "grey50") +
  scale_x_continuous(breaks = x_vals_ysd, limits = x_limits_ysd, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max(4, max(smr_by_ysd$smr, na.rm = TRUE) * 1.1)),
                  clip = "off") +
  labs(x = "Years since enrolment", y = "Standardised Mortality Ratio",
       title = "SMR by Years Since Enrolment with Study Comparisons") +
  theme_minimal(base_size = 13)

# Strip table below the plot
graph_table_ysd <- smr_by_ysd %>%
  transmute(years_since_enrolment,
            Population        = n_pop,
            `Observed deaths` = observed,
            `Male %`          = P_male,
            `Mean age`        = round(mean_age, 1)) %>%
  pivot_longer(-years_since_enrolment, names_to = "Variable", values_to = "Value")

p_table_ysd <- ggplot(graph_table_ysd,
                      aes(x = years_since_enrolment, y = Variable, label = Value)) +
  geom_text(size = 3.6, hjust = 0.5) +
  scale_x_continuous(breaks = x_vals_ysd, limits = x_limits_ysd, expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.y  = element_text(face = "bold"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 2, r = 6, b = 2, l = 6))

ggsave("Output/Figure/SMR_years_since_enrolment.png",
       p_main_ysd / p_table_ysd + plot_layout(heights = c(3, 1)),
       width = 20, height = 6, dpi = 300)


## 10.2 Figure 2: SMR by calendar year ----
# Aim:    Show how the SMR changes across calendar years, with study milestones
#         (end of recruitment, intervention and monitoring) marked. A strip
#         table below shows new recruits, population, deaths, sex mix and age.
# Output: Output/Figure/SMR_calendar_year.png
smr_by_cal <- ukpds_joined %>%
  group_by(year) %>%
  summarise(
    n_pop    = n_distinct(id),
    observed = sum(death,    na.rm = TRUE),
    expected = sum(expected, na.rm = TRUE),
    mean_age = mean(age,     na.rm = TRUE),
    P_male   = round(mean(sex == "Male", na.rm = TRUE) * 100, 1),
    .groups  = "drop"
  ) %>%
  mutate(
    smr   = if_else(expected > 0, observed / expected, NA_real_),
    lower = if_else(observed == 0, 0,
                    0.5 * qchisq(0.025, df = 2 * observed)) / expected,
    upper =          0.5 * qchisq(0.975, df = 2 * (observed + 1)) / expected
  ) %>%
  arrange(year)

# Patients newly enrolled in each calendar year (used only for the strip table).
new_recruits_by_year <- ukpds_joined %>%
  distinct(id, enrolment_year) %>%
  dplyr::count(enrolment_year, name = "n_new") %>%
  rename(year = enrolment_year)

smr_by_cal <- smr_by_cal %>%
  left_join(new_recruits_by_year, by = "year") %>%
  mutate(n_new = coalesce(n_new, 0L))

x_vals_cal   <- sort(unique(smr_by_cal$year))
x_limits_cal <- c(min(x_vals_cal) - 0.5, 2021 + 0.5)

p_main_cal <- ggplot(smr_by_cal, aes(x = year, y = smr)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.39, ymax = 1.94,
           fill = "#1b9e77", alpha = 0.15) +
  geom_hline(yintercept = 1.64, linetype = "dashed",
             colour = "#1b9e77", linewidth = 0.9) +
  annotate("text", x = min(x_vals_cal) + 1, y = 1.8,
           label = "Meta-analysis SMR: 1.64 (95% CI 1.39–1.94)",
           hjust = 0, vjust = 0, size = 3.5,
           colour = "#1b9e77", fontface = "italic") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_text(aes(label = sprintf("%.2f", smr)),
            vjust = -1.5, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 1, linewidth = 0.8, colour = "grey30") +
  annotate("text", x = 2021 - 0.2, y = 0.6,
           label = "Reference (SMR = 1)",
           hjust = 1, vjust = -0.6, size = 3.2, colour = "grey30") +
  geom_hline(yintercept = smr_overall$smr, linetype = "dotted",
             linewidth = 0.7, colour = "grey50") +
  annotate("text", x = 2006, y = smr_overall$smr - 0.4,
           label = sprintf("Overall SMR %.2f (95%% CI %.2f–%.2f)",
                           smr_overall$smr, smr_overall$lower, smr_overall$upper),
           hjust = 1, vjust = -0.6, size = 3.4, colour = "grey50") +
  # Study milestones
  geom_vline(xintercept = 1991, linetype = "dotted",
             colour = "#444444", linewidth = 0.7) +
  geom_vline(xintercept = 1997, linetype = "dotted",
             colour = "#444444", linewidth = 0.7) +
  geom_vline(xintercept = 2007, linetype = "dotted",
             colour = "#444444", linewidth = 0.7) +
  annotate("text", x = 1991, y = 4.8,
           label = "End of recruitment (1991)",
           angle = 90, hjust = 1, vjust = -0.4, size = 3,
           colour = "#444444") +
  annotate("text", x = 1997, y = 4.8,
           label = "End of intervention (1997)",
           angle = 90, hjust = 1, vjust = -0.4, size = 3,
           colour = "#444444") +
  annotate("text", x = 2007, y = 4.8,
           label = "End of monitoring (2007)",
           angle = 90, hjust = 1, vjust = -0.4, size = 3,
           colour = "#444444") +
  scale_x_continuous(breaks = x_vals_cal, limits = x_limits_cal, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 5), clip = "off") +
  labs(x = "Calendar year", y = "Standardised Mortality Ratio",
       title = "SMR by Calendar Year with Study Comparisons") +
  theme_minimal(base_size = 13)

graph_table_cal <- smr_by_cal %>%
  transmute(year,
            `New recruits`    = n_new,
            Population        = n_pop,
            `Observed deaths` = observed,
            `Male %`          = P_male,
            `Mean age`        = round(mean_age, 1)) %>%
  pivot_longer(-year, names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable,
                           levels = c("New recruits", "Population",
                                      "Observed deaths", "Male %", "Mean age")))

p_table_cal <- ggplot(graph_table_cal, aes(x = year, y = Variable, label = Value)) +
  geom_text(size = 3.6, hjust = 0.5) +
  scale_x_continuous(breaks = x_vals_cal, limits = x_limits_cal, expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.y  = element_text(face = "bold"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 2, r = 6, b = 2, l = 6))

ggsave("Output/Figure/SMR_calendar_year.png",
       p_main_cal / p_table_cal + plot_layout(heights = c(3, 1)),
       width = 20, height = 6, dpi = 300)


## 10.3 Figure 3: SMR by subgroup, 3-year bands (6-panel grid) ----
# Aim:    Show SMR across 3-year follow-up bands within each baseline subgroup
#         (sex, age, smoking, HbA1c, BMI, ethnicity), arranged in one 3x2
#         figure with shared axis labels.
# Output: Output/Figure/SMR_subgroups_combined_3yr.png

# make_subgroup_panel_3yr()
# Aim:    Build one subgroup panel showing SMR across 3-year follow-up bands
#         for one variable: a ribbon for the 95% CI, a dashed line and bold
#         label for that subgroup's overall SMR, and a solid grey reference
#         line at SMR = 1.
# Input:  group_col (column name in ukpds_analysis), a panel title, and an
#         optional y-axis cap (values above the cap are clipped for display).
# Output: a ggplot object.
make_subgroup_panel_3yr <- function(group_col, panel_title, y_max = 5) {
  band_breaks <- c(-Inf, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, Inf)
  band_labels <- c("0-3y", "3-6y", "6-9y", "9-12y",
                   "12-15y", "15-18y", "18-21y", "21-24y",
                   "24-27y", "27-30y", "30-33y", "33-36y",
                   "36-39y", "39y+")
  
  d <- ukpds_analysis %>%
    rename(grp = all_of(group_col)) %>%
    filter(!is.na(grp)) %>%
    mutate(fu_band = cut(years_since_enrolment, breaks = band_breaks,
                         labels = band_labels, right = FALSE))
  
  smr_t <- d %>%
    group_by(fu_band, grp) %>%
    summarise(observed = sum(death,    na.rm = TRUE),
              expected = sum(expected, na.rm = TRUE),
              .groups  = "drop") %>%
    filter(!is.na(fu_band)) %>%
    mutate(
      smr     = if_else(expected > 0, observed / expected, NA_real_),
      lower   = if_else(observed == 0, 0,
                        0.5 * qchisq(0.025, df = 2 * observed)) / expected,
      upper   =          0.5 * qchisq(0.975, df = 2 * (observed + 1)) / expected,
      smr_p   = pmin(smr,            y_max),
      lower_p = pmin(pmax(lower, 0), y_max),
      upper_p = pmin(upper,          y_max)
    )
  
  smr_o <- d %>%
    group_by(grp) %>%
    summarise(observed = sum(death,    na.rm = TRUE),
              expected = sum(expected, na.rm = TRUE),
              .groups  = "drop") %>%
    mutate(
      smr   = observed / expected,
      lower = if_else(observed == 0, 0,
                      0.5 * qchisq(0.025, df = 2 * observed)) / expected,
      upper =          0.5 * qchisq(0.975, df = 2 * (observed + 1)) / expected,
      label = sprintf("%.2f (%.2f-%.2f)", smr, lower, upper),
      smr_p = pmin(smr, y_max)
    )
  
  ggplot(smr_t, aes(x = fu_band, y = smr_p,
                    colour = grp, fill = grp, group = grp)) +
    geom_ribbon(aes(ymin = lower_p, ymax = upper_p), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(data = smr_o, aes(yintercept = smr_p, colour = grp),
               linetype = "dashed", alpha = 0.7) +
    geom_text(data = smr_o,
              aes(x = 1, y = smr_p, label = label, colour = grp),
              hjust = 0, vjust = -0.5, size = 2.5,   # value labels ~7 pt
              fontface = "bold", show.legend = FALSE) +
    geom_hline(yintercept = 1, colour = "grey30", linewidth = 0.5) +
    # Shorten tick labels for display only (unit is in the shared x title);
    # does not touch the SMR tables, which keep their own "0-3y" labels.
    scale_x_discrete(labels = function(b) gsub("y", "", b)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    coord_cartesian(clip = "on") +
    labs(title = panel_title, x = NULL, y = NULL, colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          legend.margin   = margin(t = -4, b = 0),
          plot.title      = element_text(size = 11, face = "bold"),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y     = element_text(size = 8),
          plot.margin     = margin(t = 4, r = 6, b = 2, l = 4))
}

# Shared y / x / figure titles around the 2x3 grid (unchanged approach)
add_shared_labels <- function(grid_plot, fig_title,
                              left = 0.05, right = 0.995, bottom = 0.04, top = 0.955) {
  wrap_elements(full = grid::grobTree(
    grid::textGrob("Standardised Mortality Ratio", rot = 90,
                   x = unit(0.012, "npc"), y = unit(0.5, "npc"),
                   gp = gpar(fontsize = 12)),
    grid::textGrob("Years since enrolment (3-year bands)",
                   x = unit(0.52, "npc"), y = unit(0.012, "npc"),
                   gp = gpar(fontsize = 12)),
    grid::textGrob(fig_title,
                   x = unit(0.52, "npc"), y = unit(0.99, "npc"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )) +
    inset_element(grid_plot, left = left, right = right, bottom = bottom, top = top)
}

panels <- list(
  make_subgroup_panel_3yr("sex",           "A  Sex"),
  make_subgroup_panel_3yr("age_cat",       "B  Baseline Age"),
  make_subgroup_panel_3yr("smoking_cat",   "C  Smoking Status"),
  make_subgroup_panel_3yr("hba1c_cat",     "D  Baseline HbA1c"),
  make_subgroup_panel_3yr("bmi_cat",       "E  Baseline BMI"),
  make_subgroup_panel_3yr("ethnicity_cat", "F  Ethnicity")
)

p_all <- add_shared_labels(
  wrap_plots(panels, ncol = 2),
  "Standardised Mortality Ratios by Subgroup in UKPDS (3-year bands)"
)

ggsave("Output/Figure/SMR_subgroups_3yr.png", p_all,
       width = 8.5, height = 11.2, units = "in", dpi = 300)


# 11. SMR TABLE BY BASELINE CHARACTERISTICS x FOLLOW-UP PERIOD =================
# Aim:    For the whole cohort and each subgroup level, report N, observed
#         deaths, the overall SMR, and the SMR within each follow-up band.
#         Two banding schemes are built: 3-year bands and single-year bands.
# Input:  ukpds_analysis, smr_fmt().
# Output: table_smr_3yr and table_smr_1yr (written as sheets in Section 12).

## 11.1 Generic table builder ----
# Aim:    Build the subgroup SMR table for any set of follow-up bands, so the
#         same code produces both the 3-year and the single-year versions.
# Input:  data (ukpds_analysis), the cut breaks and the matching band labels.
# Output: a data frame with Subgroup, N, Observed Deaths, overall SMR, and one
#         SMR (95% CI) column per band.
build_smr_subgroup_table <- function(data, period_breaks, period_labels) {
  
  # Bin follow-up time into the requested bands.
  d_all <- data %>%
    mutate(fu_period = cut(years_since_enrolment,
                           breaks = period_breaks,
                           labels = period_labels,
                           right  = FALSE))
  
  # One table row: subgroup label + overall + one cell per band.
  make_row <- function(d, label) {
    O_t <- sum(d$death,    na.rm = TRUE)
    E_t <- sum(d$expected, na.rm = TRUE)
    N_t <- n_distinct(d$id)
    row <- data.frame(
      Subgroup        = label,
      N               = format(N_t, big.mark = ","),
      Observed_Deaths = format(O_t, big.mark = ","),
      SMR_overall     = smr_fmt(O_t, E_t),
      stringsAsFactors = FALSE
    )
    for (p in period_labels) {
      dp <- d %>% filter(fu_period == p)
      row[[paste0("SMR_", p)]] <- smr_fmt(sum(dp$death,    na.rm = TRUE),
                                          sum(dp$expected, na.rm = TRUE))
    }
    row
  }
  
  # A group-name header row (all band cells left empty).
  make_header <- function(label) {
    row <- data.frame(Subgroup = label, N = "", Observed_Deaths = "",
                      SMR_overall = "", stringsAsFactors = FALSE)
    for (p in period_labels) row[[paste0("SMR_", p)]] <- ""
    row
  }
  
  # Append a subgroup block: header + one row per level + a Missing row
  # (only when that variable has any missing baseline values).
  add_subgroup_block <- function(rows, var_name, header_label) {
    rows[[length(rows) + 1]] <- make_header(header_label)
    
    lvls <- levels(d_all[[var_name]])
    if (is.null(lvls)) lvls <- sort(unique(stats::na.omit(d_all[[var_name]])))
    for (lvl in lvls) {
      rows[[length(rows) + 1]] <- make_row(
        d_all %>% filter(.data[[var_name]] == lvl),
        paste0("  ", lvl))
    }
    
    n_missing_persons <- d_all %>%
      filter(is.na(.data[[var_name]])) %>% pull(id) %>% n_distinct()
    if (n_missing_persons > 0) {
      rows[[length(rows) + 1]] <- make_row(
        d_all %>% filter(is.na(.data[[var_name]])),
        "  Missing")
    }
    rows
  }
  
  # Assemble: an overall row, then one block per subgroup variable.
  rows <- list()
  rows[[length(rows) + 1]] <- make_row(d_all, "Overall")
  rows <- add_subgroup_block(rows, "sex",           "Sex")
  rows <- add_subgroup_block(rows, "age_cat",       "Baseline age")
  rows <- add_subgroup_block(rows, "smoking_cat",   "Smoking status")
  rows <- add_subgroup_block(rows, "hba1c_cat",     "Baseline HbA1c")
  rows <- add_subgroup_block(rows, "bmi_cat",       "Baseline BMI (kg/m2)")
  rows <- add_subgroup_block(rows, "ethnicity_cat", "Ethnicity")
  
  out <- bind_rows(rows)
  names(out) <- c("Subgroup", "N", "Observed Deaths", "SMR (95% CI)",
                  paste0("SMR ", period_labels))
  out
}

## 11.2 Three-year-band table ----
breaks_3yr <- c(-Inf, seq(3, 39, by = 3), Inf)
labels_3yr <- c(sprintf("%d-%dy", seq(0, 36, by = 3), seq(3, 39, by = 3)), "39y+")

table_smr_3yr <- build_smr_subgroup_table(ukpds_analysis, breaks_3yr, labels_3yr)

cat("\n========================================================\n")
cat("SMR by Baseline Characteristics: 3-year bands (since enrolment)\n")
cat("========================================================\n\n")
print(knitr::kable(table_smr_3yr, format = "pipe",
                   align = c("l", rep("r", ncol(table_smr_3yr) - 1))))

## 11.3 Single-year-band table ----
# Same layout as the 3-year table, but each follow-up year is its own column:
# "0-1y", "1-2y", ..., "38-39y", plus a final "39y+" catch-all.
breaks_1yr <- c(-Inf, 1:39, Inf)
labels_1yr <- c(sprintf("%d-%dy", 0:38, 1:39), "39y+")

table_smr_1yr <- build_smr_subgroup_table(ukpds_analysis, breaks_1yr, labels_1yr)

cat(sprintf("\nBuilt single-year SMR table: %d rows x %d columns ",
            nrow(table_smr_1yr), ncol(table_smr_1yr)))
cat("(too wide to print here; see the Excel workbook in Section 12).\n")


# 12. WRITE RESULTS TO A SINGLE EXCEL WORKBOOK =================================
# Aim:    Save the result tables as separate sheets in one .xlsx file.
# Input:  table1_descriptive (Section 4), cohort_summary (Section 8),
#         table_smr_3yr and table_smr_1yr (Section 11).
# Output: Output/UKPDS_SMR_tables.xlsx with one sheet per table.

wb <- createWorkbook()

# A reusable header style: bold text on a light grey background.
header_style <- createStyle(textDecoration = "bold",
                            fgFill = "#F2F2F2", border = "bottom")

# add_sheet(): add one sheet, write a data frame, style the header, auto-fit
#              column widths, and freeze the header row and first column.
add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = header_style)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(data)), widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
}

add_sheet(wb, "Baseline characteristics", table1_descriptive)
add_sheet(wb, "Cohort summary",           cohort_summary)
add_sheet(wb, "SMR by subgroup (3-year)", table_smr_3yr)
add_sheet(wb, "SMR by subgroup (1-year)", table_smr_1yr)

saveWorkbook(wb, "Output/UKPDS_SMR_tables_5102.xlsx", overwrite = TRUE)

cat("\nDone.\n")
cat("Tables : Output/UKPDS_SMR_tables_5102.xlsx (4 sheets)\n")
cat("Figures: Output/Figure/SMR_years_since_enrolment.png\n")
cat("         Output/Figure/SMR_calendar_year.png\n")
cat("         Output/Figure/SMR_subgroups_combined_3yr.png\n")