# UKPDS Standardised Mortality Ratio (SMR) analysis - randomised cohort (n = 4,209) ====
#
# Author : Xinyu Li (xinyu.li@ndph.ox.ac.uk)
# Updated: 2026-06-20
#
# WHAT THIS SCRIPT DOES
#   For the randomised UKPDS patients (those with a glycaemic randomisation
#   date, GL_date), it compares the number of deaths observed with the number
#   expected if the cohort had died at UK general-population rates (HMD 1x1
#   life tables):
#                    SMR = observed deaths / expected deaths
#
#   SMR is reported on two time scales and several breakdowns:
#       - time since enrolment      (clock starts at the enrolment date)
#       - time since randomisation  (clock starts at GL_date, i.e. post run-in)
#     within these: overall, by calendar year, by baseline subgroup,
#     and by treatment arm.
#
# OUTPUTS
#   Figures  -> Output/Figure/
#     SMR_enrolment_by_followup_years_4209.png       Overall SMR vs years since enrolment
#     SMR_enrolment_by_calendar_year_4209.png        Overall SMR vs calendar year
#     SMR_enrolment_subgroups_3yr_4209.png           Subgroup SMR, 3-year bands (enrolment)
#     SMR_randomization_by_followup_years_4209.png   Overall SMR vs years since randomisation
#     SMR_randomization_subgroups_3yr_4209.png       Subgroup SMR, 3-year bands (randomisation)
#     SMR_randomization_arm_mainstudy_3yr_4209.png   Treatment-arm SMR, UKPDS 33 (randomisation)
#     SMR_randomization_arm_substudy_3yr_4209.png    Treatment-arm SMR, UKPDS 34 (randomisation)
#
#   Table    -> Output/
#     UKPDS_SMR_tables_4209.xlsx
#         tab "Enrolment"      : SMR by baseline characteristics, time since enrolment
#         tab "Randomization"  : SMR by baseline characteristics + treatment arm,
#                                time since randomisation


# 1. SETUP ====================================================================
#    Aim   : load packages, set the working directory, list input files, and
#            make sure the output folders exist.
#    Input : none.
#    Output: an R session ready to run the analysis.

setwd("K:/HERC_UKPDS/Xinyu Li Risk equations/Xinyu UKPDS SMR project")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(readr)
  library(haven)       # read Stata .dta files
  library(stringr)
  library(patchwork)   # stack the SMR curve over its data strip-table
  library(grid)        # shared axis labels on the combined subgroup figures
  library(openxlsx)    # write the two-tab Excel table
})

# Input files
path_extended_stata <- "K:/QNAP/UKPDS/DTU data/XL Death Data.dta"             # extra death/censoring dates
path_stata          <- "K:/HERC_UKPDS/UKPDS-OM2 model/Data/locf_anal1Aug11.dta"           # main longitudinal dataset
path_lt_female      <- "Data/Life tables/UK fltper_1x1.txt"
path_lt_male        <- "Data/Life tables/UK mltper_1x1.txt"

# Output folders (created if they do not already exist)
dir.create("Output",        showWarnings = FALSE, recursive = TRUE)
dir.create("Output/Figure", showWarnings = FALSE, recursive = TRUE)


# 2. LOAD AND ASSEMBLE THE ANALYSIS COHORT ====================================
#    Aim   : merge the two source datasets and reduce them to one row per
#            randomised patient, with the dates the analysis needs.
#    Input : the two .dta files above.
#    Output: UKPDS - one row per person (n = 4,209) with id, sex, dob,
#            enrolment date, follow-up end date, death indicator, and
#            randomisation date.

## 2.1 Main longitudinal dataset (many rows per person) ----
locf <- read_dta(path_stata)

## 2.2 Add death dates, censoring dates, randomisation dates and arm flags ----
extended_UKPDS <- read_dta(path_extended_stata)

locf <- merge(locf,
              extended_UKPDS[, c("ukpdsno", "anydeath_event", "anydeath_dt",
                                 "cens_dt2", "GL_date", "BP_date", "GL_study",
                                 "BP_study", "metf", "metdate", "SFUIns",
                                 "SFUinsdate")],
              by = "ukpdsno")

## 2.3 Collapse to one row per randomised patient ----
#   end   = death date if the patient died, otherwise the censoring date.
#   Keeping only patients with a randomisation date (GL_date) gives n = 4,209.
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
  filter(!is.na(dob), !is.na(enrolment), !is.na(end), end >= enrolment) %>%
  filter(!is.na(randomization))   # 4,209 randomised patients


# 3. PERSON-TIME SPLITTER =====================================================
#    Aim   : cut each person's follow-up into calendar-year slices so that
#            person-years and attained age can be matched to single-year life
#            tables. The same function serves both time scales; it is called
#            once with the enrolment date and once with the randomisation date.
#    Input : a one-row-per-person table (UKPDS) and the name of the entry-date
#            column to use.
#    Output: a long table with one row per person per calendar year, carrying
#            person-years (py), attained age, and a death flag for that slice.

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

# Enrolment time scale: split person-time from the enrolment date.
ukpds_long <- split_person_time(UKPDS, entry_col = "enrolment") %>%
  mutate(enrolment_year = year(enrolment))


# 4. LIFE TABLES (EXPECTED MORTALITY) =========================================
#    Aim   : read the HMD single-year-of-age, single-calendar-year UK mortality
#            rates and keep the age/year range present in the cohort.
#    Input : the male and female HMD 1x1 text files.
#    Output: life_1x1 - sex / Year / Age / mx (central death rate), used to turn
#            person-years into expected deaths.

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

life_1x1 <- bind_rows(
  read_hmd_1x1(path_lt_female, "Female"),
  read_hmd_1x1(path_lt_male,   "Male")
) %>%
  filter(between(Year, min(ukpds_long$year, na.rm = TRUE),
                 max(ukpds_long$year, na.rm = TRUE)),
         between(Age,  min(ukpds_long$age,  na.rm = TRUE),
                 max(ukpds_long$age,  na.rm = TRUE))) %>%
  arrange(sex, Year, Age)


# 5. HELPER FUNCTIONS =========================================================
#    Aim   : reusable building blocks for the figures, the Table 1 rows, and
#            the Excel writer. Defined once, used by Sections 7-10.

## 5.1 SMR cell formatter (for the Table 1 cells) ----
#   Input : observed (O) and expected (E) deaths.
#   Output: a string "SMR (lower-upper)" using exact Poisson 95% limits, or a
#           dash when there is no exposure.
smr_fmt <- function(O, E) {
  if (is.na(E) || E == 0) return("—")
  smr   <- O / E
  lower <- if (O == 0) 0 else 0.5 * qchisq(0.025, df = 2 * O)
  upper <-                  0.5 * qchisq(0.975, df = 2 * (O + 1))
  sprintf("%.2f (%.2f-%.2f)", smr, lower / E, upper / E)
}

## 5.2 Subgroup SMR panel, 3-year bands ----
#   Input : person-time data with death/expected, a grouping column, and a
#           3-year follow-up band column (fu_period); a panel title.
#   Output: one ggplot panel showing the SMR trajectory per subgroup level,
#           with each level's overall SMR drawn as a dashed reference line.
make_subgroup_panel_3yr <- function(data, group_col, period_col, panel_title,
                                    y_max = 5) {
  d <- data %>%
    rename(grp   = all_of(group_col),
           xband = all_of(period_col)) %>%
    filter(!is.na(grp), !is.na(xband))
  
  smr_t <- d %>%
    group_by(xband, grp) %>%
    summarise(observed = sum(death,    na.rm = TRUE),
              expected = sum(expected, na.rm = TRUE),
              .groups  = "drop") %>%
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
  
  ggplot(smr_t, aes(x = xband, y = smr_p,
                    colour = grp, fill = grp, group = grp)) +
    geom_ribbon(aes(ymin = lower_p, ymax = upper_p), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(data = smr_o, aes(yintercept = smr_p, colour = grp),
               linetype = "dashed", alpha = 0.7) +
    geom_text(data = smr_o,
              aes(x = 1, y = smr_p, label = label, colour = grp),
              hjust = 0, vjust = -0.7, size = 2.5,   # value labels ~7 pt
              fontface = "bold", show.legend = FALSE) +
    geom_hline(yintercept = 1, colour = "grey30", linewidth = 0.5) +
    scale_x_discrete(labels = function(b) gsub("y", "", b)) +   # display only
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

add_shared_labels <- function(grid_plot, fig_title, x_lab,
                              left = 0.05, right = 0.995, bottom = 0.04, top = 0.955) {
  wrap_elements(full = grid::grobTree(
    grid::textGrob("Standardised Mortality Ratio", rot = 90,
                   x = unit(0.012, "npc"), y = unit(0.5, "npc"),
                   gp = gpar(fontsize = 12)),
    grid::textGrob(x_lab,
                   x = unit(0.52, "npc"), y = unit(0.012, "npc"),
                   gp = gpar(fontsize = 12)),
    grid::textGrob(fig_title,
                   x = unit(0.52, "npc"), y = unit(0.99, "npc"),
                   gp = gpar(fontsize = 14, fontface = "bold"))
  )) +
    inset_element(grid_plot, left = left, right = right, bottom = bottom, top = top)
}

## 5.3 Treatment-arm SMR panel, 3-year bands ----
#   Input : post-randomisation person-time with death/expected, a 2-level arm
#           factor, and the 3-year band column (fu_period); a panel title.
#   Output: a combined figure - the SMR curve per arm on top, with a strip
#           table of population and observed deaths per band underneath.
make_arm_panel_3yr <- function(d, panel_title,
                               x_label = "Years since randomization",
                               y_max   = 5) {
  smr_t <- d %>%
    group_by(fu_period, arm) %>%
    summarise(n_pop    = n_distinct(id),
              observed = sum(death,    na.rm = TRUE),
              expected = sum(expected, na.rm = TRUE),
              .groups  = "drop") %>%
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
    group_by(arm) %>%
    summarise(observed = sum(death,    na.rm = TRUE),
              expected = sum(expected, na.rm = TRUE),
              .groups  = "drop") %>%
    mutate(
      smr   = observed / expected,
      lower = if_else(observed == 0, 0,
                      0.5 * qchisq(0.025, df = 2 * observed)) / expected,
      upper =          0.5 * qchisq(0.975, df = 2 * (observed + 1)) / expected,
      label = sprintf("%.2f (%.2f–%.2f)", smr, lower, upper),
      smr_p = pmin(smr, y_max)
    )
  
  arm_levels <- levels(d$arm)
  arm_colors <- setNames(c("#E41A1C", "#377EB8"), arm_levels)
  
  p_main <- ggplot(smr_t,
                   aes(x = fu_period, y = smr_p,
                       colour = arm, fill = arm, group = arm)) +
    geom_ribbon(aes(ymin = lower_p, ymax = upper_p),
                alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    geom_hline(data = smr_o, aes(yintercept = smr_p, colour = arm),
               linetype = "dashed", alpha = 0.7) +
    geom_text(data = smr_o,
              aes(x = 1, y = smr_p, label = label, colour = arm),
              hjust = 0.5, vjust = -0.7, size = 3.0,
              fontface = "bold", show.legend = FALSE) +
    # Per-band SMR value above each point
    geom_text(aes(label = sprintf("%.2f", smr)),
              position = position_dodge(width = 0),
              vjust = -1.2, size = 3.0, fontface = "bold",
              show.legend = FALSE) +
    geom_hline(yintercept = 1, colour = "grey30", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    scale_colour_manual(values = arm_colors) +
    scale_fill_manual(values   = arm_colors) +
    coord_cartesian(clip = "off") +
    labs(title = panel_title, x = x_label, y = "Standardised Mortality Ratio",
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title      = element_text(size = 12, face = "bold"),
          plot.margin     = margin(t = 5, r = 6, b = 2, l = 6))
  
  graph_table <- smr_t %>%
    transmute(fu_period, arm,
              Population        = n_pop,
              `Observed deaths` = observed) %>%
    pivot_longer(c(Population, `Observed deaths`),
                 names_to = "Variable", values_to = "Value") %>%
    mutate(row_label = paste0(arm, " — ", Variable))
  
  row_order <- c(
    paste0(arm_levels[1], " — Population"),
    paste0(arm_levels[1], " — Observed deaths"),
    paste0(arm_levels[2], " — Population"),
    paste0(arm_levels[2], " — Observed deaths")
  )
  graph_table <- graph_table %>%
    mutate(row_label = factor(row_label, levels = row_order))
  
  p_table <- ggplot(graph_table,
                    aes(x = fu_period, y = row_label,
                        label = Value, colour = arm)) +
    geom_text(size = 3.2, hjust = 0.5, show.legend = FALSE) +
    scale_y_discrete(limits = rev) +
    scale_colour_manual(values = arm_colors) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid   = element_blank(),
          axis.text.y  = element_text(face = "bold"),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin  = margin(t = 2, r = 6, b = 2, l = 6))
  
  p_main / p_table + plot_layout(heights = c(3, 1.3))
}

## 5.4 Table 1 row builders ----
#   periods             : the 3-year follow-up bands used as table columns.
#   make_row()          : one subgroup row - N, observed deaths, overall SMR,
#                         then an SMR cell for each 3-year band.
#   make_header()       : a blank labelled row that introduces a block.
#   add_baseline_block(): a header plus one row per level of a baseline
#                         variable, adding a "Missing" row when any are NA.
#   add_arm_block()     : a header (with the subset N) plus one row per arm;
#                         no Missing row, because arms are subsets not partitions.
#   add_fu_period()     : adds the fu_period band column to a person-time table.
periods <- c("0-3y", "3-6y", "6-9y", "9-12y",
             "12-15y", "15-18y", "18-21y", "21-24y",
             "24-27y", "27-30y", "30-33y", "33-36y",
             "36-39y", "39y+")

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
  for (p in periods) {
    dp <- d %>% filter(fu_period == p)
    row[[paste0("SMR_", p)]] <- smr_fmt(sum(dp$death,    na.rm = TRUE),
                                        sum(dp$expected, na.rm = TRUE))
  }
  row
}

make_header <- function(label) {
  row <- data.frame(Subgroup = label, N = "", Observed_Deaths = "",
                    SMR_overall = "", stringsAsFactors = FALSE)
  for (p in periods) row[[paste0("SMR_", p)]] <- ""
  row
}

add_baseline_block <- function(rows, data, var_name, header_label) {
  rows[[length(rows) + 1]] <- make_header(header_label)
  lvls <- levels(data[[var_name]])
  if (is.null(lvls)) lvls <- sort(unique(stats::na.omit(data[[var_name]])))
  for (lvl in lvls) {
    rows[[length(rows) + 1]] <- make_row(
      data %>% filter(.data[[var_name]] == lvl),
      paste0("  ", lvl))
  }
  n_missing_persons <- data %>%
    filter(is.na(.data[[var_name]])) %>% pull(id) %>% n_distinct()
  if (n_missing_persons > 0) {
    rows[[length(rows) + 1]] <- make_row(
      data %>% filter(is.na(.data[[var_name]])),
      "  Missing")
  }
  rows
}

add_arm_block <- function(rows, arm_data, header_label) {
  N_total <- n_distinct(arm_data$id)
  hdr <- sprintf("%s (n = %s)", header_label, format(N_total, big.mark = ","))
  rows[[length(rows) + 1]] <- make_header(hdr)
  for (lvl in levels(arm_data$arm)) {
    rows[[length(rows) + 1]] <- make_row(
      arm_data %>% filter(arm == lvl),
      paste0("  ", lvl))
  }
  rows
}

add_fu_period <- function(d, time_col) {
  d %>% mutate(fu_period = cut(.data[[time_col]],
                               breaks = c(-Inf, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, Inf),
                               labels = c("0-3y", "3-6y", "6-9y", "9-12y",
                                          "12-15y", "15-18y", "18-21y", "21-24y",
                                          "24-27y", "27-30y", "30-33y", "33-36y",
                                          "36-39y", "39y+"),
                               right = FALSE))
}

## 5.5 Two-tab Excel writer ----
#   Input : an output path and two named data frames.
#   Output: one .xlsx workbook with each data frame on its own sheet
#           (bold header row, auto column widths, frozen header + first column).
write_two_tab_excel <- function(path,
                                tab1_name, tab1_data,
                                tab2_name, tab2_data) {
  wb <- createWorkbook()
  header_style <- createStyle(textDecoration = "bold",
                              halign = "center", valign = "center",
                              fgFill = "#D9E1F2",
                              border = "TopBottom", borderColour = "#9DB0CE")
  
  add_tab <- function(sheet, data) {
    addWorksheet(wb, sheet)
    writeData(wb, sheet, data, headerStyle = header_style)
    setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
    freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 2)
  }
  
  add_tab(tab1_name, tab1_data)
  add_tab(tab2_name, tab2_data)
  saveWorkbook(wb, path, overwrite = TRUE)
}


# 6. BASELINE CHARACTERISTICS (ONE ROW PER PERSON) ============================
#    Aim   : take each person's first record and bin the baseline variables
#            into the categories used for subgroup figures and Table 1.
#    Input : the merged longitudinal data (locf).
#    Output: ukpds_baseline - id plus categorical age, smoking, HbA1c, BMI
#            and ethnicity. Joined onto the person-time tables by id.

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


# 7. ENROLMENT TIME SCALE  (clock starts at enrolment) ========================
#    Aim   : compute expected deaths and the overall SMR, then produce three
#            figures and build Table 1A (written to Excel in Section 10).
#    Input : ukpds_long, life_1x1, ukpds_baseline.
#    Output: 3 figures (Output/Figure/) and the table1_A object held in memory.

## 7.1 Attach expected deaths (person-years x population mortality rate) ----
ukpds_joined <- ukpds_long %>%
  left_join(life_1x1, by = c("sex", "year" = "Year", "age" = "Age")) %>%
  mutate(expected = mx * py)

## 7.2 Overall SMR (enrolment time scale) ----
O_all <- sum(ukpds_joined$death,    na.rm = TRUE)
E_all <- sum(ukpds_joined$expected, na.rm = TRUE)

smr_overall <- list(
  smr   = O_all / E_all,
  lower = ifelse(O_all == 0, 0, 0.5 * qchisq(0.025, df = 2 * O_all))      / E_all,
  upper =                       0.5 * qchisq(0.975, df = 2 * (O_all + 1)) / E_all
)

cat(sprintf("Enrolment time scale (n = 4,209): Overall SMR = %.2f (%.2f-%.2f)\n",
            smr_overall$smr, smr_overall$lower, smr_overall$upper))

## 7.3 Figure - overall SMR by years since enrolment ----
#   File: SMR_enrolment_by_followup_years_4209.png
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
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.39, ymax = 1.94,
           fill = "#1b9e77", alpha = 0.15) +
  geom_hline(yintercept = 1.64, linetype = "dashed",
             colour = "#1b9e77", linewidth = 0.9) +
  annotate("text", x = min(x_vals_ysd) + 1, y = 1.8,
           label = "Meta-analysis SMR: 1.64 (95% CI 1.39–1.94)",
           hjust = 0, vjust = 0, size = 3.5,
           colour = "#1b9e77", fontface = "italic") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_text(aes(label = sprintf("%.2f", smr)),
            vjust = -1.5, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 1, linewidth = 0.7, colour = "grey30") +
  annotate("text", x = max(x_vals_ysd) - 0.5, y = 0.7,
           label = "Reference (SMR = 1)",
           hjust = 1, vjust = -0.6, size = 3, colour = "grey30") +
  geom_hline(yintercept = smr_overall$smr, linetype = "dotted",
             linewidth = 0.9, colour = "grey50") +
  annotate("text", x = 20, y = smr_overall$smr - 0.25,
           label = sprintf("Overall UKPDS SMR %.2f (95%% CI %.2f–%.2f)",
                           smr_overall$smr, smr_overall$lower, smr_overall$upper),
           hjust = 1, vjust = -0.6, size = 3, colour = "grey50") +
  scale_x_continuous(breaks = x_vals_ysd, limits = x_limits_ysd, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max(4, max(smr_by_ysd$smr, na.rm = TRUE) * 1.1)),
                  clip = "off") +
  labs(x = "Years since enrolment", y = "Standardised Mortality Ratio",
       title = "SMR by Years Since Enrolment with Study Comparisons (4209 cohort)") +
  theme_minimal(base_size = 13)

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

ggsave("Output/Figure/SMR_enrolment_by_followup_years_4209.png",
       p_main_ysd / p_table_ysd + plot_layout(heights = c(3, 1)),
       width = 20, height = 6, dpi = 300)

## 7.4 Figure - overall SMR by calendar year ----
#   File: SMR_enrolment_by_calendar_year_4209.png
#   (Person-time spans the full follow-up from enrolment, so this is the
#    enrolment-anchored calendar trend.)
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
       title = "SMR by Calendar Year with Study Comparisons (4209 cohort)") +
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

ggsave("Output/Figure/SMR_enrolment_by_calendar_year_4209.png",
       p_main_cal / p_table_cal + plot_layout(heights = c(3, 1)),
       width = 20, height = 6, dpi = 300)

## 7.5 Analysis frame with baseline subgroups (enrolment time scale) ----
#   Built here because both the subgroup figure (7.6) and Table 1A (7.7) use it.
ukpds_analysis <- ukpds_joined %>%
  left_join(ukpds_baseline, by = "id") %>%
  mutate(years_since_enrolment = year - enrolment_year)

## 7.6 Figure - subgroup SMR, 3-year bands (enrolment time scale) ----
#   File: SMR_enrolment_subgroups_3yr_4209.png
ukpds_analysis_3yr <- ukpds_analysis %>% add_fu_period("years_since_enrolment")

p_enr <- add_shared_labels(
  wrap_plots(
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "sex",           "fu_period", "A  Sex"),
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "age_cat",       "fu_period", "B  Baseline Age"),
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "smoking_cat",   "fu_period", "C  Smoking Status"),
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "hba1c_cat",     "fu_period", "D  Baseline HbA1c"),
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "bmi_cat",       "fu_period", "E  Baseline BMI"),
    make_subgroup_panel_3yr(ukpds_analysis_3yr, "ethnicity_cat", "fu_period", "F  Ethnicity"),
    ncol = 2),
  "SMR by Subgroup (4,209 randomised; time since enrolment; 3-year bands)",
  x_lab = "Years since enrolment (3-year bands)"
)

ggsave("Output/Figure/SMR_enrolment_subgroups_3yr_4209.png", p_enr,
       width = 8.5, height = 11.2, units = "in", dpi = 300)

## 7.7 Build Table 1A - SMR by baseline characteristics (enrolment) ----
#   Held in memory; written to the "Enrolment" tab of the Excel in Section 10.
ukpds_t1 <- ukpds_analysis %>% add_fu_period("years_since_enrolment")

rows <- list()
rows[[length(rows) + 1]] <- make_row(ukpds_t1, "Overall")
rows <- add_baseline_block(rows, ukpds_t1, "sex",           "Sex")
rows <- add_baseline_block(rows, ukpds_t1, "age_cat",       "Baseline age")
rows <- add_baseline_block(rows, ukpds_t1, "smoking_cat",   "Smoking status")
rows <- add_baseline_block(rows, ukpds_t1, "hba1c_cat",     "Baseline HbA1c")
rows <- add_baseline_block(rows, ukpds_t1, "bmi_cat",       "Baseline BMI (kg/m2)")
rows <- add_baseline_block(rows, ukpds_t1, "ethnicity_cat", "Ethnicity")

table1_A <- bind_rows(rows)
names(table1_A) <- c("Subgroup", "N", "Observed Deaths", "SMR (95% CI)",
                     "SMR 0-3y", "SMR 3-6y", "SMR 6-9y", "SMR 9-12y",
                     "SMR 12-15y", "SMR 15-18y", "SMR 18-21y", "SMR 21-24y",
                     "SMR 24-27y", "SMR 27-30y", "SMR 30-33y", "SMR 33-36y",
                     "SMR 36-39y", "SMR 39y+")


# 8. RANDOMISATION TIME SCALE  (clock starts at GL_date, post run-in) =========
#    Aim   : repeat the expected-deaths step and the overall SMR on the
#            randomisation clock, then produce two figures.
#    Input : UKPDS, life_1x1, ukpds_baseline.
#    Output: 2 figures (Output/Figure/) and the ukpds_analysis_rand frame
#            reused for Table 1B in Section 10.

## 8.1 Person-time split from randomisation + expected deaths ----
ukpds_long_rand <- split_person_time(UKPDS, entry_col = "randomization") %>%
  mutate(randomization_year = year(randomization))

ukpds_joined_rand <- ukpds_long_rand %>%
  left_join(life_1x1, by = c("sex", "year" = "Year", "age" = "Age")) %>%
  mutate(expected = mx * py)

## 8.2 Overall SMR (randomisation time scale) ----
O_rand <- sum(ukpds_joined_rand$death,    na.rm = TRUE)
E_rand <- sum(ukpds_joined_rand$expected, na.rm = TRUE)

smr_overall_rand <- list(
  smr   = O_rand / E_rand,
  lower = ifelse(O_rand == 0, 0, 0.5 * qchisq(0.025, df = 2 * O_rand))      / E_rand,
  upper =                        0.5 * qchisq(0.975, df = 2 * (O_rand + 1)) / E_rand
)

cat(sprintf("Randomisation time scale (n = 4,209): Overall SMR = %.2f (%.2f-%.2f)\n",
            smr_overall_rand$smr, smr_overall_rand$lower, smr_overall_rand$upper))

## 8.3 Figure - overall SMR by years since randomisation ----
#   File: SMR_randomization_by_followup_years_4209.png
smr_by_ysr <- ukpds_joined_rand %>%
  mutate(years_since_randomization = year - randomization_year) %>%
  group_by(years_since_randomization) %>%
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
  arrange(years_since_randomization)

x_vals_ysr   <- sort(unique(smr_by_ysr$years_since_randomization))
x_limits_ysr <- c(min(x_vals_ysr) - 0.5, 39.5)

p_main_ysr <- ggplot(smr_by_ysr, aes(x = years_since_randomization, y = smr)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.39, ymax = 1.94,
           fill = "#1b9e77", alpha = 0.15) +
  geom_hline(yintercept = 1.64, linetype = "dashed",
             colour = "#1b9e77", linewidth = 0.9) +
  annotate("text", x = min(x_vals_ysr) + 1, y = 1.8,
           label = "Meta-analysis SMR: 1.64 (95% CI 1.39–1.94)",
           hjust = 0, vjust = 0, size = 3.5,
           colour = "#1b9e77", fontface = "italic") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_text(aes(label = sprintf("%.2f", smr)),
            vjust = -1.5, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 1, linewidth = 0.7, colour = "grey30") +
  annotate("text", x = max(x_vals_ysr) - 0.5, y = 0.7,
           label = "Reference (SMR = 1)",
           hjust = 1, vjust = -0.6, size = 3, colour = "grey30") +
  geom_hline(yintercept = smr_overall_rand$smr, linetype = "dotted",
             linewidth = 0.9, colour = "grey50") +
  annotate("text", x = 20, y = smr_overall_rand$smr - 0.25,
           label = sprintf("Overall SMR (post-randomization) %.2f (95%% CI %.2f–%.2f)",
                           smr_overall_rand$smr, smr_overall_rand$lower,
                           smr_overall_rand$upper),
           hjust = 1, vjust = -0.6, size = 3, colour = "grey50") +
  scale_x_continuous(breaks = x_vals_ysr, limits = x_limits_ysr, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 4), clip = "off") +
  labs(x = "Years since randomization", y = "Standardised Mortality Ratio",
       title = "SMR by Years Since Randomization with Study Comparisons (4209 cohort)") +
  theme_minimal(base_size = 13)

graph_table_ysr <- smr_by_ysr %>%
  transmute(years_since_randomization,
            Population        = n_pop,
            `Observed deaths` = observed,
            `Male %`          = P_male,
            `Mean age`        = round(mean_age, 1)) %>%
  pivot_longer(-years_since_randomization, names_to = "Variable", values_to = "Value")

p_table_ysr <- ggplot(graph_table_ysr,
                      aes(x = years_since_randomization, y = Variable, label = Value)) +
  geom_text(size = 3.6, hjust = 0.5) +
  scale_x_continuous(breaks = x_vals_ysr, limits = x_limits_ysr, expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.y  = element_text(face = "bold"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 2, r = 6, b = 2, l = 6))

ggsave("Output/Figure/SMR_randomization_by_followup_years_4209.png",
       p_main_ysr / p_table_ysr + plot_layout(heights = c(3, 1)),
       width = 20, height = 6, dpi = 300)

## 8.4 Analysis frame with baseline subgroups (randomisation time scale) ----
#   Reused by the subgroup figure (8.5) and Table 1B (Section 10).
ukpds_analysis_rand <- ukpds_joined_rand %>%
  left_join(ukpds_baseline, by = "id") %>%
  mutate(years_since_randomization = year - randomization_year)

## 8.5 Figure - subgroup SMR, 3-year bands (randomisation time scale) ----
#   File: SMR_randomization_subgroups_3yr_4209.png
ukpds_analysis_rand_3yr <- ukpds_analysis_rand %>%
  add_fu_period("years_since_randomization")

p_rand <- add_shared_labels(
  wrap_plots(
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "sex",           "fu_period", "A  Sex"),
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "age_cat",       "fu_period", "B  Baseline Age"),
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "smoking_cat",   "fu_period", "C  Smoking Status"),
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "hba1c_cat",     "fu_period", "D  Baseline HbA1c"),
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "bmi_cat",       "fu_period", "E  Baseline BMI"),
    make_subgroup_panel_3yr(ukpds_analysis_rand_3yr, "ethnicity_cat", "fu_period", "F  Ethnicity"),
    ncol = 2),
  "SMR by Subgroup (4,209 randomised; time since randomisation; 3-year bands)",
  x_lab = "Years since randomisation (3-year bands)"
)

ggsave("Output/Figure/SMR_randomization_subgroups_3yr_4209.png", p_rand,
       width = 8.5, height = 11.2, units = "in", dpi = 300)


# 9. TREATMENT-ARM FIGURES  (post-randomisation, 3-year bands) ================
#    Aim   : reproduce the original UKPDS trial design and plot SMR by arm.
#    Input : extended_UKPDS (treatment flags) and ukpds_joined_rand.
#    Output: 2 figures (Output/Figure/) and the two arm-attached frames
#            (ukpds_arm_substudy, ukpds_arm_main) reused for Table 1B.
#
#    Trial design:
#      Substudy (UKPDS 34): overweight only
#        Conventional (substudy)  n =   411
#        Metformin                n =   342
#      Main study (UKPDS 33): all randomised
#        Conventional             n = 1,138   (includes the 411 substudy controls)
#        Intensive SU/insulin     n = 2,729
#    The 411 conventional-substudy patients appear in BOTH analyses by design:
#    they are the correct controls for the metformin comparison.

## 9.1 Restrict the arm lookup to randomised patients ----
arms_lookup <- extended_UKPDS %>%
  filter(!is.na(GL_date)) %>%             # drop 893 run-in exclusions
  mutate(id = as.character(ukpdsno))      # match id type in ukpds_joined_rand

## 9.2 Substudy arms (UKPDS 34) ----
substudy_arms <- arms_lookup %>%
  filter(metf == 0 | metf == 1) %>%       # excludes 951 SU/insulin in substudy
  transmute(
    id,
    arm = factor(
      case_when(
        SFUIns == 0 & metf == 0 ~ "Conventional",
        metf == 1               ~ "Metformin"
      ),
      levels = c("Conventional", "Metformin")
    )
  ) %>%
  filter(!is.na(arm))

## 9.3 Main-study arms (UKPDS 33) ----
main_arms <- arms_lookup %>%
  filter(SFUIns == 0 | SFUIns == 1) %>%   # excludes 342 metformin-only
  transmute(
    id,
    arm = factor(
      case_when(
        SFUIns == 0 ~ "Conventional",
        SFUIns == 1 ~ "Intensive SU/insulin"
      ),
      levels = c("Conventional", "Intensive SU/insulin")
    )
  ) %>%
  filter(!is.na(arm))

# Sanity checks - should print 411 + 342 and 1,138 + 2,729
cat("\nSubstudy n by arm (expect 411 + 342 = 753):\n")
print(substudy_arms %>% dplyr::count(arm))
cat("\nMain study n by arm (expect 1,138 + 2,729 = 3,867):\n")
print(main_arms %>% dplyr::count(arm))

## 9.4 Attach arm to the post-randomisation person-time data ----
ukpds_arm_substudy <- ukpds_joined_rand %>%
  inner_join(substudy_arms, by = "id") %>%
  mutate(years_since_randomization = year - randomization_year) %>%
  add_fu_period("years_since_randomization")

ukpds_arm_main <- ukpds_joined_rand %>%
  inner_join(main_arms, by = "id") %>%
  mutate(years_since_randomization = year - randomization_year) %>%
  add_fu_period("years_since_randomization")

## 9.5 Figure - main study arms, 3-year bands ----
#   File: SMR_randomization_arm_mainstudy_3yr_4209.png
p_arm_main_3yr <- make_arm_panel_3yr(
  ukpds_arm_main,
  panel_title = "SMR by Treatment Arm - Glycaemic study (3-year bands)",
  x_label     = "Years since randomization (3-year bands)",
  y_max       = 5
)

ggsave("Output/Figure/SMR_randomization_arm_mainstudy_3yr_4209.png",
       p_arm_main_3yr, width = 15, height = 7, dpi = 300)

## 9.6 Figure - overweight substudy arms, 3-year bands ----
#   File: SMR_randomization_arm_substudy_3yr_4209.png
p_arm_substudy_3yr <- make_arm_panel_3yr(
  ukpds_arm_substudy,
  panel_title = "SMR by Treatment Arm - Metformin study in overweight participants (3-year bands)",
  x_label     = "Years since randomization (3-year bands)",
  y_max       = 5
)

ggsave("Output/Figure/SMR_randomization_arm_substudy_3yr_4209.png",
       p_arm_substudy_3yr, width = 15, height = 7, dpi = 300)


# 10. TABLES TO ONE EXCEL WORKBOOK ============================================
#     Aim   : build Table 1B (randomisation time scale, baseline + arms) and
#             write both Table 1A and Table 1B into a single .xlsx file, each
#             on its own tab. This replaces the two separate CSV files.
#     Input : ukpds_analysis_rand, ukpds_arm_substudy, ukpds_arm_main, table1_A.
#     Output: Output/UKPDS_SMR_tables_4209.xlsx
#               tab "Enrolment"     = table1_A (built in Section 7)
#               tab "Randomization" = table1_B (built here)

## 10.1 Build Table 1B - baseline characteristics + treatment arms ----
ukpds_t1_rand <- ukpds_analysis_rand %>%
  add_fu_period("years_since_randomization")

rows_rand <- list()
rows_rand[[length(rows_rand) + 1]] <- make_row(ukpds_t1_rand, "Overall")

# Baseline characteristics (with a Missing row wherever any are NA)
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "sex",           "Sex")
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "age_cat",       "Baseline age")
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "smoking_cat",   "Smoking status")
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "hba1c_cat",     "Baseline HbA1c")
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "bmi_cat",       "Baseline BMI (kg/m2)")
rows_rand <- add_baseline_block(rows_rand, ukpds_t1_rand, "ethnicity_cat", "Ethnicity")

# Treatment-arm blocks (subsets of the cohort, not partitions)
rows_rand <- add_arm_block(rows_rand, ukpds_arm_substudy,
                           "Metformin study in overweight participants (UKPDS 34)")
rows_rand <- add_arm_block(rows_rand, ukpds_arm_main,
                           "Glycaemic study (UKPDS 33)")

table1_B <- bind_rows(rows_rand)
names(table1_B) <- c("Subgroup", "N", "Observed Deaths", "SMR (95% CI)",
                     "SMR 0-3y", "SMR 3-6y", "SMR 6-9y", "SMR 9-12y",
                     "SMR 12-15y", "SMR 15-18y", "SMR 18-21y", "SMR 21-24y",
                     "SMR 24-27y", "SMR 27-30y", "SMR 30-33y", "SMR 33-36y",
                     "SMR 36-39y", "SMR 39y+")

## 10.2 Optional console preview ----
cat("\n========================================================\n")
cat("TABLE 1B. SMR by Baseline Characteristics and Treatment Arm\n")
cat("          (4209 cohort, time since randomisation)\n")
cat("========================================================\n\n")
print(knitr::kable(table1_B, format = "pipe", align = "lrrrrrrrrrr"))

## 10.3 Write both tables to one workbook ----
write_two_tab_excel(
  path      = "Output/UKPDS_SMR_tables_4209.xlsx",
  tab1_name = "Enrolment",
  tab1_data = table1_A,
  tab2_name = "Randomization",
  tab2_data = table1_B
)

cat("\nDone. Figures written to Output/Figure/ and the table workbook to Output/.\n")