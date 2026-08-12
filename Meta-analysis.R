# SMR meta-analysis for mortality in adults with type 2 diabetes
# =============================================================================
#
# Reads:  SMR_T2DM_input_final.xlsx
# Writes: SMR_T2DM_results.xlsx, figures_combined.pdf,
#         figures_combined_period_sorted.pdf, figures/*.png
#
# Main analysis:
#   - Each SMR is analysed on the log scale.
#   - Primary model: multilevel random-effects meta-analysis with random
#     intercepts for cohort / study / estimate.
#   - Rows from the same cohort are clustered together.
#   - The headline confidence interval is cluster-robust by cohort.
# =============================================================================



# -- 1. Setup ----------------------------------------------------------------

## 1.1 Libraries ----
library(metafor)
library(dplyr)
library(openxlsx)
library(clubSandwich)

R.version.string
sapply(c("metafor", "dplyr", "clubSandwich", "openxlsx"),
       function(p) as.character(packageVersion(p)))

## 1.2 File paths ----
ROOT     <- "J:/3rd project_SMR/Systematic review/Final review/05_meta-analysis_R"
XLSX_IN  <- file.path(ROOT, "data_in/SMR_T2DM_review.xlsx")
XLSX_OUT <- file.path(ROOT, "outputs/SMR_T2DM_results.xlsx")
FIG_DIR  <- file.path(ROOT, "outputs/figures")
OUT_DIR  <- dirname(XLSX_OUT)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

## 1.3 Constants ----
UKPDS_START <- 1977
UKPDS_END   <- 2021

## 1.4 XML-entity decoder ----
# Some spreadsheet cells contain XML character references (e.g. "&amp;",
# "&#241;", "&gt;=") that openxlsx leaves undecoded on non-UTF-8 system locales
# (e.g. GBK); a few cells are also double-encoded at source. Decode them so
# every label renders as real text ("&", "n-tilde", ">=", ...) regardless of
# the machine's locale.
decode_xml_entities <- function(x) {
  if (!is.character(x)) return(x)
  # Run the named entities twice to undo double encoding (&amp;amp; -> &amp; -> &).
  for (pass in 1:2) {
    x <- gsub("&lt;",   "<",  x, fixed = TRUE)
    x <- gsub("&gt;",   ">",  x, fixed = TRUE)
    x <- gsub("&quot;", '"',  x, fixed = TRUE)
    x <- gsub("&apos;", "'",  x, fixed = TRUE)
    x <- gsub("&amp;",  "&",  x, fixed = TRUE)
  }
  m <- gregexpr("&#[xX]?[0-9A-Fa-f]+;", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(codes) {
    vapply(codes, function(cc) {
      body <- substr(cc, 3L, nchar(cc) - 1L)
      n <- if (grepl("^[xX]", body)) strtoi(substring(body, 2L), 16L)
      else suppressWarnings(as.integer(body))
      if (is.na(n)) cc else intToUtf8(n)
    }, character(1), USE.NAMES = FALSE)
  })
  x
}

## 1.5 ASCII fold for display labels ----
# Fold accented Latin letters (and en/em dashes, sharp s) to plain ASCII so
# forest labels render identically on every machine, even if the input is
# regenerated with accented characters. Built from code points, so this source
# stays pure ASCII.
ascii_fold <- function(x) {
  if (!is.character(x)) return(x)
  acc <- intToUtf8(c(0xe0,0xe1,0xe2,0xe3,0xe4,0xe5,0xe7,0xe8,0xe9,0xea,0xeb,0xec,0xed,0xee,0xef,
                     0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf8,0xf9,0xfa,0xfb,0xfc,0xfd,
                     0xc0,0xc1,0xc2,0xc3,0xc4,0xc5,0xc7,0xc8,0xc9,0xca,0xcb,0xcc,0xcd,0xce,0xcf,
                     0xd1,0xd2,0xd3,0xd4,0xd5,0xd6,0xd8,0xd9,0xda,0xdb,0xdc,0xdd))
  pln <- "aaaaaaceeeeiiiinoooooouuuuyAAAAAACEEEEIIIINOOOOOOUUUUY"
  x <- chartr(acc, pln, x)
  x <- gsub(intToUtf8(0x2013), "-",  x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2014), "-",  x, fixed = TRUE)
  x <- gsub(intToUtf8(0x00df), "ss", x, fixed = TRUE)
  x
}

## 1.6 Year extraction from free-text period labels ----
first_year <- function(x) {
  hits <- regmatches(x, gregexpr("[0-9]{4}", x))
  vapply(hits, function(h) if (length(h)) as.numeric(h[1]) else NA_real_, numeric(1))
}
last_year <- function(x) {
  hits <- regmatches(x, gregexpr("[0-9]{4}", x))
  vapply(hits, function(h) if (length(h)) as.numeric(h[length(h)]) else NA_real_, numeric(1))
}

# -- 2. Read input data ------------------------------------------------------

paper_overall <- openxlsx::read.xlsx(XLSX_IN, sheet = "1_Paper_overall", startRow = 4)
sub_period    <- openxlsx::read.xlsx(XLSX_IN, sheet = "2_Sub_period",    startRow = 4)

# Clean any XML entities openxlsx may have left undecoded (locale-dependent).
paper_overall[] <- lapply(paper_overall, decode_xml_entities)
sub_period[]    <- lapply(sub_period,    decode_xml_entities)

cat("Paper_overall:", nrow(paper_overall), "rows x", ncol(paper_overall), "columns\n")
cat("Sub_period   :", nrow(sub_period),    "rows x", ncol(sub_period),    "columns\n")

# -- 3. Build analysis frame -------------------------------------------------

## 3.1 Paper-level rows (one all-persons estimate per paper) ----
paper_long <- paper_overall |>
  transmute(
    estimate_id = paste0(StudyID, "_overall"),
    source      = "paper-overall",
    StudyID     = as.character(StudyID),
    First_Author_Lastname,
    cohort_id,
    Country,
    WHO_region,
    WB_HIC,
    OECD,
    Setting_Type,
    Year = suppressWarnings(as.numeric(Year)),
    period_label = ifelse(
      !is.na(Study_Period_Start_Year) & !is.na(Study_Period_End_Year),
      paste0(Study_Period_Start_Year, "-", Study_Period_End_Year),
      ifelse(!is.na(Study_Period_Raw) & Study_Period_Raw != "NR",
             as.character(Study_Period_Raw), "period unknown")
    ),
    period_start = suppressWarnings(as.numeric(Study_Period_Start_Year)),
    period_end   = suppressWarnings(as.numeric(Study_Period_End_Year)),
    cause        = "all-cause",
    SMR          = suppressWarnings(as.numeric(SMR_Overall)),
    SMR_LCI      = suppressWarnings(as.numeric(SMR_Overall_LCI)),
    SMR_UCI      = suppressWarnings(as.numeric(SMR_Overall_UCI)),
    N_deaths     = suppressWarnings(as.numeric(N_Deaths_Total)),
    N_subjects   = suppressWarnings(as.numeric(N_Subjects_Total))
  )

## 3.2 Sub-period rows ----
# Add the calendar sub-period estimates from the second workbook tab.
# Rows labelled as "overall" / "whole period" / "full follow-up" are skipped so
# they do not duplicate the paper-level overall estimate from tab 1.
paper_meta <- paper_overall |>
  mutate(StudyID = as.character(StudyID)) |>
  select(StudyID, OECD, Setting_Type, Year, N_Subjects_Total) |>
  mutate(Year = suppressWarnings(as.numeric(Year)))

# Belt-and-braces: if the tab carries a sex column, keep only all-persons rows.
if ("sex" %in% names(sub_period)) {
  keep_all_persons <- is.na(sub_period$sex) |
    tolower(trimws(as.character(sub_period$sex))) %in%
    c("both", "all", "total", "both sexes", "")
  sub_period <- sub_period[keep_all_persons, , drop = FALSE]
}

sub_long <- sub_period |>
  mutate(StudyID = as.character(StudyID)) |>
  filter(
    !grepl("overall|whole period|whole-period|full follow",
           tolower(as.character(period_label))),
    # De-duplication only. Keep all-ages, all-persons calendar sub-periods and
    # drop age- or sex-restricted sub-strata (e.g. Sasaki's age<65 "lt65" rows):
    # these overlap the all-ages estimate for the same period and would
    # double-count that person-time. The marker lives in estimate_id, not
    # period_label.
    !grepl("lt[0-9]|ge[0-9]|_age|_male|_female|_men|_women",
           tolower(as.character(estimate_id)))
  ) |>
  mutate(
    SMR          = suppressWarnings(as.numeric(SMR)),
    SMR_LCI      = suppressWarnings(as.numeric(SMR_LCI)),
    SMR_UCI      = suppressWarnings(as.numeric(SMR_UCI)),
    N_deaths     = suppressWarnings(as.numeric(N_deaths)),
    period_start = first_year(as.character(period_label)),
    period_end   = last_year(as.character(period_label))
  ) |>
  left_join(paper_meta, by = "StudyID") |>
  mutate(N_subjects = suppressWarnings(as.numeric(N_Subjects_Total))) |>
  select(
    estimate_id, source, StudyID, First_Author_Lastname, cohort_id, Country,
    WHO_region, WB_HIC, OECD, Setting_Type, Year,
    period_label, period_start, period_end,
    cause, SMR, SMR_LCI, SMR_UCI, N_deaths, N_subjects
  )

## 3.3 Cohort abbreviations (for forest plots) ----
cohort_abbr <- c(
  "Verona Diabetes Study (Italy)"                                            = "Verona-DS",
  "Osaka Centre NIDDM Cohort (Japan)"                                        = "Osaka-NIDDM",
  "Canterbury T2D Clinic Survey 1989 (New Zealand)"                          = "Canterbury",
  "Vantaa Public PHC T2D Cohort (Finland)"                                   = "Vantaa-PHC",
  "ZODIAC T2D Cohort (Netherlands)"                                          = "ZODIAC",
  "ZODIAC-1 Background Population (Netherlands)"                             = "ZODIAC-1",
  "Oxford Community T2D Cohort (UK)"                                         = "Oxford",
  "Lithuanian National T2D Cohort"                                           = "LT-National",
  "Korean National Health Insurance T2D Cohort"                              = "KR-NHI",
  "Saudi National Diabetes Registry (SNDR)"                                  = "SNDR",
  "Ayrshire and Arran Diabetes Cohort (Scotland)"                            = "Ayrshire",
  "WHO Multinational Study of Vascular Disease in Diabetes (Switzerland arm)" = "WHO-MSVDD",
  "Casale Monferrato Diabetes Survey (Italy)"                                = "Casale",
  "Diabetes Incidence Study in Sweden (DISS, age 15-34)"                     = "DISS",
  "CPRD T2D Cohort (England)"                                                = "CPRD",
  "Tri-Service General Hospital T2D Cohort (Taiwan)"                         = "Tri-Service",
  "Alcaniz-Teruel T2D Cohort (Spain)"                                        = "Alcaniz",
  "Kronoberg DIK T2D Cohort (Sweden)"                                        = "Kronoberg-DIK",
  "Rio de Janeiro Hospital T2D Cohort (Brazil)"                              = "Rio-T2D",
  "German DMP National Mortality Analysis"                                   = "DE-DMP"
)

## 3.4 Combine + compute log-SMR and within-study variance ----
# Variance:
#   Preferred: vi = ((log(UCI) - log(LCI)) / (2 * 1.96))^2
#   Fallback (no CI): vi = 1 / observed deaths (Poisson approximation)
long <- bind_rows(paper_long, sub_long) |>
  mutate(
    Country = case_when(
      Country == "Taiwan"                            ~ "China",
      grepl("scotland", Country, ignore.case = TRUE) ~ "Scotland",
      grepl("wales",    Country, ignore.case = TRUE) ~ "Wales",
      grepl("england",  Country, ignore.case = TRUE) ~ "England",
      Country == "United Kingdom"                    ~ "England",
      Country == "Korea"                             ~ "South Korea",
      TRUE                                           ~ Country
    ),
    cohort_abbr = ifelse(cohort_id %in% names(cohort_abbr),
                         cohort_abbr[cohort_id], cohort_id),
    yi = log(SMR),
    vi_from_ci = ifelse(
      !is.na(SMR_LCI) & !is.na(SMR_UCI) & SMR_LCI > 0 & SMR_UCI > SMR_LCI,
      ((log(SMR_UCI) - log(SMR_LCI)) / (2 * 1.96))^2,
      NA_real_
    ),
    vi_from_deaths = ifelse(!is.na(N_deaths) & N_deaths > 0, 1 / N_deaths, NA_real_),
    vi_use = coalesce(vi_from_ci, vi_from_deaths),
    variance_source = ifelse(!is.na(vi_from_ci), "reported 95% CI",
                             ifelse(!is.na(vi_from_deaths), "1 / observed deaths",
                                    NA_character_))
  ) |>
  filter(
    cause == "all-cause" | is.na(cause),
    !is.na(SMR), SMR > 0,
    !is.na(yi), is.finite(yi),
    !is.na(vi_use), is.finite(vi_use), vi_use > 0
  )

## 3.5 Group flags (HIC, UKPDS-era overlap) ----
# Added BEFORE any pool is fit, so all res_* objects carry these columns and
# the forest plots can display them as extra columns.
long <- long |>
  mutate(
    Country               = ascii_fold(Country),
    cohort_abbr           = ascii_fold(cohort_abbr),
    First_Author_Lastname = ascii_fold(First_Author_Lastname),
    HIC_grp = case_when(
      toupper(trimws(as.character(WB_HIC))) %in%
        c("TRUE", "T", "Y", "1", "YES", "HIC") ~ "HIC",
      toupper(trimws(as.character(WB_HIC))) %in%
        c("FALSE", "F", "N", "0", "NO", "NON-HIC", "LMIC", "UMIC", "LIC") ~ "Non-HIC",
      TRUE ~ "Unknown"
    ),
    ukpds_overlap_grp = case_when(
      !is.na(period_start) & !is.na(period_end) &
        period_start <= UKPDS_END & period_end >= UKPDS_START ~ "Yes",
      !is.na(period_start) & !is.na(period_end) ~ "No",
      TRUE ~ NA_character_
    ),
    hic_ukpds_overlap_grp = case_when(
      is.na(ukpds_overlap_grp) | HIC_grp == "Unknown" ~ NA_character_,
      HIC_grp == "HIC" & ukpds_overlap_grp == "Yes"   ~ "Yes",
      TRUE ~ "No"
    )
  )

cat("Built long-format data:",
    sum(long$source == "paper-overall"), "paper-overall +",
    sum(long$source == "sub-period"),    "sub-period =", nrow(long), "rows\n")
cat("Distinct studies:", length(unique(long$StudyID)),
    " Distinct cohorts:", length(unique(long$cohort_id)), "\n")

# -- 4. Helper functions -----------------------------------------------------

## 4.1 fit_pool() - multilevel model with cluster-robust CI ----
fit_pool <- function(data, stratum_label) {
  data <- data |>
    filter(!is.na(yi), !is.na(vi_use), is.finite(yi), is.finite(vi_use), vi_use > 0)
  
  empty_summary <- data.frame(
    Stratum          = stratum_label,
    k_estimates      = nrow(data),
    k_studies        = length(unique(data$StudyID)),
    k_cohorts        = length(unique(data$cohort_id)),
    SMR              = NA_real_,
    SMR_LCI_naive    = NA_real_,
    SMR_UCI_naive    = NA_real_,
    SMR_LCI_robust   = NA_real_,
    SMR_UCI_robust   = NA_real_,
    sigma2_cohort    = NA_real_,
    sigma2_study     = NA_real_,
    sigma2_estimate  = NA_real_,
    PI_LCI           = NA_real_,
    PI_UCI           = NA_real_,
    robust_ci_source = NA_character_
  )
  
  if (nrow(data) < 2 || length(unique(data$cohort_id)) < 2) {
    return(list(model = NULL, robust = NULL, data = data, summary = empty_summary))
  }
  
  model <- tryCatch(
    metafor::rma.mv(
      yi = yi, V = vi_use,
      random = ~ 1 | cohort_id / StudyID / estimate_id,
      data = data, method = "REML", slab = estimate_id
    ),
    error = function(e) NULL
  )
  
  # Some small strata fail with the default optimizer even when estimable.
  # Try two standard optimizers before giving up.
  if (is.null(model)) {
    model <- tryCatch(
      metafor::rma.mv(
        yi = yi, V = vi_use,
        random = ~ 1 | cohort_id / StudyID / estimate_id,
        data = data, method = "REML", slab = estimate_id,
        control = list(optimizer = "optim", optmethod = "BFGS")
      ),
      error = function(e) NULL
    )
  }
  if (is.null(model)) {
    model <- tryCatch(
      metafor::rma.mv(
        yi = yi, V = vi_use,
        random = ~ 1 | cohort_id / StudyID / estimate_id,
        data = data, method = "REML", slab = estimate_id,
        control = list(optimizer = "optim", optmethod = "Nelder-Mead")
      ),
      error = function(e) {
        message("Model failed for ", stratum_label, ": ", conditionMessage(e))
        NULL
      }
    )
  }
  
  if (is.null(model)) {
    return(list(model = NULL, robust = NULL, data = data, summary = empty_summary))
  }
  
  robust_ok <- TRUE
  robust_model <- tryCatch(
    metafor::robust(model, cluster = data$cohort_id, clubSandwich = TRUE),
    error = function(e) { robust_ok <<- FALSE; model }
  )
  
  pred <- tryCatch(predict(model, transf = exp),
                   error = function(e) list(pi.lb = NA_real_, pi.ub = NA_real_))
  
  summary_row <- data.frame(
    Stratum          = stratum_label,
    k_estimates      = model$k,
    k_studies        = length(unique(data$StudyID)),
    k_cohorts        = length(unique(data$cohort_id)),
    SMR              = as.numeric(exp(model$b[1])),
    SMR_LCI_naive    = as.numeric(exp(model$ci.lb[1])),
    SMR_UCI_naive    = as.numeric(exp(model$ci.ub[1])),
    SMR_LCI_robust   = as.numeric(exp(robust_model$ci.lb[1])),
    SMR_UCI_robust   = as.numeric(exp(robust_model$ci.ub[1])),
    sigma2_cohort    = model$sigma2[1],
    sigma2_study     = model$sigma2[2],
    sigma2_estimate  = model$sigma2[3],
    PI_LCI           = as.numeric(pred$pi.lb[1]),
    PI_UCI           = as.numeric(pred$pi.ub[1]),
    robust_ci_source = ifelse(robust_ok, "cluster robust by cohort",
                              "model-based fallback")
  )
  
  list(model = model, robust = robust_model, data = data, summary = summary_row)
}

## 4.2 draw_forest() - study-level forest plot, written straight to PNG ----
# Parameters:
#   sort_by_period - TRUE  = sort rows by period_start / period_end (default)
#                    FALSE = sort by Year / Author
#   show_hic       - add a "HIC" Yes/No column
#   show_ukpds     - add a "UKPDS-era" Yes/No column
draw_forest <- function(pool, title_text, file_stub,
                        sort_by_period = TRUE,
                        show_hic = FALSE, show_ukpds = FALSE) {
  if (is.null(pool$model) || nrow(pool$data) == 0) return(invisible(NULL))
  
  data         <- pool$data
  model        <- pool$model
  robust_model <- pool$robust
  
  # ---- Row sort ----
  row_order <- if (isTRUE(sort_by_period)) {
    order(data$period_start, data$period_end, data$Year,
          data$First_Author_Lastname, data$StudyID, data$period_label)
  } else {
    order(data$Year, data$First_Author_Lastname, data$StudyID, data$period_label)
  }
  data <- data[row_order, ]
  
  # ---- Recompute display CIs from yi / vi_use ----
  data$ci_lci_plot <- exp(data$yi - 1.96 * sqrt(data$vi_use))
  data$ci_uci_plot <- exp(data$yi + 1.96 * sqrt(data$vi_use))
  
  # ---- Cell text helpers ----
  fmt_count <- function(x) {
    ifelse(is.na(x), "NR",
           format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
  }
  fmt_p <- function(p) {
    ifelse(is.na(p), "NA",
           ifelse(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p))))
  }
  shorten_text <- function(x, max_chars) {
    ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
  }
  yn_from_hic   <- function(x) ifelse(is.na(x) | x == "Unknown", "-",
                                      ifelse(x == "HIC", "Yes", "No"))
  yn_passthrough <- function(x) ifelse(is.na(x), "-", x)
  
  # ---- Cell values ----
  author_year <- paste0(shorten_text(data$First_Author_Lastname, 12), ", ", data$Year)
  death_n <- ifelse(
    is.na(data$N_deaths) | is.na(data$N_subjects), "NR",
    paste0(fmt_count(data$N_deaths), " / ", fmt_count(data$N_subjects))
  )
  hic_cells   <- if (show_hic)   yn_from_hic(data$HIC_grp)            else NULL
  ukpds_cells <- if (show_ukpds) yn_passthrough(data$ukpds_overlap_grp) else NULL
  
  # Deaths / N total: count each study once. Calendar sub-period rows carry the
  # study-level N_subjects (joined from tab 1), so a plain row-sum would
  # multiply-count any cohort contributing more than one row.
  study_tot <- data |>
    group_by(StudyID) |>
    summarise(
      N_subjects = suppressWarnings(max(N_subjects, na.rm = TRUE)),
      N_deaths   = suppressWarnings(max(N_deaths,   na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(
      N_subjects = ifelse(is.finite(N_subjects), N_subjects, NA_real_),
      N_deaths   = ifelse(is.finite(N_deaths),   N_deaths,   NA_real_)
    )
  total_deaths   <- sum(study_tot$N_deaths,   na.rm = TRUE)
  total_subjects <- sum(study_tot$N_subjects, na.rm = TRUE)
  total_death_n  <- paste0(fmt_count(total_deaths), " / ", fmt_count(total_subjects))
  
  pooled_smr  <- as.numeric(exp(model$b[1]))
  pooled_lci  <- as.numeric(exp(robust_model$ci.lb[1]))
  pooled_uci  <- as.numeric(exp(robust_model$ci.ub[1]))
  pooled_text <- sprintf("%.2f (%.2f-%.2f)", pooled_smr, pooled_lci, pooled_uci)
  
  # ---- X-axis range ----
  x_min <- max(0.05, min(c(data$ci_lci_plot, pooled_lci), na.rm = TRUE) * 0.85)
  x_max <- min(100,  max(c(data$ci_uci_plot, pooled_uci), na.rm = TRUE) * 1.15)
  x_min <- min(x_min, 0.5)
  x_max <- max(x_max, 4)
  axis_ticks <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8, 12, 20, 50, 100)
  axis_ticks <- axis_ticks[axis_ticks >= x_min & axis_ticks <= x_max]
  
  # ---- Layout coordinates ----
  n_rows           <- nrow(data)
  y_rows           <- rev(seq_len(n_rows)) + 1
  y_total          <- 1
  y_heterogeneity  <- 0.15
  y_overall_effect <- -0.55
  y_header         <- n_rows + 2.2
  y_limit          <- c(-1.0, n_rows + 2.8)
  
  # ---- Column positions inside the left panel (0..1) ----
  # Layout depends on how many extra columns we add.
  if (!show_hic && !show_ukpds) {
    col_x <- list(author = 0.00, period = 0.19, cohort = 0.33,
                  country = 0.50, deaths = 0.73,
                  hic = NA_real_, ukpds = NA_real_)
    table_w <- 0.48
  } else if (show_hic && !show_ukpds) {
    col_x <- list(author = 0.00, period = 0.17, cohort = 0.30,
                  country = 0.46, deaths = 0.68,
                  hic = 0.88, ukpds = NA_real_)
    table_w <- 0.52
  } else if (!show_hic && show_ukpds) {
    col_x <- list(author = 0.00, period = 0.17, cohort = 0.30,
                  country = 0.46, deaths = 0.68,
                  hic = NA_real_, ukpds = 0.86)
    table_w <- 0.52
  } else {
    col_x <- list(author = 0.00, period = 0.15, cohort = 0.27,
                  country = 0.41, deaths = 0.59,
                  hic = 0.79, ukpds = 0.90)
    table_w <- 0.56
  }
  layout_widths <- c(table_w, 1 - table_w - 0.16, 0.16)
  
  # ---- Strings for inferential tests ----
  heterogeneity_text <- sprintf(
    "Test for heterogeneity: QE = %.2f, p %s", model$QE, fmt_p(model$QEp)
  )
  overall_effect_text <- sprintf(
    "Test for overall effect: z = %.2f, p %s",
    as.numeric(robust_model$zval[1]), fmt_p(as.numeric(robust_model$pval[1]))
  )
  
  header_cex  <- 1.00
  row_cex     <- 0.88
  country_cex <- 0.82
  total_cex   <- 0.92
  axis_cex    <- 1.00
  stats_cex   <- 0.86
  
  # ---- Draw straight into a PNG device ----
  png_path <- file.path(FIG_DIR, paste0(file_stub, ".png"))
  png(png_path,
      width  = 3000,
      height = max(1050, 42 * n_rows + 420),
      res    = 220)
  on.exit(dev.off(), add = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = layout_widths)
  par(oma = c(1.2, 0, 2.2, 0))
  
  # Left panel - study details
  par(mar = c(4, 0.6, 0.6, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = y_limit)
  
  text(col_x$author,  y_header, "Author, Year", adj = 0, font = 2, cex = header_cex)
  text(col_x$period,  y_header, "Period",       adj = 0, font = 2, cex = header_cex)
  text(col_x$cohort,  y_header, "Cohort",       adj = 0, font = 2, cex = header_cex)
  text(col_x$country, y_header, "Country",      adj = 0, font = 2, cex = header_cex)
  text(col_x$deaths,  y_header, "Deaths / N",   adj = 0, font = 2, cex = header_cex)
  if (show_hic)   text(col_x$hic,   y_header, "HIC",       adj = 0, font = 2, cex = header_cex)
  if (show_ukpds) text(col_x$ukpds, y_header, "UKPDS-era", adj = 0, font = 2, cex = header_cex)
  
  segments(0, n_rows + 1.6, 1, n_rows + 1.6, col = "grey70")
  segments(0, 1.5,          1, 1.5,          col = "grey70")
  
  text(col_x$author,  y_rows, author_year,       adj = 0, cex = row_cex)
  text(col_x$period,  y_rows, data$period_label, adj = 0, cex = row_cex)
  text(col_x$cohort,  y_rows, shorten_text(data$cohort_abbr, 18), adj = 0, cex = row_cex)
  text(col_x$country, y_rows, shorten_text(data$Country, 18),     adj = 0, cex = country_cex)
  text(col_x$deaths,  y_rows, death_n,           adj = 0, cex = row_cex)
  if (show_hic)   text(col_x$hic,   y_rows, hic_cells,   adj = 0, cex = row_cex)
  if (show_ukpds) text(col_x$ukpds, y_rows, ukpds_cells, adj = 0, cex = row_cex)
  
  text(col_x$author, y_total, "Total",       adj = 0, font = 2, cex = total_cex)
  text(col_x$deaths, y_total, total_death_n, adj = 0, font = 2, cex = total_cex)
  text(col_x$author, y_heterogeneity,  heterogeneity_text,  adj = 0, cex = stats_cex)
  text(col_x$author, y_overall_effect, overall_effect_text, adj = 0, cex = stats_cex)
  
  # Middle panel - forest plot
  par(mar = c(4, 0.2, 0.6, 0.2))
  plot(NA, NA,
       xlim = c(x_min, x_max), ylim = y_limit, log = "x",
       xaxt = "n", yaxt = "n",
       xlab = "SMR (log scale)", ylab = "", bty = "n",
       cex.lab = axis_cex)
  axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = axis_cex)
  abline(v = 1,          lty = 3, col = "grey45")
  abline(v = pooled_smr, lty = 2, col = "#0B2C4D")
  segments(x_min, n_rows + 1.6, x_max, n_rows + 1.6, col = "grey70")
  segments(x_min, 1.5,          x_max, 1.5,          col = "grey70")
  
  segments(data$ci_lci_plot, y_rows, data$ci_uci_plot, y_rows,
           lwd = 2.0, col = "#0B2C4D")
  points(data$SMR, y_rows, pch = 15, cex = 0.95, col = "#C4373A")
  
  polygon(
    x = c(pooled_lci, pooled_smr, pooled_uci, pooled_smr),
    y = c(y_total, y_total + 0.30, y_total, y_total - 0.30),
    col = "#0B2C4D", border = "#0B2C4D"
  )
  
  # Right panel - SMR (95% CI) text
  par(mar = c(4, 0.2, 0.6, 0.6))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = y_limit)
  text(0.00, y_header, "SMR (95% CI)", adj = 0, font = 2, cex = header_cex)
  segments(0, n_rows + 1.6, 1, n_rows + 1.6, col = "grey70")
  segments(0, 1.5,          1, 1.5,          col = "grey70")
  text(0.00, y_rows,
       sprintf("%.2f (%.2f-%.2f)", data$SMR, data$ci_lci_plot, data$ci_uci_plot),
       adj = 0, cex = row_cex)
  text(0.00, y_total, pooled_text, adj = 0, font = 2, cex = total_cex)
  
  mtext(title_text, side = 3, outer = TRUE, line = 0.7, font = 2, cex = 1.25)
  
  layout(1)
  suppressWarnings(par(old_par))
  
  cat("PNG written:", png_path, "\n")
  invisible(png_path)
}

## 4.3 draw_combined_subgroup_plot() - subgroup summary forest ----
# This figure is a summary plot: each row is a subgroup-specific pooled SMR,
# not an independent study. Rows with only one cohort show "not pooled".
draw_combined_subgroup_plot <- function(subgroup_all,
                                        file_stub = "fig_subgroups_all") {
  if (nrow(subgroup_all) == 0) return(invisible(NULL))
  
  subgroup_display <- data.frame()
  
  for (this_subgroup in unique(subgroup_all$Subgroup)) {
    one_block <- subgroup_all |> filter(Subgroup == this_subgroup)
    this_p    <- unique(one_block$between_group_p)[1]
    header_label <- if (all(is.na(one_block$between_group_p))) {
      this_subgroup
    } else {
      this_p_label <- ifelse(is.na(this_p), "p = NA",
                             paste0("p = ", sprintf("%.3f", this_p)))
      paste0(this_subgroup, " (moderator ", this_p_label, ")")
    }
    
    subgroup_display <- bind_rows(
      subgroup_display,
      data.frame(
        row_type = "header", Subgroup = this_subgroup, Group = header_label,
        k_estimates = NA_integer_, k_cohorts = NA_integer_,
        SMR = NA_real_, SMR_LCI = NA_real_, SMR_UCI = NA_real_
      ),
      one_block |>
        transmute(row_type = "data", Subgroup, Group,
                  k_estimates, k_cohorts, SMR, SMR_LCI, SMR_UCI)
    )
  }
  
  row_gap <- 0.62
  subgroup_display$y <- rev(seq_len(nrow(subgroup_display))) * row_gap
  y_header <- max(subgroup_display$y) + 0.38
  y_limit  <- c(0.20, max(subgroup_display$y) + 0.60)
  
  x_min        <- 0.45
  x_max        <- 8
  x_ticks      <- c(0.5, 1, 1.5, 2, 3, 4, 6, 8)
  x_left_text  <- x_min / 1.04
  x_right_text <- x_max * 1.15
  
  png_path <- file.path(FIG_DIR, paste0(file_stub, ".png"))
  png(png_path,
      width  = 3200,
      height = max(1400, 68 * nrow(subgroup_display) + 420),
      res    = 220)
  on.exit(dev.off(), add = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5, 15, 3.5, 7), xpd = NA)
  
  plot(NA, NA,
       xlim = c(x_min, x_max), ylim = y_limit, log = "x",
       xaxt = "n", yaxt = "n",
       xlab = "SMR (log scale)", ylab = "", bty = "n")
  axis(1, at = x_ticks, labels = x_ticks)
  grid(nx = NA, ny = NULL, col = "grey90")
  par(xpd = FALSE)
  segments(x0 = 1, x1 = 1, y0 = y_limit[1], y1 = y_limit[2],
           lty = 2, col = "grey45")
  par(xpd = NA)
  
  text(x = x_left_text,  y = y_header, "Subgroup level", adj = 1, font = 2, cex = 0.75)
  text(x = x_right_text, y = y_header, "SMR (95% CI)",   adj = 0, font = 2, cex = 0.75)
  
  for (j in seq_len(nrow(subgroup_display))) {
    this_row <- subgroup_display[j, ]
    
    if (this_row$row_type == "header") {
      text(x = x_left_text, y = this_row$y, labels = this_row$Group,
           adj = 1, font = 2, cex = 0.66)
      segments(x0 = x_min, x1 = x_max,
               y0 = this_row$y - 0.26, y1 = this_row$y - 0.26, col = "grey85")
    } else {
      left_label <- sprintf("  %s (k=%d, cohorts=%d)",
                            this_row$Group, this_row$k_estimates, this_row$k_cohorts)
      text(x = x_left_text, y = this_row$y, labels = left_label, adj = 1, cex = 0.62)
      
      if (is.na(this_row$SMR)) {
        text(x = x_right_text, y = this_row$y, labels = "not pooled",
             adj = 0, cex = 0.62, col = "grey35")
      } else {
        lci_clipped <- max(this_row$SMR_LCI, x_min)
        uci_clipped <- min(this_row$SMR_UCI, x_max)
        segments(lci_clipped, this_row$y, uci_clipped, this_row$y,
                 lwd = 2, col = "#0B2C4D")
        if (this_row$SMR_LCI < x_min) {
          arrows(x0 = x_min * 1.15, y0 = this_row$y, x1 = x_min, y1 = this_row$y,
                 length = 0.05, lwd = 2, col = "#0B2C4D")
        }
        if (this_row$SMR_UCI > x_max) {
          arrows(x0 = x_max / 1.15, y0 = this_row$y, x1 = x_max, y1 = this_row$y,
                 length = 0.05, lwd = 2, col = "#0B2C4D")
        }
        points(this_row$SMR, this_row$y, pch = 18, cex = 1.05, col = "#C4373A")
        text(x = x_right_text, y = this_row$y,
             labels = sprintf("%.2f (%.2f-%.2f)",
                              this_row$SMR, this_row$SMR_LCI, this_row$SMR_UCI),
             adj = 0, cex = 0.62)
      }
    }
  }
  
  title("Combined subgroup summary forest plot")
  suppressWarnings(par(old_par))
  
  cat("PNG written:", png_path, "\n")
  invisible(png_path)
}

# -- 5. Main pooled analyses -------------------------------------------------

## 5.1 Overall pooled SMR ----
d_overall   <- long |> filter(cause == "all-cause")
res_overall <- fit_pool(d_overall, "Overall (all-cause)")
main_tab    <- bind_rows(res_overall$summary)

## 5.2 HIC vs Non-HIC ----
res_hic <- fit_pool(
  long |> filter(HIC_grp == "HIC", cause == "all-cause"),
  "HIC: overall"
)
res_nonhic <- fit_pool(
  long |> filter(HIC_grp == "Non-HIC", cause == "all-cause"),
  "Non-HIC: overall"
)
hic_tab <- bind_rows(res_hic$summary, res_nonhic$summary)

## 5.3 UKPDS-era overlap ----
# UKPDS ran 1977-2021 (recruitment + post-trial monitoring). Include any
# estimate whose follow-up period overlaps that window:
#   period_start <= 2021 AND period_end >= 1977
ukpds_overlap <- long |>
  filter(!is.na(period_start), !is.na(period_end),
         period_start <= UKPDS_END, period_end >= UKPDS_START)

res_ukpds <- fit_pool(
  ukpds_overlap |> filter(cause == "all-cause"),
  "UKPDS-era overlap: overall"
)
res_hic_ukpds <- fit_pool(
  ukpds_overlap |> filter(HIC_grp == "HIC", cause == "all-cause"),
  "HIC and UKPDS-era overlap: overall"
)
ukpds_tab <- bind_rows(res_ukpds$summary, res_hic_ukpds$summary)

ukpds_membership <- ukpds_overlap |>
  arrange(period_start, period_end, StudyID) |>
  transmute(StudyID, First_Author_Lastname, Year, cohort_id, cohort_abbr,
            Country, period_label, period_start, period_end,
            SMR, SMR_LCI, SMR_UCI, variance_source)

# -- 6. Subgroup analyses ----------------------------------------------------

## 6.1 Helpers ----
present_group_label <- function(x) {
  x_chr  <- as.character(x)
  x_norm <- toupper(trimws(x_chr))
  case_when(
    x_norm %in% c("TRUE",  "T", "Y", "YES", "1") ~ "Yes",
    x_norm %in% c("FALSE", "F", "N", "NO",  "0") ~ "No",
    TRUE ~ x_chr
  )
}

build_subgroup_rows <- function(data, group_col, group_label) {
  subgroup_data <- data |>
    mutate(grp = as.character(.data[[group_col]])) |>
    filter(!is.na(grp), grp != "", !tolower(grp) %in% c("nr", "unknown")) |>
    mutate(grp = factor(grp, levels = sort(unique(grp))))
  
  if (length(unique(subgroup_data$grp)) < 2) return(data.frame())
  
  combined_model <- tryCatch(
    metafor::rma.mv(
      yi = yi, V = vi_use, mods = ~ grp,
      random = ~ 1 | cohort_id / StudyID / estimate_id,
      data = subgroup_data, method = "REML"
    ),
    error = function(e) NULL
  )
  between_group_p <- if (!is.null(combined_model)) combined_model$QMp else NA_real_
  
  subgroup_rows <- data.frame()
  for (this_group in levels(subgroup_data$grp)) {
    one_group_data <- subgroup_data |> filter(grp == this_group)
    one_group_pool <- fit_pool(one_group_data, paste(group_label, this_group))
    subgroup_rows <- bind_rows(
      subgroup_rows,
      data.frame(
        Subgroup        = group_label,
        Group           = present_group_label(this_group),
        k_estimates     = one_group_pool$summary$k_estimates,
        k_studies       = one_group_pool$summary$k_studies,
        k_cohorts       = one_group_pool$summary$k_cohorts,
        SMR             = one_group_pool$summary$SMR,
        SMR_LCI         = one_group_pool$summary$SMR_LCI_robust,
        SMR_UCI         = one_group_pool$summary$SMR_UCI_robust,
        between_group_p = between_group_p
      )
    )
  }
  subgroup_rows
}

## 6.2 Build subgroups (region / income / OECD / UKPDS overlap) ----
subgroup_variables <- data.frame(
  variable = c("WHO_region", "WB_HIC", "OECD"),
  label    = c("WHO region", "World Bank HIC", "OECD")
)

subgroup_all <- data.frame()
for (i in seq_len(nrow(subgroup_variables))) {
  subgroup_all <- bind_rows(
    subgroup_all,
    build_subgroup_rows(d_overall,
                        subgroup_variables$variable[i],
                        subgroup_variables$label[i])
  )
}

subgroup_all <- bind_rows(
  subgroup_all,
  build_subgroup_rows(
    long |> filter(cause == "all-cause"),
    "ukpds_overlap_grp",
    paste0("UKPDS-era overlap (", UKPDS_START, "-", UKPDS_END, ")")
  ),
  build_subgroup_rows(
    long |> filter(cause == "all-cause"),
    "hic_ukpds_overlap_grp",
    paste0("HIC and UKPDS-era overlap (", UKPDS_START, "-", UKPDS_END, ")")
  )
)

# -- 7. Results workbook -----------------------------------------------------

wb <- createWorkbook()
header_style <- createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  fgFill = "#D9EAF7", border = "Bottom"
)

addWorksheet(wb, "README")
readme <- data.frame(
  Item = c("Input", "Analysis", "Population", "Primary estimate", "Primary CI",
           "Dependency handling", "UKPDS-era window", "Figures"),
  Value = c(
    basename(XLSX_IN),
    "Multilevel random-effects meta-analysis on log-SMR with random intercepts for cohort / study / estimate.",
    "All-persons (both-sex) all-cause SMRs only; no sex-specific analyses.",
    sprintf("%.2f from %d estimates, %d studies, %d cohorts",
            res_overall$summary$SMR, res_overall$summary$k_estimates,
            res_overall$summary$k_studies, res_overall$summary$k_cohorts),
    sprintf("Cluster-robust 95%% CI %.2f-%.2f by cohort",
            res_overall$summary$SMR_LCI_robust, res_overall$summary$SMR_UCI_robust),
    "Rows from the same cohort are not treated as independent.",
    paste0(UKPDS_START, "-", UKPDS_END, " overlap rule: period_start <= ",
           UKPDS_END, " and period_end >= ", UKPDS_START, "."),
    "PNG files in the outputs/figures folder; all forest plots sorted by study period."
  )
)
writeData(wb, "README", readme)
addStyle(wb, "README", header_style, rows = 1, cols = 1:ncol(readme), gridExpand = TRUE)
setColWidths(wb, "README", cols = 1:ncol(readme), widths = "auto")

result_tables <- list(
  "01_Main_pooled"       = main_tab,
  "02_Subgroups"         = subgroup_all,
  "03_HIC_pooled"        = hic_tab,
  "04_UKPDS_era_pool"    = ukpds_tab,
  "05_UKPDS_era_members" = ukpds_membership,
  "06_Analysis_rows"     = long
)
for (sheet_name in names(result_tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, result_tables[[sheet_name]])
  addStyle(wb, sheet_name, header_style, rows = 1,
           cols = 1:ncol(result_tables[[sheet_name]]), gridExpand = TRUE)
  freezePane(wb, sheet_name, firstActiveRow = 2)
  setColWidths(wb, sheet_name, cols = 1:ncol(result_tables[[sheet_name]]),
               widths = "auto")
}
saveWorkbook(wb, XLSX_OUT, overwrite = TRUE)
cat("Results workbook written:", XLSX_OUT, "\n")

# -- 7.1 Figures (PNG, all sorted by study period) ---------------------------

# Overall: show BOTH extra columns
draw_forest(res_overall, "Overall pooled SMR (sorted by study period)",
            "fig_forest_overall",
            sort_by_period = TRUE, show_hic = TRUE, show_ukpds = TRUE)

# HIC-only: show UKPDS-era overlap only (HIC is constant)
draw_forest(res_hic, "HIC-only pooled SMR (sorted by study period)",
            "fig_hic_forest",
            sort_by_period = TRUE, show_hic = FALSE, show_ukpds = TRUE)

# UKPDS-era overlap: show HIC only (UKPDS overlap is constant)
draw_forest(res_ukpds,
            paste0("UKPDS-era overlap (", UKPDS_START, "-", UKPDS_END,
                   ") pooled SMR (sorted by study period)"),
            "fig_ukpds_forest",
            sort_by_period = TRUE, show_hic = TRUE, show_ukpds = FALSE)

# HIC + UKPDS-era overlap: neither extra column (both would be constant "Yes")
draw_forest(res_hic_ukpds,
            paste0("HIC and UKPDS-era overlap (", UKPDS_START, "-", UKPDS_END,
                   ") pooled SMR (sorted by study period)"),
            "fig_hic_ukpds_forest",
            sort_by_period = TRUE, show_hic = FALSE, show_ukpds = FALSE)

# Subgroup summary forest
draw_combined_subgroup_plot(subgroup_all)

# -- 8. Console summary ------------------------------------------------------

cat("\n============================================================\n")
cat("Output folder:   ", ROOT, "\n")
cat("Results workbook:", basename(XLSX_OUT), "\n")
cat("PNG figures:     ", length(list.files(FIG_DIR, pattern = "\\.png$")), "\n")

cat("\nPrimary overall pooled SMR:\n")
cat("  SMR: ", sprintf("%.2f", res_overall$summary$SMR), "\n", sep = "")
cat("  Cluster-robust 95% CI: ",
    sprintf("%.2f-%.2f", res_overall$summary$SMR_LCI_robust,
            res_overall$summary$SMR_UCI_robust), "\n", sep = "")
cat("  Based on ", res_overall$summary$k_estimates, " estimates / ",
    res_overall$summary$k_studies, " studies / ",
    res_overall$summary$k_cohorts, " cohorts\n", sep = "")

cat("\nHIC and UKPDS-era overlap pooled SMR (main manuscript comparator):\n")
cat("  SMR: ", sprintf("%.2f", res_hic_ukpds$summary$SMR), "\n", sep = "")
cat("  Cluster-robust 95% CI: ",
    sprintf("%.2f-%.2f", res_hic_ukpds$summary$SMR_LCI_robust,
            res_hic_ukpds$summary$SMR_UCI_robust), "\n", sep = "")
cat("  Based on ", res_hic_ukpds$summary$k_estimates, " estimates / ",
    res_hic_ukpds$summary$k_studies, " studies / ",
    res_hic_ukpds$summary$k_cohorts, " cohorts\n", sep = "")
cat("============================================================\n")
