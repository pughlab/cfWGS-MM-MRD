################################################################################
##  Survival Analysis and Relapse Detection Sensitivity
##  
##  Purpose: 
##    Generate Kaplan-Meier survival curves stratified by MRD status at key
##    clinical timepoints, and calculate sensitivity of each assay for detecting
##    future relapse among frontline-treated multiple myeloma patients.
##    Includes head-to-head comparison of cfWGS, clinical assays (MFC, clonoSEQ),
##    and EasyM proteomic MRD.
##  
##  Main Analyses:
##    1. Progression-free survival (PFS) curves stratified by MRD status
##       at landmark timepoints (Post-ASCT, 1yr Maintenance, etc.)
##    2. Sensitivity calculation: % of patients who relapsed that were MRD+
##       at each timepoint, for each assay independently
##    3. Comparative sensitivity barplots (BM vs Blood-derived cfWGS models)
##    4. Optional: Time-window prediction analysis in non-frontline cohort
##  
##  Input Data Sources:
##    - all_patients_with_BM_and_blood_calls_updated6.rds
##      (from 3_1_Optimize_cfWGS_thresholds.R)
##    - EasyM_all_samples_with_optimized_calls.csv 
##      (from 3_1_A_Process_and_optimize_EasyM.R)
##    - Censor_dates_per_patient_for_PFS_updated.rds 
##      (clinical outcome tracking)
##  
##  Output Locations:
##    - Output_tables_2025/detection_progression_updated6/
##      └─ Kaplan-Meier curves (PNG, organized by timepoint)
##      └─ Sensitivity tables (CSV)
##      └─ Comparative barplots (Supp_6A, Supp_8A)
##      └─ prospective_timewindow_qc/
##         └─ prospective_Supplementary_Table_9_timewindow_results.xlsx
##            (the 16-row-per-sheet result used in Supplementary Table 9)
##  
##  Scripts this Script Depends On:
##    1. 3_1_Optimize_cfWGS_thresholds.R - cfWGS model optimization/thresholds
##    2. 3_1_A_Process_and_optimize_EasyM.R - EasyM processing and
##       isotype-specific reference-threshold calls
##    3. 2_0_Assemble_Table_With_All_Features.R - Feature integration
##
##  How to run:
##    Rscript Scripts_2025/Final_Scripts/4_1_Survival_Analysis.R
##
##  Manuscript outputs created/updated:
##    - Figure 3F: BM cfWGS time-to-event/survival panel.
##    - Figure 4E: blood cfWGS time-to-event/survival panel.
##    - Extended Data Figure 6A-K: BM cfWGS survival and relapse-detection
##      sensitivity panels.
##    - Extended Data Figure 8A-F: blood cfWGS survival,
##      relapse-detection sensitivity, and longitudinal panels.
##    - Supplementary Table 9: BM/blood time-window detection results.
##    - Paired cfWGS/EasyM landmark Cox-model source table used for the
##      manuscript Results hazard ratios.
##
##  Pipeline role:
##    Survival analyses are downstream of frozen cfWGS calls. The script tests
##    whether MRD status at clinically meaningful landmarks is associated with
##    progression-free survival and estimates how often relapsing patients were
##    detected by each assay before progression.
##
##  Analysis populations and units:
##    - Primary landmark KM curves: one earliest evaluable sample per frontline
##      patient at each landmark and assay.
##    - Landmark sensitivity summaries: one earliest evaluable sample per
##      relapsing frontline patient; denominators are assay-specific.
##    - Test-cohort time-window summaries: one row per evaluable sample and
##      prediction window, with patient counts reported separately.
##    - Longitudinal panels: one row per evaluable sample; a patient can
##      contribute multiple timepoints.
##  
##  Author: Dory Abelman
##  Updated: February 2026
##  
################################################################################
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

## ── 0. SETUP: Load Packages and Configure Paths ──────────────────────────────
##
##  This section:
##    - Loads all required R packages for survival analysis & plotting
##    - Defines input/output file paths
##    - Creates output directories if they don't exist
##
## ────────────────────────────────────────────────────────────────────────────

# Required packages for survival analysis, visualization, and data wrangling
library(tidyverse)       # dplyr, ggplot2, tidyr
library(lubridate)       # Date/time operations
library(survival)        # Survival objects, survfit()
library(survminer)       # ggsurvplot() for KM curves
library(broom)           # Tidy model outputs
library(patchwork)       # Combine multiple plots
library(tableone)        # Create summary tables (optional)
library(timeROC)         # Time-dependent ROC analysis (optional)
library(scales)          # Axis/percentage label formatting
library(glue)            # Inline text summaries used to check reported values
library(writexl)         # Simple multi-sheet Excel exports

# Shared manuscript-output helpers.
# These functions copy or save the exact figure/table files produced below into
# final_manuscript_objects/, organized by final manuscript figure/table label.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

.endpoint_helper <- file.path("Scripts_2025", "Final_Scripts", "next_event_endpoint_helpers.R")
if (!file.exists(.endpoint_helper)) {
  .endpoint_helper <- "next_event_endpoint_helpers.R"
}
source(.endpoint_helper)
rm(.endpoint_helper)

## ─────────────────────────────────────────────────────────────────────────────
## INPUT FILES: Clinical outcomes and cfWGS results
## ─────────────────────────────────────────────────────────────────────────────

# File 1: Clinical metadata - censor dates, relapse status per patient
#         Used to compute time-to-event and determine event status for survival analysis
final_tbl_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"

# File 1b: Patient-level last-known follow-up dates.
#          Used only for prospective non-frontline time-window evaluability.
#          This is intentionally separate from `final_tbl_rds`: for relapsed
#          patients, the PFS censor/event date is the first progression date,
#          while the last-known follow-up date may be much later.
patient_followup_rds <- "Exported_data_tables_clinical/patient_followup_dates_updated.rds"
patient_followup_csv <- "Exported_data_tables_clinical/patient_followup_dates_updated.csv"
latest_dates_csv     <- "Exported_data_tables_clinical/latest_dates_per_patient.csv"
latest_dates_updated_csv <- "Exported_data_tables_clinical/latest_dates_per_patient_updated.csv"
relapse_dates_full_rds <- "Exported_data_tables_clinical/Relapse_dates_full_updated.rds"

# File 2: Main data table with all cfWGS model calls and clinical MRD results
#         Each row = one sample (patient + timepoint)
#         Contains: MFC calls, clonoSEQ calls, cfWGS calls (multiple models)
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds"

## ─────────────────────────────────────────────────────────────────────────────
## OUTPUT DIRECTORY: All results saved here
## ─────────────────────────────────────────────────────────────────────────────

outdir <- "Output_tables_2025/detection_progression_updated6"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# SOURCE DATA DIRECTORY: Stores the row-level data used in figure panels
outdir_source_data <- "Output_tables_2025/Source_data"
dir.create(outdir_source_data, showWarnings = FALSE, recursive = TRUE)

# ═════════════════════════════════════════════════════════════════════════════
# FILE VERSIONING: Date-tagged working outputs retain separate runs
# ═════════════════════════════════════════════════════════════════════════════
# All output files include today's date in their names:
#   KM_assay_timepoint_updated_no_CI_2026-02-24.png
#   frontline_postASCT_sensitivity_2026-02-24.csv
# Date-tagged files retain separate runs. Fixed manuscript-output copies written
# later in the script are updated when the script is rerun.
# ═════════════════════════════════════════════════════════════════════════════

date_tag <- format(Sys.Date(), "%Y-%m-%d")

# ═════════════════════════════════════════════════════════════════════════════
# SOURCE DATA EXPORTS: Figure source data for manuscript
# ═════════════════════════════════════════════════════════════════════════════
# All figure source data saved to: Output_tables_2025/Source_data/
#
# KEY EXPORTS (exported immediately after data frame creation):
#   • Supp_6A_BM_sensitivity_barplot_source_data_YYYY-MM-DD.csv
#     └─ Source: sens_df_bm (BM subset sensitivity by assay × timepoint)
#
#   • Supp_8A_blood_sensitivity_barplot_source_data_YYYY-MM-DD.csv
#     └─ Source: sens_df_blood (Blood subset sensitivity by assay × timepoint)
#
#   • SuppFig8B_blood_HR_plot_source_data_YYYY-MM-DD.csv
#     └─ Source: hr_plot_df_blood (Blood-derived HR by landmark × assay)
#
#   • Supp_Figure_6B_BM_HR_plot_source_data_YYYY-MM-DD.csv
#     └─ Source: hr_plot_df_bm (BM-derived HR by landmark × assay)
#
# ADDITIONAL EXPORTS (analytical summaries):
#   • frontline_followup_summary_YYYY-MM-DD.csv (see line ~1046)
#   • frontline_postASCT_sensitivity_YYYY-MM-DD.csv (see line ~1175)
#   • frontline_1yr_sensitivity_YYYY-MM-DD.csv
#   • All progression metrics CSVs (see lines ~1900, 2330, 2763)
#
# NOTE: KM curve source data = filtered survival_df for each assay/timepoint
#       (generated dynamically in loop; not exported separately)
# ═════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("═", 80), "\n")
cat("SURVIVAL ANALYSIS: Frontline Cohort\n")
cat("Output directory:", outdir, "\n")
cat("Date tag for versioning:", date_tag, "\n")
cat(strrep("═", 80), "\n\n")

## ── 1. LOAD AND TIDY CORE TABLES ─────────────────────────────────────────────
##
##  This section:
##    - Loads clinical outcome metadata (follow-up dates, relapse status)
##    - Loads main cfWGS data table with all model predictions
##    - Filters to frontline cohort only (primary analysis population)
##    - Creates additional computed columns as needed
##
## ────────────────────────────────────────────────────────────────────────────

cat("1. Loading and preparing clinical data tables...\n")

# CLINICAL OUTCOMES TABLE
# Standardizes column names to lowercase and selects key columns:
#   - Patient: unique patient identifier
#   - baseline_date: treatment start date
#   - censor_date: last known follow-up date (event or censoring)
#   - relapsed: binary indicator (1=relapsed, 0=censored without relapse)
final_tbl <- readRDS(final_tbl_rds) %>%
  rename_with(
    tolower,
    any_of(c("Baseline_Date", "Censor_date", "Relapsed"))
  ) %>%
  transmute(
    Patient       = as.character(Patient),
    baseline_date = as.Date(baseline_date),
    censor_date   = as.Date(censor_date),
    relapsed      = as.integer(relapsed)
  )

cat(sprintf("  ✓ Clinical outcomes: %d patients loaded\n", n_distinct(final_tbl$Patient)))

load_patient_followup_dates <- function(rds_path, csv_path, fallback_latest_path) {
  if (file.exists(rds_path)) {
    followup <- readRDS(rds_path)
  } else if (file.exists(csv_path)) {
    followup <- readr::read_csv(csv_path, show_col_types = FALSE)
  } else if (file.exists(fallback_latest_path)) {
    followup <- readr::read_csv(fallback_latest_path, show_col_types = FALSE) %>%
      dplyr::rename(followup_end_date = latest_date) %>%
      dplyr::mutate(followup_source = fallback_latest_path)
  } else {
    message("  ! No patient-level follow-up table found; prospective QC will fall back to PFS event/censor dates.")
    return(tibble::tibble(
      Patient = character(),
      followup_end_date = as.Date(character()),
      followup_source = character()
    ))
  }

  if (!"followup_source" %in% names(followup)) {
    followup$followup_source <- "patient-level follow-up table"
  }

  followup %>%
    dplyr::transmute(
      Patient = as.character(Patient),
      followup_end_date = as.Date(followup_end_date),
      followup_source = as.character(followup_source)
    ) %>%
    dplyr::filter(!is.na(Patient), !is.na(followup_end_date)) %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      followup_end_date = max(followup_end_date),
      followup_source = paste(sort(unique(followup_source)), collapse = "; "),
      .groups = "drop"
    )
}

patient_followup_dates <- load_patient_followup_dates(
  rds_path = patient_followup_rds,
  csv_path = patient_followup_csv,
  fallback_latest_path = latest_dates_csv
)

cat(sprintf("  ✓ Patient follow-up dates: %d patients loaded for prospective QC\n",
            n_distinct(patient_followup_dates$Patient)))

# MAIN cfWGS + CLINICAL ASSAYS DATA TABLE
# Each row represents one sample (patient + timepoint)
# Contains predictions from all MRD models and timepoint classifications
dat_all <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    Cohort         = as.character(Cohort),
    timepoint_info = tolower(timepoint_info)  # standardize timepoint names for consistent grouping
  )

cat(sprintf("  ✓ Sample data: %d samples from %d patients loaded\n", 
            nrow(dat_all), n_distinct(dat_all$Patient)))

# FILTER TO FRONTLINE COHORT
# Primary landmark survival analyses use only patients in the frontline/
# induction-to-transplant cohort. The full all-cohort table is reloaded later
# for the separate non-frontline/test-cohort time-window analysis.
dat <- dat_all %>% filter(Cohort == "Frontline")

cat(sprintf("  ✓ Filtered to Frontline cohort: %d samples from %d patients\n", 
            nrow(dat), n_distinct(dat$Patient)))

# ---------------------------------------------------------------------------- #
# 2026-07-29 audit fix: REMOVED the "high sensitivity screen" decision rule.
#
# This step previously created
#   BM_zscore_only_detection_rate_screen_call =
#     as.integer(BM_zscore_only_detection_rate_prob >= 0.350)
# on both `dat` and `dat_all`, described in a comment as a "threshold optimized
# for sensitivity (detect 95% of relapsers)".
#
# Three problems. First, no derivation of 0.350 existed anywhere in the
# codebase. Second, 0.350 is not the training 95%-sensitivity threshold for this
# model, which is ~0.20-0.22 in cfWGS_model_metrics_fixed_95sens*.csv; it is
# combo_BM_prob = 0.35063, the 95%-sensitivity operating point of a different,
# legacy model specification, so a threshold derived for one model was being
# applied to another. Third, the comment implies the cutoff was tuned against
# the relapse outcome in this cohort, which would make every survival estimate
# computed at it circular.
#
# It is also obsolete. The screen was introduced because the blood models scored
# negative across the older dilution series; the current patient-derived
# dilution series are positive at the relevant levels, and the screen had
# already been dropped from the figures.
#
# Removing it changes no reported value: all manuscript survival results use the
# frozen Youden threshold (BM_zscore_only_detection_rate = 0.4216 in
# all_model_thresholds_v2_with_fragmentomics_restricted_cohorts.csv).
# ---------------------------------------------------------------------------- #

cat("  ✓ Data preparation complete\n\n")

## ── 1A. LOAD EasyM PROTEOMIC MRD DATA ─────────────────────────────────────────
##
##  This section:
##    - Loads EasyM M-protein measurements and prespecified isotype-specific calls
##    - Merges EasyM data by patient and timepoint
##    - Gracefully handles missing EasyM data with placeholder columns
##
##  EasyM Data Sources (from script 3_1_A):
##    - EasyM_value: continuous M-protein measure (%)
##    - EasyM_reference_threshold_binary: binary call (1=positive, 0=negative)
##                            using isotype-specific Rapid Novor reference thresholds
##    - EasyM_reference_threshold_call: character label of call for reporting
##
## ────────────────────────────────────────────────────────────────────────────

cat("1A. Loading EasyM proteomic MRD data...\n")

OUTPUT_DIR_EASYM <- "Output_EasyM_MRD_analysis_2025"
EasyM_file <- file.path(OUTPUT_DIR_EASYM, "EasyM_all_samples_with_optimized_calls.csv")

if (file.exists(EasyM_file)) {
  # Load EasyM predictions from script 3_1_A output
  EasyM_data <- readr::read_csv(EasyM_file, show_col_types = FALSE)
  cat(sprintf("  ✓ Loaded EasyM data: %d samples\n", nrow(EasyM_data)))
  
  # Merge EasyM data into main dataset
  # Join key: Patient (character) + Timepoint (character)
  # relationship = "many-to-one": multiple samples per patient, but one EasyM value per timepoint
  easym_join_tbl <- EasyM_data %>%
    select(Patient, Timepoint, EasyM_value, EasyM_reference_threshold_binary, EasyM_reference_threshold_call) %>%
    mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint))

  dat <- dat %>%
    mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)) %>%
    left_join(
      easym_join_tbl,
      by = c("Patient", "Timepoint"),
      relationship = "many-to-one"
    )

  dat_all <- dat_all %>%
    mutate(Patient = as.character(Patient), Timepoint = as.character(Timepoint)) %>%
    left_join(
      easym_join_tbl,
      by = c("Patient", "Timepoint"),
      relationship = "many-to-one"
    )
  
  n_easym_matched <- sum(!is.na(dat$EasyM_reference_threshold_binary))
  n_easym_matched_all <- sum(!is.na(dat_all$EasyM_reference_threshold_binary))
  cat(sprintf("  ✓ Merged EasyM data: %d frontline samples and %d all-cohort samples with EasyM calls\n",
              n_easym_matched, n_easym_matched_all))
  
} else {
  # If EasyM file not found, create placeholder columns (survival analysis will skip EasyM with NA filter)
  cat(sprintf("  ⚠ EasyM file not found at: %s\n", EasyM_file))
  cat("    Creating placeholder columns (EasyM will be skipped in downstream analyses)\n")
  
  dat <- dat %>%
    mutate(
      EasyM_value = NA_real_,
      EasyM_reference_threshold_binary = NA_integer_,
      EasyM_reference_threshold_call = NA_character_
    )
  dat_all <- dat_all %>%
    mutate(
      EasyM_value = NA_real_,
      EasyM_reference_threshold_binary = NA_integer_,
      EasyM_reference_threshold_call = NA_character_
    )
}

cat("  ✓ EasyM data preparation complete\n\n")

## ── 2. BUILD SURVIVAL DATA TABLE (sample-level) ───────────────────────────────
##
##  This section:
##    - Merges clinical outcomes (relapse, follow-up dates) with sample-level data
##    - Computes time-to-event: days from sample draw to relapse or censoring
##    - Creates binary event indicator (1=relapsed, 0=censored)
##    - Selects only columns needed for downstream survival analysis
##
##  Key computations:
##    - Time_to_event = censor_date - sample_date 
##      (how long from THIS sample until event/censoring)
##    - Relapsed_Binary = 1 if patient relapsed (regardless of when), 0 if censored
##
##  Output: survival_df
##    - One row per sample (patient × timepoint)
##    - Ready for stratified KM analysis by MRD status
##
## ────────────────────────────────────────────────────────────────────────────

cat("2. Building survival analysis table...\n")

next_event_endpoint_resources <- load_next_event_endpoint_resources(
  pfs_path = final_tbl_rds,
  relapse_dates_path = relapse_dates_full_rds,
  followup_rds_path = patient_followup_rds,
  followup_csv_path = patient_followup_csv,
  latest_dates_paths = c(latest_dates_updated_csv, latest_dates_csv),
  sample_data = dat_all,
  patient_col = "Patient",
  sample_date_col = "Date"
)

build_sample_survival_df <- function(sample_tbl) {
  sample_tbl %>%
    mutate(
      Patient = as.character(Patient),
      sample_date = as.Date(Date)
    ) %>%
    add_next_event_endpoint(
      endpoint_resources = next_event_endpoint_resources,
      sample_date_col = "sample_date",
      event_grace_days = 30L
    ) %>%
    mutate(
      censor_date = endpoint_date,
      relapsed = endpoint_status,
      Time_to_event = endpoint_days_from_sample,
      Relapsed_Binary = as.integer(endpoint_status)
    ) %>%
    select(
      Patient, Cohort, Timepoint, sample_date, censor_date, timepoint_info,
      Time_to_event, Relapsed_Binary,
      baseline_date, first_pfs_date, first_pfs_event, first_pfs_days_from_sample,
      next_progression_date, next_progression_days_raw, next_progression_days,
      endpoint_date, endpoint_status, endpoint_type, endpoint_source,
      endpoint_days_from_sample, n_prior_progressions_before_sample,
      latest_prior_progression_date, endpoint_uses_later_progression_after_first_pfs,
      endpoint_ignores_prior_progression, force_relapse_sample_day0,
      # Clinical MRD assays
      Flow_Binary, Adaptive_Binary, Rapid_Novor_Binary,
      Flow_pct_cells, Adaptive_Frequency,
      # PET imaging (if needed)
      PET_Binary,
      # cfWGS BM-derived models
      BM_zscore_only_detection_rate_call, BM_zscore_only_detection_rate_prob,
      # BM_zscore_only_detection_rate_screen_call removed 2026-07-29 (see audit
      # note at the data-preparation step): unprovenanced 0.350 operating point.
      # cfWGS blood-derived models (multiple variants)
      Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob,
      Blood_rate_only_call, Blood_rate_only_prob,
      Blood_zscore_only_detection_rate_call, Blood_zscore_only_detection_rate_prob,
      Blood_base_prob, Blood_base_call,
      Blood_plus_fragment_prob, Blood_plus_fragment_call,
      Blood_plus_fragment_min_prob, Blood_plus_fragment_min_call,
      # Fragmentomics models
      Fragmentomics_mean_coverage_only_prob, Fragmentomics_mean_coverage_only_call,
      # EasyM proteomic MRD
      EasyM_reference_threshold_binary, EasyM_value
    )
}

survival_df <- build_sample_survival_df(dat)

cat(sprintf("  ✓ Survival table created: %d samples from %d patients\n", 
            nrow(survival_df), n_distinct(survival_df$Patient)))

# QUICK DATA VALIDATION
# Verify the survival data looks reasonable
cat("\n  Data validation checks:\n")
cat(sprintf("    - Follow-up time: %.1f to %.1f days (median: %.1f days)\n",
            min(survival_df$Time_to_event, na.rm=TRUE),
            max(survival_df$Time_to_event, na.rm=TRUE),
            median(survival_df$Time_to_event, na.rm=TRUE)))

relapse_counts <- table(survival_df$Relapsed_Binary, useNA="ifany")
cat(sprintf("    - Relapse events: %d / %d (%.1f%%)\n",
            relapse_counts["1"], sum(!is.na(survival_df$Relapsed_Binary)),
            relapse_counts["1"] / sum(!is.na(survival_df$Relapsed_Binary)) * 100))

cat(sprintf("    - Timepoints represented: %s\n",
            paste(unique(survival_df$timepoint_info), collapse=", ")))

cat("  ✓ Survival table validation complete\n\n")

# Additive train+test cohort table for count-expansion sensitivity figures.
# The primary manuscript KM and longitudinal panels above/below retain the
# original frontline-only scope. This companion table keeps frontline/training
# plus non-frontline/test samples with interpretable outcome data, and excludes
# rows whose sample was collected well after the PFS event/censor anchor.
train_test_cohorts <- c("Frontline", "Non-frontline")

survival_df_train_test <- build_sample_survival_df(dat_all) %>%
  dplyr::filter(Cohort %in% train_test_cohorts) %>%
  dplyr::filter(
    !is.na(sample_date),
    !is.na(censor_date),
    !is.na(Relapsed_Binary),
    !is.na(Time_to_event),
    Time_to_event >= 0
  )

next_event_endpoint_audit <- survival_df_train_test %>%
  dplyr::filter(endpoint_ignores_prior_progression |
                  endpoint_uses_later_progression_after_first_pfs |
                  force_relapse_sample_day0) %>%
  dplyr::select(
    Patient, Cohort, Timepoint, sample_date, timepoint_info,
    first_pfs_date, first_pfs_days_from_sample,
    latest_prior_progression_date, n_prior_progressions_before_sample,
    next_progression_date, endpoint_date, endpoint_status,
    endpoint_type, endpoint_source, endpoint_days_from_sample,
    endpoint_uses_later_progression_after_first_pfs, force_relapse_sample_day0
  ) %>%
  dplyr::arrange(Patient, sample_date)

readr::write_csv(
  next_event_endpoint_audit,
  file.path(outdir_source_data, paste0("All_train_test_next_event_endpoint_audit_", date_tag, ".csv"))
)

train_test_evaluable_summary <- survival_df_train_test %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_patients = dplyr::n_distinct(Patient),
    n_frontline_samples = sum(Cohort == "Frontline", na.rm = TRUE),
    n_non_frontline_samples = sum(Cohort == "Non-frontline", na.rm = TRUE),
    n_relapsed_patients = dplyr::n_distinct(Patient[Relapsed_Binary == 1]),
    n_nonrelapsed_patients = dplyr::n_distinct(Patient[Relapsed_Binary == 0]),
    n_samples_ignoring_prior_progression = sum(endpoint_ignores_prior_progression, na.rm = TRUE),
    n_samples_using_later_progression_after_first_pfs =
      sum(endpoint_uses_later_progression_after_first_pfs, na.rm = TRUE),
    n_samples_censored_after_prior_progression =
      sum(endpoint_ignores_prior_progression & Relapsed_Binary == 0, na.rm = TRUE)
  )

readr::write_csv(
  train_test_evaluable_summary,
  file.path(outdir_source_data, paste0("All_train_test_outcome_available_summary_", date_tag, ".csv"))
)

train_test_assay_availability <- survival_df_train_test %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_patients = dplyr::n_distinct(Patient),
    n_bm_prob = sum(!is.na(BM_zscore_only_detection_rate_prob)),
    n_bm_call = sum(!is.na(BM_zscore_only_detection_rate_call)),
    n_blood_sites_prob = sum(!is.na(Blood_zscore_only_sites_prob)),
    n_blood_sites_call = sum(!is.na(Blood_zscore_only_sites_call)),
    n_blood_combined_call = sum(!is.na(Blood_plus_fragment_call)),
    n_flow_call = sum(!is.na(Flow_Binary)),
    n_clonoseq_call = sum(!is.na(Adaptive_Binary)),
    n_easym_call = sum(!is.na(EasyM_reference_threshold_binary)),
    .groups = "drop"
  )

readr::write_csv(
  train_test_assay_availability,
  file.path(outdir_source_data, paste0("All_train_test_assay_availability_", date_tag, ".csv"))
)

if (any(train_test_assay_availability$Cohort == "Non-frontline")) {
  nonfront_assay_total <- train_test_assay_availability %>%
    dplyr::filter(Cohort == "Non-frontline") %>%
    dplyr::select(
      dplyr::starts_with("n_bm_"),
      dplyr::starts_with("n_blood_"),
      n_flow_call,
      n_clonoseq_call,
      n_easym_call
    ) %>%
    as.matrix() %>%
    sum(na.rm = TRUE)
  if (identical(nonfront_assay_total, 0L) || identical(nonfront_assay_total, 0)) {
    message(
      "  ! Non-frontline/test rows have outcome data but no non-missing assay calls ",
      "for the KM/probability panels; additive outputs will not increase plotted counts."
    )
  }
}

cat(sprintf(
  "  ✓ Additive train+test outcome table: %d samples from %d patients (%d frontline, %d non-frontline samples)\n\n",
  train_test_evaluable_summary$n_samples,
  train_test_evaluable_summary$n_patients,
  train_test_evaluable_summary$n_frontline_samples,
  train_test_evaluable_summary$n_non_frontline_samples
))

## ── 3. KAPLAN-MEIER SURVIVAL ANALYSIS CONFIGURATION ──────────────────────────
##
##  This section configures settings for survival curve generation:
##    - Define which assays/models to include in KM analyses
##    - Specify timepoints to analyze
##    - Set visualization parameters (colors, DPI, minimum group size)
##    - Map human-readable labels for each timepoint
##
## ────────────────────────────────────────────────────────────────────────────

cat("3. Configuring Kaplan-Meier survival analyses...\n")

# LIST OF MRD TECHNOLOGIES/MODELS TO ANALYZE
# Each entry: column_name = "Display Label for Plots"
# Includes: clinical assays (MFC, clonoSEQ), cfWGS models, and EasyM
techs <- c(
  # Clinical MRD assays
  Flow_Binary        = "MFC",
  Adaptive_Binary    = "clonoSEQ",
  EasyM_reference_threshold_binary = "EasyM",
  # Bone Marrow-derived WGS mutations
  BM_zscore_only_detection_rate_call    = "cfWGS of BM-Derived Mutations (cVAF Model)",
  # "High Sensitivity" screen call removed 2026-07-29 (unprovenanced 0.350
  # threshold; see audit note at the data-preparation step).
  # Blood Plasma-derived WGS mutations
  Blood_zscore_only_sites_call = "cfWGS of cfDNA-Derived Mutations (Sites Model)",
  Blood_rate_only_call = "cfWGS of cfDNA-Derived Mutations (cVAF Model)",
  Blood_zscore_only_detection_rate_call = "cfWGS of cfDNA-Derived Mutations (cVAF Z-score Model)",
  Blood_plus_fragment_call = "cfWGS of cfDNA-Derived Mutations (Combined Model)"
  # NOTE: Fragmentomics models could be added here if desired
)

cat(sprintf("  - Analyzing %d MRD technologies/models\n", length(techs)))
cat(sprintf("    %s\n", paste("   ", names(techs), sep=" ✓ ", collapse="\n    "))[1:3])
cat("    ...\n")

# CLINICAL TIMEPOINTS TO ANALYZE
# Extract all unique timepoints from the data
# These represent key clinical decision points (diagnosis, post-transplant, maintenance, etc.)
tps <- unique(survival_df$timepoint_info)
cat(sprintf("  - Timepoints to analyze: %s\n", paste(tps, collapse=", ")))

# VISUALIZATION PARAMETERS
dpi_target <- 500  # Resolution for PNG output of KM curves

# GROUP SIZE THRESHOLDS
# Skip a KM curve when the complete landmark table has fewer than this many
# patients. The code separately requires both MRD groups to be present, but it
# does not require five patients in each group.
min_n <- 5

# COLOR PALETTE FOR KM CURVES
# MRD negative = black, MRD positive = red
pal_2 <- c("black", "red")

safe_logrank_pval_display <- function(surv_obj, data, group_col = "Group") {
  tryCatch(
    {
      group_values <- stats::na.omit(data[[group_col]])
      if (length(unique(group_values)) < 2) return(FALSE)

      lr <- survival::survdiff(
        stats::as.formula(paste("surv_obj ~", group_col)),
        data = data
      )
      if (!is.finite(lr$chisq) || length(lr$n) < 2 || any(lr$n == 0)) {
        return(FALSE)
      }

      p_val <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
      if (!is.finite(p_val)) return(FALSE)

      p_text <- dplyr::case_when(
        p_val < 1e-4 ~ "p < 1e-4",
        p_val < 1e-3 ~ "p < 0.001",
        p_val < 1e-2 ~ "p < 0.01",
        TRUE ~ sprintf("p = %.2f", p_val)
      )
      p_text
    },
    error = function(e) FALSE
  )
}

# HUMAN-READABLE LABELS FOR TIMEPOINTS
# Maps internal names (as stored in timepoint_info column) to plot labels
tp_labels <- c(
  `diagnosis`          = "Diagnosis",
  `post_transplant`    = "Post-ASCT",
  `1yr maintenance`    = "One-Year Maintenance", 
  `post_induction`     = "Post‑Induction",
  `post_asct`          = "Post-ASCT",
  `maintenance`        = "Maintenance",
  `1yr maint`          = "One-Year Maintenance"
)

# PREPARE BASELINE DIAGNOSIS DATES
# Compute earliest sample date per patient (used for relative time calculations if needed)
dx_tbl <- survival_df %>%
  group_by(Patient) %>%
  summarise(
    diagnosis_date = suppressWarnings(min(sample_date[timepoint_info == "diagnosis"], na.rm = TRUE)),
    .groups = "drop"
  )

cat("  ✓ Configuration complete\n\n")

# MANUSCRIPT PANEL MAP FOR KAPLAN-MEIER CURVES
# The KM loop below creates many exploratory and QC curves. Only the
# combinations listed here are active manuscript panels. Mapping by the
# internal timepoint and assay column is deliberate: it is more stable than
# matching on display text or filenames, which have changed across revisions.
km_manuscript_artifacts <- tibble::tribble(
  ~timepoint_info, ~assay_variable, ~artifact_id, ~role, ~description,
  "1yr maintenance", "BM_zscore_only_detection_rate_call", "FIG3F", "figure_panel_png", "Figure 3F: one-year maintenance PFS by BM-derived cfWGS MRD status.",
  "1yr maintenance", "Blood_zscore_only_sites_call", "FIG4E", "figure_panel_png", "Figure 4E: one-year maintenance PFS by cfDNA-derived cfWGS MRD status.",
  "1yr maintenance", "Flow_Binary", "EDFIG6C", "figure_panel_png", "Extended Data Figure 6C: one-year maintenance PFS by MFC MRD status.",
  "1yr maintenance", "Adaptive_Binary", "EDFIG6D", "figure_panel_png", "Extended Data Figure 6D: one-year maintenance PFS by clonoSEQ MRD status.",
  "post_transplant", "BM_zscore_only_detection_rate_call", "EDFIG6E", "figure_panel_png", "Extended Data Figure 6E: post-ASCT PFS by BM-derived cfWGS MRD status.",
  "post_transplant", "Flow_Binary", "EDFIG6F", "figure_panel_png", "Extended Data Figure 6F: post-ASCT PFS by MFC MRD status.",
  "1yr maintenance", "EasyM_reference_threshold_binary", "EDFIG6G", "figure_panel_png", "Extended Data Figure 6G: one-year maintenance PFS by EasyM MRD status.",
  "post_transplant", "EasyM_reference_threshold_binary", "EDFIG6H", "figure_panel_png", "Extended Data Figure 6H: post-ASCT PFS by EasyM MRD status.",
  "post_transplant", "Blood_zscore_only_sites_call", "EDFIG8C", "figure_panel_png", "Extended Data Figure 8C: post-ASCT PFS by the cfDNA sites model.",
  "post_transplant", "Blood_rate_only_call", "EDFIG8C", "alternate_cvaf_model_figure_panel_png", "Extended Data Figure 8C comparison: post-ASCT PFS by the cfDNA raw-cVAF model.",
  "post_transplant", "Blood_zscore_only_detection_rate_call", "EDFIG8C", "alternate_cvaf_zscore_model_figure_panel_png", "Extended Data Figure 8C comparison: post-ASCT PFS by the cfDNA cVAF Z-score model.",
  "post_transplant", "Blood_plus_fragment_call", "EDFIG8C", "alternate_combined_model_figure_panel_png", "Extended Data Figure 8C comparison: post-ASCT PFS by the combined cfDNA/fragmentomics model."
)


## ── 4. GENERATE KAPLAN-MEIER CURVES: Timepoint × Technology Analysis ────────
##
##  This section contains nested loops that:
##    1. Iterate through each clinical timepoint
##    2. For each timepoint, iterate through each MRD technology
##    3. For each combination, generate a stratified KM curve:
##       - Subjects split into MRD+ vs MRD- groups based on assay result
##       - Curve shows PFS probability over follow-up time
##       - Risk table shows number at risk per group at each time
##    4. Saves the PNG at 500 dpi
##
##  Nested Loop Structure:
##    for (timepoint in all timepoints)
##      for (assay in all technologies)
##        Generate KM curve for that timepoint + assay combo
##        Save to: outdir/timepoint/KM_assay_timepoint.png
##
##  Output Directory Organization:
##    detection_progression_updated6/
##    ├── diagnosis/
##    │   ├── KM_MFC_Diagnosis_updated_no_CI.png
##    │   ├── KM_clonoSEQ_Diagnosis_updated_no_CI.png
##    │   ├── KM_EasyM_Diagnosis_updated_no_CI.png
##    │   └── ...
##    ├── post_transplant/
##    │   ├── KM_MFC_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_clonoSEQ_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_EasyM_Post‑ASCT_updated_no_CI.png
##    │   ├── KM_cfWGS of BM-Derived Mutations (cVAF Model)_Post‑ASCT_updated_no_CI.png
##    │   └── ...
##    └── 1yr_maintenance/
##        └── ...
##
##  Notes on Curve Generation:
##    - Skips a timepoint/assay combination with fewer than min_n (5) patients
##      in total or with only one observed MRD group
##    - Uses Kaplan-Meier non-parametric estimator
##    - Includes log-rank p-value testing MRD+ vs MRD- groups
##    - Risk table shows N at risk below x-axis at selected timepoints
##    - X-axis = months from MRD assessment, Y-axis = PFS probability
##
## ────────────────────────────────────────────────────────────────────────────

cat("4. Generating Kaplan-Meier survival curves...\n")
cat("   (This may take a minute - creating curves for all timepoint×assay combinations)\n\n")

# Loop 1: Iterate through each timepoint
for(tp in tps) {
  # Get nice label for this timepoint
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  # Create subdirectory for this timepoint's curves
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Loop 2: Iterate through each MRD technology/model
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_no_CI_", date_tag, ".png"))
    fname_manuscript <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_no_CI.png"))
    
    # ─────────────────────────────────────────────────────────────────────────
    # PREPARE DATA FOR THIS TIMEPOINT × ASSAY
    # Filter criteria:
    #   - timepoint_info == tp: Only samples from this timepoint
    #   - !is.na(Time_to_event): Remove if follow-up time missing
    #   - !is.na(Relapsed_Binary): Remove if relapse status unknown
    #   - !is.na(.data[[var]]): Remove if MRD assay result missing
    # arrange() + group_by() + slice(1) ensures we keep only the FIRST sample per patient
    # at that timepoint (in case multiple draws on same date)
    # ─────────────────────────────────────────────────────────────────────────
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      # Create stratification group: MRD+ (assay value = 1) vs MRD- (assay value = 0)
      # This is the key predictor variable for the KM curves
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    # Convert time-to-event from days to months (multiply by 12, divide by 365)
    # Using 30.44 = average days per month (365.25 / 12) for accuracy
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # ─────────────────────────────────────────────────────────────────────────
    # VALIDATION: Skip assays with insufficient data
    # ─────────────────────────────────────────────────────────────────────────
    
    # Skip if total sample size is too small (min_n = 5)
    if(nrow(df_sub) < min_n) next
    
    # Skip if we don't have both MRD+ and MRD- groups (need both for comparison)
    if(n_distinct(df_sub$Group) < 2) next
    
    # ─────────────────────────────────────────────────────────────────────────
    # FIT KAPLAN-MEIER MODEL AND GENERATE SURVIVAL PLOT
    # ─────────────────────────────────────────────────────────────────────────
    
    # Surv() creates survival object with time and event indicator
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    
    # survfit() fits KM curves stratified by MRD status (one curve per group)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    # ggsurvplot generates the KM plot with:
    #   - Log-rank p-value comparing groups
    #   - Risk table showing number of subjects at risk over time
    #   - Customized colors, labels, and formatting
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = safe_logrank_pval_display(surv_obj, df_sub),
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = FALSE,
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )

    # MANUSCRIPT OUTPUT: stable no-date copy of active KM panels
    # The dated file above is retained as the run archive. The no-date file
    # below is the stable manuscript source file used for Figure 3F, Figure 4E,
    # and Extended Data Figures 6C-H/8C. Keeping both makes reruns auditable
    # while giving the manuscript pipeline stable filenames.
    ggsave(
      filename = fname_manuscript,
      plot     = combined,
      width    = 7,
      height   = 7,
      dpi      = dpi_target
    )

    km_artifact <- km_manuscript_artifacts %>%
      filter(timepoint_info == tp, assay_variable == var)

    if (nrow(km_artifact) == 1) {
      ms_copy_artifact(
        source_path = fname_manuscript,
        artifact_id = km_artifact$artifact_id,
        role = km_artifact$role,
        description = km_artifact$description,
        script_name = "4_1_Survival_Analysis.R"
      )
    }
  }
}

## ── 4S. SUSTAINED-MRD SENSITIVITY: POST-ASCT/MAINTENANCE ASSESSMENTS ────────
#
# The one-year landmark remains the primary Figure 3F/4E analysis. This
# additive sensitivity analysis uses post-ASCT as the first call when available;
# otherwise it uses the first evaluable maintenance call. The second call is
# the first later evaluable maintenance assessment within two years. Follow-up
# starts at the second assessment so future information is never used to define
# a group before time zero.
sustained_mrd_dir <- file.path(outdir, "sustained_mrd_first_two_maintenance")
dir.create(sustained_mrd_dir, recursive = TRUE, showWarnings = FALSE)

sustained_mrd_specs <- tibble::tribble(
  ~assay_variable, ~assay_label, ~artifact_id,
  "BM_zscore_only_detection_rate_call", "cfWGS using baseline BM-derived mutations", "FIG3F",
  "Blood_zscore_only_sites_call", "cfWGS using baseline cfDNA-derived mutations", "FIG4E"
)

build_sustained_mrd_df <- function(data, assay_variable, max_gap_days = 730L) {
  required <- c(
    "Patient", "Cohort", "timepoint_info", "sample_date", "Time_to_event",
    "Relapsed_Binary", assay_variable
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Sustained-MRD input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  evaluable <- data %>%
    filter(
      Cohort == "Frontline",
      str_to_lower(timepoint_info) == "post_transplant" |
        str_detect(str_to_lower(timepoint_info), "maint"),
      !is.na(.data[[assay_variable]]),
      .data[[assay_variable]] %in% c(0, 1),
      !is.na(sample_date)
    ) %>%
    arrange(Patient, sample_date) %>%
    distinct(Patient, sample_date, .keep_all = TRUE) %>%
    group_by(Patient) %>%
    group_modify(~ {
      post_asct <- .x %>%
        filter(str_to_lower(timepoint_info) == "post_transplant") %>%
        slice_head(n = 1)
      anchor <- if (nrow(post_asct) == 1) post_asct else slice_head(.x, n = 1)
      later_maintenance <- .x %>%
        filter(
          sample_date > anchor$sample_date[[1]],
          str_detect(str_to_lower(timepoint_info), "maint")
        ) %>%
        slice_head(n = 1)
      bind_rows(anchor, later_maintenance) %>%
        distinct(sample_date, .keep_all = TRUE)
    }) %>%
    mutate(assessment_number = row_number()) %>%
    ungroup()

  patient_summary <- evaluable %>%
    group_by(Patient) %>%
    summarise(
      n_evaluable_assessments = n(),
      first_date = first(sample_date),
      second_date = nth(sample_date, 2),
      first_timepoint = first(timepoint_info),
      second_timepoint = nth(timepoint_info, 2),
      first_call = first(.data[[assay_variable]]),
      second_call = nth(.data[[assay_variable]], 2),
      gap_days = as.integer(second_date - first_date),
      Time_to_event = nth(Time_to_event, 2),
      Relapsed_Binary = nth(Relapsed_Binary, 2),
      .groups = "drop"
    ) %>%
    mutate(
      exclusion_reason = case_when(
        n_evaluable_assessments < 2 ~ "fewer than two evaluable post-ASCT/maintenance assessments",
        is.na(gap_days) | gap_days < 0 ~ "invalid assessment-date ordering",
        gap_days > max_gap_days ~ "second evaluable maintenance assessment occurs after two years",
        is.na(Time_to_event) | Time_to_event < 0 | is.na(Relapsed_Binary) ~ "outcome unavailable at second assessment",
        TRUE ~ NA_character_
      ),
      trajectory = case_when(
        first_call == 0 & second_call == 0 ~ "Sustained MRD-negative",
        first_call == 0 & second_call == 1 ~ "MRD negative-to-positive",
        first_call == 1 & second_call == 1 ~ "Persistently MRD-positive",
        first_call == 1 & second_call == 0 ~ "MRD positive-to-negative",
        TRUE ~ NA_character_
      )
    )

  patient_summary
}

for (i in seq_len(nrow(sustained_mrd_specs))) {
  spec <- sustained_mrd_specs[i, ]
  sustained_all <- build_sustained_mrd_df(survival_df, spec$assay_variable)
  safe_assay <- gsub("[^A-Za-z0-9]+", "_", spec$assay_variable)
  write_csv(
    sustained_all,
    file.path(sustained_mrd_dir, paste0(safe_assay, "_patient_classification_audit.csv"))
  )

  sustained_plot_df <- sustained_all %>%
    filter(
      is.na(exclusion_reason),
      trajectory %in% c(
        "Sustained MRD-negative",
        "MRD negative-to-positive",
        "Persistently MRD-positive"
      )
    ) %>%
    mutate(
      trajectory = factor(
        trajectory,
        levels = c(
          "Sustained MRD-negative",
          "MRD negative-to-positive",
          "Persistently MRD-positive"
        )
      ),
      Group = trajectory,
      time_months = Time_to_event / 30.44
    )

  write_csv(
    sustained_plot_df,
    file.path(sustained_mrd_dir, paste0(safe_assay, "_KM_source_data.csv"))
  )

  group_counts <- sustained_plot_df %>% count(trajectory, name = "n")
  if (nrow(group_counts) < 2 || any(group_counts$n < 2)) {
    warning("Skipping sustained-MRD KM for ", spec$assay_variable,
            ": fewer than two groups or a group has <2 patients.")
    next
  }

  sustained_surv <- Surv(sustained_plot_df$time_months, sustained_plot_df$Relapsed_Binary)
  sustained_fit <- survfit(sustained_surv ~ Group, data = sustained_plot_df)
  sustained_km <- ggsurvplot(
    sustained_fit,
    data = sustained_plot_df,
    pval = safe_logrank_pval_display(sustained_surv, sustained_plot_df, "Group"),
    conf.int = FALSE,
    risk.table = TRUE,
    break.time.by = 12,
    palette = c("black", "#E69F00", "#D55E00"),
    legend.title = "Post-ASCT/maintenance calls",
    legend.labs = c(
      "Sustained MRD-negative",
      "MRD negative-to-positive",
      "Persistently MRD-positive"
    ),
    xlab = "Time since second evaluable assessment (months)",
    ylab = "Progression-free survival",
    title = str_wrap(paste0("PFS by sustained MRD status: ", spec$assay_label), 55),
    ggtheme = theme_classic(base_size = 12) +
      theme(
        legend.position = "top",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)
      )
  )
  sustained_combined <- ggarrange(
    sustained_km$plot + theme(axis.title.x = element_blank()),
    sustained_km$table + theme(axis.title.y = element_blank()),
    ncol = 1,
    heights = c(3, 1)
  )
  sustained_path <- file.path(
    sustained_mrd_dir,
    paste0("KM_", safe_assay, "_sustained_first_two_maintenance.png")
  )
  ggsave(
    sustained_path,
    sustained_combined,
    width = 8.5,
    height = 7,
    dpi = dpi_target,
    bg = "white"
  )
  ms_copy_artifact(
    source_path = sustained_path,
    artifact_id = spec$artifact_id,
    role = "sustained_first_two_maintenance_sensitivity_png",
    description = paste0(
      spec$artifact_id,
      " sensitivity analysis using post-ASCT when available followed by the first later evaluable maintenance call within two years; follow-up begins at the second assessment."
    ),
    script_name = "4_1_Survival_Analysis.R"
  )
}

if (identical(Sys.getenv("CFWGS_SUSTAINED_ONLY", unset = "0"), "1")) {
  message("CFWGS_SUSTAINED_ONLY=1: stopping after sustained-MRD outputs.")
  quit(save = "no", status = 0)
}

## ── 4A. ADDITIVE KM CURVES: Frontline/Training + Non-Frontline/Test ─────────
##
##  Goal:
##    Create separate exploratory/supporting KM versions that use all samples in
##    the training/frontline and test/non-frontline cohorts when patient outcome
##    information and assay calls are available.
##
##  Scope and safeguards:
##    - Original manuscript KM panels above are unchanged.
##    - Cohort scope is explicitly limited to Frontline + Non-frontline rows.
##    - Rows without outcome dates/status, assay calls, or interpretable
##      non-negative follow-up from sample to event/censor are excluded.
##    - The first sample per patient per landmark is retained, matching the
##      original KM landmark logic and avoiding repeated patient contributions.
##
##  Outputs:
##    - PNGs in Output_tables_2025/detection_progression_updated6/
##      all_train_test_outcome_available/<timepoint>/
##    - Source-data CSVs beside manuscript source data.
##
cat("4A. Generating additive train+test outcome-available KM curves...\n\n")

outdir_train_test_km <- file.path(outdir, "all_train_test_outcome_available")
dir.create(outdir_train_test_km, recursive = TRUE, showWarnings = FALSE)

km_train_test_summary <- list()
tps_train_test <- unique(survival_df_train_test$timepoint_info)
techs_train_test <- techs[names(techs) %in% names(survival_df_train_test)]

for (tp in tps_train_test) {
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)

  tp_dir <- file.path(outdir_train_test_km, gsub("\\s+", "_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)

  for (var in names(techs_train_test)) {
    assay_lab <- techs_train_test[[var]]
    file_stub <- paste0(
      "KM_", assay_lab, "_", nice_tp,
      "_all_train_test_outcome_available"
    )
    fname <- file.path(tp_dir, paste0(file_stub, "_", date_tag, ".png"))
    fname_stable <- file.path(tp_dir, paste0(file_stub, ".png"))

    df_sub <- survival_df_train_test %>%
      dplyr::filter(
        timepoint_info == tp,
        !is.na(Time_to_event),
        Time_to_event >= 0,
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      dplyr::arrange(Patient, sample_date) %>%
      dplyr::group_by(Patient) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Time_to_event_months = pmax(Time_to_event, 0) / 30.44,
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative", "Positive")
        )
      )

    km_train_test_summary[[length(km_train_test_summary) + 1]] <- df_sub %>%
      dplyr::count(Cohort, Group, name = "n_samples") %>%
      dplyr::mutate(
        timepoint_info = tp,
        assay_variable = var,
        assay_label = assay_lab,
        n_patients_total = dplyr::n_distinct(df_sub$Patient),
        n_events_total = sum(df_sub$Relapsed_Binary == 1, na.rm = TRUE)
      )

    if (nrow(df_sub) < min_n || dplyr::n_distinct(df_sub$Group) < 2) {
      next
    }

    readr::write_csv(
      df_sub,
      file.path(
        outdir_source_data,
        paste0(
          "All_train_test_KM_source_data_",
          gsub("[^A-Za-z0-9]+", "_", tp), "_",
          gsub("[^A-Za-z0-9]+", "_", var), "_",
          date_tag, ".csv"
        )
      )
    )

    surv_obj <- survival::Surv(df_sub$Time_to_event_months, df_sub$Relapsed_Binary)
    fit <- survival::survfit(surv_obj ~ Group, data = df_sub)
    pval_display <- tryCatch(
      {
        lr <- survival::survdiff(surv_obj ~ Group, data = df_sub)
        if (!is.finite(lr$chisq) || length(lr$n) < 2) {
          FALSE
        } else {
          p_val <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
          p_text <- if (is.finite(p_val) && p_val < 0.001) {
            "<0.001"
          } else if (is.finite(p_val)) {
            sprintf("%.3f", p_val)
          } else {
            NA_character_
          }
          if (!is.na(p_text)) paste0("Log-rank p = ", p_text) else FALSE
        }
      },
      error = function(e) FALSE
    )

    km <- survminer::ggsurvplot(
      fit,
      data = df_sub,
      pval = pval_display,
      break.time.by = 12,
      conf.int = FALSE,
      risk.table = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),
      palette = pal_2,
      legend.title = "MRD status",
      legend.labs = c("MRD-", "MRD+"),
      xlab = "Time since MRD assessment (months)",
      ylab = "Progression-free survival",
      title = stringr::str_wrap(
        paste0(
          "Next-event-free survival stratified by ", assay_lab, " at ", nice_tp,
          "\nTraining Cohort + Test Cohort outcome-available samples"
        ),
        width = 54
      ),
      risk.table.height = 0.25,
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          legend.position = "top",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 14)
        )
    )

    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0, face = "plain")
      )

    km$plot <- km$plot +
      theme(axis.title.x = element_blank())

    combined <- ggarrange(
      km$plot,
      km$table,
      ncol = 1,
      heights = c(3, 1)
    )

    ggsave(fname, combined, width = 7, height = 7, dpi = dpi_target)
    ggsave(fname_stable, combined, width = 7, height = 7, dpi = dpi_target)

    km_artifact <- km_manuscript_artifacts %>%
      dplyr::filter(timepoint_info == tp, assay_variable == var)

    if (nrow(km_artifact) == 1) {
      ms_copy_artifact(
        source_path = fname_stable,
        artifact_id = km_artifact$artifact_id,
        role = "all_samples_figure_panel_png",
  description = paste0(
    "All-evaluable training/test outcome-available version of ",
    km_artifact$description,
    " Uses first known progression after the MRD assessment, or censoring if no later progression is available."
  ),
        script_name = "4_1_Survival_Analysis.R"
      )
    }
  }
}

km_train_test_summary <- dplyr::bind_rows(km_train_test_summary)
readr::write_csv(
  km_train_test_summary,
  file.path(outdir_source_data, paste0("All_train_test_KM_evaluable_counts_", date_tag, ".csv"))
)

cat(sprintf(
  "  ✓ Additive train+test KM summary written for %d timepoint-assay combinations\n\n",
  dplyr::n_distinct(paste(km_train_test_summary$timepoint_info, km_train_test_summary$assay_variable))
))


## ── 4B. PAIRED cfWGS/EasyM LANDMARK HAZARD RATIOS ─────────────────────────────
##
## Goal:
##   Summarize the patient-matched cfWGS and EasyM landmark cohorts used for
##   the manuscript Results hazard ratios. The patient-level inputs are the
##   EasyM source tables written immediately above by this script.
##
## Unit of analysis and models:
##   - One patient per landmark (18 patients/8 events at one-year maintenance;
##     16 patients/7 events post-ASCT).
##   - Separate univariable Cox models for the binary cfWGS call, the
##     isotype-specific EasyM call, cfWGS probability, and EasyM residual value.
##   - Continuous predictors are standardized within landmark and reported per
##     one standard-deviation increase; ties use the Efron approximation.
##
## Safeguards:
##   Expected denominators, event counts, and the four manuscript binary-call
##   hazard ratios are checked before the source table is written. These checks
##   prevent an upstream cohort change from silently altering reported results.

paired_easym_landmarks <- tibble::tribble(
  ~landmark, ~source_file, ~expected_n, ~expected_events,
  "One-year maintenance",
  paste0(
    "All_train_test_KM_source_data_1yr_maintenance_",
    "EasyM_reference_threshold_binary_", date_tag, ".csv"
  ),
  18L, 8L,
  "Post-ASCT",
  paste0(
    "All_train_test_KM_source_data_post_transplant_",
    "EasyM_reference_threshold_binary_", date_tag, ".csv"
  ),
  16L, 7L
)

paired_easym_predictors <- tibble::tribble(
  ~predictor, ~column, ~scale_continuous, ~definition,
  "BM-informed cfWGS call",
  "BM_zscore_only_detection_rate_call",
  FALSE,
  "Locked BM-informed cVAF-model call",
  "EasyM reference-threshold call",
  "EasyM_reference_threshold_binary",
  FALSE,
  paste(
    "Isotype-specific Rapid Novor reference threshold:",
    "negative if IgG <=1% or IgA/light-chain <=0.05% of baseline"
  ),
  "BM-informed cfWGS probability per 1 SD",
  "BM_zscore_only_detection_rate_prob",
  TRUE,
  "Locked BM-informed cVAF-model probability, standardized within landmark",
  "EasyM residual value per 1 SD",
  "EasyM_value",
  TRUE,
  "Residual monoclonal-protein percentage, standardized within landmark"
)

fit_paired_easym_landmark_model <- function(data,
                                            predictor_row,
                                            landmark,
                                            source_file) {
  column <- predictor_row$column[[1]]
  required <- c("Time_to_event", "Relapsed_Binary", column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        source_file,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  analysis <- data %>%
    dplyr::select(dplyr::all_of(required)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x)))

  if (nrow(analysis) < 5L || sum(analysis$Relapsed_Binary) < 2L) {
    stop(
      sprintf("Insufficient evaluable data for %s at %s", column, landmark),
      call. = FALSE
    )
  }
  if (dplyr::n_distinct(analysis[[column]]) < 2L) {
    stop(
      sprintf("No predictor variation for %s at %s", column, landmark),
      call. = FALSE
    )
  }

  model_column <- column
  if (isTRUE(predictor_row$scale_continuous[[1]])) {
    model_column <- paste0(column, "_z")
    analysis[[model_column]] <- as.numeric(scale(analysis[[column]]))
  }

  fit <- survival::coxph(
    stats::as.formula(
      sprintf("survival::Surv(Time_to_event, Relapsed_Binary) ~ %s", model_column)
    ),
    data = analysis,
    ties = "efron"
  )
  result <- summary(fit)

  tibble::tibble(
    Landmark = landmark,
    Predictor = predictor_row$predictor[[1]],
    Source_column = column,
    N = nrow(analysis),
    Events = sum(analysis$Relapsed_Binary),
    HR = unname(exp(stats::coef(fit))[[1]]),
    CI_low = unname(result$conf.int[1, "lower .95"]),
    CI_high = unname(result$conf.int[1, "upper .95"]),
    P_value = unname(result$coefficients[1, "Pr(>|z|)"]),
    Definition = predictor_row$definition[[1]],
    Source_file = source_file
  )
}

paired_easym_results <- list()
for (landmark_index in seq_len(nrow(paired_easym_landmarks))) {
  landmark_row <- paired_easym_landmarks[landmark_index, ]
  source_path <- file.path(outdir_source_data, landmark_row$source_file[[1]])
  if (!file.exists(source_path)) {
    stop(sprintf("Missing current EasyM landmark source: %s", source_path), call. = FALSE)
  }

  landmark_data <- readr::read_csv(source_path, show_col_types = FALSE)
  if (
    nrow(landmark_data) != landmark_row$expected_n[[1]] ||
      sum(landmark_data$Relapsed_Binary, na.rm = TRUE) != landmark_row$expected_events[[1]]
  ) {
    stop(
      sprintf(
        "%s denominator drift: observed %d rows/%d events; expected %d/%d",
        landmark_row$landmark[[1]],
        nrow(landmark_data),
        sum(landmark_data$Relapsed_Binary, na.rm = TRUE),
        landmark_row$expected_n[[1]],
        landmark_row$expected_events[[1]]
      ),
      call. = FALSE
    )
  }

  for (predictor_index in seq_len(nrow(paired_easym_predictors))) {
    paired_easym_results[[length(paired_easym_results) + 1L]] <-
      fit_paired_easym_landmark_model(
        landmark_data,
        paired_easym_predictors[predictor_index, ],
        landmark_row$landmark[[1]],
        landmark_row$source_file[[1]]
      )
  }
}

paired_easym_hr_table <- dplyr::bind_rows(paired_easym_results)

paired_easym_expected_rounded <- tibble::tribble(
  ~Landmark, ~Predictor, ~HR_expected,
  "One-year maintenance", "BM-informed cfWGS call", 6.08,
  "One-year maintenance", "EasyM reference-threshold call", 4.63,
  "Post-ASCT", "BM-informed cfWGS call", 3.51,
  "Post-ASCT", "EasyM reference-threshold call", 1.29
)
paired_easym_verification <- paired_easym_hr_table %>%
  dplyr::inner_join(
    paired_easym_expected_rounded,
    by = c("Landmark", "Predictor")
  ) %>%
  dplyr::mutate(matches = round(HR, 2) == HR_expected)

if (
  nrow(paired_easym_verification) != nrow(paired_easym_expected_rounded) ||
    any(!paired_easym_verification$matches)
) {
  stop("Paired EasyM landmark hazard-ratio verification failed.", call. = FALSE)
}

paired_easym_hr_output_path <- file.path(
  outdir_source_data,
  paste0("Paired_EasyM_landmark_HR_source_data_", date_tag, ".csv")
)
readr::write_csv(paired_easym_hr_table, paired_easym_hr_output_path, na = "NA")
message("Paired cfWGS/EasyM landmark HR source written to: ", paired_easym_hr_output_path)


## Optional KM confidence-interval sensitivity exports
##
## The final manuscript KM panels use the no-CI exports above. The original
## working script also generated CI variants while evaluating plot style. Those
## files are useful for audit/QC but are not mapped to final manuscript panels,
## so they are disabled during the standard command-line manuscript run.
export_optional_km_ci_panels <- FALSE

if (isTRUE(export_optional_km_ci_panels)) {
# 5) Loop
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = safe_logrank_pval_display(surv_obj, df_sub),
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = FALSE,
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}


### Optional refit using 90% confidence limits in the survfit object
# The current ggsurvplot call sets conf.int = FALSE, so this block does not
# display a confidence band despite calculating 90% limits.
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_with_CI_90.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%           # keep just the first draw per patient if multiple at that timepoint
      ungroup() %>%
      mutate(
        Group = factor(
          ifelse(.data[[var]] == 1, "Positive", "Negative"),
          levels = c("Negative","Positive")
        )
      )
    
    df_sub <- df_sub %>% 
      mutate(Time_to_event = Time_to_event/30.44) # divide days by 30.44 to get months
    
    # skip if too few pts
    if(nrow(df_sub) < min_n) next
    
    # skip if only one group present
    if(n_distinct(df_sub$Group) < 2) next
    
    surv_obj <- Surv(df_sub$Time_to_event, df_sub$Relapsed_Binary)
    fit      <- survfit(surv_obj ~ Group, data = df_sub)
    
    # 90% CI from survfit; "log" (Greenwood on log scale) is common, "log-log" is also fine
    fit <- survfit(
      surv_obj ~ Group,
      data      = df_sub,
      conf.int  = 0.90     # ← 90% CI
    )
    
    km <- ggsurvplot(
      fit, data       = df_sub,
      pval            = safe_logrank_pval_display(surv_obj, df_sub),
      break.time.by   = 12,        # put ticks every 12 “units” (i.e. every 12 months)
      conf.int        = FALSE,
      conf.int.alpha  = 0.1,    
      risk.table      = TRUE,
      risk.table.title = "Number at risk",
      risk.table.title.theme = element_text(hjust = 0),  # ← left‑align
      palette         = pal_2,
      # legend.title    = paste0(assay_lab, " MRD"),
      legend.title    = "MRD status",
      # now we know two groups are present, so these two labels fit
      legend.labs     = c("MRD–","MRD+"),
      xlab            = "Time since MRD assessment (months)",
      ylab            = "Progression-free survival",
      title = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25, 
      ## Added theme 
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position = "top",
          axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),          # no grid
          panel.grid.minor = element_blank(),
          #  Make the tick‑labels (the numbers) larger:
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y      = element_text(size = 15),
          axis.title.x      = element_text(size = 14)
        ),
    )
    
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title      = element_text(hjust = 0, face = "plain"),
      )
    
    # # Convert survfit object to a data.frame
    # fit_df <- broom::tidy(fit)  # gives time, n.risk, n.event, surv, std.err, conf.low, conf.high, strata
    # 
    # # Overlay ribbons for CIs
    # km$plot <- km$plot +
    #   geom_ribbon(
    #     data = fit_df,
    #     aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata, group = strata),
    #     inherit.aes = FALSE, alpha = 0.2
    #   )
    
    km$plot <- km$plot +
      theme(
        axis.title.x = element_blank()
      )
    
    combined <- ggarrange(
      km$plot, km$table,
      ncol    = 1,
      heights = c(3,1)
    )
    
    ggsave(
      filename = fname,
      plot     = combined,
      width    = 7, 
      height   = 7,
      dpi      = dpi_target
    )
  }
}
}

 



#### Optional delayed-entry/from-diagnosis KM sensitivity export
### This block evaluates an alternative time scale: months since diagnosis
### rather than months since MRD assessment. It uses delayed entry
### (left truncation) so patients do not contribute risk time before their MRD
### assessment. The final manuscript panels use the simpler and more directly
### interpretable "time since MRD assessment" scale above, so these exports are
### disabled by default.
export_optional_km_from_diagnosis_panels <- FALSE

if (isTRUE(export_optional_km_from_diagnosis_panels)) {
# pretty p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 1e-4")
  if (p < 1e-3) return("p < 0.001")
  if (p < 1e-2) return("p < 0.01")
  sprintf("p = %.2f", p)
}

# risk table for delayed-entry (start–stop) data
make_risktable <- function(df, breaks) {
  tmp <- df %>%
    mutate(Group_label = dplyr::recode(as.character(Group),
                                       "Negative" = "MRD–",
                                       "Positive" = "MRD+"))
  purrr::map_dfr(breaks, function(ti) {
    tmp %>%
      dplyr::group_by(Group_label) %>%
      dplyr::summarise(n = sum(entry_m <= ti & exit_m > ti), .groups = "drop") %>%
      dplyr::mutate(time = ti)
  }) %>%
    tidyr::pivot_wider(names_from = time, values_from = n) %>%
    dplyr::arrange(factor(Group_label, levels = c("MRD–","MRD+"))) %>%
    dplyr::rename(`MRD status` = Group_label)
}


for (tp in tps) {
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  tp_dir <- file.path(outdir, gsub("\\s+","_", tp))
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_from_diagnosis.png"))
    
    df_sub <- survival_df %>%
      filter(
        timepoint_info  == tp,
        !is.na(Time_to_event),
        !is.na(Relapsed_Binary),
        !is.na(.data[[var]])
      ) %>%
      arrange(Patient, sample_date) %>%
      group_by(Patient) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(dx_tbl, by = "Patient") %>%
      mutate(
        Group   = factor(ifelse(.data[[var]] == 1, "Positive", "Negative"),
                         levels = c("Negative","Positive")),
        # entry = months from diagnosis to MRD test
        entry_m = as.numeric(sample_date - diagnosis_date) / 30.44,
        # exit = entry + observed time after MRD test
        exit_m  = entry_m + (Time_to_event / 30.44)
      ) %>%
      filter(!is.na(entry_m), !is.na(exit_m), exit_m >= entry_m)
    
    if (nrow(df_sub) < min_n) next
    if (n_distinct(df_sub$Group) < 2) next
    
    # Cox with delayed entry
    surv_obj <- Surv(time = df_sub$entry_m, time2 = df_sub$exit_m, event = df_sub$Relapsed_Binary)
    fit     <- survfit(surv_obj ~ Group, data = df_sub)
    
    cox_fit <- coxph(surv_obj ~ Group, data = df_sub, ties = "breslow")
    s <- summary(cox_fit)
    
    # Pull p-values (log-rank / score / wald) safely
    p_lrt   <- suppressWarnings(as.numeric(s$logtest[3]))  # likelihood-ratio p
    p_score <- suppressWarnings(as.numeric(s$sctest[3]))   # score test p
    p_wald  <- suppressWarnings(as.numeric(s$waldtest[3])) # Wald p
    
    # Fallback if all above are NA (e.g., separation): grab coefficient p if present
    coef_p <- NA_real_
    cs <- try(coef(summary(cox_fit)), silent = TRUE)
    if (!inherits(cs, "try-error")) {
      # find the row for Group (works even if level names change)
      r <- grep("^Group", rownames(cs), value = FALSE)
      if (length(r) >= 1) coef_p <- cs[r[1], "Pr(>|z|)"]
    }
    
    # Choose the best available p
    p_any <- if (!is.na(p_lrt)) p_lrt else if (!is.na(p_score)) p_score else if (!is.na(p_wald)) p_wald else coef_p
    
    # Format: show thresholds instead of tiny decimals
    fmt_p <- function(p) {
      if (is.na(p)) return("p = NA")
      if (p < 1e-4) return("p < 1e-4")
      if (p < 1e-3) return("p < 0.001")
      if (p < 1e-2) return("p < 0.01")
      sprintf("p = %.2f", p)
    }
    pval_str <- fmt_p(p_any)
    
    km <- ggsurvplot(
      fit, data = df_sub,
      pval              = pval_str,
      break.time.by     = 12,
      conf.int          = FALSE,
      risk.table        = TRUE,
      risk.table.title  = "Number at risk",                    # ← add
      risk.table.title.theme = element_text(hjust = 0),        # ← add (left-align)
      palette           = pal_2,
      legend.title      = "MRD status",
      legend.labs       = c("MRD-","MRD+"),  # use ASCII hyphen to avoid file/device issues
      xlab              = "Months since diagnosis",
      ylab              = "Progression-free survival",
      title             = str_wrap(paste0("PFS Stratified by ", assay_lab, " at ", nice_tp), width = 45),
      risk.table.height = 0.25,
      ggtheme = theme_classic(base_size = 12) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 17),
          legend.position  = "top",
          axis.line        = element_line(colour = "black"),   # ← add
          panel.grid.major = element_blank(),                  # ← add
          panel.grid.minor = element_blank(),                  # ← add
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12),
          axis.title.y     = element_text(size = 15),
          axis.title.x     = element_text(size = 14)
        )
    )
    
    # match the old post-processing of table/plot
    km$table <- km$table +
      theme(
        axis.title.y = element_blank(),
        plot.title   = element_text(hjust = 0, face = "plain")
      )
    
    km$plot <- km$plot +
      theme(axis.title.x = element_blank())
    
    combined <- ggarrange(km$plot, km$table, ncol = 1, heights = c(3,1))
    ggsave(fname, plot = combined, width = 7, height = 7, dpi = dpi_target)
  }
}
}




## ── 5. SENSITIVITY ANALYSIS: Frontline Cohort MRD Performance ───────────────
##
##  This section calculates sensitivity and specificity metrics for each MRD assay
##  in the Frontline cohort at two critical clinical timepoints:
##    1. Post-ASCT (post-autologous stem cell transplant)
##    2. Maintenance-1yr (1-year maintenance therapy)
##
##  For each assay and timepoint:
##    - Filter to Frontline cohort only
##    - Stratify by MRD status (positive/negative)
##    - Calculate sensitivity = % of relapsed with MRD+
##    - Calculate specificity = % of non-relapsed with MRD-
##    - Calculate PPV/NPV = predictive values
##
##  Output: Summary table with MRD performance metrics
##
## ────────────────────────────────────────────────────────────────────────────

cat("5. Calculating MRD sensitivity/specificity metrics (Frontline cohort)...\n\n")

#### Now get stats for results 
### First on frontline 
# 1) Define frontline cohort for analysis
front_patients <- dat %>%
  filter(Cohort == "Frontline") %>%
  distinct(Patient)

# ─────────────────────────────────────────────────────────────────────────────
# 2) COMPUTE FOLLOW-UP TIME AND RELAPSE STATUS FOR FRONTLINE COHORT
# Time = number of days from baseline date (diagnosis) to censor/relapse date
# This is the time-to-event used for KM curves and outcome calculations
# ─────────────────────────────────────────────────────────────────────────────

pfs_front <- final_tbl %>%
  filter(Patient %in% front_patients$Patient) %>%
  # compute time from baseline to censor/relapse
  mutate(time_days = as.numeric(censor_date - baseline_date)) 

## check one row per patient (sanity check to ensure data integrity)
pfs_front %>% 
  count(Patient) %>% 
  filter(n > 1) -> dups
if(nrow(dups)) stop("Duplicate patients found: ", paste(dups$Patient, collapse = ", "))

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY STATISTICS: Frontline cohort follow-up and relapse
# ─────────────────────────────────────────────────────────────────────────────

median_fu_mo <- median(pfs_front$time_days / 30.44, na.rm = TRUE)
n_front      <- nrow(pfs_front)
n_rel        <- sum(pfs_front$relapsed)
pct_rel      <- n_rel / n_front * 100

cat(sprintf(
  "Frontline cohort summary:\n  Median follow-up: %.1f months\n  Relapses: %d/%d (%.0f%%) patients\n\n",
  median_fu_mo, n_rel, n_front, pct_rel
))

# ─────────────────────────────────────────────────────────────────────────────
# 3) DETAILED FOLLOW-UP STATISTICS
# Comprehensive summary of follow-up times in days and months
# Includes: min/Q1/median/Q3/max/mean/sd for all subjects
# ─────────────────────────────────────────────────────────────────────────────

followup_stats <- pfs_front %>%
  summarise(
    N_patients      = dplyr::n(),
    N_relapses      = sum(relapsed),
    Relapse_rate    = N_relapses / N_patients * 100,
    
    # Follow-up time in days
    min_days        = min(time_days, na.rm = TRUE),
    q1_days         = quantile(time_days, 0.25, na.rm = TRUE),
    median_days     = median(time_days, na.rm = TRUE),
    q3_days         = quantile(time_days, 0.75, na.rm = TRUE),
    max_days        = max(time_days, na.rm = TRUE),
    mean_days       = mean(time_days, na.rm = TRUE),
    sd_days         = sd(time_days, na.rm = TRUE),
    
    # Follow-up time in months (for clinical interpretation)
    min_months      = min(time_days, na.rm = TRUE) / 30.44,
    q1_months       = quantile(time_days / 30.44, 0.25, na.rm = TRUE),
    median_months   = median(time_days / 30.44, na.rm = TRUE),
    q3_months       = quantile(time_days / 30.44, 0.75, na.rm = TRUE),
    max_months      = max(time_days, na.rm = TRUE) / 30.44,
    mean_months     = mean(time_days, na.rm = TRUE) / 30.44,
    sd_months       = sd(time_days, na.rm = TRUE) / 30.44
  )

cat("Follow-up time summary (Frontline cohort):\n")
print(followup_stats)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4) OPTIONAL: Visualization of follow-up time distribution
# Histogram showing how follow-up times are distributed across cohort
# ─────────────────────────────────────────────────────────────────────────────

followup_hist <- pfs_front %>%
  mutate(followup_months = time_days / 30.44) %>%
  ggplot(aes(x = followup_months)) +
  geom_histogram(binwidth = 3, boundary = 0) +
  labs(
    x = "Follow-up time (months)",
    y = "Number of patients",
    title = "Distribution of follow-up times in frontline cohort"
  )

# If you want to export the stats to CSV
write_csv(followup_stats, file.path(outdir, paste0("frontline_followup_summary_", date_tag, ".csv")))
## Can use 1A for this as well 

# 4.  Assays & timepoint definitions ------------------------------------------
assays <- c(
  EasyM    = "EasyM_reference_threshold_binary",
  clonoSEQ = "Adaptive_Binary",
  Flow     = "Flow_Binary",
  cfWGS_BM    = "BM_zscore_only_detection_rate_call",
  # cfWGS_BM_screen removed 2026-07-29 (unprovenanced 0.350 threshold; see audit
  # note at the data-preparation step).
  cfWGS_Blood_Sites    = "Blood_zscore_only_sites_call",
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
)

post_labels   <- c("post_transplant")
one_year_labels <- c("1yr maintenance")

# Sensitivity helper functions ------------------------------------------------
#
# Manuscript role:
#   These helpers generate the source tables used by Extended Data Figure 6A
#   (BM-subset assay sensitivity) and Extended Data Figure 8A (blood-subset
#   assay sensitivity). Sensitivity is calculated among patients who eventually
#   relapsed: N assay-positive / N tested at the landmark.
#
# Scientific assumptions:
#   - One row per patient per landmark is retained, using the earliest sample if
#     more than one sample exists at that landmark.
#   - Sensitivity denominators are assay-specific because not all assays were
#     available for every patient/timepoint.
#   - BM-subset and blood-subset tables restrict to patients with the relevant
#     cfWGS assay available, so head-to-head comparisons use matched cfWGS-tested
#     subsets.
compute_sens <- function(df, col) {
  df2      <- df %>% filter(!is.na(.data[[col]]))  # Remove rows with missing assay result
  n_tested <- nrow(df2)                            # Total relapsed patients with assay result
  n_pos    <- sum(df2[[col]] == 1, na.rm = TRUE)   # Number of those who were MRD+
  tibble(
    N_tested      = n_tested,
    N_positive    = n_pos,
    Sensitivity   = n_pos / n_tested              # Sensitivity = % MRD+ among relapsed
  )
}

build_relapsed_landmark_df <- function(dat, final_tbl, front_patients, labels, assays) {
  dat %>%
    filter(
      Patient %in% front_patients$Patient,
      str_detect(timepoint_info, paste(labels, collapse = "|"))
    ) %>%
    arrange(Patient, sample_date) %>%
    group_by(Patient) %>%
    slice(1) %>%
    ungroup() %>%
    select(Patient, one_of(assays)) %>%
    left_join(final_tbl %>% select(Patient, relapsed), by = "Patient") %>%
    filter(relapsed == 1)
}

summarise_sensitivity_by_assay <- function(df, assays) {
  map_dfr(names(assays), function(assay_label) {
    compute_sens(df, assays[[assay_label]]) %>%
      mutate(Assay = assay_label)
  }) %>%
    select(Assay, everything())
}

build_landmark_sensitivity_tables <- function(post_df, year_df, assays, subset_col = NULL) {
  if (!is.null(subset_col)) {
    post_df <- post_df %>% filter(!is.na(.data[[subset_col]]))
    year_df <- year_df %>% filter(!is.na(.data[[subset_col]]))
  }
  list(
    post = summarise_sensitivity_by_assay(post_df, assays),
    year = summarise_sensitivity_by_assay(year_df, assays)
  )
}

print_sensitivity_table <- function(tbl, title) {
  cat("\n")
  cat("═════════════════════════════════════════════════════════════════════════\n")
  cat(title, "\n", sep = "")
  cat("═════════════════════════════════════════════════════════════════════════\n")
  print(tbl)
}

write_sensitivity_table <- function(tbl, filename) {
  path <- file.path(outdir, paste0(filename, date_tag, ".csv"))
  write_csv(tbl, path)
  message("Sensitivity table written: ", path)
  invisible(path)
}

# Build the two landmark tables once, then reuse them for overall, BM-tested,
# and blood-tested sensitivity summaries.
post_df <- build_relapsed_landmark_df(dat, final_tbl, front_patients, post_labels, assays)
year_df <- build_relapsed_landmark_df(dat, final_tbl, front_patients, one_year_labels, assays)

overall_sens <- build_landmark_sensitivity_tables(post_df, year_df, assays)
post_stats <- overall_sens$post
year_stats <- overall_sens$year

bm_col <- assays[["cfWGS_BM"]]
blood_col <- assays[["cfWGS_Blood_Sites"]]

bm_subset_sens <- build_landmark_sensitivity_tables(post_df, year_df, assays, subset_col = bm_col)
post_stats_BM <- bm_subset_sens$post
year_stats_BM <- bm_subset_sens$year

blood_subset_sens <- build_landmark_sensitivity_tables(post_df, year_df, assays, subset_col = blood_col)
post_stats_blood <- blood_subset_sens$post
year_stats_blood <- blood_subset_sens$year

print_sensitivity_table(post_stats, "POST-ASCT SENSITIVITY: % of relapsed patients with MRD+ at post-ASCT")
print_sensitivity_table(year_stats, "1-YEAR MAINTENANCE SENSITIVITY: % of relapsed patients with MRD+ at 1yr")
print_sensitivity_table(post_stats_BM, "BM-SUBSET: Post-ASCT sensitivities among patients with BM-cfWGS")
print_sensitivity_table(year_stats_BM, "BM-SUBSET: 1-year sensitivities among patients with BM-cfWGS")
print_sensitivity_table(post_stats_blood, "BLOOD-SUBSET: Post-ASCT sensitivities among patients with blood-cfWGS")
print_sensitivity_table(year_stats_blood, "BLOOD-SUBSET: 1-year sensitivities among patients with blood-cfWGS")

cat("\nSaving sensitivity results...\n")
write_sensitivity_table(post_stats, "frontline_postASCT_sensitivity_")
write_sensitivity_table(year_stats, "frontline_1yr_sensitivity_")
write_sensitivity_table(post_stats_BM, "frontline_postASCT_sens_BMcfWGS_")
write_sensitivity_table(year_stats_BM, "frontline_1yr_sens_BMcfWGS_")
write_sensitivity_table(post_stats_blood, "frontline_postASCT_sens_bloodcfWGS_")
write_sensitivity_table(year_stats_blood, "frontline_1yr_sens_bloodcfWGS_")



# ─────────────────────────────────────────────────────────────────────────────
# 12. SENSITIVITY BARPLOT FOR MANUSCRIPT SUPPLEMENT
# Comparative visualization of sensitivity across assays and timepoints
# ─────────────────────────────────────────────────────────────────────────────

# Combine BM-subset results (post-ASCT and maintenance) into long format
# Filter to exclude blood-only cfWGS columns
# Recodes assay names to short manuscript labels
# Enforces desired display order: clinical assays → cfWGS

sens_df_bm <- bind_rows(
  post_stats_BM  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_BM  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  # Exclude blood-derived assays and screening variant from BM-subset comparison
  filter(Assay != "cfWGS_Blood_Sites") %>%
  filter(Assay != "cfWGS_Blood_Combined") %>%
  filter(Assay != "cfWGS_BM_screen") %>%
  mutate(
    # Convert to percentage for labeling
    Sens_pct   = Sensitivity * 100,
    # Manuscript-friendly assay names
    Assay      = recode(Assay,
                        EasyM          = "EasyM",
                        clonoSEQ       = "clonoSEQ",
                        Flow           = "MFC",
                        cfWGS_BM       = "cfWGS"),
    # Enforce display order: cfWGS, clonoSEQ, MFC, then EasyM rightmost
    Assay = factor(Assay, levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM")),
    # Enforce timepoint order for legend and grouping
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_6A (BM sensitivity barplot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  sens_df_bm,
  file.path(outdir_source_data, paste0("Supp_6A_BM_sensitivity_barplot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_6A (BM sensitivity)\n")

# ─────────────────────────────────────────────────────────────────────────────
# 13. CONFIGURE BARPLOT COLORS AND THEME
# ─────────────────────────────────────────────────────────────────────────────
# Color scheme:
#   - Post-ASCT: Deep teal (#31688E) - early therapeutic assessment
#   - Maintenance-1yr: Bright green (#35B779) - longer-term surveillance
# Theme minimizes clutter while maintaining publication quality

custom_cols <- c(
  "Post-ASCT"       = "#31688E",  # deep teal
  "Maintenance-1yr" = "#35B779"   # bright green
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 14. BUILD GROUPED BARPLOT: Sensitivity by Assay × Timepoint
# ─────────────────────────────────────────────────────────────────────────────
# Grouped bars allow visual comparison of:
#   - Sensitivity across assays (x-axis: clonoSEQ, MFC, EasyM, cfWGS)
#   - Timepoint effect (bars grouped by color: Post-ASCT vs Maintenance-1yr)
# Bar height = sensitivity %; text labels show exact percentages

p_sens <- ggplot(sens_df_bm,
                 aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  # position_dodge separates bars for same assay; width controls bar width
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",        # black outline for clarity
           size     = 0.3) +
  # Add percentage labels on top of bars
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,           # position above bar
            size     = 3.5) +
  # Use custom color mapping with explicit order
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  # Y-axis: 0-100% scale
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

print(p_sens)

# SAVE BM-SUBSET SENSITIVITY BARPLOT (500 DPI for manuscripts)
# Saves to: Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training3_{date}.png

edfig6a_path <- paste0("Final Tables and Figures/Supp_6A_Fig_sensitivity_by_tech_training3_", date_tag, ".png")
ggsave(
  filename = edfig6a_path,
  plot     = p_sens,
  width    = 6,
  height   = 4,
  dpi      = 500
)

# MANUSCRIPT OUTPUT: Extended Data Figure 6A
# BM-subset relapse-detection sensitivity by MRD assay and landmark timepoint.
ms_copy_artifact(
  source_path = edfig6a_path,
  artifact_id = "EDFIG6A",
  role = "figure_panel_png",
  description = "Extended Data Figure 6A: BM-subset sensitivity by MRD assay and landmark timepoint.",
  script_name = "4_1_Survival_Analysis.R"
)

cat("  ✓ Saved BM-subset sensitivity barplot\n\n")

# BLOOD-SUBSET SENSITIVITY BARPLOT - Same structure as BM analysis
# Restricted to blood-derived cfWGS samples for head-to-head comparison


## Now for blood
sens_df_blood <- bind_rows(
  post_stats_blood  %>% mutate(Timepoint = "Post-ASCT"),
  year_stats_blood  %>% mutate(Timepoint = "Maintenance-1yr")
) %>%
  filter(Assay != "cfWGS_BM") %>%
  filter(Assay != "cfWGS_BM_screen") %>%
  mutate(
    # Convert to percentage for labeling
    Sens_pct   = Sensitivity * 100,
    # Manuscript-friendly assay names (with line breaks for blood cfWGS variants)
    Assay      = recode(Assay,
                        EasyM                 = "EasyM",
                        clonoSEQ              = "clonoSEQ",
                        Flow                  = "MFC",
                        cfWGS_Blood_Sites     = "cfWGS\n(Sites Model)",
                        cfWGS_Blood_Combined  = "cfWGS\n(Combined Model)"),
    # Enforce timepoint order: Post-ASCT → Maintenance-1yr
    Timepoint = factor(Timepoint, levels = c("Post-ASCT", "Maintenance-1yr"))
  ) 

# Enforce assay ordering: cfWGS variants (blood) first, then clonoSEQ, MFC, EasyM rightmost.
# Note: Blood models use "\n" for line break in x-axis labels; EasyM is rightmost for consistency.
sens_df_blood <- sens_df_blood %>%
  mutate(Assay = factor(Assay,
                        levels = c("cfWGS\n(Sites Model)",
                                   "cfWGS\n(Combined Model)", 
                                   "clonoSEQ", "MFC", "EasyM")))

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_8A (Blood sensitivity barplot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  sens_df_blood,
  file.path(outdir_source_data, paste0("Supp_8A_blood_sensitivity_barplot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_8A (Blood sensitivity)\n")

# Build blood-subset barplot (same structure as BM version)
p_sens_blood <- ggplot(sens_df_blood,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Timepoint",
    values = custom_cols[c("Post-ASCT", "Maintenance-1yr")],  # enforce mapping
    limits = c("Post-ASCT", "Maintenance-1yr")                # enforce order
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of cfDNA and Clinical MRD Assays in Relapsing Patients",
    x     = "Technology",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

# Display the blood-subset barplot
p_sens_blood

# SAVE BLOOD-SUBSET SENSITIVITY BARPLOT (500 DPI for manuscripts)
# Saves to: Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood2_{date}.png

edfig8a_path <- paste0("Final Tables and Figures/Supp_8A_Fig_sensitivity_by_tech_training_blood2_", date_tag, ".png")
ggsave(
  filename = edfig8a_path,
  plot     = p_sens_blood,
  width    = 6,
  height   = 4,
  dpi      = 500
)

# MANUSCRIPT OUTPUT: Extended Data Figure 8A
# Blood-subset relapse-detection sensitivity by MRD assay and landmark timepoint.
ms_copy_artifact(
  source_path = edfig8a_path,
  artifact_id = "EDFIG8A",
  role = "figure_panel_png",
  description = "Extended Data Figure 8A: blood-subset sensitivity by MRD assay and landmark timepoint.",
  script_name = "4_1_Survival_Analysis.R"
)

cat("  ✓ Saved blood-subset sensitivity barplot\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 17. FRONTLINE LANDMARK SURVIVAL SUMMARIES FOR ED6/ED8
# ─────────────────────────────────────────────────────────────────────────────
#
# Role in manuscript:
#   The data frame `survival_df` was filtered to the Frontline cohort near the
#   top of the script. The blocks below compute descriptive 24-month RFS,
#   median RFS, Cox HRs, rank correlations, and power diagnostics at the
#   one-year-maintenance and post-ASCT landmarks. These summaries feed the ED6B
#   and ED8B hazard-ratio source tables/plots generated later in this script.
#   The non-frontline/test-cohort time-window analysis begins in the separate
#   "Time-window prediction performance in Non-frontline cohort" section below.

# Refactored landmark-summary helpers -------------------------------------------------
#
# The original working script calculated the same 24-month RFS, median RFS, Cox
# HR, Spearman correlation, and power-diagnostic values separately for each
# model/timepoint. These helpers keep the same calculations but centralize the
# repeated mechanics so the active command-line path is easier to audit.


safe_spearman <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3 || dplyr::n_distinct(x[ok]) < 2 || dplyr::n_distinct(y[ok]) < 2) {
    return(list(rho = NA_real_, p = NA_real_))
  }
  out <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(out$estimate), p = out$p.value)
}

summarise_binary_survival <- function(df, assay_col, t_days = 24 * 30.44) {
  df_assay <- df %>%
    filter(!is.na(.data[[assay_col]])) %>%
    mutate(
      .assay_numeric = as.integer(.data[[assay_col]]),
      .assay_factor = factor(.assay_numeric, levels = c(0L, 1L))
    )

  n_patients <- dplyr::n_distinct(df_assay$Patient)
  if (nrow(df_assay) == 0 || dplyr::n_distinct(df_assay$.assay_numeric) < 2) {
    return(tibble(
      n = n_patients,
      RFS24_neg = NA_real_,
      RFS24_pos = NA_real_,
      MedRFS_neg = NA_real_,
      MedRFS_pos = NA_real_,
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_
    ))
  }

  fit <- survival::survfit(
    survival::Surv(Time_to_event, Relapsed_Binary) ~ .assay_factor,
    data = df_assay
  )

  surv_at <- summary(fit, times = t_days)
  strata <- as.character(surv_at$strata)
  rfs_neg <- surv_at$surv[match(".assay_factor=0", strata)] * 100
  rfs_pos <- surv_at$surv[match(".assay_factor=1", strata)] * 100

  med_tbl <- survminer::surv_median(fit)
  med_strata <- as.character(med_tbl$strata)
  med_neg <- med_tbl$median[match(".assay_factor=0", med_strata)] / 30.44
  med_pos <- med_tbl$median[match(".assay_factor=1", med_strata)] / 30.44

  cox_tbl <- broom::tidy(
    survival::coxph(
      survival::Surv(Time_to_event, Relapsed_Binary) ~ .assay_numeric,
      data = df_assay
    ),
    exponentiate = TRUE,
    conf.int = TRUE
  )

  tibble(
    n = n_patients,
    RFS24_neg = rfs_neg,
    RFS24_pos = rfs_pos,
    MedRFS_neg = med_neg,
    MedRFS_pos = med_pos,
    HR = cox_tbl$estimate[1],
    CI_low = cox_tbl$conf.low[1],
    CI_high = cox_tbl$conf.high[1]
  )
}

power_diagnostics <- function(df, primary_call_col) {
  d <- sum(df$Relapsed_Binary, na.rm = TRUE)
  prop_p <- mean(df[[primary_call_col]] == 1, na.rm = TRUE)

  if (is.na(prop_p) || prop_p <= 0 || prop_p >= 1 || d <= 0) {
    return(tibble(Events = d, Patients = nrow(df), HR_80pct = NA_real_, Power_HR2_pct = NA_real_))
  }

  z_alpha <- qnorm(1 - 0.05 / 2)
  z_beta_80 <- qnorm(0.80)
  hr80 <- exp(2 * (z_alpha + z_beta_80) / sqrt(d * prop_p * (1 - prop_p)))

  hr_target <- 2
  z_beta <- sqrt(d * prop_p * (1 - prop_p)) * log(hr_target) - z_alpha
  pw2 <- pnorm(z_beta) * 100

  tibble(Events = d, Patients = nrow(df), HR_80pct = hr80, Power_HR2_pct = pw2)
}

build_landmark_progression_row <- function(survival_df,
                                           landmark,
                                           primary_call_col,
                                           primary_prob_col,
                                           include_easym = TRUE) {
  df_km <- survival_df %>%
    filter(timepoint_info == landmark, !is.na(.data[[primary_call_col]]))

  assay_summaries <- list(
    cf = summarise_binary_survival(df_km, primary_call_col),
    fl = summarise_binary_survival(df_km, "Flow_Binary"),
    seq = summarise_binary_survival(df_km, "Adaptive_Binary")
  )
  if (isTRUE(include_easym)) {
    assay_summaries$em <- summarise_binary_survival(df_km, "EasyM_reference_threshold_binary")
  }

  prob_cor <- safe_spearman(df_km[[primary_prob_col]], df_km$Time_to_event)
  flow_cor <- safe_spearman(df_km$Flow_pct_cells, df_km$Time_to_event)
  power_tbl <- power_diagnostics(df_km, primary_call_col)

  row <- tibble(
    Landmark = ifelse(landmark == "1yr maintenance", "1yr_maintenance", landmark),
    RFS24_cf_neg = assay_summaries$cf$RFS24_neg,
    RFS24_cf_pos = assay_summaries$cf$RFS24_pos,
    MedRFS_cf_neg = assay_summaries$cf$MedRFS_neg,
    MedRFS_cf_pos = assay_summaries$cf$MedRFS_pos,
    RFS24_fl_neg = assay_summaries$fl$RFS24_neg,
    RFS24_fl_pos = assay_summaries$fl$RFS24_pos,
    MedRFS_fl_neg = assay_summaries$fl$MedRFS_neg,
    MedRFS_fl_pos = assay_summaries$fl$MedRFS_pos,
    RFS24_seq_neg = assay_summaries$seq$RFS24_neg,
    RFS24_seq_pos = assay_summaries$seq$RFS24_pos,
    MedRFS_seq_neg = assay_summaries$seq$MedRFS_neg,
    MedRFS_seq_pos = assay_summaries$seq$MedRFS_pos,
    HR_seq = assay_summaries$seq$HR,
    CI_low_seq = assay_summaries$seq$CI_low,
    CI_high_seq = assay_summaries$seq$CI_high,
    HR_cf = assay_summaries$cf$HR,
    CI_low_cf = assay_summaries$cf$CI_low,
    CI_high_cf = assay_summaries$cf$CI_high,
    HR_fl = assay_summaries$fl$HR,
    CI_low_fl = assay_summaries$fl$CI_low,
    CI_high_fl = assay_summaries$fl$CI_high,
    Spearman_prob = prob_cor$rho,
    Spearman_flow = flow_cor$rho,
    Events = power_tbl$Events,
    Patients = power_tbl$Patients,
    HR_80pct = power_tbl$HR_80pct,
    Power_HR2_pct = power_tbl$Power_HR2_pct,
    N_cfWGS = assay_summaries$cf$n,
    N_MFC = assay_summaries$fl$n,
    N_clonoSEQ = assay_summaries$seq$n
  )

  if (isTRUE(include_easym)) {
    row <- row %>%
      mutate(
        RFS24_em_neg = assay_summaries$em$RFS24_neg,
        RFS24_em_pos = assay_summaries$em$RFS24_pos,
        MedRFS_em_neg = assay_summaries$em$MedRFS_neg,
        MedRFS_em_pos = assay_summaries$em$MedRFS_pos,
        HR_em = assay_summaries$em$HR,
        CI_low_em = assay_summaries$em$CI_low,
        CI_high_em = assay_summaries$em$CI_high,
        N_easym = assay_summaries$em$n
      )
  }

  row
}

build_landmark_progression_table <- function(survival_df,
                                             primary_call_col,
                                             primary_prob_col,
                                             include_easym = TRUE) {
  bind_rows(
    build_landmark_progression_row(
      survival_df = survival_df,
      landmark = "post_transplant",
      primary_call_col = primary_call_col,
      primary_prob_col = primary_prob_col,
      include_easym = include_easym
    ),
    build_landmark_progression_row(
      survival_df = survival_df,
      landmark = "1yr maintenance",
      primary_call_col = primary_call_col,
      primary_prob_col = primary_prob_col,
      include_easym = include_easym
    )
  )
}


export_landmark_progression_table <- function(tbl, csv_name, rds_name) {
  write_csv(tbl, file.path(outdir, paste0(csv_name, date_tag, ".csv")))
  saveRDS(tbl, file.path(outdir, paste0(rds_name, date_tag, ".rds")))
}

validate_landmark_progression_table <- function(tbl, table_name, require_easym = TRUE) {
  expected_landmarks <- c("post_transplant", "1yr_maintenance")
  required_cols <- c(
    "Landmark",
    "RFS24_cf_neg", "RFS24_cf_pos", "MedRFS_cf_neg", "MedRFS_cf_pos",
    "RFS24_fl_neg", "RFS24_fl_pos", "MedRFS_fl_neg", "MedRFS_fl_pos",
    "RFS24_seq_neg", "RFS24_seq_pos", "MedRFS_seq_neg", "MedRFS_seq_pos",
    "HR_cf", "CI_low_cf", "CI_high_cf",
    "HR_fl", "CI_low_fl", "CI_high_fl",
    "HR_seq", "CI_low_seq", "CI_high_seq",
    "Spearman_prob", "Spearman_flow",
    "Events", "Patients", "HR_80pct", "Power_HR2_pct",
    "N_cfWGS", "N_MFC", "N_clonoSEQ"
  )
  if (isTRUE(require_easym)) {
    required_cols <- c(
      required_cols,
      "RFS24_em_neg", "RFS24_em_pos", "MedRFS_em_neg", "MedRFS_em_pos",
      "HR_em", "CI_low_em", "CI_high_em", "N_easym"
    )
  }

  missing_cols <- setdiff(required_cols, names(tbl))
  if (length(missing_cols) > 0) {
    stop(
      table_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(tbl) != length(expected_landmarks)) {
    stop(
      table_name,
      " should contain one row per frontline landmark (expected ",
      length(expected_landmarks),
      ", observed ",
      nrow(tbl),
      ").",
      call. = FALSE
    )
  }
  if (!identical(sort(as.character(tbl$Landmark)), sort(expected_landmarks))) {
    stop(
      table_name,
      " has unexpected Landmark values: ",
      paste(as.character(tbl$Landmark), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(tbl)
}

use_refactored_landmark_summaries <- TRUE

if (isTRUE(use_refactored_landmark_summaries)) {
  progression_metrics <- build_landmark_progression_table(
    survival_df = survival_df,
    primary_call_col = "BM_zscore_only_detection_rate_call",
    primary_prob_col = "BM_zscore_only_detection_rate_prob",
    include_easym = TRUE
  )

  progression_metrics_blood <- build_landmark_progression_table(
    survival_df = survival_df,
    primary_call_col = "Blood_zscore_only_sites_call",
    primary_prob_col = "Blood_zscore_only_sites_prob",
    include_easym = TRUE
  )

  progression_metrics_blood_combined <- build_landmark_progression_table(
    survival_df = survival_df,
    primary_call_col = "Blood_plus_fragment_call",
    primary_prob_col = "Blood_plus_fragment_prob",
    include_easym = FALSE
  )

  validate_landmark_progression_table(progression_metrics, "BM landmark progression table", require_easym = TRUE)
  validate_landmark_progression_table(progression_metrics_blood, "Blood landmark progression table", require_easym = TRUE)
  validate_landmark_progression_table(
    progression_metrics_blood_combined,
    "Blood combined-model landmark progression table",
    require_easym = FALSE
  )


  export_landmark_progression_table(
    progression_metrics,
    csv_name = "cfWGS_vs_flow_progression_summary_",
    rds_name = "cfWGS_vs_flow_progression_summary_updated_"
  )
  export_landmark_progression_table(
    progression_metrics_blood,
    csv_name = "cfWGS_vs_flow_progression_summary_blood_muts_",
    rds_name = "cfWGS_vs_flow_progression_summaryy_blood_muts_updated_"
  )
  export_landmark_progression_table(
    progression_metrics_blood_combined,
    csv_name = "cfWGS_vs_flow_progression_summary_blood_muts_combined_model_",
    rds_name = "cfWGS_vs_flow_progression_summary_blood_muts_updated_combined_model_"
  )
}

# Legacy comparison path removed:
#   The previous manual landmark-summary implementation repeated the same
#   survival, Cox, Spearman, power, and export logic for BM, blood,
#   and combined blood models. The active helper implementation above now
#   creates the same progression_metrics objects and exports used below.
#   Historical copies remain in local or Git history if detailed audit is needed.

hr_plot_df_blood <- progression_metrics_blood %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq,
         HR_em,   CI_low_em,   CI_high_em) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq|em)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS (Sites Model)",
                   fl = "MFC",
                   seq = "clonoSEQ",
                   em = "EasyM"),
    Assay = factor(Assay,
                   levels = c("cfWGS (Sites Model)", "clonoSEQ", "MFC", "EasyM")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post-ASCT", "Maintenance-1yr"))
  )

# reshape combined model (cfWGS only)
hr_plot_df_combined <- progression_metrics_blood_combined %>%
  select(Landmark,
         HR_cf, CI_low_cf, CI_high_cf) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf)"
  ) %>%
  mutate(
    Assay = "cfWGS (Combined Model)",
    Assay = factor(Assay,
                   levels = c("cfWGS (Sites Model)", "cfWGS (Combined Model)", "clonoSEQ", "MFC", "EasyM")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post-ASCT", "Maintenance-1yr"))
  )


# bind together Sites + Combined models
hr_plot_df_blood <- bind_rows(hr_plot_df_blood, hr_plot_df_combined) %>%
  mutate(Assay = factor(Assay,
                        levels = c("cfWGS (Sites Model)", "cfWGS (Combined Model)", "clonoSEQ", "MFC", "EasyM")))

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: SuppFig8B (Blood HR plot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  hr_plot_df_blood,
  file.path(outdir_source_data, paste0("SuppFig8B_blood_HR_plot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: SuppFig8B (Blood HR)\n")

p_hr <- ggplot(hr_plot_df_blood,
               aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_continuous(
    "Hazard ratio (log scale)",
    trans        = "log10",
    limits       = c(0.04, 100),
    breaks       = c(0.1, 0.5, 1, 2, 5, 10, 100),
    minor_breaks = c(
      0.05, 0.06, 0.08,           # between 0.05 & 0.1
      0.15, 0.2, 0.3, 0.4,        # between 0.1 & 0.5
      0.6, 0.8,                   # between 0.5 & 1
      1.5, 3,                      # between 1 & 5
      6, 8,                        # between 5 & 10
      15, 20, 30, 50, 60, 80       # between 10 & 100
    ),
    labels = label_number(accuracy = .1)
  )+
  annotation_logticks(
    sides  = "b",
    short  = unit(2, "pt"),
    mid    = unit(4, "pt"),
    long   = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS (Sites Model)" = "#35608DFF",
               "cfWGS (Combined Model)" = "#440154FF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF",   # orange for clonoSEQ
               "EasyM" = "#D81B60FF")  # magenta for EasyM
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #   subtitle = "cfWGS vs. MFC (95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

edfig8b_path <- paste0("Final Tables and Figures/SuppFig8B_cfWGS_blood_HR_updated3_", date_tag, ".png")
ggsave(edfig8b_path, p_hr, width = 6, height = 4, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 8B
# Hazard ratios for relapse by blood-derived MRD assay and landmark timepoint.
ms_copy_artifact(
  source_path = edfig8b_path,
  artifact_id = "EDFIG8B",
  role = "figure_panel_png",
  description = "Extended Data Figure 8B: hazard ratios by blood-derived MRD assay and landmark timepoint.",
  script_name = "4_1_Survival_Analysis.R"
)



### Now for BM derived muts
hr_plot_df_bm <- progression_metrics %>%
  select(Landmark,
         HR_cf,   CI_low_cf,   CI_high_cf,
         HR_fl,   CI_low_fl,   CI_high_fl,
         HR_seq,   CI_low_seq,   CI_high_seq,
         HR_em,   CI_low_em,   CI_high_em) %>%
  pivot_longer(
    cols      = -Landmark,
    names_to  = c(".value", "Assay"),
    names_pattern = "(HR|CI_low|CI_high)_(cf|fl|seq|em)"
  ) %>%
  mutate(
    Assay = recode(Assay,
                   cf = "cfWGS",
                   fl = "MFC",
                   seq = "clonoSEQ",
                   em = "EasyM"),
    Assay = factor(Assay,
                   levels = c("cfWGS", "clonoSEQ", "MFC", "EasyM")),
    Landmark = factor(Landmark,
                      levels = c("post_transplant", "1yr_maintenance"),
                      labels = c("Post-ASCT", "Maintenance-1yr"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# SOURCE DATA EXPORT: Supp_Figure_6B (BM HR plot)
# ─────────────────────────────────────────────────────────────────────────────
write_csv(
  hr_plot_df_bm,
  file.path(outdir_source_data, paste0("Supp_Figure_6B_BM_HR_plot_source_data_", date_tag, ".csv"))
)
cat("  ✓ Exported source data: Supp_Figure_6B (BM HR)\n")

bm_hr_axis_limits <- c(0.04, 300)
bm_hr_ci_display_limits <- c(0.05, 250)
bm_hr_breaks <- c(0.05, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 250)
bm_hr_minor_breaks <- c(
  0.06, 0.08, 0.1, 0.15, 0.25, 0.3, 0.4, 0.6, 0.8,
  1.5, 3, 4, 6, 8, 15, 30, 40, 60, 80, 150, 200
)

hr_plot_df_bm_display <- hr_plot_df_bm %>%
  mutate(
    # Display-only caps keep zero/infinite separation intervals visible on a
    # finite log axis. The uncapped HR/CI values remain in the source-data CSV.
    CI_low_plot = if_else(
      is.finite(CI_low) & CI_low > 0,
      pmax(CI_low, bm_hr_ci_display_limits[1]),
      bm_hr_ci_display_limits[1]
    ),
    CI_high_plot = if_else(
      is.finite(CI_high) & CI_high > 0,
      pmin(CI_high, bm_hr_ci_display_limits[2]),
      bm_hr_ci_display_limits[2]
    )
  )

p_hr_bm <- ggplot(hr_plot_df_bm_display,
                  aes(x = HR, y = fct_rev(Landmark), colour = Assay)) +
  # reference line
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # 1) horizontal CIs
  geom_errorbar(
    aes(xmin = CI_low_plot, xmax = CI_high_plot),
    position = position_dodge(width = 0.6),
    orientation = "y",
    linewidth = 0.5
  ) +
  
  # 2) dots at the HR
  geom_point(
    position = position_dodge(width = 0.6),
    size     = 3
  ) +
  
  # log scale axis
  scale_x_log10(
    "Hazard ratio (log scale)",
    limits       = bm_hr_axis_limits,
    breaks       = bm_hr_breaks,
    minor_breaks = bm_hr_minor_breaks,
    labels = function(x) {
      sapply(x, function(xx) {
        if (xx > 1) {
          sprintf("%.0f", xx)
        } else {
          sprintf("%.2f", xx)
        }
      })
    }
  ) +
  annotation_logticks(
    sides = "b",
    short = unit(2, "pt"),
    mid   = unit(4, "pt"),
    long  = unit(6, "pt")
  ) +
  # colours
  scale_colour_manual(
    name   = NULL,
    values = c("cfWGS" = "#35608DFF",
               "MFC"   = "#43BF71FF",
               "clonoSEQ"= "#E69F00FF",   # orange for clonoSEQ
               "EasyM" = "#D81B60FF")    # magenta/pink for EasyM
  ) +
  
  labs(
    y        = NULL,
    title    = "Relapse hazard ratios stratified by MRD assay\nand landmark timepoint",
    #  subtitle = "(95% CI)"
  ) +
  
  # classic theme with no gridlines
  theme_classic(base_size = 11) +
  theme(
    panel.grid         = element_blank(),   # no grid at all
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold",
                                      hjust = 0.5),  # bold + centered
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  )

p_hr_bm

edfig6b_path <- paste0("Final Tables and Figures/Supp_Figure_6B_cfWGS_BM_HR_updated3_", date_tag, ".png")
ggsave(edfig6b_path, p_hr_bm, width = 6, height = 4, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 6B
# Hazard ratios for relapse by BM-derived MRD assay and landmark timepoint.
ms_copy_artifact(
  source_path = edfig6b_path,
  artifact_id = "EDFIG6B",
  role = "figure_panel_png",
  description = "Extended Data Figure 6B: hazard ratios by BM-derived MRD assay and landmark timepoint.",
  script_name = "4_1_Survival_Analysis.R"
)






### Now make time to relapse figure 
df <- survival_df %>%                           # <- the tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(BM_zscore_only_detection_rate_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      BM_zscore_only_detection_rate_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.4215524
#youden_thresh2 <- 0.35

max_mo <- max(df$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, BM_zscore_only_detection_rate_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
   # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
  shape  = 21,
  colour = "black",
  size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("cVAF Model Probability",
                     limits = c(0,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing BM-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

## Add the samples 
# Make sure the outcome labels match the scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
 #   values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )

print(p_prob2)

# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5_label.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 6J
# BM-derived cfWGS probability versus time before relapse/censoring.
ms_copy_artifact(
  source_path = "Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated5_label.png",
  artifact_id = "EDFIG6J",
  role = "figure_panel_png",
  description = "Extended Data Figure 6J: BM-derived cfWGS probability over time before relapse or censoring.",
  script_name = "4_1_Survival_Analysis.R"
)



### Instead change the scale to match what Trevor was thinking 
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = BM_zscore_only_detection_rate_prob,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to Inf
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 180   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric. Some optional component metrics are not
# present in every regenerated object, so skip those descriptively instead of
# stopping the survival pipeline.
safe_spearman <- function(df, x_col, y_col) {
  if (!all(c(x_col, y_col) %in% names(df))) return(NULL)
  x <- suppressWarnings(as.numeric(df[[x_col]]))
  y <- suppressWarnings(as.numeric(df[[y_col]]))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || dplyr::n_distinct(x[keep]) < 2 || dplyr::n_distinct(y[keep]) < 2) {
    return(NULL)
  }
  cor.test(x[keep], y[keep], method = "spearman")
}

print_spearman <- function(label, test_result) {
  cat("\n--- ", label, " ---\n", sep = "")
  if (is.null(test_result)) {
    cat("Not available: column missing, nonnumeric, or insufficient variation.\n")
  } else {
    print(test_result)
  }
}

spearman_prob <- safe_spearman(rel_df, "BM_zscore_only_detection_rate_prob", "days_before_event")
spearman_zscore <- safe_spearman(rel_df, "zscore_BM", "days_before_event")
spearman_detect_rate <- safe_spearman(rel_df, "detect_rate_BM", "days_before_event")

# 3) print them
print_spearman("cfWGS model probability vs days", spearman_prob)
print_spearman("zscore_BM vs days", spearman_zscore)
print_spearman("detect_rate_BM vs days", spearman_detect_rate)


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -30 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(BM_zscore_only_detection_rate_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nrho=%.2f, p=%s\nPre-relapse only:\nrho=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = BM_zscore_only_detection_rate_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # Inf‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or Inf for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "Inf")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0,1),
    breaks = seq(0,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.01,      # left margin
  #          y = 1650,
  #          label = sprintf("rho = %.2f\np = %s", rho, pval_str),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.01,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +

  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4D_time_to_relapse_infinity_no_reverse2_BM_muts_updated3.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


## Manuscript layout for Extended Data Figure 6K
library(cowplot)    # get_legend()
library(patchwork)  # easy assembly

pal_vals <- c("No relapse" = "black", "Relapse" = "red")
legend_df <- tibble::tibble(
  x = 1,
  y = 1,
  progress_status = factor(names(pal_vals), levels = names(pal_vals))
)

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <- ggplot(plot_df2,
                 aes(x = BM_zscore_only_detection_rate_prob, y = days_plot)) +
  geom_vline(xintercept = 0.35, linetype = "dotted", colour = "grey40") +
  geom_point(aes(colour = progress_status, fill = progress_status),
             shape = 21, size = 3, colour = "black") +
  scale_y_continuous(
    "Days until relapse (or Inf for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "Inf")
  ) +
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0.1,1),
    breaks = seq(0.1,1,by=0.2),
    labels = scales::percent_format(accuracy=1)
  )  +# colours
scale_colour_manual(
  "Patient outcome",
  values = c("No relapse" = "black", "Relapse" = "red")
) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  guides(fill = "none") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none"        # <- hide here; we'll place it below
  ) +
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using BM-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nrho=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nrho=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_BM_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 6K
# BM-derived cfWGS probability by time-to-relapse window with footer summaries.
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_BM_muts.png",
  artifact_id = "EDFIG6K",
  role = "figure_panel_png",
  description = "Extended Data Figure 6K: BM-derived cfWGS probability versus days before relapse with summary footer.",
  script_name = "4_1_Survival_Analysis.R"
)



plot_df2 %>%
  filter(
    is.na(BM_zscore_only_detection_rate_prob) | is.na(days_plot) |
      BM_zscore_only_detection_rate_prob < 0 | BM_zscore_only_detection_rate_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )

time_to_relapse_BM <- plot_df2



### Extended Data Figure 8E/8F blood-derived longitudinal and association panels
df <- survival_df %>%                           # <- the tibble
  # keep samples beyond baseline / diagnosis
  filter(!str_detect(timepoint_info, regex("Diagnosis|Baseline", TRUE))) %>%
  
  # drop rows with missing probability or time
  filter(!is.na(Blood_zscore_only_sites_prob),
         !is.na(Time_to_event)) %>%
  
  # enforce non-negative time to event for relapse event visits a few days off from CMRG date
  mutate(
    days_before_event = pmax(Time_to_event, 0), # set negative values to 0 
    mrd_status      = factor(
      Blood_zscore_only_sites_call,
      levels = c(0, 1),
      labels = c("MRD-", "MRD+")
    ),
    progress_status = factor(
      Relapsed_Binary,
      levels = c(0, 1),
      labels = c("No relapse", "Relapse")
    ),
    
    # time *before* the anchor (positive value) → plot reversed
    days_before_event = Time_to_event,    # keep positive for clarity
    months_before_event = days_before_event/30.44
  )

## Only multiple points 
df_slim <- df %>%
  group_by(Patient) %>% 
  filter(dplyr::n() > 1) %>%   # keep only patients with >1 row
  ungroup()

# ────────────────────────────────────────────────────────────────
# 2.  Plot  ──────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
youden_thresh <- 0.5166693
youden_thr <- youden_thresh # for consistency
max_mo <- max(df_slim$months_before_event, na.rm = TRUE)  

p_prob <- ggplot(df, aes(months_before_event, Blood_zscore_only_sites_prob, group = Patient)) +
  
  # 1) Youden line
  geom_hline(yintercept = youden_thresh,
             linetype = "dotted", colour = "gray40") +
  
  # 2) trajectories coloured by relapse
  geom_line(aes(colour = progress_status),
            size = 0.4, alpha = 0.4) +
  
  # 3) points: fill by relapse, stroke by MRD call, border black
  geom_point(aes(
    fill   = progress_status),
    shape  = 21,
    colour = "black",
    size   = 2
  ) +
  
  # 4) event line
  #geom_vline(xintercept = 0, linetype = "dotted", colour = "gray40") +
  
  # 5) axes
  scale_x_reverse(
    name         = "Months before event or censor",
    breaks       = seq(0, max_mo, by = 12),  # every 12 months
    minor_breaks = seq(0, max_mo, by = 6)    # every 6 months
  ) +
  scale_y_continuous("Sites Model Probability",
                     limits = c(0.3,1),
                     labels = scales::percent_format(1)) +
  
  # 6) colour for relapse status
  # scale_colour_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  # scale_fill_manual(
  #   name   = "Patient outcome",
  #   values = c("No relapse" = "#35608DFF",
  #              "Relapse"    = "#43BF71FF")
  # ) +
  
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "black",
               "Relapse"    = "red")
  ) +
  # # 7) stroke scale for MRD call
  # scale_discrete_manual(
  #   aesthetics = "stroke",
  #   values     = c("MRD-" = 0, "MRD+" = 1),
  #   guide      = guide_legend(
  #     title = "MRD call",
  #     override.aes = list(
  #       shape  = 21,
  #       fill   = "white",   # white interior in legend, makes stroke obvious
  #       size   = 4,
  #       colour = "black",
  #       stroke = c(0, 1)
  #     )
  #   )
  # ) +
  
  # 8) clean up legends
  guides(
    colour = guide_legend(order = 1),
    fill   = FALSE   # only show stroke legend for MRD
  ) +
  
  labs(
    title    = "Longitudinal cfWGS MRD Probability by Patient Outcome\nUsing cfDNA-Derived Mutation Lists"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "bottom", 
    legend.title    = element_text(size = 11),   # 
    legend.text     = element_text(size = 10)    # even smaller
  )

print(p_prob)

### Add label 
# Make sure the outcome labels match the scales
df <- df %>% mutate(progress_status = factor(progress_status,
                                             levels = c("No relapse","Relapse")))

n_patients_by <- df %>%
  distinct(Patient, progress_status) %>%
  count(progress_status, name = "n_patients")

n_timepoints_by <- df %>%
  count(progress_status, name = "n_timepoints")

# pull counts (0 if a group is absent)
get_n <- function(tbl, lvl, col) {
  val <- tbl %>% filter(progress_status == lvl) %>% pull({{col}})
  if (length(val) == 0) 0 else val
}

n_pat_nr  <- get_n(n_patients_by,  "No relapse", n_patients)
n_time_nr <- get_n(n_timepoints_by, "No relapse", n_timepoints)
n_pat_rl  <- get_n(n_patients_by,  "Relapse",    n_patients)
n_time_rl <- get_n(n_timepoints_by, "Relapse",   n_timepoints)

# text to print
lab_nr <- paste0("No relapse: n=", n_pat_nr, " patients; ", n_time_nr, " samples")
lab_rl <- paste0("Relapse: n=", n_pat_rl, " patients; ", n_time_rl, " samples")

# place labels at bottom-left (remember: x is reversed, so 'left' == large x)
x_left <- max_mo - 0.02 * max_mo  # a small inset from the left border

p_prob2 <- p_prob +
  scale_colour_manual(
    name = "Patient outcome",
#    values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
    values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  ) +
  scale_fill_manual(
    name = "Patient outcome",
   # values = c("No relapse" = "#35608DFF", "Relapse" = "#43BF71FF"),
   values = c("No relapse" = "black", "Relapse" = "red"),
    labels = c(
      paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
      paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
    )
  )



print(p_prob2)
# ────────────────────────────────────────────────────────────────
# 3.  Export  ────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood4.png",
       p_prob, width = 6, height = 4.5, dpi = 600)

ggsave("Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood2_labelled3.png",
       p_prob2, width = 6, height = 4.5, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 8E
# Blood-derived cfWGS probability versus time before relapse/censoring.
ms_copy_artifact(
  source_path = "Final Tables and Figures/F4C_cfWGS_prob_vs_time_updated3_blood2_labelled3.png",
  artifact_id = "EDFIG8E",
  role = "figure_panel_png",
  description = "Extended Data Figure 8E: longitudinal cfDNA-derived cfWGS probability over time before relapse or censoring.",
  script_name = "4_1_Survival_Analysis.R"
)


## Intermediate blood/cfDNA time-to-relapse scale check
### This plot is retained as an intermediate check. The version used in
### Extended Data Figure 8F is exported from the infinity-axis
### layout below.
df_plot <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%   # months → days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                       # keep all progressors
      row_number() == which.min(days_before_event)       # keep *latest* censor
  ) %>%
  ungroup()

max_days <- ceiling(max(df_plot$days_before_event, na.rm = TRUE) / 180) * 180

df_plot <- df_plot %>%
  mutate(Time_to_event = if_else(Time_to_event >= -30 & Time_to_event < 0, 0, Time_to_event))

## ─────────────────────────────────────────────────────────────
## 1)  Build the scatter plot                                  
## ─────────────────────────────────────────────────────────────
p_time <- ggplot(df_plot,
                 aes(x = Blood_zscore_only_sites_prob,
                     y = days_before_event)) +
  
  # ① dashed horizontal line at “event” (0 days)
  #  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  
  # ② points – colour = outcome, stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status,
                 stroke = mrd_status),
             shape  = 21,
             size   = 3,
             colour = "black") +
  
  # ③ axes
  scale_x_continuous(
    "cfWGS MRD probability",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_y_reverse(
    "Days until relapse (or censor)",
    limits = c(max_days, 0),
    breaks = seq(0, max_days, by = 180),      # every ~6 months
    minor_breaks = seq(0, max_days, by = 90)  # every 3 months
  ) +
  
  # ④ colours for relapse status
  scale_colour_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = c("No relapse" = "#35608DFF",
               "Relapse"    = "#43BF71FF")
  ) +
  
  # ⑤ stroke scale for MRD call
  scale_discrete_manual(
    aesthetics = "stroke",
    values     = c("MRD-" = 0, "MRD+" = 1.1),   # ring only if MRD+
    guide = guide_legend(
      title          = "MRD call",
      override.aes   = list(shape = 21,
                            size  = 4,
                            colour = "black",
                            fill   = "white",
                            stroke = c(0, 1.1))
    )
  ) +
  
  # ⑥ theme / labels
  labs(
    title = "Time to relapse vs. cfWGS MRD probability"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 8)
  )

print(p_time)

# save
ggsave(file.path(outdir, "Fig_time_to_relapse_vs_prob2.png"),
       p_time, width = 6, height = 4, dpi = 600)


### Show non-relapsers as infinity
# 1) compute days & a “days_plot” that sends non-relapsers to Inf
plot_df2 <- df %>%
  mutate(days_before_event = months_before_event * 30.44) %>%        # months→days
  group_by(Patient) %>%
  filter(
    progress_status == "Relapse" |                                 # keep all relapsers
      row_number() == which.min(days_before_event)                   # for non-relapsers, keep their last sample
  ) %>%
  ungroup()

# define axis maximum and “infinity” sentinel
max_days <- ceiling(max(plot_df2$days_before_event, na.rm = TRUE) / 180) * 180
overflow <- max_days + 180   # a little beyond the longest follow‑up

plot_df2 <- plot_df2 %>%
  mutate(
    days_plot = if_else(progress_status == "No relapse", overflow, days_before_event)
  )


## Check corrs 
# 1) subset to relapsers
rel_df <- plot_df2 %>%
  filter(progress_status == "Relapse")

# 2) run cor.test for each metric
spearman_prob <- cor.test(
  rel_df$Blood_zscore_only_sites_prob,
  rel_df$days_before_event,
  method = "spearman"
)

# 3) print them
cat("\n--- cfWGS model probability vs days ---\n")
print(spearman_prob)


# Edit days
plot_df2 <- plot_df2 %>%
  mutate(days_plot = if_else(days_plot >= -35 & days_plot <= 0, 0, days_plot))

# 2) compute Spearman rho on the relapsers
spearman_res <- with(
  filter(plot_df2, progress_status == "Relapse"),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

## Check other metrics
# 2A) pull out estimate + p‑value
# A) including relapse samples (the current spearman_res)
rho_all  <- spearman_res$estimate
pval_all <- spearman_res$p.value

# B) excluding relapse samples (days_plot > 0 only)
spearman_pre <- with(
  filter(plot_df2, progress_status == "Relapse", days_plot > 0),
  cor.test(Blood_zscore_only_sites_prob,
           days_before_event,
           method = "spearman")
)

rho_pre  <- spearman_pre$estimate
pval_pre <- spearman_pre$p.value

pval_all_str <- ifelse(pval_all < 0.001, "<0.001", sprintf("%.3f", pval_all))
pval_pre_str <- ifelse(pval_pre < 0.001, "<0.001", sprintf("%.3f", pval_pre))

annot_text <- sprintf("All relapse samples:\nrho=%.2f, p=%s\nPre-relapse only:\nrho=%.2f, p=%s",
                      rho_all, pval_all_str,
                      rho_pre, pval_pre_str)

## Sanity check probability range before applying retained axis limits
plot_df2 %>%
  summarise(
    min   = min(Blood_zscore_only_sites_prob, na.rm = TRUE),
    max   = max(Blood_zscore_only_sites_prob, na.rm = TRUE),
    range = max(Blood_zscore_only_sites_prob, na.rm = TRUE) -
      min(Blood_zscore_only_sites_prob, na.rm = TRUE)
  )

# 3) make the scatter
p_time_inf <- ggplot(plot_df2,
                     aes(x = Blood_zscore_only_sites_prob,
                         y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # Inf‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or Inf for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "Inf")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # annotate Spearman
  # annotate("text",
  #          x = 0.90,      # left margin
  #          y = 1650,
  #          label = sprintf("rho = %.2f\np = %s", rho, pval_fmt),
  #          hjust = 0,
  #          size = 3.5) +
  annotate("text",
           x = 0.8,
           y = 1450,
           label = annot_text,
           hjust = 0,
           size = 3.5) +
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

print(p_time_inf)

# 4) render / save
ggsave(file.path("Final Tables and Figures/Fig_4_D_time_to_relapse_infinity_no_reverse2_blood_muts4_small2.png"),
       p_time_inf, width = 5.5, height = 4.5, dpi = 600)


time_to_relapse_blood <- plot_df2

## Sanity check rows outside retained plotting limits
plot_df2 %>%
  filter(
    is.na(Blood_zscore_only_sites_prob) | is.na(days_plot) |
      Blood_zscore_only_sites_prob < 0.3 | Blood_zscore_only_sites_prob > 1 |
      days_plot < 0 | days_plot > overflow
  )



### Manuscript layout for Extended Data Figure 8F

# --- 1) build the main plot *without* a legend (legend handled below) ---
p_main <-  ggplot(plot_df2,
                  aes(x = Blood_zscore_only_sites_prob,
                      y = days_plot)) +
  
  # Youden threshold (if you still want it)
  geom_vline(xintercept = youden_thr, linetype = "dotted", colour = "grey40") +
  
  # points coloured by relapse; stroke = MRD call
  geom_point(aes(colour = progress_status,
                 fill   = progress_status),
             shape = 21, size = 3, colour = "black") +
  
  # Inf‐aware y‐axis
  scale_y_continuous(
    "Days until relapse (or Inf for censor)",
    limits = c(0, overflow),
    breaks = c(seq(0, 1620, by = 180), overflow),
    labels = c(seq(0, 1620, by = 180), "Inf")
  ) +
  
  # x‐axis as percent
  scale_x_continuous(
    "cfWGS MRD probability (%)",
    limits = c(0.3,1),
    breaks = seq(0,1,by=0.1),
    labels = scales::percent_format(accuracy=1)
  ) +
  
  # colours
  scale_colour_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  scale_fill_manual(
    "Patient outcome",
    values = c("No relapse" = "black", "Relapse" = "red")
  ) +
  
  # clean theme
  theme_classic(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  ) +
  
  labs(title = "Association Between cfWGS MRD Probability\nand Time to Relapse Using cfDNA-Derived Mutations")

# --- 2) LEGEND-ONLY PLOT (do NOT inherit guides(fill = 'none')) ---
# -----build a dummy legend that always shows both levels -----
p_legend_only <- ggplot(legend_df, aes(x, y, fill = progress_status)) +
  # make the plotting layer invisible *in the panel* …
  geom_point(shape = 21, size = 0, colour = "black", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Patient outcome",
    values = pal_vals,
    breaks = names(pal_vals)
  ) +
  guides(fill = guide_legend(
    ncol = 1,                 # <- stack items vertically
    byrow = TRUE,
    title.position = "top",   # title on its own line
    label.hjust = 0,          # left-align labels
    override.aes = list(shape = 21, size = 3, alpha = 1, colour = "black")
  )) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "vertical",         # <- vertical legend
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(6, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin       = margin(0, 0, 0, 0)
  )

# --- 3) make two small “text boxes” for the right columns ---
txt_all <- sprintf("All relapse samples:\nrho=%.2f, p=%s", rho_all, pval_all_str)
txt_pre <- sprintf("Pre-relapse only:\nrho=%.2f, p=%s", rho_pre, pval_pre_str)

mini_box <- function(s) {
  ggplot() +
    annotate("label", x = 0, y = 1, label = s,
             hjust = 0, vjust = 1, size = 3.5,
             label.size = 0, fill = scales::alpha("white", 0.7)) +
    xlim(0,1) + ylim(0,1) +
    theme_void()
}

col2 <- mini_box(txt_all)
col3 <- mini_box(txt_pre)

# --- 4) assemble: plot on top; 3 columns underneath ---
bottom_row <- p_legend_only | col2 | col3
final_plot <- p_main / bottom_row + plot_layout(heights = c(1, 0.22))

# show and save
print(final_plot)
ggsave("Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_blood_muts.png",
       final_plot, width = 5.5, height = 5.5, dpi = 600)

# MANUSCRIPT OUTPUT: Extended Data Figure 8F
# Blood-derived cfWGS probability by time-to-relapse window with footer summaries.
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_blood_muts.png",
  artifact_id = "EDFIG8F",
  role = "figure_panel_png",
  description = "Extended Data Figure 8F: association between cfDNA-derived cfWGS probability and days before relapse with summary footer.",
  script_name = "4_1_Survival_Analysis.R"
)

## ── Additive train+test longitudinal probability figures ───────────────────
##
## These figures use the preserved model probabilities and thresholds, but
## expand the plotted rows to all Frontline + Non-frontline samples with usable
## patient outcome data. They are intentionally saved as separate outputs so the
## submitted frontline-only manuscript panels remain unchanged.

format_p_for_label <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p_value)
}

safe_spearman_label <- function(df, x_col, y_col) {
  if (!all(c(x_col, y_col) %in% names(df))) {
    return("rho=NA, p=NA")
  }
  x <- suppressWarnings(as.numeric(df[[x_col]]))
  y <- suppressWarnings(as.numeric(df[[y_col]]))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || dplyr::n_distinct(x[keep]) < 2 || dplyr::n_distinct(y[keep]) < 2) {
    return("rho=NA, p=NA")
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
  sprintf("rho=%.2f, p=%s", unname(test$estimate), format_p_for_label(test$p.value))
}

plot_train_test_longitudinal_probability <- function(survival_data,
                                                      prob_col,
                                                      call_col,
                                                      threshold,
                                                      y_label,
                                                      title_suffix,
                                                      output_stem,
                                                      y_limits = c(0, 1),
                                                      title_text = NULL) {
  required_cols <- c(
    "Patient", "Cohort", "timepoint_info", "sample_date", "Time_to_event",
    "Relapsed_Binary", prob_col, call_col
  )
  missing_cols <- setdiff(required_cols, names(survival_data))
  if (length(missing_cols) > 0) {
    stop("Train+test longitudinal plot is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  plot_df <- survival_data %>%
    dplyr::filter(
      Cohort %in% train_test_cohorts,
      !stringr::str_detect(timepoint_info, stringr::regex("Diagnosis|Baseline", TRUE)),
      !is.na(.data[[prob_col]]),
      !is.na(.data[[call_col]]),
      !is.na(Time_to_event),
      !is.na(Relapsed_Binary),
      Time_to_event >= 0
    ) %>%
    dplyr::mutate(
      days_before_event = pmax(Time_to_event, 0),
      months_before_event = days_before_event / 30.44,
      mrd_status = factor(
        .data[[call_col]],
        levels = c(0, 1),
        labels = c("MRD-", "MRD+")
      ),
      progress_status = factor(
        Relapsed_Binary,
        levels = c(0, 1),
        labels = c("No relapse", "Relapse")
      ),
      cohort_label = dplyr::recode(
        Cohort,
        "Frontline" = "Training Cohort",
        "Non-frontline" = "Test Cohort",
        .default = Cohort
      )
    )

  if (nrow(plot_df) == 0) {
    warning("No rows available for ", output_stem, call. = FALSE)
    return(invisible(NULL))
  }

  cohort_levels_present <- sort(unique(plot_df$cohort_label))
  has_multiple_plot_cohorts <- length(cohort_levels_present) > 1
  cohort_title_note <- if (has_multiple_plot_cohorts) {
    "Training Cohort + Test Cohort assay-available rows"
  } else {
    "Assay-available rows after Training Cohort + Test Cohort screening: Training Cohort only"
  }

  max_mo <- ceiling(max(plot_df$months_before_event, na.rm = TRUE) / 12) * 12
  if (!is.finite(max_mo) || max_mo <= 0) max_mo <- 12

  n_patients_by <- plot_df %>%
    dplyr::distinct(Patient, progress_status) %>%
    dplyr::count(progress_status, name = "n_patients")
  n_timepoints_by <- plot_df %>%
    dplyr::count(progress_status, name = "n_timepoints")

  get_count <- function(tbl, lvl, col) {
    val <- tbl %>% dplyr::filter(progress_status == lvl) %>% dplyr::pull({{ col }})
    if (length(val) == 0) 0L else val
  }

  n_pat_nr <- get_count(n_patients_by, "No relapse", n_patients)
  n_time_nr <- get_count(n_timepoints_by, "No relapse", n_timepoints)
  n_pat_rl <- get_count(n_patients_by, "Relapse", n_patients)
  n_time_rl <- get_count(n_timepoints_by, "Relapse", n_timepoints)

  base_plot <- ggplot(
    plot_df,
    aes(x = months_before_event, y = .data[[prob_col]], group = Patient)
  ) +
    geom_hline(yintercept = threshold, linetype = "dotted", colour = "gray40") +
    geom_line(aes(colour = progress_status), size = 0.4, alpha = 0.4) +
    geom_point(
      aes(fill = progress_status),
      shape = 21,
      colour = "black",
      size = 2
    ) +
    scale_x_reverse(
      name = "Months before next event or censor",
      breaks = seq(0, max_mo, by = 12),
      minor_breaks = seq(0, max_mo, by = 6)
    ) +
    scale_y_continuous(
      y_label,
      limits = y_limits,
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_colour_manual(
      name = "Patient outcome",
      values = c("No relapse" = "black", "Relapse" = "red"),
      labels = c(
        paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
        paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
      )
    ) +
    scale_fill_manual(
      name = "Patient outcome",
      values = c("No relapse" = "black", "Relapse" = "red"),
      labels = c(
        paste0("No relapse\n(n=", n_pat_nr, " patients; ", n_time_nr, " samples)"),
        paste0("Relapse\n(n=", n_pat_rl, " patients; ", n_time_rl, " samples)")
      )
    ) +
    labs(
      title = if (is.null(title_text)) {
        paste0(
          "Longitudinal cfWGS MRD Probability by Next Patient Outcome\n",
          title_suffix, "\n",
          cohort_title_note
        )
      } else {
        title_text
      }
    ) +
    theme_classic(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  source_path <- file.path(
    outdir_source_data,
    paste0(output_stem, "_source_data_", date_tag, ".csv")
  )
  readr::write_csv(plot_df, source_path)

  ggsave(
    file.path("Final Tables and Figures", paste0(output_stem, ".png")),
    base_plot,
    width = 6.4,
    height = 4.8,
    dpi = 600
  )

  by_cohort_summary <- plot_df %>%
    dplyr::distinct(Patient, Cohort, progress_status) %>%
    dplyr::count(Cohort, progress_status, name = "n_patients") %>%
    dplyr::full_join(
      plot_df %>% dplyr::count(Cohort, progress_status, name = "n_samples"),
      by = c("Cohort", "progress_status")
    )

  readr::write_csv(
    by_cohort_summary,
    file.path(outdir_source_data, paste0(output_stem, "_counts_", date_tag, ".csv"))
  )

  invisible(plot_df)
}

plot_train_test_time_to_relapse <- function(plot_df,
                                            prob_col,
                                            threshold,
                                            x_label,
                                            title_suffix,
                                            output_stem,
                                            x_limits = c(0, 1),
                                            title_text = NULL,
                                            footer_statistics = FALSE,
                                            infinity_label = "infinity") {
  if (is.null(plot_df) || nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  relapse_plot_df <- plot_df %>%
    dplyr::arrange(Patient, sample_date) %>%
    dplyr::group_by(Patient) %>%
    dplyr::filter(
      progress_status == "Relapse" |
        dplyr::row_number() == which.min(days_before_event)
    ) %>%
    dplyr::ungroup()

  if (nrow(relapse_plot_df) == 0) {
    warning("No rows available for ", output_stem, call. = FALSE)
    return(invisible(NULL))
  }

  max_days <- ceiling(max(relapse_plot_df$days_before_event, na.rm = TRUE) / 180) * 180
  if (!is.finite(max_days) || max_days <= 0) max_days <- 180
  overflow <- max_days + 180

  relapse_plot_df <- relapse_plot_df %>%
    dplyr::mutate(
      days_plot = dplyr::if_else(progress_status == "No relapse", overflow, days_before_event)
    )

  cohort_levels_present <- sort(unique(relapse_plot_df$cohort_label))
  has_multiple_plot_cohorts <- length(cohort_levels_present) > 1
  cohort_title_note <- if (has_multiple_plot_cohorts) {
    "frontline + non-frontline assay-available rows"
  } else {
    "frontline assay-available rows only"
  }

  rel_all <- relapse_plot_df %>% dplyr::filter(progress_status == "Relapse")
  rel_pre <- rel_all %>% dplyr::filter(days_before_event > 0)
  annot_text <- paste0(
    "All relapse samples:\n",
    safe_spearman_label(rel_all, prob_col, "days_before_event"),
    "\nPre-relapse only:\n",
    safe_spearman_label(rel_pre, prob_col, "days_before_event")
  )

  readr::write_csv(
    relapse_plot_df,
    file.path(outdir_source_data, paste0(output_stem, "_source_data_", date_tag, ".csv"))
  )

  plot_title <- if (is.null(title_text)) {
    paste0(
      "cfWGS MRD Probability vs Time to Next Relapse/Progression\n",
      title_suffix, "\n",
      cohort_title_note
    )
  } else {
    title_text
  }

  p_time_main <- ggplot(
    relapse_plot_df,
    aes(x = .data[[prob_col]], y = days_plot)
  ) +
    geom_vline(xintercept = threshold, linetype = "dotted", colour = "grey40") +
    geom_point(
      aes(fill = progress_status),
      shape = 21,
      size = 3,
      colour = "black"
    ) +
    scale_y_continuous(
      "Days until next relapse/progression (or infinity for censor)",
      limits = c(0, overflow),
      breaks = c(seq(0, max_days, by = 180), overflow),
      labels = c(seq(0, max_days, by = 180), infinity_label)
    ) +
    scale_x_continuous(
      x_label,
      limits = x_limits,
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_fill_manual(
      "Patient outcome",
      values = c("No relapse" = "black", "Relapse" = "red")
    ) +
    labs(title = plot_title) +
    theme_classic(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = if (isTRUE(footer_statistics)) 13 else 12
      ),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  if (isTRUE(footer_statistics)) {
    p_time_main <- p_time_main +
      theme(legend.position = "none")

    p_legend_only <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 1.08,
        label = "Patient outcome",
        hjust = 0,
        vjust = 1,
        size = 3.5
      ) +
      annotate(
        "point",
        x = 0.08,
        y = 0.63,
        shape = 21,
        fill = "black",
        colour = "black",
        size = 3
      ) +
      annotate(
        "text",
        x = 0.18,
        y = 0.63,
        label = "No relapse",
        hjust = 0,
        size = 3.2
      ) +
      annotate(
        "point",
        x = 0.08,
        y = 0.23,
        shape = 21,
        fill = "red",
        colour = "black",
        size = 3
      ) +
      annotate(
        "text",
        x = 0.18,
        y = 0.23,
        label = "Relapse",
        hjust = 0,
        size = 3.2
      ) +
      xlim(0, 1) +
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      theme_void()

    mini_box <- function(label_text) {
      ggplot() +
        annotate(
          "text",
          x = 0,
          y = 1,
          label = label_text,
          hjust = 0,
          vjust = 1,
          size = 3.5
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void()
    }

    txt_all <- paste0(
      "All relapse samples:\n",
      safe_spearman_label(rel_all, prob_col, "days_before_event")
    )
    txt_pre <- paste0(
      "Pre-relapse only:\n",
      safe_spearman_label(rel_pre, prob_col, "days_before_event")
    )
    footer <- patchwork::wrap_plots(
      p_legend_only,
      mini_box(txt_all),
      mini_box(txt_pre),
      nrow = 1,
      widths = c(1, 1, 1)
    )
    p_time <- patchwork::wrap_plots(
      p_time_main,
      footer,
      ncol = 1,
      heights = c(1, 0.22)
    )
  } else {
    p_time <- p_time_main +
      annotate(
        "label",
        x = x_limits[1] + 0.58 * diff(x_limits),
        y = max_days,
        label = annot_text,
        hjust = 0,
        size = 3.2,
        fill = scales::alpha("white", 0.75)
      ) +
      theme(legend.position = "bottom")
  }

  ggsave(
    file.path("Final Tables and Figures", paste0(output_stem, ".png")),
    p_time,
    width = if (isTRUE(footer_statistics)) 6.5 else 5.8,
    height = if (isTRUE(footer_statistics)) 5.5 else 4.9,
    dpi = 600
  )

  invisible(relapse_plot_df)
}

bm_train_test_longitudinal_df <- plot_train_test_longitudinal_probability(
  survival_data = survival_df_train_test,
  prob_col = "BM_zscore_only_detection_rate_prob",
  call_col = "BM_zscore_only_detection_rate_call",
  threshold = 0.4215524,
  y_label = "cVAF Model Probability",
  title_suffix = "Using BM-Derived Mutation Lists",
  output_stem = "F4C_cfWGS_prob_vs_time_all_train_test_BM_outcome_available",
  y_limits = c(0, 1),
  title_text = paste0(
    "Longitudinal cfWGS MRD Probability by Patient Outcome\n",
    "Using BM-Derived Mutation Lists"
  )
)

ms_copy_artifact(
  source_path = "Final Tables and Figures/F4C_cfWGS_prob_vs_time_all_train_test_BM_outcome_available.png",
  artifact_id = "EDFIG6J",
  role = "all_samples_figure_panel_png",
  description = "All-evaluable-sample training/test combined version of Extended Data Figure 6J using the next progression/censor endpoint after each sample.",
  script_name = "4_1_Survival_Analysis.R"
)

plot_train_test_time_to_relapse(
  plot_df = bm_train_test_longitudinal_df,
  prob_col = "BM_zscore_only_detection_rate_prob",
  threshold = 0.4215524,
  x_label = "cfWGS MRD probability (%)",
  title_suffix = "Using BM-Derived Mutations",
  output_stem = "Fig_4D_time_to_relapse_footer3cols_BM_muts_all_train_test_outcome_available",
  x_limits = c(0, 1),
  title_text = paste0(
    "Association Between cfWGS MRD Probability\n",
    "and Time to Relapse Using BM-Derived Mutations"
  ),
  footer_statistics = TRUE,
  infinity_label = "Inf"
)

ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_BM_muts_all_train_test_outcome_available.png",
  artifact_id = "EDFIG6K",
  role = "all_samples_figure_panel_png",
  description = "All-evaluable-sample training/test combined version of Extended Data Figure 6K using the next progression/censor endpoint after each sample.",
  script_name = "4_1_Survival_Analysis.R"
)

blood_train_test_longitudinal_df <- plot_train_test_longitudinal_probability(
  survival_data = survival_df_train_test,
  prob_col = "Blood_zscore_only_sites_prob",
  call_col = "Blood_zscore_only_sites_call",
  threshold = 0.5166693,
  y_label = "Sites Model Probability",
  title_suffix = "Using cfDNA-Derived Mutation Lists",
  output_stem = "F4C_cfWGS_prob_vs_time_all_train_test_blood_outcome_available",
  y_limits = c(0.3, 1),
  title_text = paste0(
    "Longitudinal cfWGS MRD Probability by Patient Outcome\n",
    "Using cfDNA-Derived Mutation Lists"
  )
)

ms_copy_artifact(
  source_path = "Final Tables and Figures/F4C_cfWGS_prob_vs_time_all_train_test_blood_outcome_available.png",
  artifact_id = "EDFIG8E",
  role = "all_samples_figure_panel_png",
  description = "All-evaluable-sample training/test combined version of Extended Data Figure 8E using the next progression/censor endpoint after each sample.",
  script_name = "4_1_Survival_Analysis.R"
)

plot_train_test_time_to_relapse(
  plot_df = blood_train_test_longitudinal_df,
  prob_col = "Blood_zscore_only_sites_prob",
  threshold = 0.5166693,
  x_label = "cfWGS MRD probability (%)",
  title_suffix = "Using cfDNA-Derived Mutations",
  output_stem = "Fig_4D_time_to_relapse_footer3cols_blood_muts_all_train_test_outcome_available",
  x_limits = c(0.3, 1),
  title_text = paste0(
    "Association Between cfWGS MRD Probability\n",
    "and Time to Relapse Using cfDNA-Derived Mutations"
  ),
  footer_statistics = TRUE,
  infinity_label = "Inf"
)

ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig_4D_time_to_relapse_footer3cols_blood_muts_all_train_test_outcome_available.png",
  artifact_id = "EDFIG8F",
  role = "all_samples_figure_panel_png",
  description = "All-evaluable-sample training/test combined version of Extended Data Figure 8F using the next progression/censor endpoint after each sample.",
  script_name = "4_1_Survival_Analysis.R"
)













### Non-frontline/test-cohort time-window evaluation
## Reintroduce the full all-cohort table because the previous sections limited
## `dat` and `survival_df` to Cohort == "Frontline".
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


## Do rescored
# 2026-07-29 audit fix: the "screen column" re-creation was removed here too.
# See the audit note at the data-preparation step: the 0.350 cutoff was an
# unprovenanced operating point borrowed from a legacy combo_BM specification.
# The non-frontline time-window analyses below use the frozen Youden calls.

################################################################################
##  Time-window prediction performance in Non-frontline cohort
################################################################################

# 1) The assays vector
assays <- c(
  EasyM        = "EasyM_reference_threshold_binary",
  Flow         = "Flow_Binary",
  cfWGS_BM     = "BM_zscore_only_detection_rate_call",
  cfWGS_Blood  = "Blood_zscore_only_sites_call", 
  cfWGS_Blood_Combined = "Blood_plus_fragment_call"
)

## Get additional dates
Relapse_dates_full <- read_csv(
  "Exported_data_tables_clinical/Relapse dates cfWGS updated2.csv",
  col_types = cols(
    Patient          = col_character(),
    Progression_date = col_date(format = "%Y-%m-%d")
  )
)

# 2) Build df_sf: non-frontline samples + assays + relapse info
df_sf <- dat %>%
  filter(tolower(Cohort) == "non-frontline") %>%
  transmute(
    Patient,
    sample_date = as.Date(Date),
    Flow_Binary,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_call,
    Blood_plus_fragment_call
  ) %>%
  # keep rows with at least one available assay result. Some assay columns,
  # such as EasyM, are not present in every regenerated table.
  {
    assays_present <- intersect(unname(assays), names(.))
    if (length(assays_present) == 0) {
      stop("No non-frontline assay columns were available for time-window performance analysis.")
    }
    filter(., if_any(all_of(assays_present), ~ !is.na(.x)))
  } %>%
  # join in per‐patient relapse_date + relapsed flag
  left_join(
    final_tbl %>% 
      transmute(
        Patient,
        relapse_date = as.Date(censor_date),
        relapsed     = as.integer(relapsed)
      ),
    by = "Patient"
  )

## Edit if progression date earlier than sample collected, since some patients had multiple
df_sf2 <- df_sf %>%
  select(-relapse_date, -relapsed) %>%           # ① drop the old pair
  left_join(Relapse_dates_full, by = "Patient") %>%
  filter(Progression_date >= sample_date) %>%
  group_by(Patient, sample_date,
           Flow_Binary,
           BM_zscore_only_detection_rate_call,
           Blood_zscore_only_sites_call) %>%
  slice_min(Progression_date, with_ties = FALSE) %>%
  ungroup() %>%
  rename(relapse_date = Progression_date) %>%    # ② now safe to rename
  mutate(relapsed = as.integer(!is.na(relapse_date)))

assays_available <- assays[unname(assays) %in% names(df_sf2)]
missing_assays <- setdiff(names(assays), names(assays_available))
if (length(missing_assays) > 0) {
  message("Skipping absent optional assay columns in non-frontline window analysis: ",
          paste(missing_assays, collapse = ", "))
}

# 3) Define windows (days) to evaluate
windows <- c(90, 180, 365, 730)

# 4) Helper to compute metrics for one assay + one window
calc_metrics <- function(df, assay_label, col_name, win_d) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[col_name]])) %>%
    dplyr::mutate(
      test_pos        = (.data[[col_name]] == 1),
      event_in_window = relapsed == 1 &
        !is.na(relapse_date) &
        relapse_date >= sample_date &
        relapse_date <= sample_date + lubridate::days(win_d)   # << here
    )
  
  tab <- table(
    factor(df2$test_pos,        levels = c(FALSE, TRUE)),
    factor(df2$event_in_window, levels = c(FALSE, TRUE))
  )
  
  tp <- tab["TRUE","TRUE"]; fn <- tab["FALSE","TRUE"]
  fp <- tab["TRUE","FALSE"]; tn <- tab["FALSE","FALSE"]
  
  tibble::tibble(
    Window_days = win_d,
    Assay       = assay_label,
    N_samples   = nrow(df2),
    N_patients  = dplyr::n_distinct(df2$Patient),
    TP = tp, FN = fn, FP = fp, TN = tn,
    Sensitivity = if((tp+fn)>0) tp/(tp+fn) else NA_real_,
    Specificity = if((tn+fp)>0) tn/(tn+fp) else NA_real_,
    PPV         = if((tp+fp)>0) tp/(tp+fp) else NA_real_,
    NPV         = if((tn+fn)>0) tn/(tn+fn) else NA_real_
  )
}


# Inspect all-assay, all-sample time-window metrics using the curated future
# progression-date table. This is retained as a QC/sensitivity view because it
# excludes samples where the updated progression-date table places progression
# before the sample draw.
results <- imap_dfr(assays_available,
                    # .x = column name, .y = assay label
                    .f = function(col_name, assay_label) {
                      map_dfr(windows, function(win_d) {
                        calc_metrics(df_sf2, assay_label, col_name, win_d)
                      })
                    }
)

results %>%
  arrange(Window_days, desc(Sensitivity), desc(Specificity)) %>%
  print(n = Inf)

build_timewindow_metrics <- function(df, assays_available, windows) {
  purrr::map_dfr(windows, function(w) {
    purrr::map_dfr(names(assays_available), function(a) {
      calc_metrics(df, a, assays_available[[a]], w)
    })
  }) %>%
    arrange(Window_days, desc(Sensitivity), desc(Specificity))
}

# Earlier matched cfWGS-tested subsets retained for comparison. The final
# manuscript block below is overwritten with the prospective current labels,
# because revision samples can occur after an earlier patient-level PFS event.
bm_col <- assays[["cfWGS_BM"]]
blood_col <- assays[["cfWGS_Blood"]]
manuscript_assays_available <- assays[unname(assays) %in% names(df_sf)]
df_sf_BM <- df_sf %>% filter(!is.na(.data[[bm_col]]))
df_sf_blood <- df_sf %>% filter(!is.na(.data[[blood_col]]))

legacy_results_BM <- build_timewindow_metrics(df_sf_BM, manuscript_assays_available, windows)
legacy_results_blood <- build_timewindow_metrics(df_sf_blood, manuscript_assays_available, windows)

message("=== BM-cfWGS subset ===")
print(legacy_results_BM)

message("=== Blood-cfWGS subset ===")
print(legacy_results_blood)

## Prospective robust time-window analysis for future test-cohort updates.
##
## Scientific rule:
##   For each non-frontline sample and each prediction window, a sample is
##   evaluable only if one of the following is true:
##     1. a curated progression date occurs after the sample draw and within
##        the prediction window;
##     2. a curated progression date occurs after the sample draw but outside
##        the prediction window; or
##     3. no later curated progression is observed, but clinical follow-up
##        extends at least through the end of the prediction window.
##
## This avoids counting patients with insufficient follow-up as true negatives.
## It also avoids dropping all samples without a future progression date, which
## was too conservative for this small test cohort. These prospective outputs
## provide the calculation used for Extended Data Figures 6I and 8D and
## Supplementary Table 9.
max_date_pair <- function(x, y) {
  out <- pmax(as.numeric(as.Date(x)), as.numeric(as.Date(y)), na.rm = TRUE)
  out[is.infinite(out)] <- NA_real_
  as.Date(out, origin = "1970-01-01")
}

build_prospective_timewindow_labels <- function(sample_df,
                                                final_tbl,
                                                progression_dates,
                                                windows,
                                                followup_dates = NULL) {
  required_cols <- c("Patient", "sample_date")
  missing_cols <- setdiff(required_cols, names(sample_df))
  if (length(missing_cols) > 0) {
    stop("Prospective time-window input is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  patient_outcomes <- final_tbl %>%
    dplyr::transmute(
      Patient = as.character(Patient),
      pfs_event_or_censor_date = as.Date(censor_date),
      patient_relapsed = as.integer(relapsed)
    )

  if (!is.null(followup_dates) && nrow(followup_dates) > 0) {
    patient_outcomes <- patient_outcomes %>%
      dplyr::left_join(
        followup_dates %>%
          dplyr::transmute(
            Patient = as.character(Patient),
            explicit_followup_end_date = as.Date(followup_end_date),
            followup_source = as.character(followup_source)
          ),
        by = "Patient"
      )
  } else {
    patient_outcomes <- patient_outcomes %>%
      dplyr::mutate(
        explicit_followup_end_date = as.Date(NA),
        followup_source = "PFS event/censor date fallback"
      )
  }

  patient_outcomes <- patient_outcomes %>%
    dplyr::mutate(
      followup_end_date = max_date_pair(pfs_event_or_censor_date, explicit_followup_end_date),
      followup_source = dplyr::case_when(
        !is.na(explicit_followup_end_date) & explicit_followup_end_date >= pfs_event_or_censor_date ~ followup_source,
        !is.na(explicit_followup_end_date) & is.na(pfs_event_or_censor_date) ~ followup_source,
        TRUE ~ "PFS event/censor date fallback"
      )
    )

  sample_base <- sample_df %>%
    dplyr::mutate(
      sample_row_id = dplyr::row_number(),
      Patient = as.character(Patient),
      sample_date = as.Date(sample_date)
    ) %>%
    dplyr::left_join(patient_outcomes, by = "Patient")

  future_progression <- sample_base %>%
    dplyr::select(sample_row_id, Patient, sample_date) %>%
    dplyr::left_join(
      progression_dates %>%
        dplyr::transmute(
          Patient = as.character(Patient),
          Progression_date = as.Date(Progression_date)
        ),
      by = "Patient",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(is.na(Progression_date) | Progression_date >= sample_date) %>%
    dplyr::group_by(sample_row_id) %>%
    dplyr::summarise(
      future_progression_date = if (all(is.na(Progression_date))) {
        as.Date(NA)
      } else {
        min(Progression_date, na.rm = TRUE)
      },
      n_future_progression_dates = sum(!is.na(Progression_date)),
      .groups = "drop"
    )

  labelled_samples <- sample_base %>%
    dplyr::left_join(future_progression, by = "sample_row_id") %>%
    dplyr::mutate(
      days_followup_after_sample = as.numeric(followup_end_date - sample_date),
      days_to_future_progression = as.numeric(future_progression_date - sample_date),
      no_future_progression_with_adequate_followup = !is.na(followup_end_date)
    )

  purrr::map_dfr(windows, function(win_d) {
    labelled_samples %>%
      dplyr::mutate(
        Window_days = win_d,
        window_end_date = sample_date + lubridate::days(win_d),
        event_in_window = !is.na(future_progression_date) &
          future_progression_date <= window_end_date,
        evaluable_in_window = event_in_window |
          (!is.na(future_progression_date) & future_progression_date > window_end_date) |
          (is.na(future_progression_date) & !is.na(followup_end_date) & followup_end_date >= window_end_date),
        nonevaluable_reason = dplyr::case_when(
          evaluable_in_window ~ NA_character_,
          is.na(followup_end_date) ~ "missing follow-up/censor date",
          followup_end_date < sample_date ~ "follow-up/censor date before sample date",
          followup_end_date < window_end_date ~ "insufficient follow-up for window",
          TRUE ~ "not evaluable for unknown reason"
        )
      )
  })
}

calc_prospective_metrics <- function(labelled_df, assay_label, col_name, win_d) {
  if (!col_name %in% names(labelled_df)) {
    return(tibble::tibble())
  }

  df2 <- labelled_df %>%
    dplyr::filter(
      Window_days == win_d,
      evaluable_in_window,
      !is.na(.data[[col_name]])
    ) %>%
    dplyr::mutate(test_pos = .data[[col_name]] == 1)

  tab <- table(
    factor(df2$test_pos, levels = c(FALSE, TRUE)),
    factor(df2$event_in_window, levels = c(FALSE, TRUE))
  )

  tp <- tab["TRUE", "TRUE"]
  fn <- tab["FALSE", "TRUE"]
  fp <- tab["TRUE", "FALSE"]
  tn <- tab["FALSE", "FALSE"]

  tibble::tibble(
    Window_days = win_d,
    Assay = assay_label,
    N_samples = nrow(df2),
    N_patients = dplyr::n_distinct(df2$Patient),
    N_nonevaluable_with_assay = labelled_df %>%
      dplyr::filter(Window_days == win_d, !evaluable_in_window, !is.na(.data[[col_name]])) %>%
      nrow(),
    TP = tp, FN = fn, FP = fp, TN = tn,
    Sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    Specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    PPV = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    NPV = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  )
}

build_prospective_metrics <- function(labelled_df, assay_lookup, windows) {
  purrr::map_dfr(windows, function(w) {
    purrr::map_dfr(names(assay_lookup), function(a) {
      calc_prospective_metrics(labelled_df, a, assay_lookup[[a]], w)
    })
  }) %>%
    dplyr::arrange(Window_days, dplyr::desc(Sensitivity), dplyr::desc(Specificity))
}

prospective_timewindow_labels <- build_prospective_timewindow_labels(
  sample_df = df_sf %>% dplyr::select(-relapse_date, -relapsed),
  final_tbl = final_tbl,
  progression_dates = Relapse_dates_full,
  windows = windows,
  followup_dates = patient_followup_dates
)

prospective_timewindow_dir <- file.path(outdir, "prospective_timewindow_qc")
dir.create(prospective_timewindow_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(
  prospective_timewindow_labels,
  file.path(prospective_timewindow_dir, "prospective_timewindow_sample_level_labels.csv")
)

prospective_assays_available <- assays[unname(assays) %in% names(prospective_timewindow_labels)]
prospective_results_BM <- prospective_timewindow_labels %>%
  dplyr::filter(!is.na(.data[[bm_col]])) %>%
  build_prospective_metrics(prospective_assays_available, windows)
prospective_results_blood <- prospective_timewindow_labels %>%
  dplyr::filter(!is.na(.data[[blood_col]])) %>%
  build_prospective_metrics(prospective_assays_available, windows)

readr::write_csv(
  prospective_results_BM,
  file.path(prospective_timewindow_dir, "prospective_BM_cfWGS_timewindow_results.csv")
)
readr::write_csv(
  prospective_results_blood,
  file.path(prospective_timewindow_dir, "prospective_blood_cfWGS_timewindow_results.csv")
)
writexl::write_xlsx(
  list(
    BM_models = prospective_results_BM,
    Blood_models = prospective_results_blood
  ),
  path = file.path(prospective_timewindow_dir, "prospective_Supplementary_Table_9_timewindow_results.xlsx")
)

prospective_evaluability_summary <- prospective_timewindow_labels %>%
  dplyr::group_by(Window_days) %>%
  dplyr::summarise(
    Samples_total = dplyr::n_distinct(sample_row_id),
    Patients_total = dplyr::n_distinct(Patient),
    Samples_evaluable = dplyr::n_distinct(sample_row_id[evaluable_in_window]),
    Samples_not_evaluable = dplyr::n_distinct(sample_row_id[!evaluable_in_window]),
    Events_in_window = sum(event_in_window, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(
  prospective_evaluability_summary,
  file.path(prospective_timewindow_dir, "prospective_timewindow_evaluability_summary.csv")
)

message("Prospective robust time-window QC outputs written to: ", prospective_timewindow_dir)
print(prospective_evaluability_summary)

## Final manuscript source for time-window panels and Supplementary Table 9.
##
legacy_results_dir <- file.path(outdir, "legacy_patient_level_pfs_timewindow_qc")
dir.create(legacy_results_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(
  legacy_results_BM,
  file.path(legacy_results_dir, "legacy_BM_cfWGS_timewindow_results.csv")
)
readr::write_csv(
  legacy_results_blood,
  file.path(legacy_results_dir, "legacy_blood_cfWGS_timewindow_results.csv")
)

results_BM <- prospective_results_BM
results_blood <- prospective_results_blood

count_relapses_by_window <- function(df, wins = c(180, 365)) {
  purrr::map_dfr(wins, function(w) {
    df %>%
      dplyr::mutate(event_in_window = relapsed == 1 &
                      !is.na(relapse_date) &
                      relapse_date >= sample_date &
                      relapse_date <= sample_date + lubridate::days(w)) %>%
      dplyr::summarise(
        Window_days       = w,
        Samples_relapsed  = sum(event_in_window, na.rm = TRUE),
        Patients_relapsed = dplyr::n_distinct(Patient[event_in_window])
      )
  })
}

count_prospective_events_by_window <- function(labelled_df, assay_cols, wins) {
  purrr::map_dfr(wins, function(w) {
    labelled_df %>%
      dplyr::filter(
        Window_days == w,
        evaluable_in_window,
        dplyr::if_any(dplyr::all_of(assay_cols), ~ !is.na(.x))
      ) %>%
      dplyr::summarise(
        Window_days = w,
        N_patients = dplyr::n_distinct(Patient[event_in_window]),
        N_samples = sum(event_in_window, na.rm = TRUE),
        .groups = "drop"
      )
  })
}

event_counts_BM <- count_prospective_events_by_window(
  prospective_timewindow_labels,
  bm_col,
  windows
)
event_counts_blood <- count_prospective_events_by_window(
  prospective_timewindow_labels,
  blood_col,
  windows
)

message("Using prospective current time-window metrics for ED6I, ED8D, and Supplementary Table 9.")

wins <- c(180, 365)

# Narrative counts that match the exact subsets used for the plotted
# denominator tables.
summ_bm <- event_counts_BM %>%
  dplyr::filter(Window_days %in% wins) %>%
  dplyr::select(Window_days, Patients_relapsed = N_patients, Samples_relapsed = N_samples)

# B) Blood-only counts using the cfWGS_Blood evaluability denominators
summ_blood <- event_counts_blood %>%
  dplyr::filter(Window_days %in% wins) %>%
  dplyr::select(Window_days, Patients_relapsed = N_patients, Samples_relapsed = N_samples)

# C) Union of evaluable rows with BM or blood cfWGS assays
summ_any <- count_prospective_events_by_window(
  prospective_timewindow_labels,
  c(bm_col, blood_col),
  wins
)


## Get event counts 
# define the windows of interest
windows <- c(90, 180, 365, 730)

# this will count, for each window:
# - how many distinct patients relapsed within that window
# - how many samples fall into that window

print(event_counts_BM)

## Make figure 
# ────────────────────────────────────────────────────────────────────────────
# 1) Prepare the data
# ────────────────────────────────────────────────────────────────────────────
sens_BM_df <- results_BM %>%
  # ED6I is restricted to the BM-derived cfWGS model and its MFC comparator.
  filter(Assay %in% c("cfWGS_BM", "Flow")) %>%

  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow        = "MFC",
      cfWGS_BM    = "cfWGS",
    ),
    Assay = factor(Assay, levels = c("cfWGS", "MFC"))
  )
# ────────────────────────────────────────────────────────────────────────────
# 2) Colours & theme (match the existing style)
# ────────────────────────────────────────────────────────────────────────────
custom_cols <- c(
  "90 days"  = "#440154FF",
  "180 days" = "#31688EFF",
  "365 days" = "#35B779FF",
  "730 days" = "#E69F00FF"
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title      = element_text(size = 11),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.line       = element_line(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(10, 10, 30, 10)
  )

# ────────────────────────────────────────────────────────────────────────────
# 3) Build the grouped bar‑plot
# ────────────────────────────────────────────────────────────────────────────
p_sens_bm <- ggplot(sens_BM_df,
                    aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow-up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 14)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_6_Fig_sensitivity_windows_BM_test_cohort_updated2.png",
       plot = p_sens_bm,
       width = 4.75, height = 6.25, dpi = 500)

# MANUSCRIPT OUTPUT: Extended Data Figure 6I
# Test-cohort BM sensitivity across relapse follow-up windows.
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Fig_6_Fig_sensitivity_windows_BM_test_cohort_updated2.png",
  artifact_id = "EDFIG6I",
  role = "figure_panel_png",
  description = "Extended Data Figure 6I: BM test-cohort sensitivity across relapse follow-up windows.",
  script_name = "4_1_Survival_Analysis.R"
)


### Now remake for blood muts
sens_blood_df <- results_blood %>%
  # if you want to drop the blood‑only assay, uncomment:
  # filter(Assay != "cfWGS_Blood") %>%
  
  # turn Window_days into a nice factor
  mutate(
    Timepoint = factor(
      Window_days,
      levels = c(90, 180, 365, 730),
      labels = c("90 days", "180 days", "365 days", "730 days")
    ),
    Sens_pct = Sensitivity * 100,
    Assay = recode(
      Assay,
      Flow                   = "MFC",
      cfWGS_Blood            = "cfWGS (Sites Model)",
      cfWGS_Blood_Combined   = "cfWGS (Combined Model)",
    )
  )

sens_blood_df <- sens_blood_df %>%
  filter(Assay != "cfWGS_BM") %>%
  mutate(
    Assay = factor(
      Assay,
      levels = c(
        "cfWGS (Sites Model)",
        "cfWGS (Combined Model)",
        "MFC"
      )
    )
  )

p_sens_blood <- ggplot(sens_blood_df,
                       aes(x = Assay, y = Sens_pct, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           colour   = "black",
           size     = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Sens_pct)),
            position = position_dodge(width = 0.8),
            vjust    = -0.3,
            size     = 3.5) +
  scale_fill_manual(
    name   = "Window",
    values = custom_cols,
    breaks = names(custom_cols)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(scale = 1)
  ) +
  labs(
    title = "Sensitivity of MRD assays over\nfollow-up windows (Test Cohort)",
    x     = "Assay",
    y     = "Sensitivity"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ────────────────────────────────────────────────────────────────────────────
# 4) (Optional) Save
# ────────────────────────────────────────────────────────────────────────────
ggsave("Final Tables and Figures/Supp_Fig_8_Fig_sensitivity_windows_blood_test_cohort3.png",
       plot = p_sens_blood,
       width = 5, height = 5, dpi = 500)

# MANUSCRIPT OUTPUT: Extended Data Figure 8D
# Test-cohort blood sensitivity across relapse follow-up windows.
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Fig_8_Fig_sensitivity_windows_blood_test_cohort3.png",
  artifact_id = "EDFIG8D",
  role = "figure_panel_png",
  description = "Extended Data Figure 8D: blood test-cohort sensitivity across relapse follow-up windows.",
  script_name = "4_1_Survival_Analysis.R"
)





### Export this
# full results (all patients with any assay)
#write_csv(
#  results,
#  file.path(outdir, "all_assays_timewindow_results.csv")
#)

# BM‐cfWGS subset
write_csv(
  results_BM,
  file.path(outdir, "BM_cfWGS_timewindow_results2.csv")
)

# blood‐cfWGS subset
write_csv(
  results_blood,
  file.path(outdir, "blood_cfWGS_timewindow_results2.csv")
)

## export 
# --- Clean BM results ---
results_BM_clean <- results_BM %>%
  # remove blood assays
  filter(!Assay %in% c("cfWGS_Blood", "cfWGS_Blood_Combined")) %>%
  # rename Flow → MFC
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      TRUE ~ Assay
    )
  )

# --- Clean Blood results ---
results_blood_clean <- results_blood %>%
  # remove BM assays
  filter(!Assay %in% c("cfWGS_BM")) %>%
  # rename assays
  mutate(
    Assay = case_when(
      Assay == "Flow" ~ "MFC",
      Assay == "cfWGS_Blood" ~ "cfWGS_Blood (Sites Model)",
      Assay == "cfWGS_Blood_Combined" ~ "cfWGS_Blood (Combined Model)",
      TRUE ~ Assay
    )
  )

# --- Export to Excel ---
export_list <- list(
  "BM_models"    = results_BM_clean,
  "Blood_models" = results_blood_clean
)

# Write to Excel
write_xlsx(
  export_list,
  path = file.path("Final Tables and Figures/Supplementary_Table_9_timewindow_results_test_cohort.xlsx")
)

# Historical reduced export/staging block.
# `results_BM` and `results_blood` contain the prospective calculations above,
# but the two filters immediately before this block remove comparator rows. The
# retained Supplementary Table 9 instead uses the complete 16-row-per-sheet
# workbook written to `prospective_timewindow_qc/`. This block is retained
# unchanged for now; do not treat its staged workbook as the final table.
ms_copy_artifact(
  source_path = file.path(
    prospective_timewindow_dir,
    "prospective_Supplementary_Table_9_timewindow_results.xlsx"
  ),
  artifact_id = "STABLE9",
  role = "supplementary_table_xlsx",
  description = paste(
    "Supplementary Table 9: complete prospective time-window results for",
    "BM and blood models, including all comparator rows."
  ),
  script_name = "4_1_Survival_Analysis.R"
)

# event counts for BM‐cfWGS
write_csv(
  event_counts_BM,
  file.path(outdir, "BM_cfWGS_event_counts2.csv")
)

write_csv(
  event_counts_blood,
  file.path(outdir, "blood_cfWGS_event_counts2.csv")
)




### See at what pont an increase occured 
d2m <- function(days, digits = 1) round(as.numeric(days) / 30.44, digits)

# 1) choose the probability column you want to analyze:
assay_prob <- "BM_zscore_only_detection_rate_prob"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_BM %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )




##### Now redo for blood muts 
# 1) choose the probability column you want to analyze:
assay_prob <- "Blood_zscore_only_sites_call"  # or "BM_zscore_only_detection_rate_prob"
# Optional: restrict to surveillance timepoints only (post-induction / ASCT / maintenance)
surveillance_only <- FALSE
surv_regex <- "(post[_ -]?induction|pre[_ -]?(asct|transplant)|post[_ -]?(asct|transplant)|maintenance)"

# ==== BASE (relapse patients only, consistent with figure) ====
df0 <- time_to_relapse_blood %>%
  mutate(timepoint_info = tolower(timepoint_info)) %>%
  filter(progress_status == "Relapse",
         !is.na(.data[[assay_prob]]))

if (surveillance_only) {
  df0 <- df0 %>% filter(grepl(surv_regex, timepoint_info))
}

# ==== A) STRICT pre-progression monitoring set ====
# Use days_before_event > 0 to enforce strictly pre-progression (same logic as sample_date < censor_date)
df_pre <- df0 %>%
  filter(days_before_event > 0) %>%
  arrange(Patient, sample_date)

n_patients <- n_distinct(df_pre$Patient)
n_samples  <- nrow(df_pre)

sample_timing_stats <- df_pre %>%
  summarise(
    median_days_before = median(days_before_event, na.rm = TRUE),
    iqr_days_before    = IQR(days_before_event,    na.rm = TRUE),
    min_days_before    = min(days_before_event,    na.rm = TRUE),
    max_days_before    = max(days_before_event,    na.rm = TRUE),
    .groups = "drop"
  )

sample_timing_stats <- sample_timing_stats %>%
  mutate(
    median_days_before_mo = d2m(median_days_before),
    iqr_days_before_mo    = d2m(iqr_days_before),
    min_days_before_mo    = d2m(min_days_before),
    max_days_before_mo    = d2m(max_days_before)
  )

# ==== B) Per-patient nadir & first increase (pre-progression only) ====
advance_df <- df0 %>%
  # Keep rows up to relapse (+30d tolerance) so patients with slight date mismatches are still present
  dplyr::filter(sample_date <= (censor_date + 30)) %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(sample_date, .by_group = TRUE) %>%
  dplyr::group_modify(~{
    df <- .
    prog_date <- df$censor_date[1]
    
    # STRICT pre-relapse rows for nadir / first-increase logic
    df_pre <- dplyr::filter(df, sample_date < prog_date)
    
    if (nrow(df_pre) == 0L) {
      return(tibble::tibble(
        nadir_date = as.Date(NA),  nadir_prob = NA_real_,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    # Nadir = minimum probability before relapse
    nadir_row <- df_pre %>%
      dplyr::slice_min(.data[[assay_prob]], with_ties = FALSE)
    
    nadir_date <- nadir_row$sample_date
    nadir_prob <- nadir_row[[assay_prob]]
    
    # First increase after nadir (strictly later in time AND strictly higher prob)
    post_nadir <- df_pre %>%
      dplyr::filter(sample_date > nadir_date,
                    .data[[assay_prob]] > nadir_prob) %>%
      dplyr::slice_head(n = 1)
    
    if (nrow(post_nadir) == 0L) {
      return(tibble::tibble(
        nadir_date = nadir_date,  nadir_prob = nadir_prob,
        first_inc_date = as.Date(NA), first_inc_prob = NA_real_,
        days_to_first_increase  = NA_real_,
        days_before_progression = NA_real_,
        prog_date = prog_date
      ))
    }
    
    first_inc_date <- post_nadir$sample_date
    first_inc_prob <- post_nadir[[assay_prob]]
    
    tibble::tibble(
      nadir_date = nadir_date,
      nadir_prob = nadir_prob,
      first_inc_date = first_inc_date,
      first_inc_prob = first_inc_prob,
      days_to_first_increase  = as.numeric(first_inc_date - nadir_date),
      days_before_progression = as.numeric(prog_date - first_inc_date),
      prog_date = prog_date
    )
  }) %>%
  dplyr::ungroup()


# Nadir timing relative to progression
nadir_timing_stats <- advance_df %>%
  mutate(days_nadir_before = as.numeric(prog_date - nadir_date)) %>%
  summarise(
    median_nadir_days = median(days_nadir_before, na.rm = TRUE),
    iqr_nadir_days    = IQR(days_nadir_before,    na.rm = TRUE),
    min_nadir_days    = min(days_nadir_before,    na.rm = TRUE),
    max_nadir_days    = max(days_nadir_before,    na.rm = TRUE),
    .groups = "drop"
  )

nadir_timing_stats <- nadir_timing_stats %>%
  mutate(
    median_nadir_mo = d2m(median_nadir_days),
    iqr_nadir_mo    = d2m(iqr_nadir_days),
    min_nadir_mo    = d2m(min_nadir_days),
    max_nadir_mo    = d2m(max_nadir_days)
  )

# First increase timing (after nadir, and its lead time to progression)
increase_stats2 <- advance_df %>%
  summarise(
    median_after_nadir = median(days_to_first_increase,   na.rm = TRUE),
    iqr_after_nadir    = IQR(days_to_first_increase,      na.rm = TRUE),
    min_after_nadir    = min(days_to_first_increase,      na.rm = TRUE),
    max_after_nadir    = max(days_to_first_increase,      na.rm = TRUE),
    median_before_prog = median(days_before_progression,  na.rm = TRUE),
    iqr_before_prog    = IQR(days_before_progression,     na.rm = TRUE),
    min_before_prog    = min(days_before_progression,     na.rm = TRUE),
    max_before_prog    = max(days_before_progression,     na.rm = TRUE),
    .groups            = "drop"
  )

increase_stats2 <- increase_stats2 %>%
  mutate(
    median_after_nadir_mo = d2m(median_after_nadir),
    iqr_after_nadir_mo    = d2m(iqr_after_nadir),
    min_after_nadir_mo    = d2m(min_after_nadir),
    max_after_nadir_mo    = d2m(max_after_nadir),
    median_before_prog_mo = d2m(median_before_prog),
    iqr_before_prog_mo    = d2m(iqr_before_prog),
    min_before_prog_mo    = d2m(min_before_prog),
    max_before_prog_mo    = d2m(max_before_prog)
  )











### Analyst note on non-frontline progression labels
# The non-frontline/test-cohort timing summaries above use the curated
# progression dates available in `Relapse dates cfWGS updated2.csv`. They do
# not further distinguish clinical versus biochemical progression unless that
# distinction is present in the staged input file. No additional patient
# filtering is applied here.

# ═════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY: Script Output and Results
# ═════════════════════════════════════════════════════════════════════════════
#
# This comprehensive survival analysis script has produced the following:
#
# ─────────────────────────────────────────────────────────────────────────────
# PRIMARY OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────
#
# 1. KAPLAN-MEIER SURVIVAL CURVES (PNG files @ 500 DPI)
#    Location: detection_progression_updated6/[timepoint]/
#    Files: KM_[Assay]_[Timepoint]_updated_no_CI.png
#    - Diagnosis: 4 assays (EasyM, clonoSEQ, MFC, cfWGS)
#    - Post-ASCT: 7 assays (+ BM/blood variants)
#    - Maintenance-1yr: 7 assays (+ BM/blood variants)
#    Each plot shows: PFS curves, risk table, log-rank p-value
#
# 2. SENSITIVITY TABLES (CSV files)
#    Location: Output directory root
#    - frontline_postASCT_sensitivity.csv: % MRD+ among relapsed at post-ASCT
#    - frontline_1yr_sensitivity.csv: % MRD+ among relapsed at 1-year
#    - frontline_postASCT_sens_BMcfWGS.csv: Sensitivity in BM-cfWGS subset
#    - frontline_1yr_sens_BMcfWGS.csv: Sensitivity in BM-cfWGS subset
#    - frontline_postASCT_sens_bloodcfWGS.csv: Sensitivity in blood-cfWGS subset
#    - frontline_1yr_sens_bloodcfWGS.csv: Sensitivity in blood-cfWGS subset
#
# 3. SENSITIVITY BARPLOTS (PNG files @ 500 DPI)
#    Location: Final Tables and Figures/
#    - Supp_6A_Fig_sensitivity_by_tech_training3.png (BM-subset analysis)
#    - Supp_8A_Fig_sensitivity_by_tech_training_blood2.png (blood-subset analysis)
#    Each plot shows: Grouped bars comparing assay sensitivities at post-ASCT vs 1-year
#
# ─────────────────────────────────────────────────────────────────────────────
# SECONDARY ANALYSES
# ─────────────────────────────────────────────────────────────────────────────
#
# 4. FRONTLINE LANDMARK RFS/HAZARD-RATIO SUMMARIES
#    - 24-month RFS by cfWGS BM at 1-year maintenance
#    - Median RFS and hazard ratios by assay
#    - Comparable metrics for Flow, clonoSEQ, EasyM
#
# 5. NON-FRONTLINE/TEST-COHORT TIME-WINDOW ANALYSES
#    - Sensitivity, specificity, PPV, and NPV across fixed relapse windows
#    - BM cfWGS subset results used for Extended Data Figure 6I
#    - Blood cfWGS subset results used for Extended Data Figure 8D
#    - Reduced compatibility workbook exported after the complete prospective
#      Supplementary Table 9 workbook is written above
#
# 6. CONSOLE SUMMARY STATISTICS
#    - Patient/sample counts and demographics
#    - Timing analysis: Days from baseline to collection/relapse
#    - Nadir timing and relapse progression timing
#    - Formatted values used to check methods/results text
#
# ─────────────────────────────────────────────────────────────────────────────
# KEY STATISTICAL FINDINGS EXPORTED FOR VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────
#
# - PFS probability at each timepoint (KM curves)
# - Log-rank p-values testing MRD+/- difference
# - Hazard ratios with 95% confidence intervals
# - Sensitivity metrics (% detection among relapsers)
# - Time-to-event distributions (days and months)
#
# ═════════════════════════════════════════════════════════════════════════════

cat("\n")
cat(strrep("═", 80), "\n")
cat("SURVIVAL ANALYSIS COMPLETE\n")
cat("Scripts completed: All KM curves, sensitivity analyses, and barplots generated\n")
cat(strrep("═", 80), "\n\n")
