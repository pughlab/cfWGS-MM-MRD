## ============================================================
## SCRIPT PURPOSE: Process and Optimize EasyM Proteomic MRD Data
## ============================================================
# This active intermediate script processes EasyM proteomic MRD data from the
# clinical collaborator files and creates the EasyM call table required by
# downstream manuscript-output scripts.
#
# Manuscript outputs created/updated:
#   - None directly. This script does not copy a final figure/table into
#     final_manuscript_objects/ because no current manuscript artifact is mapped
#     directly to this script in docs/manuscript_artifact_source_map.tsv.
#
# Required downstream role:
#   - 3_2_Plot_optimal_cutoff_and_clinical_concordance.R uses the prespecified
#     isotype-specific EasyM reference-threshold calls for Figure 3E, Figure 4C,
#     Figure 4D, and Supplementary Table 8.
#   - 4_1_Survival_Analysis.R uses the same isotype-specific EasyM calls for the
#     survival and relapse-detection panels in Figure 3F/Figure 4E and Extended
#     Data Figures 6/8.
#
# Key outputs consumed downstream:
#   - Output_EasyM_MRD_analysis_2025/EasyM_all_samples_with_optimized_calls.csv
#   - Output_EasyM_MRD_analysis_2025/EasyM_threshold_values_by_timepoint.csv
#   - Output_EasyM_MRD_analysis_2025/landmark_analyses/
#
# Main inputs:
#   - Aimee additional data/RAPID NOVOR VALUES with values with relapse.csv
#   - Aimee additional data/RAPID NOVOR pos-neg with relapse.csv
#   - Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds
#   - Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds
#
# EasyM units and manuscript call:
#   - EasyM_value is recorded in percentage points relative to the patient's
#     baseline bone marrow immunoglobulin signal: 1.00 means 1% of baseline and
#     0.05 means 0.05% of baseline.
#   - The manuscript call is negative at <=1% for IgG and <=0.05% for
#     IgA/light-chain disease. EasyM_reference_threshold_binary stores this
#     prespecified call; the optimized timepoint-specific calls are exploratory.
#
# Support-only outputs:
#   The EasyM-only plots, exploratory landmark tables, optimized cutoff RDS
#   files, and joint risk summaries are retained for audit and interpretation.
#   They are not final manuscript artifacts unless the source map is updated to
#   assign a specific figure/table to them.
#
# Workflow:
#   0. Setup and output directories.
#   1. Helper functions for robust IO, data shaping, and survival summaries.
#   2. Load EasyM quantitative and binary source files.
#   3. Generate EasyM-only descriptive support plots.
#   4. Merge EasyM with cfWGS calls/probabilities.
#   5. Run descriptive landmark survival/comparison analyses.
#   6. Derive exploratory optimized thresholds and the prespecified
#      isotype-specific reference-threshold calls.
#   7. Export the downstream EasyM call table and threshold reference tables.
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/3_1_A_Process_and_optimize_EasyM.R
# ============================================================
# Pipeline status:
#   Active required intermediate. Keep this script in the command-line pipeline
#   even though it does not directly own a final manuscript figure/table.
#

## ============================================================
## 0) Setup
## ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(broom)
  library(scales)
})

## ---- Project paths ----
# Input: project-relative raw EasyM data from clinical collaborators.
# Output: processed EasyM call/threshold tables used by 3_2 and 4_1, plus
# support-only figures/tables for audit.
easym_input_dir <- "Aimee additional data"

# Main output directory for the EasyM analyses and downstream call table.
out_dir  <- "Output_EasyM_MRD_analysis_2025/"

# Create directory if it doesn't exist
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "landmark_analyses"), showWarnings = FALSE, recursive = TRUE)

# Print confirmation
cat(sprintf("\n✓ Output directory created: %s\n\n", normalizePath(out_dir)))

## ---- Visit ordering (keep consistent everywhere) ----
# Standardized visit labels used across all scripts
# V1 = diagnosis, V3 = post-induction, V5 = post-transplant, 
# V7 = 1-year maintenance, R = relapse
visit_levels <- c("V1","V3","V5","V7","V8","V9","V10","V11","V12","V13","R")

## ---- Landmark settings (edit to match abstract) ----
## If you only have visit labels (no dates), use a visit landmark (eg "V7").
landmark_visit <- "V7"

## Fixed horizon PFS estimate (months). Set the horizon to match the reported survival estimate.
fixed_horizon_months <- 24

## Proteomics high/low thresholding rule for “joint modality” groups
## Options implemented below: "median", "predefined"
proteomics_threshold_rule <- "median"

## If you want a predefined cutoff, set it here and switch rule to "predefined"
proteomics_predefined_cutoff <- NA_real_

## ============================================================
## 1) Helpers (robust IO + consistent shaping)
## ============================================================
# These functions ensure consistent data formatting and error checking
# throughout the analysis pipeline
stop_if_missing <- function(df, cols, df_name = "data") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      df_name, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
}

read_csv_safely <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

values_to_long <- function(df, tech_label) {
  ## Assumes first column is patient ID, remaining columns are visits
  df %>%
    rename(patient = 1) %>%
    mutate(patient = as.character(patient)) %>%
    pivot_longer(
      cols      = -patient,
      names_to  = "visit",
      values_to = "value"
    ) %>%
    mutate(
      technology = tech_label,
      visit      = factor(visit, levels = visit_levels)
    )
}

pos_to_long <- function(df, tech_label) {
  ## Assumes first column is patient ID, remaining columns are visits
  df %>%
    rename(patient = 1) %>%
    mutate(patient = as.character(patient)) %>%
    pivot_longer(
      cols      = -patient,
      names_to  = "visit",
      values_to = "code"
    ) %>%
    mutate(
      technology = tech_label,
      visit      = factor(visit, levels = visit_levels)
    )
}

code_to_mrd_status <- function(code) {
  ## MRD status encoding: 100 = MRD+, else MRD-
  ## Also supports strings like "pos"/"neg" if they show up later.
  code_chr <- as.character(code)
  
  case_when(
    is.na(code_chr) ~ NA_character_,
    code_chr %in% c("100", "1", "POS", "Pos", "pos", "MRD+", "positive", "Positive") ~ "MRD+",
    code_chr %in% c("0", "NEG", "Neg", "neg", "MRD-", "negative", "Negative") ~ "MRD-",
    TRUE ~ {
      ## fall back: numeric > 0 treated as positive
      suppressWarnings({
        x <- as.numeric(code_chr)
      })
      ifelse(!is.na(x) & x > 0, "MRD+", "MRD-")
    }
  )
}

classify_easym_reference_isotype_text <- function(details, subtype) {
  text <- str_to_lower(str_squish(paste(details, subtype, sep = " ")))

  case_when(
    is.na(text) | text == "" ~ NA_character_,
    str_detect(text, "\\bigg\\b|ig\\s*g") ~ "IgG",
    str_detect(text, "\\biga\\b|ig\\s*a") ~ "IgA",
    str_detect(
      text,
      "light\\s*chain|light-chain|free\\s+kappa|free\\s+lambda|kappa\\s+light|lambda\\s+light|lc\\s*only"
    ) ~ "Light Chain",
    TRUE ~ NA_character_
  )
}

infer_easym_reference_isotype <- function(details, subtype, iga, igg, igm) {
  text_isotype <- classify_easym_reference_isotype_text(details, subtype)
  ifelse(
    !is.na(text_isotype),
    text_isotype,
    case_when(
      !is.na(igg) & (is.na(iga) | igg >= iga) & (is.na(igm) | igg >= igm) ~ "IgG",
      !is.na(iga) & (is.na(igg) | iga > igg) & (is.na(igm) | iga >= igm) ~ "IgA",
      TRUE ~ NA_character_
    )
  )
}

easym_reference_threshold_percent <- function(isotype) {
  case_when(
    isotype == "IgG" ~ 1.00,
    isotype %in% c("IgA", "Light Chain") ~ 0.05,
    TRUE ~ NA_real_
  )
}

build_patient_easym_reference_isotypes <- function(dat) {
  required <- c("Patient", "Timepoint", "Date", "timepoint_info", "DETAILS", "Subtype", "IgA", "IgG", "IgM")
  stop_if_missing(dat, required, "EasyM/cfWGS joined clinical table")

  candidates <- dat %>%
    mutate(
      Patient = as.character(Patient),
      Timepoint = as.character(Timepoint),
      diagnosis_priority = case_when(
        Timepoint == "01" ~ 1L,
        str_to_lower(as.character(timepoint_info)) %in% c("diagnosis", "baseline") ~ 2L,
        TRUE ~ 3L
      ),
      EasyM_reference_isotype = infer_easym_reference_isotype(
        DETAILS, Subtype,
        suppressWarnings(as.numeric(IgA)),
        suppressWarnings(as.numeric(IgG)),
        suppressWarnings(as.numeric(IgM))
      ),
      EasyM_reference_isotype_source = case_when(
        !is.na(classify_easym_reference_isotype_text(DETAILS, Subtype)) ~ "clinical_text_DETAILS_or_Subtype",
        !is.na(EasyM_reference_isotype) ~ "baseline_heavy_chain_dominance",
        TRUE ~ "unclassified"
      )
    ) %>%
    filter(!is.na(Patient), !is.na(EasyM_reference_isotype)) %>%
    arrange(Patient, diagnosis_priority, Timepoint, Date)

  patient_isotypes <- candidates %>%
    group_by(Patient) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      Patient,
      EasyM_reference_isotype,
      EasyM_reference_threshold_percent = easym_reference_threshold_percent(EasyM_reference_isotype),
      EasyM_reference_isotype_source,
      EasyM_reference_source_timepoint = Timepoint,
      EasyM_reference_source_details = DETAILS,
      EasyM_reference_source_subtype = Subtype
    )

  patient_isotypes
}

compute_easym_reference_contingency <- function(df, pred_col, pred_label, cohort_label, timepoint_label) {
  required <- c(pred_col, "EasyM_reference_threshold_binary")
  stop_if_missing(df, required, "EasyM reference contingency input")

  eval_df <- df %>%
    filter(!is.na(.data[[pred_col]]), !is.na(EasyM_reference_threshold_binary))

  if (nrow(eval_df) == 0) {
    return(tibble())
  }

  pred <- as.integer(eval_df[[pred_col]])
  truth <- as.integer(eval_df$EasyM_reference_threshold_binary)
  tp <- sum(pred == 1L & truth == 1L, na.rm = TRUE)
  fp <- sum(pred == 1L & truth == 0L, na.rm = TRUE)
  fn <- sum(pred == 0L & truth == 1L, na.rm = TRUE)
  tn <- sum(pred == 0L & truth == 0L, na.rm = TRUE)
  n <- tp + fp + fn + tn

  tibble(
    Cohort = cohort_label,
    Timepoint = timepoint_label,
    Pred_Column = pred_col,
    Pred_Label = pred_label,
    Comparator = "EasyM reference threshold",
    N = n,
    Agree = tp + tn,
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Concordance = if_else(n > 0, (tp + tn) / n, NA_real_),
    Sensitivity = if_else((tp + fn) > 0, tp / (tp + fn), NA_real_),
    Specificity = if_else((tn + fp) > 0, tn / (tn + fp), NA_real_),
    PPV = if_else((tp + fp) > 0, tp / (tp + fp), NA_real_),
    NPV = if_else((tn + fn) > 0, tn / (tn + fn), NA_real_),
    Matrix = sprintf("TN=%d; FP=%d; FN=%d; TP=%d", tn, fp, fn, tp)
  )
}

theme_mrd <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.6),
    axis.ticks       = element_line(colour = "black", linewidth = 0.6),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "plain", size = 14),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

## ============================================================
## 2) EasyM ONLY: load quantitative + binary EasyM data
## ============================================================
# EasyM provides two types of data:
# 1. Quantitative: residual immunoglobulin signal as a percentage of the
#    patient's baseline bone marrow immunoglobulin signal
# 2. Binary: MRD positive (100) vs negative (0) calls
# Both are measured longitudinally across treatment visits

## EasyM source files are staged inside the repository/project so command-line
## and Code Ocean runs do not depend on user-specific sibling folders.
EasyM_quant_path <- file.path(easym_input_dir, "RAPID NOVOR VALUES with values with relapse.csv")
EasyM_bin_path   <- file.path(easym_input_dir, "RAPID NOVOR pos-neg with relapse.csv")

EasyM_values_raw <- read_csv_safely(EasyM_quant_path)
EasyM_pos_raw    <- read_csv_safely(EasyM_bin_path)

EasyM_values_long <- values_to_long(EasyM_values_raw, tech_label = "EasyM") %>%
  filter(!is.na(value)) %>%
  mutate(value = suppressWarnings(as.numeric(value)))

EasyM_pos_long <- pos_to_long(EasyM_pos_raw, tech_label = "EasyM") %>%
  mutate(mrd_status = code_to_mrd_status(code)) %>%
  filter(!is.na(mrd_status))

## Optional: join EasyM quant + binary at the same patient/visit when both exist
EasyM_joined <- EasyM_values_long %>%
  select(patient, visit, EasyM_value = value) %>%
  full_join(
    EasyM_pos_long %>% select(patient, visit, EasyM_mrd = mrd_status),
    by = c("patient", "visit")
  ) %>%
  arrange(patient, visit)

write_csv(EasyM_joined, file.path(out_dir, "EasyM_long_quant_and_binary.csv"))

## ============================================================
## 3) EasyM ONLY: plots restricted to EasyM measurements
## ============================================================
# Generate descriptive plots showing:
# - Longitudinal trajectories of M-protein levels per patient
# - Proportion of MRD+ patients at each visit

## 3A) Longitudinal quantitative trajectories
p_EasyM_traj <- ggplot(EasyM_values_long, aes(x = visit, y = value, group = patient)) +
  geom_line(linewidth = 0.7, colour = "black", alpha = 0.5) +
  geom_point(size = 2, colour = "black") +
  labs(
    title = "EasyM longitudinal trajectories",
    x = "Visit",
    y = "Residual immunoglobulin (% of baseline BM)"
  ) +
  theme_mrd

ggsave(
  filename = file.path(out_dir, "EasyM_trajectories_quantitative.png"),
  plot     = p_EasyM_traj,
  width    = 6,
  height   = 4,
  dpi      = 600
)

## 3B) MRD+ proportion over time (EasyM only)
EasyM_pos_summary <- EasyM_pos_long %>%
  group_by(visit) %>%
  summarise(
    n        = n(),
    n_pos    = sum(mrd_status == "MRD+"),
    prop_pos = n_pos / n,
    .groups  = "drop"
  )

p_EasyM_prop_pos <- ggplot(EasyM_pos_summary, aes(x = visit, y = prop_pos, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "EasyM MRD positivity over time",
    x = "Visit",
    y = "MRD+ patients (%)"
  ) +
  theme_mrd

ggsave(
  filename = file.path(out_dir, "EasyM_MRDpos_proportion_over_time.png"),
  plot     = p_EasyM_prop_pos,
  width    = 6,
  height   = 4,
  dpi      = 600
)


## ============================================================
## 4) Merge EasyM with cfWGS data
## ============================================================
# Load cfWGS MRD detection data (from prior pipeline steps)
# and merge with EasyM to enable head-to-head comparisons

### Now load in the cfWGS data 
dat_rds       <- "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6.rds"

## Add the other info that we have to this to get the actual total patient counts including baseline 
dat <- readRDS(dat_rds) %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )

## ============================================================
## Harmonize visit labels between datasets
## ============================================================
# EasyM uses "V1", "V3" etc. while cfWGS uses "01", "03" etc.
# Convert EasyM format to match cfWGS numbering

## Now get together in one table 
library(stringr)

EasyM_joined <- EasyM_joined %>%
  mutate(
    visit = as.character(visit),
    visit = case_when(
      visit == "R" ~ "R",
      str_detect(visit, "^V\\d+$") ~ {
        vnum <- as.integer(str_remove(visit, "^V"))
        if_else(vnum < 10, str_c("0", vnum), as.character(vnum))
      },
      TRUE ~ visit
    ),
    visit = factor(visit, levels = c("01","03","05","07","08","09","10","11","12","13","R"))
  )

colnames(dat)


dat_joined <- dat %>%
  mutate(
    Patient   = as.character(Patient),
    Timepoint = as.character(Timepoint)
  ) %>%
  left_join(
    EasyM_joined %>%
      mutate(
        patient = as.character(patient),
        visit   = as.character(visit)
      ) %>%
      select(
        patient, visit,
        EasyM_value, EasyM_mrd
      ),
    by = c("Patient" = "patient", "Timepoint" = "visit")
  )


## ============================================================
## Create analysis dataset: paired samples only
## ============================================================
# Restrict to samples where BOTH EasyM and cfWGS measurements exist
# Exclude patients with only baseline/diagnosis samples (need follow-up)

## Keep only samples with both assays, then drop patients who only have diagnosis (Timepoint == "01")
dat_joined <- dat_joined %>%
  filter(
    !is.na(EasyM_value),
    !is.na(BM_zscore_only_detection_rate_prob)
  ) %>%
  mutate(
    Patient   = as.character(Patient),
    Timepoint = as.character(Timepoint)
  ) %>%
  group_by(Patient) %>%
  filter(
    ## keep patients with at least one non-diagnosis timepoint
    any(Timepoint != "01", na.rm = TRUE)
  ) %>%
  ungroup()

## ============================================================
## 5) Cohort description: sample and patient counts
## ============================================================
# Summarize the longitudinal cohort for reporting in manuscript
# - Total patients with paired data
# - Sample counts at each timepoint
# - Number of clinical sites

## Discuss patients with actual longitudinal samples 
## Now get the patient, sample and timepoint counts
tp_summary <- dat_joined %>%
  mutate(
    Patient   = as.character(Patient),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info)
  ) %>%
  group_by(Timepoint, timepoint_info) %>%
  summarise(
    n_samples  = n(),
    n_patients = n_distinct(Patient),
    .groups = "drop"
  ) %>%
  arrange(Timepoint)

tp_summary

n_patients_total <- dat_joined %>%
  summarise(n = n_distinct(Patient)) %>%
  pull(n)




## 1) Count unique sites as the 2-letter prefix in Patient (e.g., "CA" from "CA-01")
n_sites <- dat_joined %>%
  transmute(site = str_extract(as.character(Patient), "^[A-Za-z]{2}")) %>%
  filter(!is.na(site)) %>%
  summarise(n_sites = n_distinct(site)) %>%
  pull(n_sites)

n_sites







## ============================================================
## 6) Load survival/relapse data and prepare for PFS analysis
## ============================================================
# Merge clinical outcomes data (PFS events and censor dates)
# with the paired MRD measurements to enable time-to-event analysis

### Now add in the info regarding relapse 
## ── 0.  SETUP ────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(tableone)      # optional – baseline table
library(timeROC)       # optional – AUC vs time

## INPUT (already saved from the pipeline) ------------------------------
final_tbl_rds <- "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds"

## OUTPUT ----------------------------------------------------------------------

## ── 1.  LOAD & TIDY CORE TABLES ──────────────────────────────────────────────
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

dat <- dat_joined %>%
  mutate(
    Patient        = as.character(Patient),
    sample_date    = as.Date(Date),
    timepoint_info = tolower(timepoint_info)
  )


## ============================================================
## 7) Build patient-level survival data frame
## ============================================================
# For each MRD measurement, compute time from that sample to PFS event
# This enables landmark survival analysis from any visit

## ── 2.  BUILD PFS TABLE (patient-level) ──────────────────────────────────────
survival_df <- dat %>%
  mutate(
    Patient     = as.character(Patient),
    sample_date = as.Date(Date)                     # the draw date for each row
  ) %>%
  # bring in censor_date + relapsed from final_tbl
  left_join(
    final_tbl %>% 
      select(Patient, censor_date, relapsed),
    by = "Patient"
  ) %>%
  # compute days from the sample to the event/censor
  mutate(
    Time_to_event   = as.numeric(censor_date - sample_date),
    Relapsed_Binary = as.integer(relapsed)
  ) %>%
  # keep only the columns the KM‐loop needs
  select(
    Patient, Timepoint, sample_date, censor_date, timepoint_info,
    Time_to_event, Relapsed_Binary,
    Flow_Binary, Adaptive_Binary, Rapid_Novor_Binary,
    Flow_pct_cells, Adaptive_Frequency,
    PET_Binary,
    BM_zscore_only_detection_rate_call, BM_zscore_only_detection_rate_prob,
    Blood_zscore_only_sites_call, Blood_zscore_only_sites_prob, Blood_base_prob, Blood_base_call,
    Blood_plus_fragment_prob, Blood_plus_fragment_call,
    Blood_plus_fragment_min_prob, Blood_plus_fragment_min_call, Fragmentomics_prop_short_only_call, Fragmentomics_prop_short_only_prob,
    Fragmentomics_min_prob, Fragmentomics_min_call, EasyM_value, EasyM_mrd
  )

# sanity checks
table(survival_df$Relapsed_Binary, useNA="ifany")
summary(survival_df$Time_to_event)
table(survival_df$timepoint_info)

## Make EasyM binary 
survival_df <- survival_df %>%
  mutate(
    EasyM_Binary = case_when(
      EasyM_mrd %in% c("MRD+", "Positive", "pos", "POS", 1, "1") ~ 1L,
      EasyM_mrd %in% c("MRD-", "Negative", "neg", "NEG", 0, "0") ~ 0L,
      TRUE ~ NA_integer_
    )
  )

## ============================================================
## 8) Generate Kaplan-Meier survival curves stratified by MRD status
## ============================================================
# For each timepoint and MRD technology, create KM plots showing
# PFS stratified by MRD positive vs negative status
#
# Technologies compared:
# - MFC (Flow cytometry)
# - clonoSEQ (NGS-based Ig sequencing)
# - EasyM (Proteomic M-protein)
# - cfWGS BM-informed (This study's primary method)

#### get function ready 
# 2) Friendly tech names
techs <- c(
  Flow_Binary        = "MFC",
  Adaptive_Binary    = "clonoSEQ",
  EasyM_Binary       = "EasyM",
  BM_zscore_only_detection_rate_call = "BM-informed cfWGS (cVAF Zscore)"
)


# 3) All timepoints to cover
tps <- unique(survival_df$timepoint_info)
#tps <- "diagnosis"
dpi_target <- 500

# 4) Minimum n per group to plot
min_n <- 2

#pal_2 <- viridis(2, option = "D", begin = 0.3, end = 0.7)   # nice mid‑range hues
# old colors: palette = c("#E7B800","#2E9FDF")
pal_2 <- c("black", "red") # updated 

tp_labels <- c(
  `diagnosis` = "Diagnosis",
  `post_transplant`    = "Post‑ASCT",
  `1yr maintenance` = "One-Year Maintenance", 
  `post_induction`    = "Post‑Induction"
)

# once before the loops
dx_tbl <- survival_df %>%
  group_by(Patient) %>%
  summarise(
    diagnosis_date = suppressWarnings(min(sample_date[timepoint_info == "diagnosis"], na.rm = TRUE)),
    .groups = "drop"
  )



### Describe some info on the samples 
## Use baseline/diagnosis rows only: one row per patient for follow-up + event
baseline_df <- survival_df %>%
  filter(Timepoint == "01") %>%              # diagnosis
  mutate(
    Patient = as.character(Patient),
    Time_to_event = as.numeric(Time_to_event)
  ) %>%
  group_by(Patient) %>%
  summarise(
    time_to_event_days = first(Time_to_event),
    progressed         = first(Relapsed_Binary),
    .groups = "drop"
  )

## Follow-up time in months (from diagnosis)
fu_mo <- baseline_df$time_to_event_days / 30.4375

## Use censored patients for follow-up summary (standard)
fu_cens_mo <- fu_mo[baseline_df$progressed == 0 & is.finite(fu_mo)]

median_followup_months <- median(fu_cens_mo, na.rm = TRUE)
min_followup_months    <- min(fu_cens_mo, na.rm = TRUE)
max_followup_months    <- max(fu_cens_mo, na.rm = TRUE)

## Number progressed
n_progressed <- sum(baseline_df$progressed == 1, na.rm = TRUE)
n_total      <- nrow(baseline_df)


## Median follow-up in months (from diagnosis)
median_followup_months <- median(baseline_df$time_to_event_days, na.rm = TRUE) / 30.4375

## Number progressed
n_progressed <- sum(baseline_df$progressed == 1, na.rm = TRUE)
n_total <- nrow(baseline_df)




# 5) Loop
for(tp in tps) {
  #  nice_tp <- tp_labels[tp] %||% tp   # fall back to tp if no mappiht
  # instead of  %||% line:
  nice_tp <- as.character(tp_labels[tp])
  if (is.na(nice_tp) || nice_tp == "") nice_tp <- as.character(tp)
  
  
  tp_dir <- file.path(out_dir, gsub("\\s+","_", tp))  # sanitize folder name
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(var in names(techs)) {
    assay_lab <- techs[[var]]
    fname     <- file.path(tp_dir, paste0("KM_", assay_lab, "_", nice_tp, "_updated_no_CI.png"))
    
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
      pval            = TRUE,
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



## ============================================================
## 9) Landmark analysis: detailed statistical comparisons
## ============================================================
# At key timepoints (e.g., post-transplant, 1-year maintenance),
# perform comprehensive head-to-head comparison of EasyM vs cfWGS:
#
# Analyses performed:
# 1. Binary agreement (kappa, percent concordance)
# 2. Continuous correlation (Spearman rho)
# 3. Cox proportional hazards models for PFS prediction
# 4. C-index comparison (discrimination ability)
# 5. Integrated models (cfWGS + EasyM together)
# 6. Optimal threshold determination via log-rank maximization

## ============================================================
## 1) Helpers
## ============================================================

days_to_months <- function(x_days) x_days / 30.4375

## Robust Cohen's kappa for 2x2 tables (manual)
kappa_2x2 <- function(tab) {
  # tab is a 2x2 matrix/table with same ordering for both variables
  if (!all(dim(tab) == c(2,2))) return(NA_real_)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (abs(1 - pe) < 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}

## Compute concordance index for a Cox model
c_index_cox <- function(fit, data, time_col = "time_months", event_col = "event") {
  lp <- predict(fit, newdata = data, type = "lp")
  cc <- survival::concordance(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ lp")),
    data = data
  )
  unname(cc$concordance)
}

## Pretty HR extraction for single coefficient models
extract_hr <- function(fit) {
  s <- summary(fit)
  out <- data.frame(
    term = rownames(s$coefficients),
    HR   = exp(s$coefficients[, "coef"]),
    lo   = exp(s$conf.int[, "lower .95"]),
    hi   = exp(s$conf.int[, "upper .95"]),
    p    = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
  out
}

## ============================================================
## 2) Clean survival_df into analysis-ready format
## ============================================================
df <- survival_df %>%
  mutate(
    Patient   = as.character(Patient),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info),
    
    ## survival
    time_months = days_to_months(as.numeric(Time_to_event)),
    event       = as.integer(Relapsed_Binary),
    
    ## categorical predictors
    EasyM_Binary = as.integer(EasyM_Binary),
    BM_call      = as.integer(BM_zscore_only_detection_rate_call),
    
    ## continuous predictors
    EasyM_value  = as.numeric(EasyM_value),
    BM_prob      = as.numeric(BM_zscore_only_detection_rate_prob),
    
    ## useful transforms (avoid log(0); EasyM sometimes hits 0.xx etc)
    EasyM_log10  = log10(pmax(EasyM_value, 1e-6)),
    BM_prob_clamped = pmin(pmax(BM_prob, 1e-12), 1 - 1e-12) # if you ever logit it later
  )

## ============================================================
## 3) Baseline cohort counts (diagnosis timepoint)
## ============================================================
baseline_df <- df %>%
  filter(Timepoint == "01") %>%
  group_by(Patient) %>%
  summarise(
    time_months = first(time_months),
    event       = first(event),
    .groups = "drop"
  )

median_fu_months <- median(baseline_df$time_months, na.rm = TRUE)
n_prog <- sum(baseline_df$event == 1, na.rm = TRUE)
n_tot  <- nrow(baseline_df)


## ============================================================
## 4) Landmark analysis function
## ============================================================
run_landmark <- function(tp,
                         require_paired = TRUE,
                         easyM_cut = c("median", "tertile_top", "threshold"),
                         easyM_thresholds = NULL) {
  
  easyM_cut <- match.arg(easyM_cut)
  
  lm_df <- df %>%
    filter(Timepoint == tp) %>%
    # keep only patients with defined survival from this landmark
    filter(!is.na(time_months), !is.na(event)) %>%
    # optionally require both assays present (continuous + binary as needed)
    { if (require_paired) filter(., !is.na(EasyM_value), !is.na(BM_prob), !is.na(BM_call)) else . } %>%
    distinct(Patient, .keep_all = TRUE)
  
  cat("============================================================\n")
  cat(sprintf("Landmark Timepoint %s (%s)\n", tp,
              paste(unique(na.omit(lm_df$timepoint_info)), collapse = ", ")))
  cat(sprintf("N (paired) = %d; events = %d\n",
              nrow(lm_df), sum(lm_df$event == 1, na.rm = TRUE)))
  
  if (nrow(lm_df) < 5) {
    cat("Too few samples at this landmark to run models reliably.\n\n")
    return(invisible(NULL))
  }
  
  ## ---- 4A) Binary agreement between EasyM_Binary and BM_call ----
  cat("\n[Binary agreement: EasyM_Binary vs BM_call]\n")
  if (dplyr::n_distinct(lm_df$EasyM_Binary[!is.na(lm_df$EasyM_Binary)]) < 2) {
    cat("EasyM_Binary is (near-)constant at this landmark; skipping agreement/kappa.\n")
  } else {
    tab <- table(BM_call = lm_df$BM_call, EasyM_Binary = lm_df$EasyM_Binary, useNA = "no")
    print(tab)
    po <- sum(diag(tab)) / sum(tab)
    kap <- kappa_2x2(tab)
    cat(sprintf("Percent agreement = %.1f%%; Cohen's kappa = %s\n",
                100 * po,
                ifelse(is.na(kap), "NA (undefined)", sprintf("%.3f", kap))))
  }
  
  ## ---- 4B) Continuous correlation (BM_prob vs EasyM_value) ----
  cat("\n[Continuous correlation: BM_prob vs EasyM_value]\n")
  cor_dat <- lm_df %>% filter(!is.na(BM_prob), !is.na(EasyM_log10))
  if (nrow(cor_dat) >= 5) {
    ct <- cor.test(cor_dat$BM_prob, cor_dat$EasyM_log10, method = "spearman", exact = FALSE)
    cat(sprintf("Spearman rho = %.3f; p = %.3g (BM_prob vs log10(EasyM_value))\n",
                unname(ct$estimate), ct$p.value))
  } else {
    cat("Too few complete pairs for correlation.\n")
  }
  
  ## ---- 4C) Cox models (binary and continuous) ----
  cat("\n[Cox models for PFS from landmark]\n")
  
  # 1) cfWGS (BM_call) binary
  fit_cf_bin <- coxph(Surv(time_months, event) ~ BM_call, data = lm_df)
  print(extract_hr(fit_cf_bin))
  c_cf <- c_index_cox(fit_cf_bin, lm_df)
  cat(sprintf("C-index (cfWGS binary) = %.3f\n", c_cf))
  
  # 2) EasyM continuous (log10 + z-score scaling)
  lm_df <- lm_df %>% mutate(EasyM_z = as.numeric(scale(EasyM_log10)),
                            BMprob_z = as.numeric(scale(BM_prob)))
  
  fit_easy_cont <- coxph(Surv(time_months, event) ~ EasyM_z, data = lm_df)
  print(extract_hr(fit_easy_cont))
  c_easy <- c_index_cox(fit_easy_cont, lm_df)
  cat(sprintf("C-index (EasyM continuous) = %.3f\n", c_easy))
  
  # 3) cfWGS continuous
  fit_cf_cont <- coxph(Surv(time_months, event) ~ BMprob_z, data = lm_df)
  print(extract_hr(fit_cf_cont))
  c_cfcont <- c_index_cox(fit_cf_cont, lm_df)
  cat(sprintf("C-index (cfWGS continuous) = %.3f\n", c_cfcont))
  
  # 4) Integrated (cfWGS binary + EasyM continuous)
  fit_int <- coxph(Surv(time_months, event) ~ BM_call + EasyM_z, data = lm_df)
  print(extract_hr(fit_int))
  c_int <- c_index_cox(fit_int, lm_df)
  cat(sprintf("C-index (Integrated: BM_call + EasyM_z) = %.3f\n", c_int))
  
  # Likelihood ratio test: does adding EasyM improve over cfWGS binary alone?
  lrt <- anova(fit_cf_bin, fit_int, test = "LRT")
  print(lrt)
  
  ## ---- 4D) Clinically interpretable joint risk groups ----
  cat("\n[Joint risk groups: cfWGS binary x EasyM burden (high/low)]\n")
  
  if (easyM_cut == "median") {
    thr <- median(lm_df$EasyM_log10, na.rm = TRUE)
    lm_df <- lm_df %>% mutate(EasyM_high = if_else(EasyM_log10 >= thr, 1L, 0L))
    cut_label <- "above-median"
    
  } else if (easyM_cut == "tertile_top") {
    thr <- quantile(lm_df$EasyM_log10, probs = 2/3, na.rm = TRUE)
    lm_df <- lm_df %>% mutate(EasyM_high = if_else(EasyM_log10 >= thr, 1L, 0L))
    cut_label <- "top tertile"
    
  } else if (easyM_cut == "threshold") {
    if (is.null(easyM_thresholds) || is.na(easyM_thresholds[tp])) {
      stop(sprintf("No EasyM threshold provided for timepoint %s", tp), call. = FALSE)
    }
    thr_raw <- unname(easyM_thresholds[tp])
    lm_df <- lm_df %>% mutate(EasyM_high = if_else(EasyM_value >= thr_raw, 1L, 0L))
    cut_label <- paste0(">= ", signif(thr_raw, 3), "%")
  }
  
  lm_df <- lm_df %>%
    mutate(
      cf_pos = if_else(BM_call == 1, "cfWGS+", "cfWGS-"),
      em_pos = if_else(EasyM_high == 1, paste0("EasyM_high (", cut_label, ")"), paste0("EasyM_low (", cut_label, ")")),
      joint_group = factor(paste(cf_pos, em_pos, sep = " / "))
    )
  
  print(table(lm_df$joint_group, useNA = "ifany"))
  
  fit_grp <- coxph(Surv(time_months, event) ~ joint_group, data = lm_df)
  print(extract_hr(fit_grp))
  
  invisible(list(
    lm_df = lm_df,
    fits = list(cf_bin = fit_cf_bin, easy_cont = fit_easy_cont, cf_cont = fit_cf_cont, integrated = fit_int, group = fit_grp),
    cindex = c(cf_bin = c_cf, easy_cont = c_easy, cf_cont = c_cfcont, integrated = c_int)
  ))
}

## ============================================================
## 5) Run landmark analyses (choose timepoints needed for this analysis)
## ============================================================
# Example: post-induction, post-transplant, 1yr maintenance, relapse
timepoints_to_run <- c("05", "07")

threshold_seed_path <- file.path(out_dir, "EasyM_threshold_values_by_timepoint.csv")
if (file.exists(threshold_seed_path)) {
  threshold_seed <- readr::read_csv(threshold_seed_path, show_col_types = FALSE)
  easyM_thresholds_calc <- stats::setNames(
    threshold_seed$Threshold_raw_percent,
    threshold_seed$Timepoint
  )
} else {
  easyM_thresholds_calc <- numeric(0)
  message(
    "No existing EasyM threshold table found for early exploratory landmark block; ",
    "using median EasyM cutoffs in that block. The primary optimized threshold ",
    "section below still computes and exports the manuscript thresholds."
  )
}

#results_list <- lapply(timepoints_to_run, run_landmark, require_paired = TRUE, easyM_cut = "median")
if (length(easyM_thresholds_calc) > 0) {
  results_list <- lapply(
    timepoints_to_run,
    run_landmark,
    require_paired = TRUE,
    easyM_cut = "threshold",
    easyM_thresholds = easyM_thresholds_calc
  )
} else {
  results_list <- lapply(timepoints_to_run, run_landmark, require_paired = TRUE, easyM_cut = "median")
}

names(results_list) <- timepoints_to_run





### Detailed landmark analyses and supporting tables

library(glue)
library(stringr)

## ----------------------------
## 1) Helpers
## ----------------------------
days_to_months <- function(x_days) x_days / 30.4375

fmt_num <- function(x, digits = 2) {
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}


extract_hr <- function(fit) {
  s  <- summary(fit)
  ci <- s$conf.int
  data.frame(
    term = rownames(ci),
    HR   = unname(ci[, "exp(coef)"]),
    lo   = unname(ci[, "lower .95"]),
    hi   = unname(ci[, "upper .95"]),
    p    = unname(s$coefficients[, "Pr(>|z|)"]),
    row.names = NULL
  )
}

c_index_cox <- function(fit, data, time_col = "time_months", event_col = "event") {
  lp <- predict(fit, newdata = data, type = "lp")
  cc <- survival::concordance(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ lp")),
    data = data
  )
  unname(cc$concordance)
}

kappa_2x2 <- function(tab) {
  # tab must be 2x2, no NA
  if (!all(dim(tab) == c(2,2))) return(NA_real_)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (abs(1 - pe) < 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}

km_surv_at <- function(time, event, group, t_months) {
  # returns named vector of survival at t_months per group (if estimable)
  sf <- survfit(Surv(time, event) ~ group)
  s  <- summary(sf, times = t_months)
  out <- tapply(s$surv, s$strata, function(v) v[1])
  out
}

## ============================================================
## KEY FUNCTION: find_opt_cut_logrank()
## ============================================================
# This function finds the optimal EasyM cutoff threshold at a given landmark
# by testing all possible cutpoints and selecting the one that maximizes
# separation between "cleared" vs "residual" groups (via log-rank test).
#
# INPUTS:
#   - time: survival time in months
#   - event: 1 = progressed, 0 = censored
#   - marker: continuous EasyM values (log10-transformed)
#   - min_group_n: minimum patients required in each group
#
# OUTPUT:
#   A list containing:
#   - cut: optimal cutoff value (log10 scale)
#   - p: log-rank p-value at this cutoff
#   - chisq: chi-square statistic
#   - n_low: patients below cutoff ("cleared")
#   - n_high: patients above cutoff ("residual")
#
# Find an "optimal" cutoff for a continuous marker at a landmark by scanning cutpoints.
# Criterion: maximize log-rank chi-square (equivalently minimize log-rank p).
# To avoid absurd splits, enforce a minimum group size per side.
find_opt_cut_logrank <- function(time, event, marker, min_group_n = 3) {
  stopifnot(length(time) == length(event), length(event) == length(marker))
  ok <- is.finite(time) & !is.na(event) & is.finite(marker)
  time   <- time[ok]
  event  <- event[ok]
  marker <- marker[ok]
  
  # Candidate cutpoints: unique marker values (excluding extremes to ensure both sides exist)
  u <- sort(unique(marker))
  if (length(u) < 3) return(NULL)
  
  # Use midpoints between consecutive unique values as candidate thresholds
  cuts <- (u[-1] + u[-length(u)]) / 2
  
  best <- list(chisq = -Inf, p = NA_real_, cut = NA_real_,
               n_low = NA_integer_, n_high = NA_integer_)
  for (cut in cuts) {
    grp <- ifelse(marker <= cut, "low_or_cleared", "high_or_residual")
    n_low  <- sum(grp == "low_or_cleared")
    n_high <- sum(grp == "high_or_residual")
    if (min(n_low, n_high) < min_group_n) next
    
    sd <- try(survdiff(Surv(time, event) ~ grp), silent = TRUE)
    if (inherits(sd, "try-error")) next
    
    chisq <- unname(sd$chisq)
    p <- pchisq(chisq, df = 1, lower.tail = FALSE)
    
    if (is.finite(chisq) && chisq > best$chisq) {
      best <- list(chisq = chisq, p = p, cut = cut, n_low = n_low, n_high = n_high)
    }
  }
  
  if (!is.finite(best$chisq)) return(NULL)
  best
}

## ----------------------------
## 2) Prep analysis df
## ----------------------------
df <- survival_df %>%
  mutate(
    Patient   = as.character(Patient),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info),
    
    time_months = days_to_months(as.numeric(Time_to_event)),
    event       = as.integer(Relapsed_Binary),
    
    ## Categorical
    EasyM_Binary = as.integer(EasyM_Binary),
    BM_call      = as.integer(BM_zscore_only_detection_rate_call),
    
    ## Continuous
    EasyM_value  = as.numeric(EasyM_value),
    BM_prob      = as.numeric(BM_zscore_only_detection_rate_prob),
    
    ## Common transforms
    EasyM_log10  = log10(pmax(EasyM_value, 1e-6))
  )

## ----------------------------
## 3) Overall follow-up summary (diagnosis rows)
## ----------------------------
baseline_df <- df %>%
  filter(Timepoint == "01") %>%
  group_by(Patient) %>%
  summarise(
    time_months = first(time_months),
    event       = first(event),
    .groups = "drop"
  )

median_fu <- median(baseline_df$time_months, na.rm = TRUE)
n_prog    <- sum(baseline_df$event == 1, na.rm = TRUE)
n_total   <- nrow(baseline_df)



## Follow-up time in months (from diagnosis)
fu_mo <- baseline_df$time_to_event_days / 30.4375

## Use censored patients for follow-up summary (standard)
fu_cens_mo <- fu_mo[baseline_df$event == 0 & is.finite(fu_mo)]

median_followup_months <- median(fu_cens_mo, na.rm = TRUE)
min_followup_months    <- min(fu_cens_mo, na.rm = TRUE)
max_followup_months    <- max(fu_cens_mo, na.rm = TRUE)

## Number progressed
n_progressed <- sum(baseline_df$event == 1, na.rm = TRUE)
n_total      <- nrow(baseline_df)


## ============================================================
## CRITICAL FUNCTION: build_landmark_results()
## ============================================================
# This is the main analysis engine that performs comprehensive landmark
# analysis at a specific timepoint, comparing EasyM vs cfWGS for PFS prediction.
#
# INPUT PARAMETERS:
#   tp: Timepoint code (e.g., "05" = post-transplant, "07" = 1-yr maintenance)
#   horizon_months: Time horizon for survival estimates (default 24 months)
#   easyM_cut: Method for dichotomizing EasyM ("median", "top_tertile", or "opt")
#   easyM_fixed_raw: Pre-specified threshold (if not using optimization)
#   opt_cut_min_group_n: Minimum patients per group for optimization
#
# KEY OUTPUTS (stored in returned list):
#   1. landmark_df: Patient-level data with MRD calls and survival
#   2. outputs$agreement: Agreement tables (binary concordance, kappa)
#   3. outputs$joint: Joint risk groups (cfWGS+/- x EasyM high/low)
#   4. outputs$opt_cut: Optimized threshold details
#   5. outputs$models: Cox model results and c-indices
#
## ----------------------------
## 4) Landmark analysis and numerical summaries
## ----------------------------
build_landmark_results <- function(tp = "07",
                                   horizon_months = 24,
                                   easyM_cut = c("median","top_tertile","opt"),
                                   easyM_fixed_raw = NULL,
                                   opt_cut_min_group_n = 3) {
  
  easyM_cut <- match.arg(easyM_cut)
  
  lm_df <- df %>%
    filter(Timepoint == tp) %>%
    filter(!is.na(time_months), !is.na(event)) %>%
    # Require paired data for comparisons
    filter(!is.na(EasyM_value), !is.na(BM_prob), !is.na(BM_call)) %>%
    distinct(Patient, .keep_all = TRUE) %>%
    mutate(
      EasyM_z  = as.numeric(scale(EasyM_log10)),
      BMprob_z = as.numeric(scale(BM_prob))
    )
  
  n_paired <- nrow(lm_df)
  n_events <- sum(lm_df$event == 1, na.rm = TRUE)
  tp_info  <- paste(unique(na.omit(lm_df$timepoint_info)), collapse = ", ")
  if (tp_info == "") tp_info <- "landmark"
  
  ## NEW: counts of positives/negatives at the landmark
  # cfWGS binary
  n_cf_pos <- sum(lm_df$BM_call == 1, na.rm = TRUE)
  n_cf_neg <- sum(lm_df$BM_call == 0, na.rm = TRUE)
  
  # EasyM binary from EasyM_Binary if present; otherwise from EasyM_mrd
  if ("EasyM_Binary" %in% colnames(lm_df) && any(!is.na(lm_df$EasyM_Binary))) {
    n_em_pos <- sum(lm_df$EasyM_Binary == 1, na.rm = TRUE)
    n_em_neg <- sum(lm_df$EasyM_Binary == 0, na.rm = TRUE)
  } else {
    n_em_pos <- sum(lm_df$EasyM_mrd == "MRD+", na.rm = TRUE)
    n_em_neg <- sum(lm_df$EasyM_mrd == "MRD-", na.rm = TRUE)
  }
  
  
  ## Binary agreement (only if both vary)
  easyM_bin_varies <- n_distinct(lm_df$EasyM_Binary[!is.na(lm_df$EasyM_Binary)]) >= 2
  bm_call_varies   <- n_distinct(lm_df$BM_call[!is.na(lm_df$BM_call)]) >= 2
  
  if (easyM_bin_varies && bm_call_varies) {
    tab <- table(BM_call = lm_df$BM_call, EasyM_Binary = lm_df$EasyM_Binary, useNA = "no")
    po  <- sum(diag(tab)) / sum(tab)
    kap <- kappa_2x2(tab)
  }
  
  ## Continuous correlation
  cor_dat <- lm_df %>% filter(!is.na(BM_prob), !is.na(EasyM_log10))
  rho <- NA_real_; cor_p <- NA_real_
  if (nrow(cor_dat) >= 5) {
    ct <- cor.test(cor_dat$BM_prob, cor_dat$EasyM_log10, method = "spearman", exact = FALSE)
    rho   <- unname(ct$estimate)
    cor_p <- ct$p.value
  }
  
  ## Cox models
  # cfWGS binary
  fit_cf_bin <- coxph(Surv(time_months, event) ~ BM_call, data = lm_df)
  hr_cf_bin  <- extract_hr(fit_cf_bin)[1, ]
  c_cf_bin   <- c_index_cox(fit_cf_bin, lm_df)
  
  # cfWGS continuous
  fit_cf_cont <- coxph(Surv(time_months, event) ~ BMprob_z, data = lm_df)
  hr_cf_cont  <- extract_hr(fit_cf_cont)[1, ]
  c_cf_cont   <- c_index_cox(fit_cf_cont, lm_df)
  
  # EasyM binary (NEW)
  c_em_bin <- NA_real_
  hr_em_bin <- NULL
  if (easyM_bin_varies) {
    fit_em_bin <- coxph(Surv(time_months, event) ~ EasyM_Binary, data = lm_df)
    hr_em_bin  <- extract_hr(fit_em_bin)[1, ]
    c_em_bin   <- c_index_cox(fit_em_bin, lm_df)
  }
  
  # EasyM continuous
  fit_em_cont <- coxph(Surv(time_months, event) ~ EasyM_z, data = lm_df)
  hr_em_cont  <- extract_hr(fit_em_cont)[1, ]
  c_em_cont   <- c_index_cox(fit_em_cont, lm_df)
  
  ## Integrated model: cfWGS binary + EasyM continuous
  fit_int <- coxph(Surv(time_months, event) ~ BM_call + EasyM_z, data = lm_df)
  c_int   <- c_index_cox(fit_int, lm_df)
  
  # LRT: cfWGS bin vs integrated
  lrt <- anova(fit_cf_bin, fit_int, test = "LRT")
  lrt_p <- lrt$`Pr(>|Chi|)`[2]
  
  
  ## NEW: “optimal clearance threshold” for EasyM (exploratory)
  opt <- find_opt_cut_logrank(
    time   = lm_df$time_months,
    event  = lm_df$event,
    marker = lm_df$EasyM_log10,
    min_group_n = opt_cut_min_group_n
  )
  
  if (!is.null(opt)) {
    # Convert cut back to raw EasyM_value scale
    opt_cut_log10 <- opt$cut
    opt_cut_raw   <- 10^opt_cut_log10
    
    # Build a binary split based on "cleared" (<= cut) vs "residual" (> cut)
    lm_df <- lm_df %>%
      mutate(
        EasyM_clear_opt = if_else(EasyM_log10 <= opt_cut_log10, 1L, 0L),
        EasyM_opt_group = if_else(EasyM_clear_opt == 1L, "EasyM_cleared_opt", "EasyM_residual_opt")
      )
    
    # Cox for this split (optional, but useful)
    fit_opt <- coxph(Surv(time_months, event) ~ EasyM_opt_group, data = lm_df)
    hr_opt  <- extract_hr(fit_opt)[1, ]
    
    ## --- Agreement: cfWGS binary vs optimized EasyM clearance threshold ---
    # Define EasyM "residual" call using the optimized cutoff:
    # cleared (<= cutoff) -> 0, residual (> cutoff) -> 1
    lm_df <- lm_df %>%
      mutate(
        EasyM_opt_binary = if_else(EasyM_log10 > opt_cut_log10, 1L, 0L)  # 1=residual, 0=cleared
      )
    
    
    # Only compute if both variables vary
    if (n_distinct(lm_df$BM_call[!is.na(lm_df$BM_call)]) >= 2 &&
        n_distinct(lm_df$EasyM_opt_binary[!is.na(lm_df$EasyM_opt_binary)]) >= 2) {
      
      tab_opt <- table(
        BM_call = lm_df$BM_call,
        EasyM_opt_binary = lm_df$EasyM_opt_binary,
        useNA = "no"
      )
      
      # tab_opt already computed above
      n_total_opt <- sum(tab_opt)
      n_agree_opt <- sum(diag(tab_opt))
      po_opt      <- n_agree_opt / n_total_opt
      kap_opt     <- kappa_2x2(tab_opt)
      
      discord_df <- lm_df %>%
        filter(!is.na(BM_call), !is.na(EasyM_opt_binary)) %>%
        mutate(
          discord_type = case_when(
            BM_call == 1L & EasyM_opt_binary == 0L ~ "cfWGS+ / EasyM_cleared",
            BM_call == 0L & EasyM_opt_binary == 1L ~ "cfWGS- / EasyM_residual",
            TRUE ~ NA_character_
          )
        ) %>%
        filter(!is.na(discord_type)) %>%
        transmute(
          Patient,
          discord_type,
          cfWGS_call = BM_call,
          EasyM_opt_residual = EasyM_opt_binary,
          relapsed = event,
          time_to_event_months = time_months
        )
      
      discord_sum <- discord_df %>%
        group_by(discord_type) %>%
        summarise(
          n = n(),
          patients = paste(Patient, collapse = ", "),
          n_relapsed = sum(relapsed == 1L, na.rm = TRUE),
          med_tte_mo = if (sum(relapsed == 1L, na.rm = TRUE) > 0)
            median(time_to_event_months[relapsed == 1L], na.rm = TRUE) else NA_real_,
          min_tte_mo = if (sum(relapsed == 1L, na.rm = TRUE) > 0)
            min(time_to_event_months[relapsed == 1L], na.rm = TRUE) else NA_real_,
          max_tte_mo = if (sum(relapsed == 1L, na.rm = TRUE) > 0)
            max(time_to_event_months[relapsed == 1L], na.rm = TRUE) else NA_real_,
          .groups = "drop"
        )
      
      
      
    }
    
  }
  
  ## Joint risk groups (cfWGS +/- x EasyM high/low by chosen rule)
  if (easyM_cut == "median") {
    thr <- median(lm_df$EasyM_log10, na.rm = TRUE)
    cut_label <- "median"
    
  } else if (easyM_cut == "top_tertile") {
    thr <- quantile(lm_df$EasyM_log10, probs = 2/3, na.rm = TRUE)
    cut_label <- "top tertile"
    
  } else if (easyM_cut == "opt") {
    if (is.null(opt)) stop("opt is NULL; cannot use easyM_cut='opt'")
    thr <- opt$cut
    cut_label <- paste0("opt=", fmt_num(10^thr, 2))
  }
  
  lm_df <- lm_df %>%
    mutate(
      EasyM_high = if_else(EasyM_log10 >= thr, 1L, 0L),
      cf_group   = if_else(BM_call == 1, "cfWGS+", "cfWGS-"),
      em_group   = if_else(EasyM_high == 1, paste0("EasyM_high (", cut_label, ")"), paste0("EasyM_low (", cut_label, ")")),
      joint_group = paste(cf_group, em_group, sep = " / ")
    )
  
  joint_counts <- lm_df %>%
    count(joint_group, name = "n") %>%
    arrange(desc(n))
  
  
  ## --- Joint risk groups summary with relapse timing ---
  joint_summary <- lm_df %>%
    group_by(joint_group) %>%
    summarise(
      n = n(),
      n_relapsed = sum(event == 1L, na.rm = TRUE),
      med_tte_mo = if (sum(event == 1L, na.rm = TRUE) > 0) median(time_months[event == 1L], na.rm = TRUE) else NA_real_,
      min_tte_mo = if (sum(event == 1L, na.rm = TRUE) > 0) min(time_months[event == 1L], na.rm = TRUE) else NA_real_,
      max_tte_mo = if (sum(event == 1L, na.rm = TRUE) > 0) max(time_months[event == 1L], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      joint_group = factor(
        joint_group,
        levels = c(
          paste0("cfWGS- / EasyM_low (", cut_label, ")"),
          paste0("cfWGS- / EasyM_high (", cut_label, ")"),
          paste0("cfWGS+ / EasyM_low (", cut_label, ")"),
          paste0("cfWGS+ / EasyM_high (", cut_label, ")")
        )
      )
    ) %>%
    arrange(joint_group)
  
  
  
  
  
  
  
  
  list(
    landmark_df = lm_df,
    summary = list(
      tp = tp,
      tp_info = tp_info,
      n = n_paired,
      events = n_events,
      n_cf_pos = n_cf_pos, n_cf_neg = n_cf_neg,
      n_em_pos = n_em_pos, n_em_neg = n_em_neg,
      hr_cf_bin = hr_cf_bin, hr_em_bin = hr_em_bin,
      hr_cf_cont = hr_cf_cont, hr_em_cont = hr_em_cont,
      c_cf_bin = c_cf_bin, c_cf_cont = c_cf_cont,
      c_em_bin = c_em_bin, c_em_cont = c_em_cont,
      c_int = c_int, lrt_p = lrt_p,
      rho = rho, cor_p = cor_p,
      opt_cut = opt
    )
  )
  
  ### Get full exports 
  
  # --- Return the landmark results and supporting tables from build_landmark_results():
  
  outputs <- list(
    # Core objects
    landmark_df   = lm_df,
    
    
    # Agreement / discordance tables
    agreement = list(
      # “tab” only exists if you computed binary agreement for EasyM_Binary vs BM_call
      tab_binary = if (exists("tab")) tab else NULL,
      
      # optimized threshold agreement table
      tab_opt    = if (exists("tab_opt")) tab_opt else NULL,
      
      # discordant cases and their summary
      discordant_df  = if (exists("discord_df")) discord_df else NULL,
      discordant_sum = if (exists("discord_sum")) discord_sum else NULL
    ),
    
    # Joint stratification summary (counts + relapse timing)
    joint = list(
      joint_counts  = if (exists("joint_counts")) joint_counts else NULL,
      joint_summary = if (exists("joint_summary")) joint_summary else NULL
    ),
    
    # Correlation stats
    correlation = list(
      n     = nrow(cor_dat),
      rho   = rho,
      p     = cor_p
    ),
    
    # Model summaries (tables you printed)
    models = list(
      hr_cf_bin  = hr_cf_bin,
      hr_cf_cont = hr_cf_cont,
      hr_em_bin  = hr_em_bin,
      hr_em_cont = hr_em_cont,
      cindex = list(
        c_cf_bin  = c_cf_bin,
        c_cf_cont = c_cf_cont,
        c_em_bin  = c_em_bin,
        c_em_cont = c_em_cont,
        c_int     = c_int
      ),
      lrt = list(
        p = lrt_p
      )
    ),
    
    # Optional: store full model objects for later re-use (forest plots, diagnostics, etc.)
    fits = list(
      fit_cf_bin   = fit_cf_bin,
      fit_cf_cont  = fit_cf_cont,
      fit_em_cont  = fit_em_cont,
      fit_int      = fit_int,
      fit_em_bin   = if (exists("fit_em_bin")) fit_em_bin else NULL,
      fit_opt      = if (exists("fit_opt")) fit_opt else NULL
    ),
    
    # Optimal cutoff details (retain the output of find_opt_cut_logrank)
    opt_cut = list(
      opt = if (exists("opt")) opt else NULL,
      # convenience copies (if you created them)
      cut_log10 = if (exists("opt_cut_log10")) opt_cut_log10 else NULL,
      cut_raw   = if (exists("opt_cut_raw")) opt_cut_raw else NULL
    )
  )
  
  return(list(
    landmark_df = lm_df,
    summary     = summary,   # keep what already built
    outputs     = outputs
  ))
  
  
}

## ============================================================
## EXPORT FUNCTION: export_landmark_outputs()
## ============================================================
# Saves all results from a landmark analysis to CSV files for review
#
# Landmark-analysis summary table:
# The file named "{prefix}_landmark_df.csv" contains:
#   - Patient ID
#   - Timepoint
#   - EasyM_opt_binary: Whether patient is positive (1) or negative (0) 
#                       using the OPTIMIZED cutoff for that timepoint
#   - EasyM_value: Raw M-protein percentage
#   - BM_call: cfWGS binary MRD call
#   - time_months: Time to progression from this landmark
#   - event: Whether patient progressed (1) or was censored (0)
#   - joint_group: Combined risk stratification (e.g., "cfWGS- / EasyM_low")
#
# Other exported tables:
#   - *_agreement_opt_table.csv: 2x2 table of cfWGS vs EasyM agreement
#   - *_discordant_cases.csv: List of patients where assays disagree
#   - *_joint_summary.csv: Relapse rates by combined risk group
#   - *_opt_cut_summary.csv: The calculated optimal threshold value
#   - *_metrics.csv: C-indices, correlation, p-values

export_landmark_outputs <- function(res, out_dir = "Output_EasyM_MRD_analysis_2025", prefix = NULL) {
  stopifnot(is.list(res), !is.null(res$outputs), !is.null(res$outputs$landmark_df))
  
  if (is.null(prefix)) {
    tp <- res$summary$tp %||% "TP"
    prefix <- paste0("landmark_", tp)
  }
  
  # Create subdirectory for this timepoint's results
  tp_subdir <- file.path(out_dir, "landmark_analyses")
  dir.create(tp_subdir, showWarnings = FALSE, recursive = TRUE)
  
  # helper: safe write_csv for NULL objects
  safe_write_csv <- function(x, path) {
    if (is.null(x)) return(invisible(FALSE))
    # tables can be matrix/table; coerce cleanly
    if (inherits(x, "table") || is.matrix(x)) {
      x <- as.data.frame.matrix(x)
      x <- tibble::rownames_to_column(x, var = "row")
    }
    readr::write_csv(as.data.frame(x), path)
    invisible(TRUE)
  }
  
  
  # 1) Landmark analysis dataframe
  lm_path <- file.path(out_dir, paste0(prefix, "_landmark_df.csv"))
  readr::write_csv(res$outputs$landmark_df, lm_path)
  
  # 2) Key tables
  safe_write_csv(res$outputs$agreement$tab_binary,
                 file.path(out_dir, paste0(prefix, "_agreement_binary_table.csv")))
  safe_write_csv(res$outputs$agreement$tab_opt,
                 file.path(out_dir, paste0(prefix, "_agreement_opt_table.csv")))
  safe_write_csv(res$outputs$agreement$discordant_df,
                 file.path(out_dir, paste0(prefix, "_discordant_cases.csv")))
  safe_write_csv(res$outputs$agreement$discordant_sum,
                 file.path(out_dir, paste0(prefix, "_discordant_summary.csv")))
  
  safe_write_csv(res$outputs$joint$joint_counts,
                 file.path(out_dir, paste0(prefix, "_joint_counts.csv")))
  safe_write_csv(res$outputs$joint$joint_summary,
                 file.path(out_dir, paste0(prefix, "_joint_summary.csv")))
  
  # 3) Model summary tables (HRs + c-index + correlation)
  safe_write_csv(res$outputs$models$hr_cf_bin,
                 file.path(out_dir, paste0(prefix, "_hr_cfWGS_binary.csv")))
  safe_write_csv(res$outputs$models$hr_cf_cont,
                 file.path(out_dir, paste0(prefix, "_hr_cfWGS_continuous.csv")))
  safe_write_csv(res$outputs$models$hr_em_bin,
                 file.path(out_dir, paste0(prefix, "_hr_EasyM_binary.csv")))
  safe_write_csv(res$outputs$models$hr_em_cont,
                 file.path(out_dir, paste0(prefix, "_hr_EasyM_continuous.csv")))
  
  readr::write_csv(
    tibble::tibble(
      metric = c("c_cf_bin","c_cf_cont","c_em_bin","c_em_cont","c_int","rho_spearman","p_spearman","lrt_p"),
      value  = c(
        res$outputs$models$cindex$c_cf_bin,
        res$outputs$models$cindex$c_cf_cont,
        res$outputs$models$cindex$c_em_bin,
        res$outputs$models$cindex$c_em_cont,
        res$outputs$models$cindex$c_int,
        res$outputs$correlation$rho,
        res$outputs$correlation$p,
        res$outputs$models$lrt$p
      )
    ),
    file.path(out_dir, paste0(prefix, "_metrics.csv"))
  )
  
  # 4) Optimal cutoff (if present)
  opt <- res$outputs$opt_cut$opt
  if (!is.null(opt)) {
    # store a human-readable one-row summary
    readr::write_csv(
      tibble::tibble(
        opt_cut_log10 = res$outputs$opt_cut$cut_log10,
        opt_cut_raw   = res$outputs$opt_cut$cut_raw,
        opt_p_logrank = opt$p %||% NA_real_,
        n_low         = opt$n_low %||% NA_integer_,
        n_high        = opt$n_high %||% NA_integer_
      ),
      file.path(out_dir, paste0(prefix, "_opt_cut_summary.csv"))
    )
    
    # store the full object (often contains the full search grid)
    saveRDS(opt, file.path(out_dir, paste0(prefix, "_opt_cut_full.rds")))
  }
  
  # 5) Save all numerical results and model objects as an RDS
  saveRDS(res, file.path(out_dir, paste0(prefix, "_FULL_RESULT.rds")))
  
  # 6) Session info for reproducibility
  writeLines(capture.output(sessionInfo()),
             file.path(out_dir, paste0(prefix, "_sessionInfo.txt")))
  
  message("Exported outputs to: ", normalizePath(out_dir))
  invisible(list(out_dir = out_dir, prefix = prefix))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


## ============================================================
## EXECUTE ANALYSES: Run landmark analyses for key timepoints
## ============================================================
## 5) Calculate Thresholds: TP05 & TP07 (ORIGINAL + PRIMARY) 
##    + Optimized for ALL OTHER Timepoints
## ============================================================
# IMPORTANT: To preserve consistency with previous analyses,
# we use the SAME METHODS that were originally used for TP05 and TP07:
#   - TP05: Uses MEDIAN threshold (as originally analyzed)
#   - TP07: Uses OPTIMIZED threshold (as originally analyzed)
# Then calculate optimized thresholds for all OTHER timepoints

cat("\n============================================================\n")
cat("Threshold Calculation Strategy:\n")
cat("  TP05: MEDIAN (original method for consistency)\n")
cat("  TP07: OPTIMIZED (original method for consistency)\n")
cat("  Others: OPTIMIZED (exploratory extension)\n")
cat("============================================================\n\n")

# Initialize storage
all_landmark_results <- list()
easyM_thresholds_calc <- numeric(0)

# ============================================================
# STEP 1: Primary timepoints (TP05 & TP07) - ORIGINAL METHODS
# ============================================================

cat("\n>>> PRIMARY ANALYSES (using original methods)...\n\n")

# TP07: 1-year maintenance (OPTIMIZED - original method)
cat("  TP07 (1-year maintenance): Optimized threshold\n")
res_07 <- build_landmark_results(tp = "07", horizon_months = 24, easyM_cut = "opt", opt_cut_min_group_n = 3)
all_landmark_results[["07"]] <- res_07
if (!is.null(res_07$outputs$opt_cut$cut_raw)) {
  easyM_thresholds_calc["07"] <- res_07$outputs$opt_cut$cut_raw
  cat(sprintf("    ✓ TP07 threshold: %.4f %% (optimized)\n", easyM_thresholds_calc["07"]))
} else {
  cat("    ⚠ Could not calculate TP07 threshold\n")
}
export_landmark_outputs(res_07, out_dir = out_dir, prefix = "TP07_1yrMaint")

# TP05: Post-transplant (MEDIAN - original method)
cat("  TP05 (post-transplant): Median threshold\n")
res_05 <- build_landmark_results(tp = "05", horizon_months = 24, easyM_cut = "median", opt_cut_min_group_n = 3)
all_landmark_results[["05"]] <- res_05
if (!is.null(res_05$outputs$opt_cut$cut_raw)) {
  easyM_thresholds_calc["05"] <- res_05$outputs$opt_cut$cut_raw
  cat(sprintf("    ✓ TP05 threshold: %.4f %% (median-based)\n", easyM_thresholds_calc["05"]))
} else {
  cat("    ⚠ Could not calculate TP05 threshold\n")
}
export_landmark_outputs(res_05, out_dir = out_dir, prefix = "TP05_Post_transplant")

cat("\n✓ Primary timepoints complete (results consistent with previous analyses)\n")

# ============================================================
# STEP 2: All OTHER timepoints - OPTIMIZED (NEW)
# ============================================================

cat("\n>>> EXPLORATORY ANALYSES (other timepoints with optimization)...\n\n")

# Get all other timepoints
all_timepoints <- df %>%
  filter(Timepoint != "01") %>%
  pull(Timepoint) %>%
  unique() %>%
  sort()

other_timepoints <- setdiff(all_timepoints, c("05", "07"))

if (length(other_timepoints) > 0) {
  cat(sprintf("Found %d additional timepoints to analyze: %s\n\n", 
              length(other_timepoints), paste(other_timepoints, collapse = ", ")))
  
  for (tp in other_timepoints) {
    tryCatch({
      cat(sprintf("  TP%s: Optimized threshold\n", tp))
      
      res <- build_landmark_results(
        tp = tp,
        horizon_months = 24,
        easyM_cut = "opt",
        opt_cut_min_group_n = 3
      )
      
      all_landmark_results[[tp]] <- res
      
      if (!is.null(res$outputs$opt_cut$opt) && !is.null(res$outputs$opt_cut$cut_raw)) {
        threshold_raw <- res$outputs$opt_cut$cut_raw
        easyM_thresholds_calc[tp] <- threshold_raw
        
        cat(sprintf("    ✓ TP%s threshold: %.4f %% (optimized)\n", tp, threshold_raw))
      } else {
        cat(sprintf("    ⚠ Could not calculate TP%s threshold (insufficient data)\n", tp))
      }
      
      # Export detailed results
      prefix <- sprintf("TP%s_landmark", tp)
      export_landmark_outputs(res, out_dir = out_dir, prefix = prefix)
      
    }, error = function(e) {
      cat(sprintf("    ✗ Error processing TP%s: %s\n", tp, e$message))
    })
  }
  
  cat("\n✓ Exploratory timepoints complete\n")
} else {
  cat("No additional timepoints beyond TP05 and TP07.\n\n")
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("FINAL THRESHOLD SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Primary Analyses (Original Methods - MAINTAINED):\n"))
cat(sprintf("  TP07: %.4f %% (optimized)\n", easyM_thresholds_calc["07"]))
cat(sprintf("  TP05: %.4f %% (median-based)\n", easyM_thresholds_calc["05"]))
if (length(other_timepoints) > 0) {
  cat(sprintf("\nExploratory Analyses (Optimized - NEW):\n"))
  for (tp in other_timepoints) {
    if (!is.na(easyM_thresholds_calc[tp])) {
      cat(sprintf("  TP%s: %.4f %% (optimized)\n", tp, easyM_thresholds_calc[tp]))
    }
  }
}
cat("============================================================\n")
cat("\nAll thresholds saved to: EasyM_threshold_values_by_timepoint.csv\n\n")

## ============================================================
## POOLED CORRELATION ANALYSIS: All post-treatment timepoints
## ============================================================
# Examine correlation between cfWGS and EasyM across all follow-up visits
# (excluding diagnosis, which has different biology)

## Get overall correlation 

# Pooled across post-treatment timepoints 03/05/07
cor_df <- survival_df %>%
  mutate(
    Timepoint   = as.character(Timepoint),
    Patient     = as.character(Patient),
    BM_prob     = as.numeric(BM_zscore_only_detection_rate_prob),
    EasyM_value = as.numeric(EasyM_value),
    EasyM_log10 = log10(pmax(EasyM_value, 1e-6))
  ) %>%
  filter(!(Timepoint %in% c("01"))) %>%
  filter(!is.na(BM_prob), !is.na(EasyM_log10)) %>%
  distinct(Patient, Timepoint, .keep_all = TRUE)

## Check normality 
## -----------------------
## 1) Quick visuals (overall)
## -----------------------
op <- par(mfrow = c(2, 2))
hist(cor_df$BM_prob, main = "BM_prob (post-treatment)", xlab = "BM_prob")
qqnorm(cor_df$BM_prob); qqline(cor_df$BM_prob)

hist(cor_df$EasyM_log10, main = "log10(EasyM_value) (post-treatment)", xlab = "log10(EasyM)")
qqnorm(cor_df$EasyM_log10); qqline(cor_df$EasyM_log10)
par(op)

## -----------------------
## 2) Formal normality tests (overall)
## Shapiro-Wilk is fine here, but can be overly sensitive; interpret with plots.
## -----------------------
sh_bm <- shapiro.test(cor_df$BM_prob)
sh_em <- shapiro.test(cor_df$EasyM_log10)

cat("\nShapiro-Wilk (overall):\n")
cat("BM_prob:     W =", round(sh_bm$statistic, 3), " p =", signif(sh_bm$p.value, 3), "\n")
cat("EasyM_log10: W =", round(sh_em$statistic, 3), " p =", signif(sh_em$p.value, 3), "\n")

## -----------------------
## 3) Normality by timepoint (03/05/07)
## (Shapiro requires at least 3 non-missing values.)
## -----------------------
shapiro_by_tp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(W = NA_real_, p = NA_real_, n = length(x)))
  out <- shapiro.test(x)
  c(W = unname(out$statistic), p = out$p.value, n = length(x))
}

by_tp <- cor_df %>%
  group_by(Timepoint) %>%
  summarise(
    BM_W = shapiro_by_tp(BM_prob)["W"],
    BM_p = shapiro_by_tp(BM_prob)["p"],
    EM_W = shapiro_by_tp(EasyM_log10)["W"],
    EM_p = shapiro_by_tp(EasyM_log10)["p"],
    n = dplyr::n(),
    n_patients = dplyr::n_distinct(Patient),
    .groups = "drop"
  )

print(by_tp)

## Optional: QQ plots by timepoint
plot_qq_by_tp <- function(df, col, label) {
  tps <- sort(unique(df$Timepoint))
  op <- par(mfrow = c(1, length(tps)))
  for (tp in tps) {
    x <- df %>% filter(Timepoint == tp) %>% pull({{col}})
    x <- x[is.finite(x)]
    qqnorm(x, main = paste0(label, " TP ", tp)); qqline(x)
  }
  par(op)
}

plot_qq_by_tp(cor_df, BM_prob, "BM_prob")
plot_qq_by_tp(cor_df, EasyM_log10, "EasyM_log10")

## -----------------------
## 4) For Pearson, what matters more is whether the *residuals* of a linear model
## are ~normal (not necessarily each variable alone).
## -----------------------
fit <- lm(BM_prob ~ EasyM_log10, data = cor_df)
res <- residuals(fit)

op <- par(mfrow = c(1, 2))
hist(res, main = "LM residuals", xlab = "Residual")
qqnorm(res); qqline(res)
par(op)

sh_res <- shapiro.test(res)
cat("\nShapiro-Wilk (LM residuals): W =", round(sh_res$statistic, 3),
    " p =", signif(sh_res$p.value, 3), "\n")

summary(fit)


# Spearman correlation
ct <- cor.test(cor_df$BM_prob, cor_df$EasyM_log10, method = "spearman", exact = FALSE)

rho <- unname(ct$estimate)
pval <- ct$p.value

# Helpful counts for reporting
n_samples  <- nrow(cor_df)
n_patients <- n_distinct(cor_df$Patient)







## ============================================================
## BASELINE CLINICAL CHARACTERISTICS
## ============================================================
# Describe the cohort's clinical features at diagnosis

### Lastly get clinical info 
### Now get info 
# 1. Define the 8 patients of interest
patients <- cor_df$Patient
patients <- patients %>% 
  unique()

# 2. Subset joined_consolidated to diagnosis timepoints for those patients
joined_consolidated <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")
clin_diag <- joined_consolidated %>%
  filter(timepoint_info == "Diagnosis", Timepoint == "01",
         Patient %in% patients)

# (3) Quick helper: simplify heavy‑/light‑chain subtype (IgG, IgA, LC‑only, Other)
clin_diag <- clin_diag %>%
  mutate(
    Ig_subtype = case_when(
      str_detect(Subtype, regex("IgG", ignore_case = TRUE)) ~ "IgG",
      str_detect(Subtype, regex("IgA", ignore_case = TRUE)) ~ "IgA",
      str_detect(Subtype, regex("light", ignore_case = TRUE)) ~ "Light‑chain",
      TRUE ~ "Other"
    )
  )

# (4) Collect counts & basic stats clinicians care about
counts <- clin_diag %>%
  summarise(
    n_total          = dplyr::n(),
    
    # Gender
    n_female         = sum(Gender == "Female", na.rm = TRUE),
    n_male           = sum(Gender == "Male",   na.rm = TRUE),
    pct_female       = 100 * n_female / n_total,
    pct_male         = 100 * n_male   / n_total,
    
    n_IgG            = sum(Ig_subtype == "IgG"),
    n_IgA            = sum(Ig_subtype == "IgA"),
    n_LC             = sum(Ig_subtype == "Light‑chain"),
    n_t11_14         = sum(T_11_14  == "Positive",          na.rm = TRUE),
    n_t4_14          = sum(T_4_14   == "Positive",          na.rm = TRUE),
    n_t14_16         = sum(T_14_16  == "Positive",          na.rm = TRUE),
    n_del17p         = sum(DEL_17P  == "Positive",          na.rm = TRUE),
    n_amp1q          = sum(AMP_1Q   == "Positive",          na.rm = TRUE),
    n_hyperdip       = sum(hyperdiploid %in% c("Yes", "Hyperdiploid"), na.rm = TRUE),
    n_high_risk      = sum(str_detect(Cytogenetic_Risk, regex("High", ignore_case = TRUE)), na.rm = TRUE),
    age_median       = median(AGE, na.rm = TRUE),
    tf_min           = min(WGS_Tumor_Fraction_Blood_plasma_cfDNA, na.rm = TRUE) * 100,
    tf_max           = max(WGS_Tumor_Fraction_Blood_plasma_cfDNA, na.rm = TRUE) * 100,
    tf_median        = median(WGS_Tumor_Fraction_Blood_plasma_cfDNA, na.rm = TRUE) * 100
  )




## ============================================================
## FINAL INTERMEDIATE EXPORT: EasyM table for 3_2 and 4_1
## ============================================================
# Create the MAIN EasyM table used by downstream manuscript-output scripts:
# 1. ALL EasyM measurements (not just landmark pairs)
# 2. Raw EasyM values
# 3. Clinician's binary calls (original)
# 4. Optimized threshold-based calls (per timepoint)
# 5. Log10 and z-score transformed values
#
# Downstream consumers:
#   - 3_2 loads this file for the Figure 3E/Figure 4C/Figure 4D EasyM panels
#     and Supplementary Table 8 source workbook.
#   - 4_1 loads this file for EasyM comparisons in the survival and
#     relapse-detection analyses.

EasyM_reference_isotypes <- build_patient_easym_reference_isotypes(dat)

readr::write_csv(
  EasyM_reference_isotypes,
  file.path(out_dir, "tables", "EasyM_reference_isotype_thresholds_by_patient.csv")
)

EasyM_comprehensive <- dat %>%
  mutate(
    Patient     = as.character(Patient),
    Timepoint   = as.character(Timepoint),
    EasyM_value = as.numeric(EasyM_value),
    EasyM_log10 = log10(pmax(EasyM_value, 1e-6))
  ) %>%
  left_join(EasyM_reference_isotypes, by = "Patient") %>%
  select(
    Patient, Timepoint, timepoint_info, Cohort, Sample_Code, Date,
    DETAILS, Subtype,
    EasyM_value, EasyM_mrd,
    BM_zscore_only_detection_rate_prob,
    BM_zscore_only_detection_rate_call,
    Blood_zscore_only_sites_prob,
    Blood_zscore_only_sites_call,
    EasyM_reference_isotype,
    EasyM_reference_threshold_percent,
    EasyM_reference_isotype_source,
    EasyM_reference_source_timepoint,
    EasyM_reference_source_details,
    EasyM_reference_source_subtype
  ) %>%
  mutate(
    # ============================================================
    # ORIGINAL CLINICIAN CALLS (for reference/comparison)
    # ============================================================
    EasyM_clinician_binary = case_when(
      EasyM_mrd == "MRD+" ~ 1L,
      EasyM_mrd == "MRD-" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # ============================================================
    # OPTIMIZED THRESHOLD-BASED CALLS (per timepoint)
    # ============================================================
    # Dynamically apply the calculated threshold for each timepoint
    EasyM_optimized_binary = case_when(
      # For each timepoint with a calculated threshold, compare volume to that threshold
      Timepoint %in% names(easyM_thresholds_calc) & !is.na(EasyM_value) ~ 
        ifelse(EasyM_value > easyM_thresholds_calc[Timepoint], 1L, 0L),
      
      # For timepoints without a calculated threshold: fall back to clinician call
      TRUE ~ case_when(
        EasyM_mrd == "MRD+" ~ 1L,
        EasyM_mrd == "MRD-" ~ 0L,
        TRUE ~ NA_integer_
      )
    ),

    # ============================================================
    # RAPID NOVOR REFERENCE THRESHOLD CALLS (isotype-specific)
    # ============================================================
    # Rapid Novor/EasyM reference negativity is defined as:
    #   IgG:          residual M-protein <= 1% of baseline
    #   IgA/LC-only:  residual M-protein <= 0.05% of baseline
    # Values above the relevant threshold are scored positive/residual.
    EasyM_reference_threshold_binary = case_when(
      !is.na(EasyM_value) & !is.na(EasyM_reference_threshold_percent) ~
        as.integer(EasyM_value > EasyM_reference_threshold_percent),
      TRUE ~ NA_integer_
    ),
    EasyM_reference_threshold_call = case_when(
      EasyM_reference_threshold_binary == 1L ~ paste0(
        "Positive/residual (>",
        sprintf("%.3f", EasyM_reference_threshold_percent),
        "%; ", EasyM_reference_isotype, ")"
      ),
      EasyM_reference_threshold_binary == 0L ~ paste0(
        "Negative/cleared (≤",
        sprintf("%.3f", EasyM_reference_threshold_percent),
        "%; ", EasyM_reference_isotype, ")"
      ),
      is.na(EasyM_reference_threshold_percent) ~ "Not scored: unsupported or unresolved isotype",
      TRUE ~ NA_character_
    ),
    EasyM_reference_threshold_method = case_when(
      !is.na(EasyM_reference_threshold_percent) ~
        "Rapid Novor reference: negative if IgG <=1%; IgA/Light Chain <=0.05% of baseline M-protein",
      TRUE ~ NA_character_
    ),
    
    # ============================================================
    # INTERPRETABLE CALL LABELS
    # ============================================================
    EasyM_optimized_call = case_when(
      # For timepoints with calculated thresholds
      Timepoint %in% names(easyM_thresholds_calc) ~ 
        ifelse(EasyM_optimized_binary == 1L,
               paste0("Residual (>", sprintf("%.3f", easyM_thresholds_calc[Timepoint]), "%)"),
               paste0("Cleared (≤", sprintf("%.3f", easyM_thresholds_calc[Timepoint]), "%)")),
      # For other timepoints, use clinician status
      TRUE ~ case_when(
        EasyM_mrd == "MRD+" ~ "Clinician: Positive",
        EasyM_mrd == "MRD-" ~ "Clinician: Negative",
        TRUE ~ NA_character_
      )
    ),
    
    # ============================================================
    # ADDITIONAL USEFUL METRICS FOR PLOTTING
    # ============================================================
    EasyM_log10 = log10(pmax(EasyM_value, 1e-6)),
    
    # Track which threshold method was used for this sample
    threshold_method = case_when(
      Timepoint %in% names(easyM_thresholds_calc) ~ "optimized",
      TRUE ~ "clinician"
    )
  ) %>%
  filter(!is.na(EasyM_value)) %>%
  arrange(Patient, Timepoint)

# Export the comprehensive table
readr::write_csv(
  EasyM_comprehensive,
  file.path(out_dir, "EasyM_all_samples_with_optimized_calls.csv")
)

cat("✓ Exported comprehensive EasyM table: EasyM_all_samples_with_optimized_calls.csv\n")
cat("  Columns: Patient, Timepoint, timepoint_info, EasyM_value, \n")
cat("           EasyM_clinician_binary, EasyM_optimized_binary, EasyM_reference_threshold_binary\n\n")

reference_status_summary <- EasyM_comprehensive %>%
  group_by(
    Timepoint,
    timepoint_info,
    EasyM_reference_isotype,
    EasyM_reference_threshold_percent,
    EasyM_reference_threshold_binary
  ) %>%
  summarise(
    n_samples = n(),
    n_patients = n_distinct(Patient),
    .groups = "drop"
  ) %>%
  mutate(
    EasyM_reference_status = case_when(
      EasyM_reference_threshold_binary == 1L ~ "Positive/residual",
      EasyM_reference_threshold_binary == 0L ~ "Negative/cleared",
      TRUE ~ "Not scored"
    )
  ) %>%
  arrange(Timepoint, EasyM_reference_isotype, EasyM_reference_status)

readr::write_csv(
  reference_status_summary,
  file.path(out_dir, "tables", "EasyM_reference_threshold_status_counts_by_timepoint.csv")
)

build_easym_reference_supplement_tables <- function(easym_df) {
  analysis_sets <- list(
    list(
      cohort = "All",
      timepoint = "All non-diagnosis",
      data = easym_df %>%
        filter(!Timepoint %in% c("01"), !str_to_lower(as.character(timepoint_info)) %in% c("diagnosis", "baseline"))
    ),
    list(
      cohort = "Frontline",
      timepoint = "Post-ASCT",
      data = easym_df %>% filter(Cohort == "Frontline", str_to_lower(as.character(timepoint_info)) == "post_transplant")
    ),
    list(
      cohort = "Frontline",
      timepoint = "Maintenance-1yr",
      data = easym_df %>% filter(Cohort == "Frontline", str_to_lower(as.character(timepoint_info)) == "1yr maintenance")
    ),
    list(
      cohort = "Non-frontline",
      timepoint = "All non-diagnosis",
      data = easym_df %>%
        filter(
          Cohort == "Non-frontline",
          !Timepoint %in% c("01"),
          !str_to_lower(as.character(timepoint_info)) %in% c("diagnosis", "baseline")
        )
    )
  )

  pred_specs <- tibble(
    pred_col = c("BM_zscore_only_detection_rate_call", "Blood_zscore_only_sites_call"),
    pred_label = c("BM cVAF Z-score cfWGS", "Blood MM-sites cfWGS")
  )

  metrics <- purrr::map_dfr(analysis_sets, function(analysis_set) {
    purrr::map2_dfr(
      pred_specs$pred_col,
      pred_specs$pred_label,
      ~ compute_easym_reference_contingency(
        df = analysis_set$data,
        pred_col = .x,
        pred_label = .y,
        cohort_label = analysis_set$cohort,
        timepoint_label = analysis_set$timepoint
      )
    )
  }) %>%
    mutate(
      Agree_str = sprintf("%d/%d", Agree, N),
      Concordance_pct = scales::percent(Concordance, accuracy = 1),
      Sensitivity_pct = scales::percent(Sensitivity, accuracy = 1),
      Specificity_pct = scales::percent(Specificity, accuracy = 1),
      PPV_pct = scales::percent(PPV, accuracy = 1),
      NPV_pct = scales::percent(NPV, accuracy = 1)
    ) %>%
    select(
      Cohort, Timepoint, Pred_Column, Pred_Label, Comparator,
      Agree_str, Concordance_pct, Sensitivity_pct, Specificity_pct, PPV_pct, NPV_pct,
      N, Agree, TP, FP, FN, TN, Concordance, Sensitivity, Specificity, PPV, NPV, Matrix
    )

  cells <- metrics %>%
    select(Cohort, Timepoint, Pred_Label, TN, FP, FN, TP) %>%
    pivot_longer(
      cols = c(TN, FP, FN, TP),
      names_to = "Cell",
      values_to = "Count"
    ) %>%
    mutate(
      Observed_EasyM_status = case_when(Cell %in% c("TN", "FP") ~ "neg", TRUE ~ "pos"),
      Predicted_cfWGS_status = case_when(Cell %in% c("TN", "FN") ~ "neg", TRUE ~ "pos")
    ) %>%
    select(Cohort, Timepoint, Pred_Label, Observed_EasyM_status, Predicted_cfWGS_status, Cell, Count)

  sample_level <- easym_df %>%
    select(
      Patient, Timepoint, timepoint_info, Cohort, Sample_Code, Date,
      EasyM_value, EasyM_reference_isotype, EasyM_reference_threshold_percent,
      EasyM_reference_threshold_binary, EasyM_reference_threshold_call,
      BM_zscore_only_detection_rate_call, Blood_zscore_only_sites_call
    ) %>%
    filter(!is.na(EasyM_reference_threshold_binary)) %>%
    arrange(Patient, Timepoint)

  list(Metrics = metrics, Contingency_Cells = cells, Sample_Level = sample_level)
}

easym_reference_supp_tables <- build_easym_reference_supplement_tables(EasyM_comprehensive)

readr::write_csv(
  easym_reference_supp_tables$Metrics,
  file.path(out_dir, "tables", "EasyM_reference_threshold_contingency_vs_cfWGS.csv")
)
readr::write_csv(
  easym_reference_supp_tables$Contingency_Cells,
  file.path(out_dir, "tables", "EasyM_reference_threshold_2x2_cells_vs_cfWGS.csv")
)
readr::write_csv(
  easym_reference_supp_tables$Sample_Level,
  file.path(out_dir, "tables", "EasyM_reference_threshold_sample_level_calls.csv")
)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    easym_reference_supp_tables,
    path = file.path(out_dir, "tables", "EasyM_reference_threshold_supplement_contingency_tables.xlsx")
  )
}

cat("✓ Exported Rapid Novor reference-threshold EasyM contingency tables\n")
cat("  - tables/EasyM_reference_threshold_contingency_vs_cfWGS.csv\n")
cat("  - tables/EasyM_reference_threshold_2x2_cells_vs_cfWGS.csv\n")
cat("  - tables/EasyM_reference_threshold_sample_level_calls.csv\n\n")

# ============================================================
# CUTOFF SUMMARY TABLE
# ============================================================
# Document all thresholds used for reference

easyM_cutoff_summary <- data.frame(
  Timepoint = names(easyM_thresholds_calc),
  Threshold_raw_percent = as.numeric(easyM_thresholds_calc),
  Threshold_log10 = log10(as.numeric(easyM_thresholds_calc)),
  Interpretation = case_when(
    names(easyM_thresholds_calc) == "01" ~ "Diagnosis",
    names(easyM_thresholds_calc) == "03" ~ "Post-Induction",
    names(easyM_thresholds_calc) == "05" ~ "Post-Transplant",
    names(easyM_thresholds_calc) == "07" ~ "1-Year Maintenance",
    names(easyM_thresholds_calc) == "08" ~ "Follow-up",
    names(easyM_thresholds_calc) == "09" ~ "Follow-up",
    names(easyM_thresholds_calc) == "10" ~ "Follow-up",
    names(easyM_thresholds_calc) == "11" ~ "Follow-up",
    names(easyM_thresholds_calc) == "12" ~ "Follow-up",
    names(easyM_thresholds_calc) == "13" ~ "Follow-up",
    TRUE ~ "Relapse"
  ),
  stringsAsFactors = FALSE
) %>%
  arrange(Timepoint)

readr::write_csv(
  easyM_cutoff_summary,
  file.path(out_dir, "EasyM_threshold_values_by_timepoint.csv")
)

cat("✓ Exported threshold summary: EasyM_threshold_values_by_timepoint.csv\n")
cat(sprintf("  Calculated thresholds for %d timepoints:\n", nrow(easyM_cutoff_summary)))
print(easyM_cutoff_summary)
cat("\n")

# ============================================================
# SAMPLE COUNTS BY STATUS
# ============================================================
# Quick sanity check on how many patients are positive/negative

status_summary <- EasyM_comprehensive %>%
  group_by(Timepoint, EasyM_optimized_binary, threshold_method) %>%
  summarise(
    n_samples = n(),
    n_patients = n_distinct(Patient),
    .groups = "drop"
  ) %>%
  mutate(
    Status = ifelse(EasyM_optimized_binary == 1L, "Positive", "Negative"),
    .keep = "unused"
  ) %>%
  arrange(Timepoint, Status)

cat("EasyM Status Distribution by Timepoint (Optimized Calls):\n")
print(status_summary, n = Inf)
cat("\n")

# ============================================================
# COMPLETION SUMMARY
# ============================================================
cat("\n", strrep("=", 70), "\n")
cat("✓ SCRIPT COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 70), "\n\n")

# Print summary of key files exported
cat("KEY OUTPUT FILES FOR DOWNSTREAM ANALYSIS:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Main export table
main_export <- file.path(out_dir, "EasyM_all_samples_with_optimized_calls.csv")
cat(sprintf("1. Comprehensive EasyM table with reference, any-detect, and optimized calls:\n   ├─ Path: %s\n\n", main_export))

# Threshold reference
threshold_export <- file.path(out_dir, "EasyM_threshold_values_by_timepoint.csv")
cat(sprintf("2. Threshold values reference:\n   ├─ Path: %s\n\n", threshold_export))

# Landmark analyses directory
landmarks_dir <- file.path(out_dir, "landmark_analyses")
cat(sprintf("3. Detailed landmark analysis results:\n   ├─ Directory: %s\n   ├─ Contains per-timepoint analysis results\n\n", landmarks_dir))

cat("THRESHOLD METHODS USED:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("✓ TP05 (post-ASCT):        Median (original method - preserved)\n")
cat("✓ TP07 (1-year maintenance): Optimized via log-rank (original method - preserved)\n")
cat("✓ Other timepoints:         Optimized via log-rank (exploratory - new)\n\n")

cat("NEXT STEPS:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("→ 3_2 loads EasyM_all_samples_with_optimized_calls.csv for EasyM panels in Figure 3E, Figure 4C, and Figure 4D\n")
cat("→ 4_1 loads EasyM_all_samples_with_optimized_calls.csv for EasyM survival and relapse-detection summaries\n")
cat("→ Use EasyM_reference_threshold_binary for EasyM binary comparisons reported in the manuscript\n")
cat("→ EasyM_optimized_binary is retained as an exploratory timepoint-specific analysis\n\n")

## End of script
