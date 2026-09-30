# =============================================================================
# Script: 2_1_Clinical_Demographics_Table.R
#
# Description:
#   This script processes the final aggregated cfWGS feature table (with clinical
#   and demographic annotations), selects one Baseline/Diagnosis record per
#   patient, attaches the existing Frontline/Non-frontline cohort assignment,
#   and builds the categorical manuscript Table 1 of baseline characteristics.
#
# Unit of analysis:
#   One patient. When multiple eligible baseline/diagnosis records exist, the
#   documented selection rules below retain exactly one record for Table 1.
#
# Steps:
#   1. Load the aggregated RDS of cfWGS features + clinical/demographics
#   2. Attach the current patient-level cohort labels
#   3. Subset to Diagnosis/Baseline and resolve duplicate baseline candidates
#   4. Fill missing patient-level fields from curated clinical sources, with audit
#   5. Clean and recode the categorical variables shown in Table 1
#   6. Export deterministic count/denominator companions and the formatted table
#   7. Stage the current DOCX/PDF manuscript artifacts
#   8. Optionally build historical continuous-variable and cohort-review helpers
#
# Inputs:
#   • RDS: Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds
#   • RDS: baseline_high_quality_patients_updated.rds
#   • CSV: Output_tables_2025/patient_cohort_assignment.csv (comparison audit)
#   • Curated cohort-specific clinical workbooks listed in
#     read_table1_patient_clinical_sources() (used only to fill missing fields)
#
# Outputs:
#   • Word: table1_clinical_categorical_updated_final_3.docx
#   • CSV: Output_tables_2025/Table_1_clinical_demographics_computed_source_counts.csv
#   • TSV: Output_tables_2025/Table_1_clinical_demographics_computed_source_qc.tsv
#   • DOCX/PDF: current Table 1 artifacts under Final Tables and Figures/ and
#     Figures_Exported/Tables_exported/
#   • Optional Word helper: baseline_characteristics_updated.docx
#   • Optional exploratory cohort/fragmentomics audit output printed to console
#
# Dependencies:
#   library(tidyverse)
#   library(gtsummary)
#   library(officer)
#   library(flextable)
#
# Usage:
#   Rscript Scripts_2025/Final_Scripts/2_1_Clinical_Demographics_Table.R
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_1_Clinical_Demographics_Table.R
#
# Manuscript outputs created/updated:
#   - Table 1: baseline clinical/demographic characteristics, exported as the
#     generated Word source and, when present, the manually edited final DOCX/PDF
#     used for the manuscript.
#
# Pipeline role:
#   Table 1 is built from the baseline/diagnosis sample per patient after
#   resolving known duplicate baseline records. The cohort assignment exported
#   by the upstream cohort workflow is consumed here so the Training/Testing
#   display groups remain consistent with downstream analyses.
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

# -----------------------------------------------------------
# 0.  PACKAGES  (install once, then keep only library() lines)
# -----------------------------------------------------------
# install.packages(c("tidyverse", "gtsummary", "gt", "officer"))   # ← run once
library(tidyverse)
library(gtsummary)
library(officer)
library(flextable)

# Shared manuscript-output helpers.
# These copy the exact final Table 1 artifacts into final_manuscript_objects/
# while retaining this script as the place where the Table 1 source is built.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

.publication_export_helper <- file.path(
  "Scripts_2025", "Final_Scripts", "publication_export_helpers.R"
)
if (!file.exists(.publication_export_helper)) {
  .publication_export_helper <- "publication_export_helpers.R"
}
source(.publication_export_helper)
rm(.publication_export_helper)


# -----------------------------------------------------------
# 1. Load baseline table inputs
# -----------------------------------------------------------
# Table 1 should reflect the current revision-inclusive aggregate table.
dat <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")


### Load the current Table 1 cohort labels
# Preserve the original Table 1 baseline-row selection logic below. The only
# denominator change is to use the current full Table 1 cohort labels from
# `baseline_high_quality_patients_updated.rds` instead of filtering by
# Output_tables_2025/patient_cohort_assignment.csv, which is a narrower
# 51-patient intermediate and excludes current Non-frontline/test cases.
export_dir <- "Output_tables_2025"
baseline_high_quality_path <- "baseline_high_quality_patients_updated.rds"
if (!file.exists(baseline_high_quality_path)) {
  stop("Missing current Table 1 cohort helper: ", baseline_high_quality_path, call. = FALSE)
}
baseline_high_quality_patients <- readRDS(baseline_high_quality_path)

required_baseline_cols <- c("Patient", "Cohort")
missing_baseline_cols <- setdiff(required_baseline_cols, names(baseline_high_quality_patients))
if (length(missing_baseline_cols)) {
  stop(
    "Current Table 1 cohort helper is missing required columns: ",
    paste(missing_baseline_cols, collapse = ", "),
    call. = FALSE
  )
}

table1_current_cohort_tbl <- baseline_high_quality_patients %>%
  mutate(
    Patient = as.character(Patient),
    Cohort = as.character(Cohort)
  ) %>%
  distinct(Patient, Cohort)

cohort_dups <- table1_current_cohort_tbl %>%
  count(Patient) %>%
  filter(n > 1)
if (nrow(cohort_dups)) {
  stop(
    "Current Table 1 cohort labels must have one cohort per patient but duplicate patients were found: ",
    paste(cohort_dups$Patient, collapse = ", "),
    call. = FALSE
  )
}

required_clinical_cols <- c(
  "Patient", "Date", "Timepoint", "Sample_Code", "timepoint_info",
  "Gender", "AGE_GROUP", "ISS_STAGE", "Cytogenetic_Risk", "Subtype",
  "ECOG_SCORE", "T_4_14", "T_14_16", "DEL_17P"
)
missing_clinical_cols <- setdiff(required_clinical_cols, names(dat))
if (length(missing_clinical_cols)) {
  stop(
    "Final aggregate Table 1 input is missing required clinical columns: ",
    paste(missing_clinical_cols, collapse = ", "),
    call. = FALSE
  )
}

### Select baseline/diagnosis rows using the original Table 1 logic.
# 1. Subset to only Diagnosis or Baseline timepoints
dat_tb <- dat %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  mutate(Date = as.Date(Date))

# 2. Remove CA-03 timepoint 02 when encoded as timepoint_info
dat_tb2 <- dat_tb %>%
  filter(!(Patient == "CA-03" & timepoint_info == "02"))

# 3. Consolidate the two CA-02 rows
resp_CA02 <- dat_tb2 %>%
  filter(Patient == "CA-02") %>%
  arrange(factor(timepoint_info, levels = c("01","02"))) %>%
  summarise(across(everything(), ~{
    vals <- .
    first_val <- vals[which(!is.na(vals))[1]]
    other <- vals[!is.na(vals) & vals != first_val][1]
    if (!is.na(first_val) &&
        first_val %in% c("Unknown","Other") &&
        !is.na(other) &&
        !other %in% c("Unknown","Other")) {
      other
    } else {
      first_val
    }
  }))

if (!nrow(resp_CA02)) {
  stop("Table 1 duplicate-resolution logic expected at least one CA-02 row.", call. = FALSE)
}

# 4. Drop all CA-02 originals and keep the intended SPORE_0009 baseline row
dat_tb_final <- dat_tb2 %>%
  filter(Patient != "CA-02") %>%
  filter_manual_baseline_row_selection("spore0009_baseline_only") %>%
  bind_rows(resp_CA02) %>%
  arrange(Patient)

dat_base <- dat_tb_final %>%
  filter(Patient %in% table1_current_cohort_tbl$Patient) %>%
  filter(!(Patient == "CA-03" & Timepoint == "02"))

baseline_duplicate_audit <- dat_base %>%
  group_by(Patient) %>%
  mutate(n_patient_baseline_candidates = n()) %>%
  ungroup() %>%
  filter(n_patient_baseline_candidates > 1) %>%
  mutate(
    patient_baseline_timepoint_rank = case_when(
      str_detect(as.character(Timepoint), regex("^T?0$", ignore_case = TRUE)) ~ 0L,
      str_detect(as.character(Timepoint), regex("^T?1$|^0?1$", ignore_case = TRUE)) ~ 1L,
      TRUE ~ 2L
    ),
    patient_baseline_has_feature_evidence = !is.na(BM_Mutation_Count) | !is.na(Blood_Mutation_Count)
  ) %>%
  group_by(Patient) %>%
  arrange(
    is.na(Date),
    Date,
    patient_baseline_timepoint_rank,
    desc(patient_baseline_has_feature_evidence),
    Timepoint,
    .by_group = TRUE
  ) %>%
  mutate(
    selected_for_patient_baseline_table = row_number() == 1L,
    patient_baseline_selection_reason = case_when(
      selected_for_patient_baseline_table ~ "selected earliest dated baseline/diagnosis candidate, preferring T0/T1 and rows with WGS mutation evidence",
      TRUE ~ "not selected for patient-level Table 1 to enforce one baseline/diagnosis row per patient"
    )
  ) %>%
  ungroup()

dir.create(file.path(export_dir, "clinical_support"), recursive = TRUE, showWarnings = FALSE)
if (nrow(baseline_duplicate_audit) > 0L) {
  readr::write_csv(
    baseline_duplicate_audit,
    file.path(export_dir, "clinical_support", "table1_patient_baseline_duplicate_selection_audit.csv")
  )
  dat_base <- dat_base %>%
    left_join(
      baseline_duplicate_audit %>%
        select(Patient, Sample_Code, Timepoint, selected_for_patient_baseline_table),
      by = c("Patient", "Sample_Code", "Timepoint")
    ) %>%
    filter(is.na(selected_for_patient_baseline_table) | selected_for_patient_baseline_table) %>%
    select(-selected_for_patient_baseline_table)
}

dups <- dat_base %>%
  count(Patient) %>%
  filter(n > 1)

if (nrow(dups)) {
  stop(
    "Still multiple Diagnosis/Baseline rows after Table 1 duplicate selection for: ",
    paste(unique(dups$Patient), collapse = ", "),
    call. = FALSE
  )
} else {
  message("All patients now have at most one Diagnosis/Baseline row.")
}

dat_base <- dat_base %>%
  left_join(table1_current_cohort_tbl, by = "Patient")

valid_cohorts <- c("Frontline", "Non-frontline")
invalid_cohorts <- setdiff(unique(na.omit(dat_base$Cohort)), valid_cohorts)
if (length(invalid_cohorts)) {
  stop(
    "Unexpected Table 1 cohort labels. Expected Frontline/Non-frontline but found: ",
    paste(invalid_cohorts, collapse = ", "),
    call. = FALSE
  )
}

### Fill missing Table 1 clinical fields from curated patient-level sources.
# This preserves the original row-selection logic above and only rescues missing
# patient-level fields for Table 1. Every filled value is audited below.
is_missing_label <- function(x) {
  x_chr <- trimws(as.character(x))
  is.na(x_chr) | !nzchar(x_chr) |
    toupper(x_chr) %in% c("NA", "N/A", "UNKNOWN", "UNKNOWN/MISSING", "MISSING")
}

first_non_missing_value <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr <- x_chr[!is_missing_label(x_chr)]
  if (length(x_chr)) x_chr[[1]] else NA_character_
}

coalesce_missing_label <- function(current, rescue) {
  dplyr::if_else(is_missing_label(current) & !is_missing_label(rescue), as.character(rescue), as.character(current))
}

normalize_gender <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    toupper(x_chr) %in% c("M", "MALE") ~ "Male",
    toupper(x_chr) %in% c("F", "FEMALE") ~ "Female",
    TRUE ~ x_chr
  )
}

age_group_from_age <- function(age) {
  age_num <- suppressWarnings(as.numeric(age))
  dplyr::case_when(
    age_num < 50 ~ "<50",
    age_num >= 50 & age_num < 60 ~ "50-59",
    age_num >= 60 & age_num < 70 ~ "60-69",
    age_num >= 70 & age_num < 80 ~ "70-79",
    age_num >= 80 ~ "80+",
    TRUE ~ NA_character_
  )
}

normalize_stage <- function(x) {
  x_chr <- stringr::str_squish(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x_chr, stringr::regex("^stage\\s*1$|^stage\\s*i$", ignore_case = TRUE)) ~ "Stage I",
    stringr::str_detect(x_chr, stringr::regex("^stage\\s*2$|^stage\\s*ii$", ignore_case = TRUE)) ~ "Stage II",
    stringr::str_detect(x_chr, stringr::regex("^stage\\s*3$|^stage\\s*iii$", ignore_case = TRUE)) ~ "Stage III",
    TRUE ~ x_chr
  )
}

normalize_marker <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    toupper(x_chr) %in% c("YES", "Y", "TRUE", "1", "POSITIVE") ~ "Positive",
    toupper(x_chr) %in% c("NO", "N", "FALSE", "0", "NEGATIVE") ~ "Negative",
    TRUE ~ NA_character_
  )
}

build_isotype_subtype <- function(heavy_chain, light_chain) {
  heavy <- trimws(as.character(heavy_chain))
  light <- trimws(as.character(light_chain))
  heavy[is_missing_label(heavy)] <- NA_character_
  light[is_missing_label(light)] <- NA_character_
  dplyr::case_when(
    is.na(heavy) & is.na(light) ~ NA_character_,
    toupper(heavy) %in% c("FLC", "LIGHT CHAIN", "LIGHT CHAIN ONLY", "NEGATIVE") & !is.na(light) ~
      paste("Light Chain", light),
    !is.na(heavy) & !is.na(light) ~ paste(heavy, light),
    !is.na(heavy) ~ paste(heavy, "(Subtype Unknown)"),
    TRUE ~ NA_character_
  )
}

normalize_img_patient_id <- function(x) {
  x_chr <- trimws(as.character(x))
  direct <- stringr::str_extract(x_chr, "IMG-[0-9]+")
  numeric_id <- stringr::str_extract(x_chr, "[0-9]+")
  dplyr::case_when(
    !is.na(direct) ~ sprintf("IMG-%03d", as.integer(stringr::str_extract(direct, "[0-9]+"))),
    !is.na(numeric_id) ~ sprintf("IMG-%03d", as.integer(numeric_id)),
    TRUE ~ NA_character_
  )
}

read_table1_patient_clinical_sources <- function() {
  clinical_sources <- list()

  imagine_curated_path <- "Clinical data/IMMAGINE/Clinical data for IMG patients at diagnosis_filled_DA_edited.xlsx"
  if (file.exists(imagine_curated_path)) {
    clinical_sources[["imagine_curated_diagnosis"]] <- readxl::read_excel(imagine_curated_path) %>%
      transmute(
        Patient = as.character(Patient_ID),
        Gender = normalize_gender(GENDER),
        AGE_GROUP = age_group_from_age(AGE_AT_CONSENT),
        ISS_STAGE = normalize_stage(ISS_STAGE),
        Cytogenetic_Risk = as.character(Cytogenetic_Risk),
        Subtype = as.character(Subtype),
        ECOG_SCORE = as.character(ECOG_score),
        T_4_14 = NA_character_,
        T_14_16 = NA_character_,
        DEL_17P = NA_character_,
        table1_clinical_rescue_source = imagine_curated_path
      )
  }

  spore_curated_path <- "Clinical data/SPORE/Spore_baseline_clinical_demographics_DA_edited.xlsx"
  if (file.exists(spore_curated_path)) {
    clinical_sources[["spore_curated_baseline"]] <- readxl::read_excel(spore_curated_path) %>%
      transmute(
        Patient = as.character(Patient_ID),
        Gender = normalize_gender(GENDER),
        AGE_GROUP = age_group_from_age(AGE_AT_CONSENT),
        ISS_STAGE = normalize_stage(ISS_STAGE),
        Cytogenetic_Risk = NA_character_,
        Subtype = as.character(Subtype),
        ECOG_SCORE = as.character(ECOG_score),
        T_4_14 = NA_character_,
        T_14_16 = NA_character_,
        DEL_17P = NA_character_,
        table1_clinical_rescue_source = spore_curated_path
      )
  }

  oicr_summary_path <- "New OICR Submissions/derived_metadata/oicr_submission_patient_clinical_summary.csv"
  if (file.exists(oicr_summary_path)) {
    clinical_sources[["oicr_patient_summary"]] <- readr::read_csv(oicr_summary_path, show_col_types = FALSE) %>%
      transmute(
        Patient = as.character(patient_img_id),
        Gender = normalize_gender(Sex),
        AGE_GROUP = age_group_from_age(`Age at Diagnosis`),
        ISS_STAGE = normalize_stage(`Risk Status at Diagnosis`),
        Cytogenetic_Risk = NA_character_,
        Subtype = build_isotype_subtype(`Heavy Chain Isotype`, `Light Chain Isotype`),
        ECOG_SCORE = NA_character_,
        T_4_14 = normalize_marker(`t(4;14)`),
        T_14_16 = normalize_marker(`t(14;16)`),
        DEL_17P = normalize_marker(del17p),
        table1_clinical_rescue_source = oicr_summary_path
      )
  }

  esther_liberate_paths <- c(
    "Clinical data/IMMAGINE/IMMAGINE  LIBERATE_ID_noMRNorDoB_6Feb2026.xlsx",
    "Clinical data/IMMAGINE/IMMAGINE  LIBERATE_ID_noMRN_6Feb2026.xlsx",
    "Clinical data/IMMAGINE/IMMAGINE  LIBERATE_MRD_withoutMRN_or_DOB_18Dec2025.xlsx"
  )
  for (esther_path in esther_liberate_paths[file.exists(esther_liberate_paths)]) {
    esther_sheets <- readxl::excel_sheets(esther_path)
    for (sheet_name in esther_sheets) {
      esther_sheet <- tryCatch(
        readxl::read_excel(esther_path, sheet = sheet_name),
        error = function(e) tibble()
      )
      if (!nrow(esther_sheet) || !"GENDER" %in% names(esther_sheet)) {
        next
      }
      id_column <- dplyr::case_when(
        "STUDY_NAME" %in% names(esther_sheet) ~ "STUDY_NAME",
        "PATIENT_ID" %in% names(esther_sheet) ~ "PATIENT_ID",
        TRUE ~ NA_character_
      )
      if (is.na(id_column)) {
        next
      }
      source_name <- paste0(esther_path, " [", sheet_name, "]")
      clinical_sources[[paste0("esther_imagine_liberate_", length(clinical_sources) + 1L)]] <-
        esther_sheet %>%
        transmute(
          Patient = normalize_img_patient_id(.data[[id_column]]),
          Gender = normalize_gender(GENDER),
          AGE_GROUP = NA_character_,
          ISS_STAGE = NA_character_,
          Cytogenetic_Risk = NA_character_,
          Subtype = NA_character_,
          ECOG_SCORE = NA_character_,
          T_4_14 = NA_character_,
          T_14_16 = NA_character_,
          DEL_17P = NA_character_,
          table1_clinical_rescue_source = source_name
        )
    }
  }

  imagine_pcm_path <- "Clinical data/IMMAGINE/IMMAGINE_PCM_Import_24Mar2023/data_clinical_patients.txt"
  if (file.exists(imagine_pcm_path)) {
    clinical_sources[["imagine_cbio_patient_table"]] <- readr::read_tsv(
      imagine_pcm_path,
      comment = "#",
      show_col_types = FALSE
    ) %>%
      transmute(
        Patient = as.character(PATIENT_ID),
        Gender = normalize_gender(SEX_AT_BIRTH),
        AGE_GROUP = age_group_from_age(AGE),
        ISS_STAGE = normalize_stage(RISK_STATUS_DIAGNOSIS),
        Cytogenetic_Risk = NA_character_,
        Subtype = build_isotype_subtype(IGH_ISOTYPE, IGL_ISOTYPE),
        ECOG_SCORE = NA_character_,
        T_4_14 = normalize_marker(TX_T4_14),
        T_14_16 = normalize_marker(TX_T14_16),
        DEL_17P = normalize_marker(CNV_DEL17P),
        table1_clinical_rescue_source = imagine_pcm_path
      )
  }

  if (!length(clinical_sources)) {
    return(tibble())
  }

  rescue_fields <- c(
    "Gender", "AGE_GROUP", "ISS_STAGE", "Cytogenetic_Risk", "Subtype",
    "ECOG_SCORE", "T_4_14", "T_14_16", "DEL_17P"
  )

  rescue_raw <- bind_rows(clinical_sources) %>%
    filter(!is_missing_label(Patient))

  rescue_values <- rescue_raw %>%
    filter(!is_missing_label(Patient)) %>%
    group_by(Patient) %>%
    summarise(
      across(
        all_of(rescue_fields),
        first_non_missing_value
      ),
      table1_clinical_rescue_source = paste(unique(table1_clinical_rescue_source), collapse = "; "),
      .groups = "drop"
    )

  rescue_field_sources <- rescue_raw %>%
    pivot_longer(
      cols = all_of(rescue_fields),
      names_to = "field",
      values_to = "value"
    ) %>%
    filter(!is_missing_label(value)) %>%
    group_by(Patient, field) %>%
    summarise(
      table1_clinical_rescue_field_source = first(table1_clinical_rescue_source),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = field,
      values_from = table1_clinical_rescue_field_source,
      names_prefix = "table1_clinical_rescue_source_"
    )

  rescue_values %>%
    left_join(rescue_field_sources, by = "Patient")
}

table1_clinical_rescue <- read_table1_patient_clinical_sources()

table1_rescue_fields <- c(
  "Gender", "AGE_GROUP", "ISS_STAGE", "Cytogenetic_Risk", "Subtype",
  "ECOG_SCORE", "T_4_14", "T_14_16", "DEL_17P"
)

if (nrow(table1_clinical_rescue)) {
  dat_before_clinical_rescue <- dat_base %>%
    select(Patient, all_of(table1_rescue_fields)) %>%
    mutate(across(all_of(table1_rescue_fields), as.character))

  dat_base_rescue_joined <- dat_base %>%
    left_join(
      table1_clinical_rescue %>%
        rename_with(~ paste0(.x, "_table1_rescue"), all_of(table1_rescue_fields)),
      by = "Patient"
    )

  for (field in table1_rescue_fields) {
    rescue_field <- paste0(field, "_table1_rescue")
    dat_base_rescue_joined[[field]] <- coalesce_missing_label(
      dat_base_rescue_joined[[field]],
      dat_base_rescue_joined[[rescue_field]]
    )
  }

  dat_base <- dat_base_rescue_joined %>%
    select(-ends_with("_table1_rescue"))

  dat_after_clinical_rescue <- dat_base %>%
    select(Patient, all_of(table1_rescue_fields)) %>%
    mutate(across(all_of(table1_rescue_fields), as.character))

  table1_clinical_rescue_field_sources_long <- table1_clinical_rescue %>%
    select(Patient, starts_with("table1_clinical_rescue_source_")) %>%
    pivot_longer(
      cols = -Patient,
      names_to = "field",
      values_to = "table1_clinical_rescue_field_source",
      names_prefix = "table1_clinical_rescue_source_"
    )

  table1_clinical_rescue_audit <- dat_before_clinical_rescue %>%
    rename_with(~ paste0(.x, "_before"), all_of(table1_rescue_fields)) %>%
    left_join(
      dat_after_clinical_rescue %>%
        rename_with(~ paste0(.x, "_after"), all_of(table1_rescue_fields)),
      by = "Patient"
    ) %>%
    pivot_longer(
      cols = -Patient,
      names_to = c("field", ".value"),
      names_pattern = "(.+)_(before|after)"
    ) %>%
    left_join(
      table1_clinical_rescue_field_sources_long,
      by = c("Patient", "field")
    ) %>%
    filter(is_missing_label(before), !is_missing_label(after)) %>%
    arrange(Patient, field)

  readr::write_csv(
    table1_clinical_rescue_audit,
    file.path(export_dir, "clinical_support", "table1_patient_clinical_rescue_audit.csv")
  )
}

table1_special_case_audit <- dat_base %>%
  filter(Patient %in% c("CA-02", "CA-03", "SPORE_0009", "IMG-142", "IMG-235")) %>%
  transmute(
    Patient,
    Date,
    Timepoint,
    Sample_Code,
    timepoint_info,
    Cohort,
    table1_special_case_expectation = case_when(
      Patient == "CA-02" ~ "selected diagnosis row CA-02-01, not duplicate 02 row",
      Patient == "CA-03" ~ "selected diagnosis row CA-03-01, not duplicate 02 row",
      Patient == "SPORE_0009" ~ "selected audited baseline row from docs/manual_baseline_row_selection.csv",
      Patient %in% c("IMG-142", "IMG-235") ~ "selected one raw T1/T0 baseline row with original Table 1 duplicate-ranking logic",
      TRUE ~ NA_character_
    )
  )

readr::write_csv(
  table1_special_case_audit,
  file.path(export_dir, "clinical_support", "table1_special_case_selection_audit.csv")
)

table1_special_case_failures <- table1_special_case_audit %>%
  filter(
    (Patient == "CA-02" & !(Sample_Code == "CA-02-01" & Timepoint == "01")) |
      (Patient == "CA-03" & !(Sample_Code == "CA-03-01" & Timepoint == "01")) |
      (Patient == "SPORE_0009" & Date != load_manual_baseline_row_selection("spore0009_baseline_only")$Date[[1]])
  )
if (nrow(table1_special_case_failures)) {
  stop(
    "Table 1 special-case baseline selection drifted for: ",
    paste(table1_special_case_failures$Patient, collapse = ", "),
    ". Inspect ",
    file.path(export_dir, "clinical_support", "table1_special_case_selection_audit.csv"),
    call. = FALSE
  )
}

patient_cohort_assignment_path <- file.path(export_dir, "patient_cohort_assignment.csv")
patient_cohort_tbl_csv <- if (file.exists(patient_cohort_assignment_path)) {
  read.csv(patient_cohort_assignment_path, stringsAsFactors = FALSE)
} else {
  data.frame(Patient = character(), Cohort = character())
}
if (!all(c("Patient", "Cohort") %in% names(patient_cohort_tbl_csv))) {
  warning("Narrow patient cohort assignment audit CSV lacks Patient/Cohort columns: ", patient_cohort_assignment_path)
  patient_cohort_tbl_csv <- data.frame(Patient = character(), Cohort = character())
}

message("Table 1 baseline cohort patients: ", nrow(dat_base))


dat_base <- dat_base %>%
  mutate(
    cohort = factor(
      publication_cohort_label(Cohort),
      levels = c("Training", "Testing")
    )
  )


# -----------------------------------------------------------
# 3. Recode variables shown in Table 1
# -----------------------------------------------------------
# ------------------------- 2.1  Convert continuous vars to numeric
# 1. Collapse Subtype
dat_base <- dat_base %>%
  mutate(Subtype = fct_collapse(
    Subtype,
    IgG              = c("IgG Kappa", "IgG Lambda", "IgG (Subtype Unknown)"),
    IgA              = c("IgA Kappa", "IgA Lambda"),
    `Light Chain Only` = c(
      "Light Chain Only (Kappa)", "Light Chain Only (Lambda)",
      "Light Chain Kappa", "Light Chain Lambda", "FLC Kappa", "FLC Lambda"
    ),
    Other            = c("Other", NA),
    .default         = "Other"
  ))

# Fix issue 
dat_base <- dat_base %>%
  mutate(Subtype = na_if(as.character(Subtype), ".default")) %>%
  mutate(Subtype = factor(Subtype))  # re-factor after conversion


# 2. Define the final var lists
vars_cat <- c("Gender", "AGE_GROUP", "ISS_STAGE",
              "Cytogenetic_Risk", "Subtype", "ECOG_SCORE")

### Most things NA, so trying best with what isn't NA for everything at least across both groups 
## Plasma cll percent NA for most

# 2. Recode “Unknown” and “Unknown/Missing” to NA and drop them
dat_base <- dat_base %>%
  # 1) Convert categorical variables to character
  mutate(across(all_of(vars_cat), as.character)) %>%
  
  # 2) Replace "Unknown" or "Unknown/Missing" with NA
  mutate(across(
    all_of(vars_cat),
    ~ if_else(.x %in% c("Unknown", "Unknown/Missing"), NA_character_, .x)
  )) %>%
  
  # 3) Back to factor and drop any empty levels
  mutate(across(all_of(vars_cat), as.factor)) %>%
  mutate(across(all_of(vars_cat), fct_drop))

### Fill missing cytogenetic risk where FISH markers are sufficient
# This does not overwrite an existing risk label. It only derives risk for rows
# where Cytogenetic_Risk is missing and the high-risk FISH markers are present.
dat_base <- dat_base %>%
  mutate(
    # 1. Calculate risk from FISH results
    Cytogenetic_Risk_calc = case_when(
      # Any high-risk marker positive → High Risk
      T_4_14  == "Positive" |
        T_14_16 == "Positive" |
        DEL_17P == "Positive"    ~ "High Risk",
      
      # All three markers explicitly NEGATIVE → Standard Risk
      !is.na(T_4_14)  & T_4_14  == "Negative" &
        !is.na(T_14_16) & T_14_16 == "Negative" &
        !is.na(DEL_17P) & DEL_17P == "Negative" ~ "Standard Risk",
      
      # Otherwise (any missing or contradictory data) → remain NA
      TRUE                      ~ NA_character_
    ),
    
    # 2. Fill only the NA entries in the original Cytogenetic_Risk
    Cytogenetic_Risk = coalesce(Cytogenetic_Risk, Cytogenetic_Risk_calc),
    
    # 3. (re)factor with the two levels
    Cytogenetic_Risk = factor(
      Cytogenetic_Risk,
      levels = c("Standard Risk", "High Risk")
    )
  ) %>%
  select(-Cytogenetic_Risk_calc)   # drop the helper column if you like

### Export deterministic Table 1 source-data companions
# These CSV/TSV files are not the formatted manuscript table. They are a
# reproducibility audit trail for the categorical counts and denominators that
# feed the gtsummary Table 1 object below. Keeping them in the
# original Table 1 script makes the final_manuscript_objects source-data folder
# reproducible without relying on the separate reproducible_workflow generator.
format_table1_count_percent <- function(n, denominator, missing_level = FALSE) {
  if (isTRUE(missing_level)) {
    return(as.character(ifelse(is.na(n), 0L, n)))
  }
  if (is.na(n) || is.na(denominator) || denominator == 0) {
    return("0 (0%)")
  }
  pct <- 100 * n / denominator
  digits <- ifelse(pct > 0 && pct < 10, 1L, 0L)
  paste0(n, " (", formatC(pct, format = "f", digits = digits), "%)")
}

table1_labels <- c(
  Gender = "Gender",
  AGE_GROUP = "Age Group",
  ISS_STAGE = "ISS Stage",
  Cytogenetic_Risk = "Cytogenetic Risk",
  Subtype = "Myeloma Ig Subtype",
  ECOG_SCORE = "ECOG Performance Status"
)
table1_cohort_levels <- c("Training", "Testing")

table1_source_rows <- list()
for (variable in vars_cat) {
  values <- as.character(dat_base[[variable]])
  display_values <- ifelse(is.na(values), "(Missing)", values)
  levels_to_show <- unique(display_values)
  levels_to_show <- levels_to_show[order(levels_to_show == "(Missing)", levels_to_show)]

  for (level in levels_to_show) {
    source_row <- data.frame(
      variable = unname(table1_labels[[variable]]),
      level = level,
      stringsAsFactors = FALSE
    )
    for (cohort_level in table1_cohort_levels) {
      in_cohort <- dat_base$cohort == cohort_level
      n_level <- sum(in_cohort & display_values == level, na.rm = TRUE)
      denominator <- sum(in_cohort & !is.na(values), na.rm = TRUE)
      source_row[[cohort_level]] <- format_table1_count_percent(
        n_level,
        denominator,
        missing_level = identical(level, "(Missing)")
      )
    }
    source_row[["Total"]] <- format_table1_count_percent(
      sum(display_values == level, na.rm = TRUE),
      sum(!is.na(values)),
      missing_level = identical(level, "(Missing)")
    )
    table1_source_rows[[length(table1_source_rows) + 1]] <- source_row
  }
}

table1_source_counts <- bind_rows(table1_source_rows)
dir.create("Output_tables_2025", recursive = TRUE, showWarnings = FALSE)
table1_source_counts_path <- file.path(
  "Output_tables_2025",
  "Table_1_clinical_demographics_computed_source_counts.csv"
)
write.csv(table1_source_counts, table1_source_counts_path, row.names = FALSE, quote = TRUE)

table1_source_qc <- data.frame(
  input_rds = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds",
  baseline_row_selection = "original Table 1 raw baseline/diagnosis selection with CA-02, CA-03, SPORE_0009, and duplicate-ranking logic retained",
  cohort_source = "current full Table 1 cohort labels from baseline_high_quality_patients_updated.rds",
  baseline_high_quality_rds = baseline_high_quality_path,
  narrow_cohort_csv_for_audit = patient_cohort_assignment_path,
  n_narrow_cohort_csv_patients = length(unique(patient_cohort_tbl_csv$Patient)),
  n_patients_added_beyond_narrow_csv = length(setdiff(dat_base$Patient, patient_cohort_tbl_csv$Patient)),
  n_baseline_table_patients = nrow(dat_base),
  n_training = sum(dat_base$cohort == "Training", na.rm = TRUE),
  n_testing = sum(dat_base$cohort == "Testing", na.rm = TRUE),
  n_missing_cohort = sum(is.na(dat_base$cohort)),
  output_csv = table1_source_counts_path,
  stringsAsFactors = FALSE
)
table1_source_qc_path <- file.path(
  "Output_tables_2025",
  "Table_1_clinical_demographics_computed_source_qc.tsv"
)
write.table(
  table1_source_qc,
  file = table1_source_qc_path,
  sep = "\t",
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

ms_copy_artifact(
  source_path = table1_source_counts_path,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_source_counts_csv",
  description = "Table 1 computed categorical count/percent source data used to build the gtsummary table.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

ms_copy_artifact(
  source_path = table1_source_qc_path,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_source_qc_tsv",
  description = "Table 1 source-data QC summary with input paths and cohort denominators.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

### Build and export categorical manuscript Table 1 source
tbl1_cat <- dat_base %>%
  select(all_of(c(vars_cat, "cohort"))) %>%
  tbl_summary(
    by           = cohort,
    missing      = "ifany",                # show NA counts
    missing_text = "(Missing)",
    statistic    = all_categorical() ~ "{n} ({p}%)",
    label        = list(
      Gender            ~ "Gender",
      AGE_GROUP         ~ "Age Group",
      ISS_STAGE         ~ "ISS Stage",
      Cytogenetic_Risk  ~ "Cytogenetic Risk",
      Subtype           ~ "Myeloma Ig Subtype",
      ECOG_SCORE        ~ "ECOG Performance Status"
    )
  ) %>%
  add_overall(last = TRUE, col_label = "Total") %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

# Export the scripted DOCX source. After Spring 2026 test-cohort integration,
# the retained DA-edited DOCX/PDF can be stale relative to the regenerated
# cohort counts, so this script also refreshes the mapped final DOCX/PDF from
# this deterministic source table before staging manuscript artifacts.
tbl1_flex <- as_flex_table(tbl1_cat)
doc <- read_docx() %>%
  body_add_flextable(tbl1_flex) %>%
  body_end_section_portrait()
table1_generated_docx <- "table1_clinical_categorical_updated_final_3.docx"
print(doc, target = table1_generated_docx)

write_table1_pdf_from_source_counts <- function(table1_source_counts, table1_source_qc, pdf_path) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    warning("Package gridExtra is not available; cannot write fallback Table 1 PDF: ", pdf_path)
    return(FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    warning("Package grid is not available; cannot write fallback Table 1 PDF: ", pdf_path)
    return(FALSE)
  }

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  pdf_device(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::grid.text(
    "Table 1. Baseline clinical and demographic characteristics",
    x = 0.03,
    y = 0.965,
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
  grid::grid.text(
    paste0(
      "Current cohort: n = ", table1_source_qc$n_baseline_table_patients,
      " (Training n = ", table1_source_qc$n_training,
      "; Testing n = ", table1_source_qc$n_testing, ")."
    ),
    x = 0.03,
    y = 0.925,
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9)
  )

  table_grob <- gridExtra::tableGrob(
    table1_source_counts,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 7,
      padding = grid::unit(c(3, 3), "pt"),
      core = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  grid::grid.draw(gridExtra::arrangeGrob(table_grob, vp = grid::viewport(width = 0.98, height = 0.84)))
  TRUE
}

render_table1_pdf <- function(docx_path, pdf_path, table1_source_counts, table1_source_qc) {
  soffice <- Sys.which("soffice")
  if (!nzchar(soffice)) {
    warning(
      "LibreOffice/soffice was not found; cannot refresh Table 1 PDF from current DOCX: ",
      pdf_path
    )
    return(write_table1_pdf_from_source_counts(table1_source_counts, table1_source_qc, pdf_path))
  }

  temp_dir <- tempfile("table1_pdf_render_")
  profile_dir <- tempfile("table1_lo_profile_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(c(temp_dir, profile_dir), recursive = TRUE, force = TRUE), add = TRUE)

  status <- system2(
    soffice,
    args = c(
      paste0("-env:UserInstallation=file://", normalizePath(profile_dir, mustWork = TRUE)),
      "--headless",
      "--convert-to",
      "pdf",
      "--outdir",
      normalizePath(temp_dir, mustWork = TRUE),
      normalizePath(docx_path, mustWork = TRUE)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  converted_pdf <- file.path(
    temp_dir,
    paste0(tools::file_path_sans_ext(basename(docx_path)), ".pdf")
  )
  if (!file.exists(converted_pdf)) {
    warning(
      "LibreOffice did not produce the expected Table 1 PDF: ",
      converted_pdf,
      "\nCommand output:\n",
      paste(status, collapse = "\n")
    )
    return(write_table1_pdf_from_source_counts(table1_source_counts, table1_source_qc, pdf_path))
  }

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(converted_pdf, pdf_path, overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy refreshed Table 1 PDF to: ", pdf_path, call. = FALSE)
  }
  TRUE
}

# MANUSCRIPT OUTPUT: Table 1
# This script regenerates the categorical Table 1 source DOCX above and promotes
# it to the mapped final manuscript artifact paths. This keeps the manuscript
# bundle synchronized with the current revision-inclusive cohort assignment
# instead of preserving a stale manual Word/PDF export from an older cohort.
table1_final_docx <- file.path(
  "Final Tables and Figures",
  "table1_clinical_categorical_updated_final_3_DA_edited.docx"
)
table1_final_pdf <- file.path(
  "Final Tables and Figures",
  "table1_clinical_categorical_updated_final_3_DA_edited.pdf"
)
table1_renamed_docx <- file.path(
  "Figures_Exported",
  "Tables_exported",
  "Renamed",
  "Table_1_Baseline_Clinical_Demographics.docx"
)
table1_renamed_pdf <- file.path(
  "Figures_Exported",
  "Tables_exported",
  "Renamed",
  "Table_1_Baseline_Clinical_Demographics.pdf"
)
table1_old_name_docx <- file.path(
  "Figures_Exported",
  "Tables_exported",
  "Old name",
  "table1_clinical_categorical_updated_final_3_DA_edited.docx"
)
table1_old_name_pdf <- file.path(
  "Figures_Exported",
  "Tables_exported",
  "Old name",
  "table1_clinical_categorical_updated_final_3_DA_edited.pdf"
)

dir.create(dirname(table1_final_docx), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(table1_generated_docx, table1_final_docx, overwrite = TRUE)) {
  stop("Failed to refresh final Table 1 DOCX from regenerated source: ", table1_final_docx, call. = FALSE)
}
render_table1_pdf(table1_final_docx, table1_final_pdf, table1_source_counts, table1_source_qc)

for (docx_path in c(table1_renamed_docx, table1_old_name_docx)) {
  dir.create(dirname(docx_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(table1_final_docx, docx_path, overwrite = TRUE)) {
    stop("Failed to refresh mapped Table 1 DOCX path: ", docx_path, call. = FALSE)
  }
}
if (file.exists(table1_final_pdf)) {
  for (pdf_path in c(table1_renamed_pdf, table1_old_name_pdf)) {
    dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(table1_final_pdf, pdf_path, overwrite = TRUE)) {
      stop("Failed to refresh mapped Table 1 PDF path: ", pdf_path, call. = FALSE)
    }
  }
}

ms_copy_artifact(
  source_path = table1_final_docx,
  artifact_id = "TABLE1_DOCX",
  role = "main_table_docx",
  description = "Table 1: clinical and demographic baseline characteristics regenerated from the current revision-inclusive cohort.",
  script_name = "2_1_Clinical_Demographics_Table.R"
)

if (file.exists(table1_final_pdf)) {
  ms_copy_artifact(
    source_path = table1_final_pdf,
    artifact_id = "TABLE1_PDF",
    role = "main_table_pdf",
    description = "Table 1: PDF rendered from the current regenerated revision-inclusive DOCX table.",
    script_name = "2_1_Clinical_Demographics_Table.R"
  )
} else {
  warning(
    "Final edited Table 1 PDF was not found. The PDF is a manuscript/export artifact and is not regenerated by this R script: ",
    table1_final_pdf
  )
}

### Optional continuous Table 1 companion
# This block is retained from the original script for traceability, but it is
# not the mapped final manuscript Table 1. Set the flag below to TRUE only when
# you intentionally want to regenerate the historical continuous-variable helper
# table (`baseline_characteristics_updated.docx`).
export_exploratory_continuous_table1 <- FALSE

if (isTRUE(export_exploratory_continuous_table1)) {
# ------------------------- 2.1  Convert continuous vars to numeric
vars_cont <- c(
  "AGE", "dFLC", "Albumin", "B2_micro", "Calcium",
  "Creatinine", "Hemoglobin", "LDH", "Plasma_pct"
)

dat_base <- dat_base %>%
  mutate(across(all_of(vars_cont), as.numeric))

# ------------------------- 2.2  Make categorical vars factors
vars_cat <- c(
  "Gender", "AGE_GROUP", "ISS_STAGE", "R_ISS_STAGE",
  "Cytogenetic_Risk", "Subtype", "ECOG_SCORE", "KPS_SCORE"
)

dat_base <- dat_base %>%
  mutate(across(all_of(vars_cat), ~ factor(.x)))


# -----------------------------------------------------------
# 4.  BUILD TABLE 1  ---------------------------------------
# -----------------------------------------------------------
tbl1 <- dat_base %>%
  select(all_of(c(vars_cat, vars_cont, "cohort"))) %>%
  tbl_summary(
    by        = cohort,
    type      = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    digits    = all_continuous() ~ 1,
    missing   = "always",                 # <‑‑ show missing counts explicitly
    label     = list(
      AGE        ~ "Age (years)",
      dFLC       ~ "dFLC (mg/L)",
      B2_micro   ~ "β₂‑microglobulin (mg/L)",
      Plasma_pct ~ "BM plasma cells (%)"
    )
  ) %>%
  add_overall(last = TRUE, col_label = "Total") %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

# -----------------------------------------------------------
# 5.  EXPORT TO WORD  --------------------------------------
# -----------------------------------------------------------
# Convert to flextable
tbl1_flex <- as_flex_table(tbl1)

# Export using flextable functions
doc <- read_docx() %>%
  body_add_flextable(tbl1_flex) %>%
  body_end_section_portrait()

print(doc, target = "baseline_characteristics_updated.docx")
}



### Exploratory fragmentomics/clinical MRD audit
# This optional audit does not create a manuscript figure or table. It lists
# possible additional non-IMG/non-SPORE patients with fragmentomics plus
# clinical MRD data who are not in the current cohort assignment file. Keep
# FALSE for the manuscript pipeline; set TRUE during cohort-expansion review.
run_fragmentomics_clinical_mrd_audit <- FALSE

if (isTRUE(run_fragmentomics_clinical_mrd_audit)) {

# Define which columns count as “clinical MRD”
clinical_mrd_cols <- c(
  "MRD_by_clinical_testing",
  "MRD_by_clinical_testing_stringent",
  "MRD_Clinical_Binary",
  "MRD_Clinical_Stringent_Binary"
)

patients_with_FS_and_clinical_MRD <- dat %>%
  # 1) require FS not missing
  filter(!is.na(FS)) %>%
  # 2) require at least one clinical‐MRD column not missing
  filter(if_any(all_of(clinical_mrd_cols), ~ !is.na(.))) %>%
  # 3) Check not at baseline 
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  # 4) pull unique patient IDs
  distinct(Patient) %>%
  pull(Patient)

# Print the result
cat("Patients with FS and clinical MRD data:\n")
print(patients_with_FS_and_clinical_MRD)
cat("Total:", length(patients_with_FS_and_clinical_MRD), "patients\n")

# Get patients that do not start with IMG or SPORE
non_img_spore_patients <- patients_with_FS_and_clinical_MRD[
  !grepl("^(IMG|SPORE)", patients_with_FS_and_clinical_MRD)
]

# Which of these are NOT in the cohort_df Patient list?
new_possible_patients <- setdiff(non_img_spore_patients, cohort_df$Patient)

# Show them
cat("Patients with FS+clinical MRD, not in cohort_df, and not IMG/SPORE:\n")
print(new_possible_patients)
cat("Total:", length(new_possible_patients), "patients\n")
}


### Archived idea: excluded from the manuscript workflow
# The original working script had considered adding more complexity after this
# point, but that path is excluded from all final manuscript figure/table outputs.
