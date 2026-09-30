# =============================================================================
# 2_1_Part2_Cohort_Swim_Plot.R
#
# Purpose:
#   Build a cohort swim plot showing treatment timelines
#   (induction, transplant, maintenance, progression/relapse) for the manuscript
#   cohort drawn from M4, SPORE, and IMMAGINE. Each row is one patient;
#   horizontal bars represent treatment lines and vertical markers indicate
#   key clinical events (e.g. ASCT, relapse). Saved as a high-resolution PNG
#   for manuscript Figure 1 / supplementary.
#
# Inputs:
#   - Clinical data/SPORE/tidy_treatments.csv
#   - M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx
#   - Clinical data/IMMAGINE/Cleaned_IMMAGINE_chemotherapy.csv
#   - combined_clinical_data_updated_April2025.csv  (for relapse/event dates)
#   - cohort_assignment_table_updated.rds            (M4/SPORE/IMMAGINE label)
#   - Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6_full.rds
#       (required to add eligible Spring 2026 longitudinal patients)
#   - id_map.rds (fixed patient-to-index mapping used before this script later
#       recreates the map; a clean run currently requires the existing file)
#
# Outputs:
#   - Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated10A.png
#       (Figure 1A plotted component)
#   - Final Tables and Figures/
#       Supp_Table_1_all_events_for_swim_plot_INDEX_DATES_privacy_protected.csv
#       (generated Supplementary Table 1 candidate)
#   - Supporting event tables and review plots described in the labelled
#       export sections below
#
# Dependencies:
#   tidyverse, readxl, lubridate, patchwork, forcats, purrr
#   This script does not depend on the Table 1 output from
#   2_1_Clinical_Demographics_Table.R; it requires the clinical, cohort,
#   MRDetect-call, and fixed ID-map inputs listed above.
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_1_Part2_Cohort_Swim_Plot.R
#
# Manuscript outputs created/updated:
#   - Figure 1A: cohort treatment/sample-timing swim plot component.
#   - Supplementary Table 1: patient-level clinical event/treatment timeline
#     table supporting the swim plot and cohort description.
#
# Pipeline role:
#   This script harmonizes treatment and event dates across M4, SPORE, and
#   IMMAGINE source files, then converts them into a patient-by-time display.
#   The output is the traceable plot component used in Figure 1A plus the
#   accompanying timeline table.
#
# Author:    Dory Abelman
# Last update: May 2025
# =============================================================================
#### Now make swim plot 
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)
library(forcats)
library(purrr)

# Support-only review outputs from this script are written here so that the
# project root stays clean and manuscript outputs remain easy to identify.
swim_support_dir <- file.path("Final Tables and Figures", "swim_plot_support")
dir.create(swim_support_dir, recursive = TRUE, showWarnings = FALSE)

# Shared helper for final manuscript-organized outputs.
# This keeps the scientific code in this script, while also copying the final
# figure/table components into Scripts_2025/Final_Scripts/final_manuscript_objects
# with manuscript labels such as Figure_1A and Supplementary_Table_1.
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

empty_swim_events <- function() {
  tibble(
    patient = character(),
    event = character(),
    start = as.Date(character()),
    end = as.Date(character()),
    details = character()
  )
}

load_swim_mrdetect_calls <- function(path = "Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6_full.rds") {
  # The swim plot only adds Spring 2026 patients when processed MRDetect calls
  # are available. Missing this table is a warning rather than a hard failure so
  # older figure builds can still run without silently adding unscored patients.
  if (!file.exists(path)) {
    warning("Missing MRDetect call table used to gate revision swim-plot patients: ", path, call. = FALSE)
    return(NULL)
  }
  readRDS(path)
}

eligible_spring2026_swim_patients <- function() {
  # ## Decide which Spring 2026 patients belong in the swim plot
  # Revision metadata can contain many submitted samples, but the manuscript swim
  # plot should only add patients who have analyzable baseline BM/blood BAMs and
  # at least one longitudinal MRDetect call. This keeps the revision set as a
  # scored test-cohort addition rather than a broad submission inventory.
  revision <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(revision)) return(character())

  calls <- load_swim_mrdetect_calls()
  if (is.null(calls)) return(character())

  require_columns(
    revision,
    c("Patient", "Sample_type", "timepoint_info", "has_bam"),
    "Spring 2026 revision metadata for swim plot"
  )
  require_columns(
    calls,
    c("Patient", "timepoint_info"),
    "MRDetect call table for swim plot"
  )

  baseline_bam_patients <- revision %>%
    filter(
      # Baseline/Diagnosis BAMs establish the patient's molecular context for the
      # longitudinal row. Other revision timepoints alone are not sufficient.
      .data$timepoint_info %in% c("Diagnosis", "Baseline"),
      .data$Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
      .data$has_bam %in% TRUE
    ) %>%
    distinct(Patient) %>%
    pull(Patient)

  call_cols <- intersect(
    # The gate accepts either BM-derived or blood-derived MRDetect call columns.
    # Column names vary across exports, so only available recognized columns are
    # tested.
    c(
      "BM_zscore_only_sites_call",
      "BM_zscore_only_detection_rate_call",
      "Blood_zscore_only_sites_call",
      "Blood_zscore_only_detection_rate_call"
    ),
    names(calls)
  )
  if (!length(call_cols)) {
    warning("MRDetect call table has no recognized BM/Blood call columns for revision swim-plot gating.", call. = FALSE)
    return(character())
  }

  longitudinal_mrdetect_patients <- calls %>%
    filter(
      .data$Patient %in% revision$Patient,
      # Require a non-baseline timepoint with at least one call value. Baseline
      # calls alone do not demonstrate longitudinal monitoring for the swim plot.
      !.data$timepoint_info %in% c("Diagnosis", "Baseline")
    ) %>%
    filter(if_any(all_of(call_cols), ~ !is.na(.x))) %>%
    distinct(Patient) %>%
    pull(Patient)

  eligible <- sort(intersect(baseline_bam_patients, longitudinal_mrdetect_patients))
  message(
    "Spring 2026 swim-plot test cohort: ",
    length(eligible), " patients with baseline BM/blood BAM and longitudinal MRDetect calls."
  )
  eligible
}

load_swim_cohort_assignment <- function(path = "cohort_assignment_table_updated.rds") {
  # Start from the historical cohort assignment and append eligible revision
  # patients as Non-frontline. Existing cohort labels are not overwritten.
  cohort_df <- readRDS(path)
  require_columns(cohort_df, c("Patient", "Cohort"), "Swim-plot cohort assignment")

  revision_patients <- eligible_spring2026_swim_patients()
  if (!length(revision_patients)) {
    return(cohort_df %>% distinct(.data$Patient, .keep_all = TRUE))
  }

  revision_cohort <- tibble(
    Patient = revision_patients,
    Cohort = "Non-frontline"
  )

  cohort_df %>%
    mutate(Patient = as.character(.data$Patient)) %>%
    bind_rows(revision_cohort %>% filter(!.data$Patient %in% cohort_df$Patient)) %>%
    distinct(.data$Patient, .keep_all = TRUE)
}

build_spring2026_revision_swim_events <- function(eligible_patients) {
  # ## Convert Spring 2026 metadata into swim-plot event rows
  # Events are intentionally conservative: baseline is anchored to the earliest
  # baseline/diagnosis BAM date, endpoint events use progression date when the
  # patient progressed otherwise censor date, and clinical MRD events come from
  # explicit mrd_test_date rows.
  revision <- load_spring2026_revision_metadata(required = FALSE)
  if (is.null(revision)) return(empty_swim_events())

  revision <- revision %>%
    filter(.data$Patient %in% eligible_patients) %>%
    mutate(
      Date_of_sample_collection = parse_date_safely(.data$Date_of_sample_collection),
      first_progression_date = parse_date_safely(.data$first_progression_date),
      relapse_or_censor_date = parse_date_safely(.data$relapse_or_censor_date),
      mrd_test_date = parse_date_safely(.data$mrd_test_date)
    )

  if (!nrow(revision)) return(empty_swim_events())

  baseline_events <- revision %>%
    filter(
      # Use BAM-backed baseline/diagnosis BM or blood rows only; submitted rows
      # without BAMs are not plotted as molecular baseline events.
      .data$timepoint_info %in% c("Diagnosis", "Baseline"),
      .data$Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
      .data$has_bam %in% TRUE,
      !is.na(.data$Date_of_sample_collection)
    ) %>%
    group_by(Patient) %>%
    summarize(start = min(.data$Date_of_sample_collection), .groups = "drop") %>%
    transmute(
      patient = .data$Patient,
      event = "Baseline",
      start = .data$start,
      end = .data$start,
      details = "Spring 2026 baseline BM/blood BAM"
    )

  endpoint_events <- revision %>%
    distinct(
      Patient,
      first_progression_date,
      relapse_or_censor_date,
      relapse_or_censor_status
    ) %>%
    mutate(
      endpoint_date = case_when(
        # If progression occurred and a first progression date is present, that
        # is the event date. Otherwise use the relapse/censor date from metadata.
        .data$relapse_or_censor_status == "relapse/progression" & !is.na(.data$first_progression_date) ~ .data$first_progression_date,
        TRUE ~ .data$relapse_or_censor_date
      ),
      event = if_else(
        .data$relapse_or_censor_status == "censored_last_followup",
        "Last follow-up",
        "Relapse"
      ),
      details = if_else(
        .data$relapse_or_censor_status == "censored_last_followup",
        "Relapsed: 0",
        "Relapsed: 1"
      )
    ) %>%
    filter(!is.na(.data$endpoint_date)) %>%
    transmute(
      patient = .data$Patient,
      event = .data$event,
      start = .data$endpoint_date,
      end = .data$endpoint_date,
      details = .data$details
    )

  clinical_mrd_events <- revision %>%
    # Clinical MRD is a separate event type from cfWGS MRDetect calls. It is
    # plotted only when the metadata provides an explicit test date.
    distinct(Patient, mrd_test_date, mrd_result_comprehensive) %>%
    filter(!is.na(.data$mrd_test_date)) %>%
    transmute(
      patient = .data$Patient,
      event = "MRD (clinical)",
      start = .data$mrd_test_date,
      end = .data$mrd_test_date,
      details = paste0("Clinical MRD result: ", .data$mrd_result_comprehensive)
    )

  bind_rows(baseline_events, endpoint_events, clinical_mrd_events) %>%
    distinct(.data$patient, .data$event, .data$start, .data$end, .data$details)
}

build_spring2026_revision_treatment_events <- function(
    eligible_patients,
    clinical_dir = file.path("New OICR Submissions", "derived_metadata"),
    support_dir = swim_support_dir
) {
  # ## Build OICR revision treatment bars for swim plots
  # This function is restricted to patients already eligible for the Spring 2026
  # swim-plot test cohort. It prefers directly dated treatment rows. Relative
  # treatment timelines are used only for patients without dated treatment rows
  # and require a usable diagnosis date as the time origin.
  eligible_patients <- sort(intersect(as.character(eligible_patients), eligible_spring2026_swim_patients()))
  if (!length(eligible_patients)) return(empty_swim_events())

  summary_path <- file.path(clinical_dir, "oicr_submission_patient_clinical_summary.csv")
  treatment_path <- file.path(clinical_dir, "oicr_submission_clinical_treatment_rows.csv")
  timeline_path <- file.path(clinical_dir, "oicr_submission_clinical_treatment_timeline_rows.csv")
  esther_treatment_path <- file.path(
    "Clinical data",
    "Additional sample sheets from Esther",
    "Chem_LiberateID_19Feb2026.xlsx"
  )

  required_paths <- c(summary_path, treatment_path, timeline_path)
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths)) {
    warning(
      "Missing OICR revision treatment source file(s): ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
    return(empty_swim_events())
  }

  clean_date <- function(x) {
    # Treat Excel-origin sentinel/invalid dates at or before 1900-01-01 as
    # missing rather than plotting them as real clinical events.
    x_chr <- as.character(x)
    parsed <- parse_date_safely(x_chr)
    excel_serial_idx <- is.na(parsed) &
      !is.na(x_chr) &
      str_detect(x_chr, "^\\d+(\\.0+)?$")
    if (any(excel_serial_idx)) {
      serial <- suppressWarnings(as.numeric(x_chr[excel_serial_idx]))
      valid_serial <- !is.na(serial) & serial >= 20000 & serial <= 60000
      excel_dates <- rep(as.Date(NA), length(serial))
      excel_dates[valid_serial] <- as.Date(serial[valid_serial], origin = "1899-12-30")
      parsed[which(excel_serial_idx)] <- excel_dates
    }
    parsed[!is.na(parsed) & parsed <= as.Date("1900-01-01")] <- NA
    parsed
  }

  is_transplant_regimen <- function(x) {
    # Transplants are plotted as point events separate from chemotherapy bars.
    str_detect(
      coalesce(as.character(x), ""),
      regex("\\bASCT\\b|\\bBMT\\b|transplant", ignore_case = TRUE)
    )
  }

  patient_summary <- read_csv(summary_path, col_types = cols(.default = "c"))
  dated_treatments <- read_csv(treatment_path, col_types = cols(.default = "c"))
  timeline_treatments <- read_csv(timeline_path, col_types = cols(.default = "c"))
  revision_metadata <- load_spring2026_revision_metadata(required = FALSE)

  require_columns(
    patient_summary,
    c("patient_img_id", "patient_numeric_id", "study_name", "date_diagnosis", "date_of_last_followup"),
    "OICR revision patient clinical summary"
  )
  require_columns(
    dated_treatments,
    c("patient_numeric_id", "study_name", "line_regimen", "line_of_treatment",
      "treatment_start_date", "treatment_end_date", "treatment_intent"),
    "OICR revision dated treatment rows"
  )
  require_columns(
    timeline_treatments,
    c("patient_img_id", "treatment_start_days", "treatment_stop_days",
      "treatment_name", "treatment_line", "treatment_category", "treatment_setting"),
    "OICR revision relative treatment timeline rows"
  )

  patient_key <- patient_summary %>%
    # Dated treatment rows use numeric patient IDs plus study name, while the
    # rest of the workflow uses patient_img_id. Build a bridge between the two.
    transmute(
      patient_img_id = as.character(.data$patient_img_id),
      patient_numeric_id = as.character(.data$patient_numeric_id),
      study_name = as.character(.data$study_name)
    ) %>%
    filter(!is.na(.data$patient_img_id)) %>%
    distinct(.data$patient_numeric_id, .data$study_name, .keep_all = TRUE)

  if (!is.null(revision_metadata)) {
    require_columns(
      revision_metadata,
      c("Patient", "patient_numeric_id", "Study"),
      "Spring 2026 revision metadata patient ID bridge"
    )
    revision_patient_key <- revision_metadata %>%
      transmute(
        patient_img_id = as.character(.data$Patient),
        patient_numeric_id = as.character(.data$patient_numeric_id),
        study_name = as.character(.data$Study)
      ) %>%
      filter(!is.na(.data$patient_img_id), !is.na(.data$patient_numeric_id))

    patient_key <- bind_rows(patient_key, revision_patient_key) %>%
      distinct(.data$patient_numeric_id, .data$study_name, .keep_all = TRUE)
  }

  patient_key_by_id <- patient_key %>%
    filter(!is.na(.data$patient_numeric_id), .data$patient_img_id %in% eligible_patients) %>%
    group_by(.data$patient_numeric_id) %>%
    filter(n_distinct(.data$patient_img_id) == 1L) %>%
    summarize(patient_img_id_by_id = first(.data$patient_img_id), .groups = "drop")

  primary_dated_patient_ids <- dated_treatments %>%
    mutate(patient_numeric_id = as.character(.data$patient_numeric_id)) %>%
    left_join(patient_key, by = c("patient_numeric_id", "study_name")) %>%
    filter(.data$patient_img_id %in% eligible_patients) %>%
    distinct(.data$patient_numeric_id) %>%
    pull(.data$patient_numeric_id)

  esther_treatments <- tibble(
    patient_numeric_id = character(),
    study_name = character(),
    line_regimen = character(),
    line_of_treatment = character(),
    treatment_start_date = character(),
    treatment_end_date = character(),
    treatment_study_drug = character(),
    treatment_intent = character(),
    best_response_date = character(),
    best_response = character(),
    progression = character(),
    progression_date = character(),
    progression_details = character(),
    treatment_source_file = character()
  )
  if (file.exists(esther_treatment_path)) {
    esther_raw <- read_excel(esther_treatment_path, col_types = "text")
    require_columns(
      esther_raw,
      c("PATIENT_ID", "REGIMEN_NAME", "LINE_OF_TREATMENT", "START_DATE",
        "END_DATE", "STUDY_DRUG", "INTENT", "BEST_RESPONSE", "PROGRESSION",
        "PROGRESSION_DATE"),
      "Esther IMMAGINE chemotherapy workbook"
    )
    esther_treatments <- esther_raw %>%
      transmute(
        patient_numeric_id = as.character(.data$PATIENT_ID),
        study_name = NA_character_,
        line_regimen = as.character(.data$REGIMEN_NAME),
        line_of_treatment = as.character(.data$LINE_OF_TREATMENT),
        treatment_start_date = as.character(.data$START_DATE),
        treatment_end_date = as.character(.data$END_DATE),
        treatment_study_drug = as.character(.data$STUDY_DRUG),
        treatment_intent = as.character(.data$INTENT),
        best_response_date = as.character(.data$BEST_RESPONSE_DATE),
        best_response = as.character(.data$BEST_RESPONSE),
        progression = as.character(.data$PROGRESSION),
        progression_date = as.character(.data$PROGRESSION_DATE),
        progression_details = NA_character_,
        treatment_source_file = "Chem_LiberateID_19Feb2026.xlsx"
      ) %>%
      filter(
        .data$patient_numeric_id %in% patient_key$patient_numeric_id,
        !.data$patient_numeric_id %in% primary_dated_patient_ids,
        !is.na(.data$line_regimen),
        nzchar(.data$line_regimen)
      )
  }

  dated_treatments <- dated_treatments %>%
    mutate(
      treatment_source_file = "oicr_submission_clinical_treatment_rows.csv",
      best_response_date = NA_character_
    ) %>%
    bind_rows(esther_treatments)

  date_index <- patient_summary %>%
    # Diagnosis anchors relative treatment timelines; last follow-up can be used
    # as an end date when a relative stop date is unavailable.
    transmute(
      patient_img_id = as.character(.data$patient_img_id),
      date_diagnosis = clean_date(.data$date_diagnosis),
      date_of_last_followup = clean_date(.data$date_of_last_followup)
    )

  if (!is.null(revision_metadata)) {
    # Add diagnosis/follow-up dates from the curated revision metadata as a
    # fallback date index. This covers patients whose clinical summary rows are
    # incomplete but whose repo-style metadata has endpoint dates.
    require_columns(
      revision_metadata,
      c("Patient", "date_diagnosis", "relapse_or_censor_date"),
      "Spring 2026 revision metadata treatment date index"
    )
    revision_dates <- revision_metadata %>%
      transmute(
        patient_img_id = as.character(.data$Patient),
        date_diagnosis = clean_date(.data$date_diagnosis),
        date_of_last_followup = clean_date(.data$relapse_or_censor_date)
      )

    date_index <- bind_rows(date_index, revision_dates)
  }

  date_index <- date_index %>%
    filter(.data$patient_img_id %in% eligible_patients) %>%
    group_by(.data$patient_img_id) %>%
    summarize(
      date_diagnosis = suppressWarnings(min(.data$date_diagnosis, na.rm = TRUE)),
      date_of_last_followup = suppressWarnings(max(.data$date_of_last_followup, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      date_diagnosis = if_else(is.infinite(.data$date_diagnosis), as.Date(NA), .data$date_diagnosis),
      date_of_last_followup = if_else(is.infinite(.data$date_of_last_followup), as.Date(NA), .data$date_of_last_followup)
    )

  dated_joined <- dated_treatments %>%
    # Directly dated treatment rows are the highest-confidence source because
    # start/end dates are already absolute calendar dates.
    mutate(
      patient_numeric_id = as.character(.data$patient_numeric_id),
      study_name = as.character(.data$study_name),
      treatment_start_date = clean_date(.data$treatment_start_date),
      treatment_end_date = clean_date(.data$treatment_end_date),
      best_response_date = clean_date(.data$best_response_date),
      progression_date = clean_date(.data$progression_date)
    ) %>%
    left_join(patient_key, by = c("patient_numeric_id", "study_name")) %>%
    left_join(patient_key_by_id, by = "patient_numeric_id") %>%
    mutate(patient_img_id = coalesce(.data$patient_img_id, .data$patient_img_id_by_id)) %>%
    filter(
      .data$patient_img_id %in% eligible_patients,
      !is.na(.data$treatment_start_date)
    ) %>%
    mutate(
      regimen = coalesce(na_if(.data$line_regimen, ""), "Unspecified OICR treatment"),
      line_label = if_else(
        is.na(.data$line_of_treatment) | !nzchar(.data$line_of_treatment),
        NA_character_,
        paste0("L", .data$line_of_treatment)
      ),
      details = str_squish(paste(coalesce(.data$line_label, ""), .data$regimen))
    )

  dated_bars <- dated_joined %>%
    filter(!is_transplant_regimen(.data$regimen)) %>%
    transmute(
      patient = .data$patient_img_id,
      event = "Chemotherapy",
      start = .data$treatment_start_date,
      end = coalesce(
        .data$treatment_end_date,
        .data$progression_date,
        .data$best_response_date,
        .data$treatment_start_date
      ),
      details = .data$details,
      provenance = .data$treatment_source_file
    )

  dated_transplants <- dated_joined %>%
    filter(is_transplant_regimen(.data$regimen)) %>%
    transmute(
      patient = .data$patient_img_id,
      event = "Transplant",
      start = .data$treatment_start_date,
      end = .data$treatment_start_date,
      details = .data$details,
      provenance = .data$treatment_source_file
    )

  dated_patients <- dated_joined %>%
    distinct(.data$patient_img_id) %>%
    pull(.data$patient_img_id)

  timeline_joined <- timeline_treatments %>%
    # Relative treatment rows are used only for patients lacking dated rows. The
    # start/stop day offsets are converted to calendar dates using diagnosis date.
    mutate(
      patient_img_id = as.character(.data$patient_img_id),
      treatment_start_days = suppressWarnings(as.numeric(.data$treatment_start_days)),
      treatment_stop_days = suppressWarnings(as.numeric(.data$treatment_stop_days)),
      treatment_name = coalesce(na_if(.data$treatment_name, ""), "Unspecified OICR treatment")
    ) %>%
    filter(
      .data$patient_img_id %in% setdiff(eligible_patients, dated_patients),
      !is.na(.data$treatment_start_days)
    ) %>%
    left_join(date_index, by = "patient_img_id") %>%
    filter(!is.na(.data$date_diagnosis)) %>%
    mutate(
      start = .data$date_diagnosis + days(.data$treatment_start_days),
      end = case_when(
        # Prefer explicit relative stop day; otherwise extend to last follow-up
        # when available, or use a point event at start as the least assumptive
        # fallback.
        !is.na(.data$treatment_stop_days) ~ .data$date_diagnosis + days(.data$treatment_stop_days),
        !is.na(.data$date_of_last_followup) ~ .data$date_of_last_followup,
        TRUE ~ .data$date_diagnosis + days(.data$treatment_start_days)
      ),
      details = str_squish(paste(
        coalesce(na_if(.data$treatment_line, ""), NA_character_),
        .data$treatment_name,
        sep = " "
      ))
    )

  timeline_bars <- timeline_joined %>%
    transmute(
      patient = .data$patient_img_id,
      event = "Chemotherapy",
      start = .data$start,
      end = pmax(.data$end, .data$start, na.rm = TRUE),
      details = .data$details,
      provenance = "oicr_submission_clinical_treatment_timeline_rows.csv"
    )

  timeline_transplants <- timeline_joined %>%
    filter(is_transplant_regimen(.data$treatment_name)) %>%
    transmute(
      patient = .data$patient_img_id,
      event = "Transplant",
      start = .data$start,
      end = .data$start,
      details = .data$details,
      provenance = "oicr_submission_clinical_treatment_timeline_rows.csv"
    )

  treatment_events <- bind_rows(
    dated_bars,
    dated_transplants,
    timeline_bars,
    timeline_transplants
  ) %>%
    filter(!is.na(.data$start)) %>%
    distinct(.data$patient, .data$event, .data$start, .data$end, .data$details, .keep_all = TRUE) %>%
    arrange(.data$patient, .data$start, .data$event)

  coverage_audit <- tibble(patient = sort(unique(eligible_patients))) %>%
    # One-row-per-patient audit showing which treatment source contributed plot
    # events and which patients still lack structured treatment bars.
    left_join(
      dated_joined %>% count(.data$patient_img_id, name = "dated_treatment_rows"),
      by = c("patient" = "patient_img_id")
    ) %>%
    left_join(
      timeline_treatments %>%
        count(.data$patient_img_id, name = "relative_timeline_rows"),
      by = c("patient" = "patient_img_id")
    ) %>%
    left_join(
      treatment_events %>% count(.data$patient, name = "plot_events_added"),
      by = "patient"
    ) %>%
    mutate(
      across(
        c("dated_treatment_rows", "relative_timeline_rows", "plot_events_added"),
        ~ replace_na(.x, 0L)
      ),
      treatment_source_used = case_when(
        .data$dated_treatment_rows > 0 ~ "dated_treatment_rows",
        .data$plot_events_added > 0 & .data$relative_timeline_rows > 0 ~ "relative_timeline_rows_with_diagnosis_anchor",
        .data$relative_timeline_rows > 0 ~ "relative_timeline_rows_without_usable_anchor",
        TRUE ~ "no_structured_treatment_rows_found"
      )
    )

  write_csv(
    treatment_events,
    file.path(support_dir, "oicr_revision_treatment_events_used.csv")
  )
  write_csv(
    coverage_audit,
    file.path(support_dir, "oicr_revision_treatment_coverage_audit.csv")
  )

  missing_plot_treatment <- coverage_audit %>%
    filter(.data$plot_events_added == 0) %>%
    pull(.data$patient)
  if (length(missing_plot_treatment)) {
    warning(
      "No structured OICR treatment rows were available for revision swim-plot patients: ",
      paste(missing_plot_treatment, collapse = ", "),
      call. = FALSE
    )
  }

  treatment_events %>%
    select("patient", "event", "start", "end", "details")
}

### Load cohort assignments up front
# The cohort table is needed both for early event-table filtering and for the
# final Figure 1A annotation tracks.
cohort_df <- load_swim_cohort_assignment()

## --------------------------
## 1. Load primary treatment-event sources
## --------------------------

# 1a. “tidy_treatments.csv” you already have
tt_raw <- read_csv("Clinical data/SPORE/tidy_treatments.csv",
                   col_types = cols(.default = "c"))   # keep all as character; we’ll parse later



# 1b. M4 cohort chemotherapy Excel
m4_raw <- read_excel("M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx",
                     col_types = "text")                # read everything as character

# 1c. IMMAGINE chemotherapy CSV
imm_raw <- read_csv("Clinical data/IMMAGINE/Cleaned_IMMAGINE_chemotherapy.csv",
                    col_types = cols(.default = "c"))

## -------------------------------------------------
## 2. Tidy each data set into (patient, event, date)
## -------------------------------------------------

## First SPORE 

# ─────────────────────────────────────────────────────────────────────────────
# A. SPORE “tidy_treatments.csv”
# ─────────────────────────────────────────────────────────────────────────────
## 2a. tidy_treatments  (already long, just rename / parse), includes stransplant
tt <- tt_raw %>%
  transmute(
    patient = patient,                        # the CSV’s patient column
    event   = str_to_sentence(line),          # “line” → e.g. “Diagnosis”, “Transplant”
    start   = ymd(start_date),                # parse yyyy-mm-dd
    end     = NA,                          # single-day events, set as same
    details = regimen                         # whatever you’d like to show here
  )


# ─────────────────────────────────────────────────────────────────────────────
# B. M4 cohort 
# ─────────────────────────────────────────────────────────────────────────────
### Next M4
## 2b. M4 chemotherapy rows  → “Chemotherapy” events
m4_chemo <- m4_raw %>%
  mutate(
    # convert Excel serial → R Date
    start = as.Date(as.numeric(START_DATE), origin = "1899-12-30"),
    end   = as.Date(as.numeric(END_DATE),   origin = "1899-12-30"),
    details = REGIMEN_NAME
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Chemotherapy",
    start,
    end,
    details
  )

m4_chemo <- m4_chemo %>% filter(!is.na(start))
write.csv(
  m4_chemo %>% filter(patient %in% cohort_df$Patient) %>% filter(is.na(end)),
  file.path(swim_support_dir, "m4_chemo_missing_end_dates.csv"),
  row.names = FALSE
)

### Add transplant info for M4 
# Read the raw sheet
m4_trans_raw <- read_excel(
  "M4_CMRG_Data/M4_COHORT_STEM_CELL_TRANSPLANT.xlsx",
  col_types = "text"
)

# Robust date parser with NA‐guards
parse_m4_date <- function(x) {
  n <- length(x)
  out <- rep(as.Date(NA), n)
  
  # only non-NA entries
  not_na <- !is.na(x)
  
  # 2a. Excel serials: purely digits
  is_serial <- not_na & str_detect(x, "^[0-9]+$")
  out[is_serial] <- as.Date(as.numeric(x[is_serial]), origin = "1899-12-30")
  
  # 2b. Remaining non-blank, non-serial text
  is_text <- not_na & !is_serial & x != ""
  out[is_text] <- ymd(x[is_text])
  
  out
}

# 2c. Transplant events
m4_transplant_events <- m4_trans_raw %>%
  mutate(
    tx_date = parse_m4_date(TRANSPLANT_DATE),
    details = str_c(
      PROCEDURE_TYPE_TEXT,
      TRANSPLANT_TYPE,
      if_else(is.na(INJECTED_CD34) | INJECTED_CD34 == "",
              "",
              str_c("CD34:", INJECTED_CD34)),
      sep = " | "
    )
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Transplant",
    start   = tx_date,
    end     = tx_date,
    details
  )

# 2d. Check for any remaining failures in parsing
bad_dates <- m4_trans_raw %>%
  filter(!is.na(TRANSPLANT_DATE) & TRANSPLANT_DATE != "" &
           is.na(parse_m4_date(TRANSPLANT_DATE))) %>%
  pull(TRANSPLANT_DATE) %>%
  unique()

if (length(bad_dates)) {
  message("These raw TRANSPLANT_DATE values still failed to parse:\n",
          paste(bad_dates, collapse = ", "))
}



## 2e. Immune‐response event
m4_immun_resp <- m4_trans_raw %>%
  mutate(IMMUN_RESPONSE_DATE = as.numeric(IMMUN_RESPONSE_DATE)) %>%
  # keep only rows with a non‐blank IMMUN_RESPONSE_DATE
  filter(!is.na(IMMUN_RESPONSE_DATE) & IMMUN_RESPONSE_DATE != "") %>%
  mutate(
    # parse with the same helper
    resp_date = parse_m4_date(IMMUN_RESPONSE_DATE)
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Immune response",
    start   = resp_date,
    end     = resp_date,
    details = IMMUN_RESPONSE_TYPE
  )

## 2f. Best‐response event
m4_best_resp <- m4_trans_raw %>%
  mutate(BEST_RESPONSE_DATE = as.numeric(BEST_RESPONSE_DATE)) %>%
  filter(!is.na(BEST_RESPONSE_DATE) & BEST_RESPONSE_DATE != "") %>%
  mutate(
    best_date = parse_m4_date(BEST_RESPONSE_DATE)
  ) %>%
  transmute(
    patient = M4_id,
    event   = "Best response",
    start   = best_date,
    end     = best_date,
    details = BEST_RESPONSE
  )




# ─────────────────────────────────────────────────────────────────────────────
#  C. IMMAGINE 
# ─────────────────────────────────────────────────────────────────────────────
### Lastly IMMAGINE
## 2d. IMMAGINE rows  (already long but rename fields)
imm_clean <- imm_raw %>%
  rename(all = `Patient,Event,Date,Details`) %>%
  separate(
    col  = all,
    into = c("patient", "event", "date", "details"),
    sep  = ",",
    fill = "right"      # in case some rows end with a comma
  )

## Fix dates 
imm_clean <- imm_clean %>%
  # Capitalize event names if you like
  mutate(event = str_to_sentence(event)) %>%
  
  # Count dashes, pad to first-of-month or first-of-year
  mutate(
    n_dash   = str_count(date, fixed("-")),
    date_full = case_when(
      n_dash == 2 ~ date,                   # “YYYY-MM-DD”
      n_dash == 1 ~ paste0(date, "-01"),    # “YYYY-MM”    → “YYYY-MM-01”
      n_dash == 0 ~ paste0(date, "-01-01"), # “YYYY”       → “YYYY-01-01”
      TRUE        ~ date
    )
  ) %>%
  
  # Parse into real Date, set end = start
  mutate(
    start = ymd(date_full),
    end   = NA
  ) %>%
  
  # Drop helper columns
  select(patient, event, start, end, details)




# ─────────────────────────────────────────────────────────────────────────────
# D. Relapse, Baseline & Last follow‐up events
# ─────────────────────────────────────────────────────────────────────────────
#### Add diagnosis dates to this and the censor dates
# a) relapse dates (one row per patient, many relapse‐cols)
relapse_dates_full <- read_csv(
  "Relapse dates cfWGS updated.csv",
  col_types = cols(.default = "c")
)

# b) sample‐collection / censor dates per patient
censor_tbl <- readRDS("Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds")
censor_tbl$Baseline_Date <- censor_tbl$baseline_date # for consistency

# 2) Build a “Relapse” events table ----------------------------------------
relapse_events <- relapse_dates_full %>%
  # if the patient ID column is called something else, rename it:
  rename(patient = Patient) %>%  
  # pivot all the relapse columns (e.g. Relapse1, Relapse2, …) into long form:
  pivot_longer(
    cols      = -patient,
    names_to  = "which_relapse",
    values_to = "date"
  ) %>%
  filter(!is.na(date) & date != "") %>%            # drop missing
  mutate(
    start   = ymd(date),                           # parse to Date
    end     = start,
    event   = "Relapse",
    details = which_relapse                        # e.g. “Relapse1”
  ) %>%
  select(patient, event, start, end, details)


# 3) Build a baseline and date of last followup events table -----------------------------
# 1. Baseline events from the Baseline_Date column
baseline_events <- censor_tbl %>%
  transmute(
    patient = Patient,
    event   = "Baseline",
    # Baseline_Date is POSIXct; convert to Date
    start   = as_date(Baseline_Date),
    end     = start,
    details = NA_character_
  )

## Export this for use in future 
write.csv(baseline_events, file = "Final Tables and Figures/Baseline dates for samples.csv")

# 2. Last follow-up events from the Censor_date column
followup_events <- censor_tbl %>%
  transmute(
    patient = Patient,
    event   = "Last follow-up",
    start   = censor_date,
    end     = start,
    details = paste0("Relapsed: ", relapsed)
  )




# ─────────────────────────────────────────────────────────────────────────────
# E. Sample collections (BM vs cfDNA)
# ─────────────────────────────────────────────────────────────────────────────
### Now add all the sample collection dates 
.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

combined_clinical_data_updated <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)
revision_swim_events <- build_spring2026_revision_swim_events(cohort_df$Patient)
revision_treatment_events <- build_spring2026_revision_treatment_events(cohort_df$Patient)

# 1. BM sample collection events
bm_events <- combined_clinical_data_updated %>%
  filter(Sample_type == "BM_cells") %>%
  transmute(
    patient = Patient,
    event   = "BM sample collection",
    start   = Date_of_sample_collection,
    end     = start,
    details = Sample_ID
  )

# 2. cfDNA (plasma) sample collection events
cfDNA_events <- combined_clinical_data_updated %>%
  # assuming plasma‐derived cfDNA are labeled "Blood_plasma"
  filter(str_detect(Sample_type, regex("plasma", ignore_case = TRUE))) %>%
  transmute(
    patient = Patient,
    event   = "cfDNA sample collection",
    start   = Date_of_sample_collection,
    end     = start,
    details = Sample_ID
  )


# ─────────────────────────────────────────────────────────────────────────────
# F. MRD test dates
# ─────────────────────────────────────────────────────────────────────────────
dat <- read.csv("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.csv")


mfc_events <- dat %>%
  filter(!is.na(Flow_Binary)) %>%
  mutate(
    start   = as_date(Date),
    end     = start,
    details = paste0("Flow result: ", Flow_Binary)
  ) %>%
  transmute(
    patient = Patient,
    event   = "MRD (MFC)",
    start,
    end,
    details
  )

# 2. clonoSEQ (adaptive) MRD test events
clonoseq_events <- dat %>%
  filter(!is.na(Adaptive_Binary)) %>%
  mutate(
    start   = as_date(Date),
    end     = start,
    details = paste0("Adaptive result: ", Adaptive_Binary)
  ) %>%
  transmute(
    patient = Patient,
    event   = "MRD (clonoSEQ)",
    start,
    end,
    details
  )

# Get everythong as date 
# e.g. if bm_events has character dates in “YYYY‑MM‑DD” form:
bm_events <- bm_events %>% 
  mutate(
    start = as_date(start),      # parse “2020‑07‑23” into Date
    end   = as_date(end)         # if end is also character
  )

# repeat for any other dfs with char dates...
cfDNA_events   <- cfDNA_events   %>% mutate(start = as_date(start), end = as_date(end))
followup_events<- followup_events%>% mutate(start = as_date(start), end = as_date(end))
revision_swim_events <- revision_swim_events %>% mutate(start = as_date(start), end = as_date(end))
revision_treatment_events <- revision_treatment_events %>% mutate(start = as_date(start), end = as_date(end))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Combine Everything & Check
# ─────────────────────────────────────────────────────────────────────────────
# Historical interactive versions sometimes had an `m4_progression` object in
# memory. Define an empty compatible table here so a clean Rscript run does not
# depend on workspace state; relapse/progression events are already represented
# in `relapse_events`.
if (!exists("m4_progression")) {
  m4_progression <- tibble(
    patient = character(),
    event = character(),
    start = as.Date(character()),
    end = as.Date(character()),
    details = character(),
    Patient = character(),
    Progression_date = as.Date(character())
  )
}

all_events <- bind_rows(
  tt,
  m4_chemo,
  m4_transplant_events,
  m4_immun_resp,
  m4_best_resp,
  m4_progression,
  imm_clean,
  relapse_events,
  baseline_events,
  followup_events,
  bm_events,
  cfDNA_events,
  revision_treatment_events,
  revision_swim_events,
  clonoseq_events, 
  mfc_events
) %>%
  arrange(patient, start)

all_events <- all_events %>% filter(!is.na(start))

# Quick parse check
bad <- all_events %>%
  filter(is.na(start)) %>%
  distinct(patient, event, start) 

if (nrow(bad)) {
  message("Warning: the following rows have NA start dates:\n")
  print(bad)
} else {
  message("All events parsed successfully!")
}


### Now edit 
all_events <- all_events %>%
  mutate(
    event = case_when(
      # 1) For any SPORE patient whose event name has a digit → Chemotherapy
      str_detect(patient, "^SPORE") & str_detect(event, "\\d+") ~ "Chemotherapy",
      # 2) For any SPORE patient whose details mention ASCT → Transplant
      str_detect(patient, "^SPORE") & 
        str_detect(details, regex("ASCT", ignore_case = TRUE)) ~ "Transplant",
      # 3) Globally, any “Induction” → Chemotherapy
      event == "Induction" ~ "Chemotherapy",
      # 4) Otherwise, keep whatever was there
      TRUE ~ event
    )
  )

all_events <- all_events %>%
  mutate(
    # Collapse both “Best response” and “Immune response” into “Response”
    event = case_when(
      event %in% c("Best response", "Immune response") ~ "Response",
      TRUE                                             ~ event
    )
  ) %>%
  unique() %>%
  # Drop any stray Excel‐origin rows that parsed to 1899-12-31
  filter(start != as.Date("1899-12-31"))

cohort_df <- load_swim_cohort_assignment()

## Filter to the manuscript cohort.
all_events <- all_events %>% filter(patient %in% cohort_df$Patient)

## add end to be start if unsure 
all_events <- all_events %>%
  mutate(
    end = if_else(is.na(end), start, end)
  )

all_events <- all_events %>% select(-any_of(c("Patient", "Progression_date")))


### Add new info from Sarah and Esther 
M4_new <- read_excel("Clinical data/M4/Updated dates from Sarah for swim plot - DA edited.xlsx")
IMG_new <- read_excel("Clinical data/IMMAGINE/Updated_data_esther.xlsx")

  # 1) Parse the M4_new dates into Date class
  M4_new2 <- M4_new %>%
  mutate(
    start = dmy(start),
    end   = dmy(end)
  )

# 2) Drop from all_events any rows that share patient+event+start with M4_new2
all_events_clean <- all_events %>%
  anti_join(
    M4_new2 %>% select(patient, event, start),
    by = c("patient", "event", "start")
  )

# 3) Bind them back together and (optionally) re‑order
all_events_updated <- bind_rows(all_events_clean, M4_new2) %>%
  arrange(patient, start)


# 1) Convert IMG_new’s POSIXct columns to Date
IMG_new2 <- IMG_new %>%
  mutate(
    start = as_date(start),
    end   = as_date(end)
  )

# 2) Drop all “Transplant” rows for IMG‑181 and IMG‑098 from the master table
all_events_clean <- all_events_updated %>%
  filter(!(patient %in% c("IMG-181","IMG-098") & event == "Transplant"))

# 3) (Optional) also guard against exact dupes on patient+event+start
all_events_clean <- all_events_clean %>%
  anti_join(
    IMG_new2 %>% filter(event=="Transplant") %>% select(patient, event, start),
    by = c("patient","event","start")
  )

# 4) Finally bind the cleaned master with the new rows
all_events <- bind_rows(all_events_clean, IMG_new2) %>%
  arrange(patient, start)

## Ensure ASCT treated as transplant 
all_events <- all_events %>%
  mutate(
    event = if_else(
      details == "ASCT",    # when TRUE…
      "Transplant",         # …use this
      event,                # when FALSE…
      missing = event       # …and when NA, also use this
    )
  )

## Now if ongoing put to the latest timepoint have on patient for chemotherapy or progression 
# define which events to fall back on if there is no next Chemo/Transplant
fallback_events <- c("Relapse", "Last follow-up")

# 2) For each patient, stretch only Chemo rows whose end == start
all_events_updated <- all_events %>%
  group_by(patient) %>%
  arrange(start) %>%
  group_modify(~ {
    df <- .
    
    for (i in seq_len(nrow(df))) {
      if (
        !is.na(df$event[i]) &&
        df$event[i] == "Chemotherapy" &&
        !is.na(df$start[i]) &&
        !is.na(df$end[i]) &&
        df$end[i] == df$start[i]
      ) {
        # 1) Next Chemo or Transplant
        next_chemo_tx <- df$start[
          df$event %in% c("Chemotherapy", "Transplant") &
            df$start > df$start[i]
        ]
        
        if (length(next_chemo_tx) > 0) {
          df$end[i] <- min(next_chemo_tx)
        } else {
          # 2) Fallback: earliest Relapse/Last followup AFTER this chemo start
          future_fallback <- df$start[
            df$event %in% fallback_events &
              df$start > df$start[i]
          ]
          if (length(future_fallback) > 0) {
            df$end[i] <- min(future_fallback)
          }
          # else leave df$end[i] as-is (or set to NA if you prefer)
        }
      }
    }
    
    df
  }) %>%
  ungroup()



all_events <- all_events_updated
# Export all_events to CSV. This is a full event-table checkpoint; the mapped
# Supplementary Table 1 export below is the privacy-protected indexed version.
write_csv(all_events %>% select(-any_of("details_2")), "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated.csv")

# Export all_events to RDS
saveRDS(all_events, "Final Tables and Figures/all_events_for_swim_plot_combined_updated2.rds")

#all_events <- readRDS("Final Tables and Figures/all_events_for_swim_plot_combined_updated2.rds")
### Export clean version with correct id 

### this is a supplementary table used in manuscript 

# 1. Load the ID map (must have columns Patient, New_ID)
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

# 2. Merge all_events with id_map and replace Patient with New_ID
all_events_updated <- all_events %>%
  rename(Patient = patient) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(Patient = coalesce(New_ID, Patient),  # use New_ID when available
  details = if_else(
    event %in% c("BM sample collection", "cfDNA sample collection"),
    NA_character_,
    details
  )) %>%
  select(-New_ID, -any_of("details_2"))            # drop temp + details_2 col if present

# 3. Write to CSV
write_csv(all_events_updated,
          "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated_with_new_patient_ID.csv")


# Support-only export for collaborator review of SPORE chemotherapy dates.
# This file is not copied to final_manuscript_objects/ and is not required for
# Figure 1A or Supplementary Table 1 regeneration.
spore_chemo_events <- all_events %>%
  filter(event == "Chemotherapy", grepl("^SPORE", patient)) %>% 
  mutate(end = NA) %>% 
  filter(details != "ASCT")

# Export to CSV
write_csv(spore_chemo_events, file.path(swim_support_dir, "spore_chemo_events.csv"))

# Support-only export for collaborator review of all SPORE timeline events.
# This file is not copied to final_manuscript_objects/ and is not required for
# Figure 1A or Supplementary Table 1 regeneration.
spore_all_events <- all_events %>%
  filter(grepl("^SPORE", patient)) %>% 
  mutate(end = NA) %>% 
  filter(details != "ASCT")

# Export to CSV
write_csv(spore_all_events, file.path(swim_support_dir, "spore_all_events.csv"))


 #### Now assemble plot 
# ──────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ──────────────────────────────────────────────────────────────────────
events   <- read_csv("Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated.csv",
                     col_types = cols(
                       patient = col_character(),
                       event   = col_character(),
                       start   = col_date(),
                       end     = col_date(),
                       details = col_character()
                     ))

events <- events %>%
  # recode any “Relapse” rows into “Progression”
  mutate(
    event = if_else(event == "Relapse", "Progression", event)
  )

cohort_df <- load_swim_cohort_assignment()

# ──────────────────────────────────────────────────────────────────────
# 2. MERGE COHORT INFO  &  KEEP PATIENTS WITH A BASELINE
# ──────────────────────────────────────────────────────────────────────
events <- events %>%
  left_join(cohort_df,  by = c("patient" = "Patient")) %>%
  mutate(cohort = if_else(Cohort == "Frontline", "Front-line cohort",
                          "Non-front-line cohort"))


# (a) Front‐line patients use the "Baseline" event date
baseline_front <- events %>%
  filter(cohort == "Front-line cohort", event == "Baseline") %>%
  select(patient, baseline_date = start)

# (b) Non‐front‐line patients use the first BM or blood draw
#     where Timepoint is “Diagnosis” or “Baseline” in the sample table
non_ids <- cohort_df %>%
  filter(Cohort != "Frontline") %>%
  pull(Patient)

# then compute baseline_non only for those non-frontline patients
baseline_non <- combined_clinical_data_updated %>%
  filter(
    Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
    timepoint_info %in% c("Diagnosis", "Baseline")
  ) %>%
  group_by(Patient) %>%
  summarize(
    baseline_date = if (all(is.na(Date_of_sample_collection))) {
      as.Date(NA)
    } else {
      min(Date_of_sample_collection, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  rename(patient = Patient) %>%
  filter(patient %in% non_ids)

# bind the two sets of baseline dates
# Preserve already-parsed Date columns. Using ymd() on a Date vector can coerce
# valid dates to NA and silently drop non-front-line timelines from the plot.
baseline_non <- baseline_non %>%
  mutate(baseline_date = as_date(baseline_date))

baseline_event_fallback <- events %>%
  filter(.data$patient %in% non_ids, .data$event == "Baseline") %>%
  group_by(.data$patient) %>%
  summarize(baseline_date_event = min(.data$start, na.rm = TRUE), .groups = "drop") %>%
  mutate(baseline_date_event = as_date(.data$baseline_date_event))

baseline_non <- baseline_non %>%
  left_join(baseline_event_fallback, by = "patient") %>%
  mutate(
    baseline_date = if_else(
      is.na(.data$baseline_date) | !is.finite(as.numeric(.data$baseline_date)),
      .data$baseline_date_event,
      .data$baseline_date
    )
  ) %>%
  select(-baseline_date_event)

# (Optional) if baseline_front ended up as character too, parse it likewise:
baseline_front <- baseline_front %>%
  mutate(baseline_date = as_date(baseline_date))

# now bind
baseline_tbl <- bind_rows(baseline_front, baseline_non)

# re‐join and drop any patients w/o a computed baseline
events <- events %>%
  left_join(baseline_tbl, by = "patient") %>%
  filter(!is.na(baseline_date))


# 3b. Create index-date version for privacy (days from baseline)
# First, we need to get baseline_date for each patient from the earlier analysis
baseline_tbl_recoded <- baseline_tbl %>%
  left_join(id_map, by = c("patient" = "Patient")) %>%
  mutate(
    patient = coalesce(New_ID, patient)  # replace when mapping exists
  ) %>%
  select(-New_ID)

all_events_indexed <- all_events_updated %>%
  left_join(baseline_tbl_recoded, by = c("Patient" = "patient")) %>%
  mutate(
    # Calculate days from baseline for each event
    start_day = as.numeric(start - baseline_date),
    end_day   = as.numeric(end - baseline_date)
  ) %>%
  # Replace absolute dates with index dates for privacy
  mutate(
    start = start_day,  # Now contains days from baseline (as numeric)
    end   = end_day
  ) %>%
  # Drop the baseline_date and day columns as they're now in start/end
  select(-baseline_date, -start_day, -end_day) %>%
  # Update column names to reflect index dates
  rename(start_day_from_baseline = start,
         end_day_from_baseline   = end)

write_csv(all_events_indexed,
          "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_INDEX_DATES_privacy_protected.csv")

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 1
#
# What this is:
#   Privacy-protected, indexed event table for the treatment/sample timing swim
#   plot. Dates are represented as days from baseline rather than absolute
#   calendar dates.
#
# Why it is here:
#   This table supports Main Figure 1A and is the code-generated counterpart to
#   final Supplementary Table 1.
#
# Current provenance note:
#   docs/manuscript_artifact_source_map.tsv records that the retained renamed
#   manuscript CSV is used until the exact protected-input
#   version of this old-name export is fully reconciled. We still copy the table
#   generated here so this script is visibly responsible for the regenerated
#   Supplementary Table 1 candidate.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_INDEX_DATES_privacy_protected.csv",
  artifact_id = "STABLE1",
  role = "regenerated_table_csv",
  description = paste(
    "Regenerated privacy-protected indexed event table for Supplementary Table 1;",
    "supports Main Figure 1A swim plot."
  ),
  script_name = "2_1_Part2_Cohort_Swim_Plot.R"
)


# ──────────────────────────────────────────────────────────────────────
# 3. DAYS-FROM-BASELINE & INTERVAL FLAG
# ──────────────────────────────────────────────────────────────────────
events <- events %>%
  mutate(
    start_day   = as.numeric(start   - baseline_date),
    end_day     = as.numeric(end     - baseline_date),
    is_interval = start_day != end_day
  )

## edit typo
events <- events %>%
  # if the details column is exactly "ASCT" (or contains ASCT), recode event
  mutate(
    event = if_else(
      str_detect(details, regex("^ASCT$", ignore_case = TRUE)),
      "Transplant",
      event
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4. EXTEND ONE-DAY CHEMO BARS → 30 d BEFORE NEXT CHEMO START
# ─────────────────────────────────────────────────────────────────────────────
#chemo_adj <- events %>%
#  filter(event == "Chemotherapy") %>%
#  arrange(patient, start_day) %>%
#  group_by(patient) %>%
#  mutate(next_start = lead(start_day)) %>%
#  mutate(end_day = if_else(end_day == start_day & !is.na(next_start) &
#                             (next_start - 30 > start_day),
#                           next_start - 30, end_day)) %>%
#  select(patient, start_day, end_day)

#events <- events %>%
#  left_join(chemo_adj, by = c("patient", "start_day"),
#            suffix = c("", "_new")) %>%
#  mutate(
#    end_day = coalesce(end_day_new, end_day),
#    end     = baseline_date + days(end_day),
#    is_interval = if_else(event == "Chemotherapy" & end_day > start_day,
#                          TRUE, is_interval)
#  ) %>%
#  select(-ends_with("_new"))

# ─────────────────────────────────────────────────────────────────────────────
# 5. REGIMEN → COLOUR GROUP
# ─────────────────────────────────────────────────────────────────────────────

## See events
all_chemo <- events %>%
  filter(event == "Chemotherapy") %>%
  distinct(details) %>%
  arrange(details)

events <- events %>%
  mutate(
    chemo_group = case_when(
      event != "Chemotherapy"                                 ~ NA_character_,
      str_detect(details, regex("CY[Bb]OR.?[DP]", ignore_case = TRUE))  
      ~ "CyBorD", 
      str_detect(details, regex("\\bV?RD\\b|Lenalidomide|Rev", TRUE))
      ~ "R/VRD ± Len",
      str_detect(details, regex("Dara", TRUE))                ~ "Dara-based",
      str_detect(details, regex("Carfilzomib|\\bKD\\b|CAR", TRUE))
      ~ "Carfilzomib-based",
      str_detect(details, regex("Ixazomib|\\bIxa\\b", TRUE))  ~ "Ixazomib-based",
      str_detect(details, regex("Elranatamab", TRUE))         ~ "Elranatamab",
      str_detect(details, regex("Iberdomide", TRUE))          ~ "Iberdomide",
      TRUE                                                    ~ "Other"
    )
  )

## Updated 
events <- events %>%
  mutate(
    chemo_group = case_when(
      event != "Chemotherapy"                                            ~ NA_character_,
      # CyBorD (any D or P variant)
      str_detect(details, regex("CY[Bb]OR.?[DP]", ignore_case=TRUE))    ~ "CyBorD",
      # VRD / RVD / Lenalidomide / Revlimid variants
      str_detect(details, regex("\\bV?RD\\b|Lenalidomide|Rev", TRUE))   ~ "Lenalidomide-based",
      # Daratumumab‐based combos (even DaraPom, DaraRVD, etc.)
      str_detect(details, regex("Dara|Daratumumab", TRUE))              ~ "Dara-based",
      # Pomalidomide‐based (PomDex, KPomD, CAR POM, ELO POM, DAR POM, PomDex)
      str_detect(details, regex("Pom(Dex|D)|POM|KPomD|CAR POM|ELO POM|DAR POM|PomDex", TRUE))
      ~ "Pomalidomide-based",
      # Carfilzomib combinations (Carfil-, KD-, but avoid “DAR CAR” which is caught above)
      str_detect(details, regex("Carfilozomib|\\bKD\\b|CAR(?! *REV)", TRUE))
      ~ "Carfilzomib-based",
      # Ixazomib combos
      str_detect(details, regex("Ixa|Ixazomib", TRUE))                  ~ "Ixazomib-based",
      # Elranatamab (including the /daratumumab trial)
      str_detect(details, regex("Elranatamab", TRUE))                   ~ "Elranatamab",
      # Iberdomide combos
      str_detect(details, regex("Iberdomide", TRUE))                    ~ "Iberdomide",
      # Cyclophosphamide + dexamethasone (“CYCLO DEX” or “CYCLONE”)
      str_detect(details, regex("CYCLO DEX|CYCLONE", TRUE))             ~ "Cyclophosphamide-based",
      # Early-phase/trial drugs
 #     str_detect(details, regex("MEDI|TAK-|VSV", TRUE))                 ~ "Clinical trial",
      # Catch-all for everything else
      TRUE                                                              ~ "Other"
    )
  )

chemo_cols <- c(
  "CyBorD"              = "#4477AA",
  "R/VRD ± Len"         = "#CC6677",
  "Dara-based"          = "#228833",
  "Carfilzomib-based"   = "#AA3377",
  "Ixazomib-based"      = "#66C2A5",
  "Elranatamab"         = "#EE7733",
  "Iberdomide"          = "#994F00",
  "Other"               = "#999999"
)


### Updated
chemo_cols <- c(
  "CyBorD"                    = "#0072B2",  # blue
  "Lenalidomide-based"        = "#CC6677",  # pink 
  "Dara-based"                = "#009E73",  # green
  "Pomalidomide-based"        = "#D55E00",  # vermillon
  "Carfilzomib-based"         = "#E69F00",  # orange
  "Ixazomib-based"            = "#56B4E9",  # sky blue
  "Elranatamab"               = "#F0E442",  # yellow
  "Iberdomide"                = "#800080",  # purple
  "Cyclophosphamide-based"    = "#994F00",  # brown-orange
#  "Clinical trial"            = "#F7C6C7",  # light pink
  "Other"                     = "#999999"   # grey
)

# ─────────────────────────────────────────────────────────────────────────────
# 6. SHAPE MAP FOR ALL POINT-EVENTS
# ─────────────────────────────────────────────────────────────────────────────
shape_map <- c(
  "BM sample collection"   = 3,   # +
  "cfDNA sample collection"= 4,   # x
  "MRD (MFC)"              = 1,   # circle
  "MRD (clonoSEQ)"         = 18,  # diamond
  "MRD (clinical)"         = 17,  # triangle
  "Transplant"             = 7,  # square with x
  "Progression"            = 16   # filled circle
)

# Treat transplants & progression as points
events <- events %>%
  mutate(is_interval = if_else(event %in% c("Transplant", "Progression"),
                               FALSE, is_interval))

# ─────────────────────────────────────────────────────────────────────────────
# 7. PATIENT ORDER (earliest baseline → top)
# ─────────────────────────────────────────────────────────────────────────────
## Change to based on tumor fraction change from ichorCNA
dat_tf <- dat %>%
  mutate(Date = as_date(Date)) %>%
  transmute(
    Patient,
    Date,
    timepoint_info,
    TF = suppressWarnings(as.numeric(WGS_Tumor_Fraction_Blood_plasma_cfDNA))
  )

# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat_tf %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(TF)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = TF,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat_tf %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(TF),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = TF,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
tf_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
tf_change <- tf_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

## Adjust for NAs 
tf_change <- tf_change %>%
  mutate(
    pct_for_plot = case_when(
      !is.na(pct_change) ~ pct_change,
      is.na(pct_change) & baseline_TF == 0 & first_post_TF > 0 ~ 200,  # big value for sorting
      TRUE ~ NA_real_
    ),
    pct_label = case_when(
      !is.na(pct_change) ~ sprintf("%.0f%%", pct_change),
      is.na(pct_change) & baseline_TF == 0 & first_post_TF > 0 ~ ">100%",
      TRUE ~ NA_character_
    )
  )


### Now try with cVAF
## Add more info to patients 
Additional_info <- read_csv("MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated Sep with all patients.csv")

add_lookup <- Additional_info %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         !is.na(sites_rate_zscore_charm)) %>%
  transmute(
    Patient,
    add_Date  = as_date(Date_of_sample_collection_Sample_ID_Bam),
    zscore_blood_from_additional = as.numeric(sites_rate_zscore_charm)
  )

# keep a row id so no row is ever lost
dat_clean <- dat %>%
  mutate(Date = as_date(Date),
         .row_id = row_number())

dat2 <- dat_clean %>%
  # many-to-many join on Patient
  left_join(add_lookup, by = "Patient", relationship = "many-to-many") %>%
  mutate(
    add_Date = as_date(add_Date),
    # if there's no candidate date, set diff = Inf so the original row survives
    day_diff = ifelse(is.na(add_Date), Inf, abs(as.integer(add_Date - Date)))
  ) %>%
  # pick the closest Additional_info date per ORIGINAL ROW (not per Patient+Date)
  group_by(.row_id) %>%
  slice_min(order_by = day_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # fill only if original was NA and the nearest candidate is within 7 days
  mutate(
    zscore_blood = if_else(
      is.na(zscore_blood) &
        is.finite(day_diff) & day_diff <= 7 &
        !is.na(zscore_blood_from_additional),
      zscore_blood_from_additional,
      zscore_blood
    )
  ) %>%
  select(-add_Date, -day_diff, -zscore_blood_from_additional, -.row_id)

# 3) Quick report of how many were filled
filled_n <- sum(is.na(dat$zscore_blood) & !is.na(dat2$zscore_blood))
total_na_before <- sum(is.na(dat$zscore_blood))
total_na_after  <- sum(is.na(dat2$zscore_blood))

message("zscore_blood backfilled: ", filled_n,
        " (NA before: ", total_na_before, ", NA after: ", total_na_after, ")")

# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat2 %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(detect_rate_BM)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = detect_rate_BM,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat2 %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(detect_rate_BM),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = detect_rate_BM,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
BM_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
BM_change <- BM_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

### Now for blood
# 1) First non-NA baseline TF where timepoint_info is Baseline/Diagnosis
baseline_df <- dat2 %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis"),
         !is.na(detect_rate_blood)) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            baseline_TF   = detect_rate_blood,
            baseline_date = Date)

# 2) First non-NA TF strictly AFTER that baseline date
followup_df <- dat2 %>%
  inner_join(baseline_df, by = "Patient") %>%
  filter(!is.na(detect_rate_blood),
         Date > baseline_date) %>%
  arrange(Patient, Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(Patient,
            first_post_TF   = detect_rate_blood,
            first_post_date = Date)

# 3) Combine + percent change (baseline → first post)
Blood_change <- baseline_df %>%
  left_join(followup_df, by = "Patient") %>%
  mutate(
    pct_change = if_else(!is.na(baseline_TF) & !is.na(first_post_TF) & baseline_TF > 0,
                         100 * (first_post_TF - baseline_TF) / baseline_TF,
                         NA_real_)
  )

## Add cohort 
Blood_change <- Blood_change %>%
  left_join(events %>% distinct(patient, cohort) %>% rename(Patient = patient),
            by = "Patient")

#  Add suffixes to every column EXCEPT Patient and cohort
blood_renamed <- Blood_change %>%
  rename_with(~ paste0(.x, "_blood"),
              .cols = -c(Patient, cohort))

bm_renamed <- BM_change %>%
  rename_with(~ paste0(.x, "_bm"),
              .cols = -c(Patient, cohort))

# Join side-by-side on Patient + cohort
change_combined <- full_join(
  blood_renamed,
  bm_renamed,
  by = c("Patient", "cohort")
)


### Identify review-only differences without printing them into a figure build.
# `change_combined` contains support rows beyond the manuscript cohort; it must
# not determine who is present in the final Figure 1A patient order.
out_of_cohort_change_patients <- setdiff(change_combined$Patient, cohort_df$Patient)
cohort_patients_missing_change <- setdiff(cohort_df$Patient, change_combined$Patient)


### Use the tumor fraction instead since too many with NAs for the zscore due to poor quality lists

## Original method, by length of tracking 
# patient_order <- events %>%
#   group_by(cohort, patient) %>%
#   summarize(first_day = min(start_day), .groups = "drop") %>%
#   arrange(cohort, first_day) %>%
#   mutate(y = rev(row_number()))
# 
# events <- events %>% left_join(patient_order, by = c("cohort", "patient"))


# ## Updated, by tumor fraciton 
# patient_order <- tf_change %>%
#   mutate(
#     cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort"))
#   ) %>%
#   arrange(cohort, is.na(pct_for_plot), pct_for_plot) %>%     # NAs last, then ascending %Δ
#   group_by(cohort) %>%
#   mutate(y = row_number()) %>%                               # 1,2,3... within cohort
#   ungroup() %>%
#   transmute(
#     cohort,
#     patient = Patient,                                       # match events$patient
#     y,
#     pct_for_plot,
#     pct_label
#   )

# patient_order_tf <- patient_order


### by the change in the cVAF
patient_order <- change_combined %>%
  mutate(
    cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort")),
    in_cohort = Patient %in% cohort_df$Patient,
    
    # Use BM if present, else blood
    primary_change = coalesce(pct_change_bm, pct_change_blood),
    
    # Record which metric was used
    source_metric = case_when(
      !is.na(pct_change_bm) ~ "bm",
      is.na(pct_change_bm) & !is.na(pct_change_blood) ~ "blood",
      TRUE ~ NA_character_
    ),
    
    # Send to end if not in cohort_df or both metrics missing
    send_to_end = (!in_cohort) | (is.na(pct_change_bm) & is.na(pct_change_blood))
  ) %>%
  arrange(
    cohort,
    send_to_end,
    desc(primary_change),
    Patient
  ) %>%
  group_by(cohort) %>%
  mutate(y = row_number()) %>%
  ungroup() %>%
  transmute(
    cohort,
    patient = Patient,
    y,
    pct_for_plot = primary_change,
    pct_label = if_else(is.na(primary_change), NA_character_, sprintf("%.1f%%", primary_change)),
    source_metric
  )

# ---- Add missing patients at bottom ----
missing_patients <- setdiff(cohort_df$Patient, patient_order$patient)

if (length(missing_patients) > 0) {
  add_rows <- cohort_df %>%
    filter(Patient %in% missing_patients) %>%
    transmute(
      cohort = factor(Cohort, levels = c("Front", "Non-frontline")),
      patient = Patient,
      y = NA_integer_,
      pct_for_plot = NA_real_,
      pct_label = NA_character_,
      source_metric = NA_character_
    )
  
  patient_order <- bind_rows(patient_order, add_rows)
}

patient_order <- patient_order %>%
  mutate(
    cohort = case_when(
      cohort == "Non-frontline" ~ "Non-front-line cohort",  # fix naming
      patient == "IMG-127" ~ "Front-line cohort",           # manual override
      TRUE ~ as.character(cohort)
    ),
    cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort"))
  )


patient_order_cVAF <- patient_order


# 2) Join onto eventpatient_order# 2) Join onto events for plotting
events <- events %>%
  left_join(patient_order, by = c("cohort", "patient"))


# ─────────────────────────────────────────────────────────────────────────────
# 8. FRONT-LINE COHORT PLOT
# ─────────────────────────────────────────────────────────────────────────────
front_data <- events %>% filter(cohort == "Front-line cohort")

p_front <- ggplot() +
  geom_segment(
    data = front_data %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y,
        colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = front_data %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = front_data %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
 # scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen") +
 # scale_shape_manual(values = shape_map, name = "Point events") +
  scale_colour_manual(
    values = chemo_cols,
    name   = "Chemotherapy regimen",
    drop   = FALSE            # <-- keep all colours in legend
  ) +
  scale_shape_manual(
    values = shape_map,
    name   = "Point events",
    drop   = FALSE            # <-- keep all shapes in legend
  ) +
  scale_x_continuous("Days from baseline", expand = expansion(mult = 0.02)) +
  scale_y_continuous(
    NULL,
    breaks = patient_order$y[patient_order$cohort == "Front-line cohort"],
    labels = patient_order$patient[patient_order$cohort == "Front-line cohort"]
  ) +
  ggtitle("Training Cohort") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 6),
    legend.position    = "bottom",
    legend.box         = "vertical",
    plot.title         = element_text(face = "bold", hjust = 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 9. NON-FRONT-LINE COHORT PLOT
# ─────────────────────────────────────────────────────────────────────────────
non_data <- events %>% filter(cohort == "Non-front-line cohort")

p_non <- ggplot() +
  geom_segment(
    data = non_data %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y,
        colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = non_data %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = non_data %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", guide = "none") +
  scale_shape_manual(values = shape_map, name = "Point events" , guide = "none")  +
  scale_x_continuous("Days from baseline", expand = expansion(mult = 0.02))+  scale_y_continuous(
    NULL,
    breaks = patient_order$y[patient_order$cohort == "Non-front-line cohort"],
    labels = patient_order$patient[patient_order$cohort == "Non-front-line cohort"]
  ) +
  ggtitle("Test Cohort") +
  guides(colour = FALSE, shape = FALSE) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 6),
    legend.position    = "none",
  #  legend.box         = "vertical",
    plot.title         = element_text(face = "bold", hjust = 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 10. COMBINE & EXPORT
# ─────────────────────────────────────────────────────────────────────────────
combined_plot <- p_front / p_non +
  plot_layout(guides = "collect", heights = c(3, 1)) &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.x = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) &
  guides(
    colour = guide_legend(
      title   = "Chemotherapy regimen",
      nrow    = 2,
      byrow   = TRUE
    ),
    shape = guide_legend(
      title   = "Point events",
      nrow    = 2,
      byrow   = TRUE
    )
  )

ggsave(file.path(swim_support_dir, "swimplot_by_cohort_v5.pdf"), combined_plot,
       width = 12, height = 10, device = cairo_pdf)

ggsave(file.path(swim_support_dir, "swimplot_by_cohort_v5.png"), combined_plot,
       width = 12, height = 10, dpi = 300)



## Separate cohort-level review plots.
# These are useful QA views but are not the final Figure 1A component.
ggsave(
  filename = file.path(swim_support_dir, "swimplot_nonfrontline_cohort_updated.png"),
  plot     = p_non,
  width    = 10,
  height   = 2,
  dpi      = 500
)

ggsave(
  filename = file.path(swim_support_dir, "swimplot_frontline_cohort_updated.png"),
  plot     = p_front,
  width    = 10,
  height   = 10,
  dpi      = 500
)




# ─────────────────────────────────────────────────────────────────────────────
# 11. SINGLE PLOT WITH BOTH COHORTS, FRONT‐LINE FIRST
# ─────────────────────────────────────────────────────────────────────────────
# A) Build a combined patient ordering: front-line patients first, then non-front-line,
#    each ordered by their first start_day.
# patient_order_combined <- events %>%
#   group_by(patient, cohort) %>%
#   summarize(first_day = min(start_day), .groups = "drop") %>%
#   mutate(
#     cohort = factor(cohort, levels = c("Front-line cohort", "Non-front-line cohort"))
#   ) %>%
#   arrange(cohort, first_day) %>%
#   mutate(y = row_number())

## To use tumor fraction 
patient_order_combined <- patient_order

# B) Join the new y positions back into `events`
events_combined <- events %>%
  select(-any_of("y")) %>%
  left_join(patient_order_combined %>% select(patient, y), by = "patient")

# C) Single swim‐plot
p_combined <- ggplot() +
  # 1) Chemo duration bars
  geom_segment(
    data = events_combined %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y, colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  # 2) Other intervals (e.g. maintenance) in black
  geom_segment(
    data = events_combined %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day, xend = end_day, y = y, yend = y),
    colour = "black", size = 5, lineend = "round"
  ) +
  # 3) Point events (samples, MRD, transplant, progression)
  geom_point(
    data = events_combined %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day, y = y, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  # 4) Scales: chemo colours & shapes, keep all keys
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map,   name = "Point events",          drop = FALSE) +
  # 5) X-axis: days from baseline
  scale_x_continuous("Days from baseline", expand = expansion(mult = c(0.02, 0.02))) +
  # 6) Y-axis: patient names in the combined order, top-down
  scale_y_reverse(
    breaks = patient_order_combined$y,
    labels = patient_order_combined$patient,
    name   = "Patient"
  ) +
  # 7) Theme & legend at bottom
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 10),
    axis.title.x        = element_text(size = 12, face = "bold", vjust = -0.5),
    axis.title.y        = element_text(size = 12, face = "bold", vjust =  1),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.margin        = margin(5, 10, 5, 10)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  )

p_combined <- p_combined +
  labs(
    title = "Swim Plot of Treatment Timelines") +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# D) Display / Save
print(p_combined)
ggsave(file.path(swim_support_dir, "swimplot_combined_cohorts_v3.png"), p_combined,
       width = 10, height = 12, dpi = 500)
ggsave(file.path(swim_support_dir, "swimplot_combined_cohorts_wide_v3.png"), p_combined,
       width = 16, height = 12, dpi = 500)



## Above is figure 1A swim plot 




#### Remove EK-09 pre-treatment and instead add star
## Add other features to make nice 
## Need to specify which patients have baseline BM and not 
### ──────────────────────────────────────────────────────────────
### 1.  Prepare the tumour-fraction table  (ord_df)
### ──────────────────────────────────────────────────────────────
# A) Standardise patient names (strip the "_Baseline")
ord_df <- read.csv("ordering_df_for_Figure_1.csv")
ord_df <- ord_df %>%
  mutate(
    patient = sub("_Baseline$", "", Sample),
    
    # map the 'Cohort' to the swim-plot cohort names
    cohort  = recode(Cohort,
                     "Training" = "Front-line cohort",
                     "Test"    = "Non-front-line cohort"),
    
    # sample-type flag for the y-axis label
    sample_type = case_when(
      Paired                     ~ "Paired",
      !Paired & is.na(TumourFraction) ~ "BM only",
      TRUE                       ~ "Blood only"
    )
  )

baseline_order_availability <- combined_clinical_data_updated %>%
  filter(
    .data$Patient %in% cohort_df$Patient,
    .data$Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
    .data$timepoint_info %in% c("Diagnosis", "Baseline")
  ) %>%
  group_by(Patient) %>%
  summarize(
    has_baseline_bm = any(.data$Sample_type == "BM_cells", na.rm = TRUE),
    has_baseline_blood = any(.data$Sample_type == "Blood_plasma_cfDNA", na.rm = TRUE),
    .groups = "drop"
  )

tumour_fraction_order_lookup <- dat %>%
  mutate(
    Date = as_date(.data$Date),
    TumourFraction = coalesce(
      suppressWarnings(as.numeric(.data$WGS_Tumor_Fraction_Blood_plasma_cfDNA)),
      suppressWarnings(as.numeric(.data$WGS_Tumor_Fraction_BM_cells))
    )
  ) %>%
  filter(
    .data$Patient %in% cohort_df$Patient,
    .data$timepoint_info %in% c("Diagnosis", "Baseline"),
    !is.na(.data$TumourFraction)
  ) %>%
  arrange(.data$Patient, .data$Date) %>%
  group_by(Patient) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Patient, TumourFraction)

missing_order_patients <- setdiff(cohort_df$Patient, ord_df$patient)
if (length(missing_order_patients)) {
  revision_order_rows <- cohort_df %>%
    filter(.data$Patient %in% missing_order_patients) %>%
    left_join(baseline_order_availability, by = "Patient") %>%
    left_join(tumour_fraction_order_lookup, by = "Patient") %>%
    mutate(
      has_baseline_bm = replace_na(.data$has_baseline_bm, FALSE),
      has_baseline_blood = replace_na(.data$has_baseline_blood, FALSE),
      Sample = paste0(.data$Patient, "_Baseline"),
      Cohort = if_else(.data$Cohort == "Frontline", "Training", "Test"),
      Paired = .data$has_baseline_bm & .data$has_baseline_blood,
      patient = .data$Patient,
      cohort = recode(
        .data$Cohort,
        "Training" = "Front-line cohort",
        "Test" = "Non-front-line cohort"
      ),
      sample_type = case_when(
        .data$has_baseline_bm & .data$has_baseline_blood ~ "Paired",
        .data$has_baseline_bm ~ "BM only",
        .data$has_baseline_blood ~ "Blood only",
        TRUE ~ "No baseline BAM"
      )
    ) %>%
    select(Sample, Cohort, TumourFraction, Paired, patient, cohort, sample_type)

  ord_df <- bind_rows(ord_df, revision_order_rows)
}

ord_df <- ord_df %>%
  distinct(.data$patient, .keep_all = TRUE) %>%
  arrange(factor(.data$Cohort, levels = c("Training", "Test")), desc(.data$TumourFraction), .data$patient)

### ──────────────────────────────────────────────────────────────
### 2.  Re-order patients by cohort → descending tumour-fraction
### ──────────────────────────────────────────────────────────────
patient_order_combined <- ord_df %>%
  mutate(
    cohort = factor(Cohort,
                    levels = c("Training", "Test"))
  ) %>%
  group_by(cohort) %>%
  mutate(y = row_number()) %>%
  ungroup() %>%
  select(patient, cohort, y, Paired, sample_type, TumourFraction)


# tack the new 'y' onto events
events_combined <- events %>%
  select(-any_of("y")) %>%
  left_join(patient_order_combined %>% select(patient, y), by = "patient")

patient_levels <- patient_order_combined$patient
patient_order_combined <- patient_order_combined %>%
  mutate(patient = factor(patient, levels = patient_levels))

events_combined <- events_combined %>%
  mutate(patient = factor(patient, levels = patient_levels))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Tumour‐fraction strip
# ─────────────────────────────────────────────────────────────────────────────
ann_tf <- ggplot(patient_order_combined,
                 aes(x = TumourFraction, y = patient, group = 1)) +
  geom_path(colour = "grey70", size = 0.4) +
  geom_point(colour = "grey20", size = 2) +
  scale_x_continuous(
    name   = "cfDNA\ntumour fraction",
    limits = c(0, max(patient_order_combined$TumourFraction, na.rm = TRUE) * 1.05),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    limits = rev(patient_levels),
    expand = c(0, 0)    # shrink vertical padding between rows
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_text(size = 10, hjust = 1, lineheight = 0.8),
    plot.margin     = margin(0, 1, 0, 1)
  )




# ─────────────────────────────────────────────────────────────────────────────
# 3) Cohort colour bar
# ─────────────────────────────────────────────────────────────────────────────
# 1) Recode the cohort labels
#patient_order_combined <- patient_order_combined %>%
#  mutate(cohort = recode(cohort,
#                         "Front-line cohort"     = "Train",
#                         "Non-front-line cohort" = "Test"))

# 2) Define the new colour mapping
cohort_cols <- c(
  "Training" = "#1f77b4",
  "Test"  = "#e6550d"
)


ann_cohort <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = cohort)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  scale_x_continuous(name   = "Cohort", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0),   limits = rev(patient_levels)) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )




# ─────────────────────────────────────────────────────────────────────────────
# 4) Paired‐status bar
# ─────────────────────────────────────────────────────────────────────────────
paired_cols <- c(
  "Paired"     = "#9467bd",  # purple
  "BM only"    = "#1f77b4",  # blue (same as Front‑line cohort)
  "Blood only" = "#CC6677",  # red (same as LEN‑based chemo)
  "No baseline BAM" = "#F7F7F7"
)


ann_paired <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = factor(sample_type))) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = paired_cols, name = "Baseline samples available") +
  scale_x_continuous(name   = "Samples\navailable", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_levels),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )



# ─────────────────────────────────────────────────────────────────────────────
# 5) Swim‐plot (rebuilt on events_combined)
# ─────────────────────────────────────────────────────────────────────────────
events_combined <- events_combined %>%
  mutate(
    start_day_plot = pmax(start_day, -500),
    end_day_plot   = pmax(end_day,   -500)
  )

xmax <- max(events_combined$end_day_plot, na.rm=TRUE) * 1.05

p_swim <- ggplot() +
  geom_segment(
    data = events_combined %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_day_plot, xend = end_day_plot, y = patient, yend = patient, colour = chemo_group),
    size = 5, lineend = "round"
  ) +
  geom_segment(
    data = events_combined %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_day_plot, xend = end_day_plot, y = patient, yend = patient),
    colour = "black", size = 5, lineend = "round"
  ) +
  geom_point(
    data = events_combined %>% filter(!is_interval, event %in% names(shape_map)),
    aes(x = start_day_plot, y = patient, shape = event),
    colour = "black", fill = "white", stroke = 0.35, size = 2.6
  ) +
  scale_colour_manual(values = chemo_cols, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map, name = "Point events", drop = FALSE) +
  scale_x_continuous(
    "Days from baseline",
    limits = c(-500, xmax),
    # pick whatever breaks you like; here's an example every 500 days
    breaks = c(-500, seq(0, ceiling(xmax/500)*500, by = 500)),
    # label -750 as "<-750", everything else as its value
    labels = function(x) ifelse(x == -500, "<-500", as.character(x)),
    expand = expansion(mult = c(0.02,0.02))
  ) +
  scale_y_discrete(labels = NULL, limits = rev(levels(events_combined$patient))) +  # reverse so first is on top
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_blank(),  # we show them in ann_paired
    axis.text.x        = element_text(size = 10),             # ← bigger tick labels
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(title = "Swim Plot of Treatment Timelines", y = NULL)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Patchwork: three left panels + swim plot
# ─────────────────────────────────────────────────────────────────────────────
final_plot <- ann_tf + ann_cohort + ann_paired + p_swim +
  plot_layout(widths = c(0.05, 0.035, 0.035, 0.88),
              guides = "collect") &
  theme(
    # position & spacing
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    # unify title/text/key sizes
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.8, "lines")
  ) &
  guides(
    fill   = guide_legend(
      nrow          = 1,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    colour = guide_legend(
      nrow          = 4,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    shape  = guide_legend(
      nrow          = 3,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8),
      override.aes  = list(size = 2.6)  # match the point size
    )
  )

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_updated.png",
       final_plot,
       width  = 15,
       height = 10,
       dpi    = 500)

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated.png",
       final_plot,
       width  = 18,
       height = 10,
       dpi    = 500)



### Edit to be simpler - less shapes, clearer
# 1) Collapse to two point types

lastfu_tbl <- all_events %>%
  filter(
    event   == "Last follow-up",
    details == "Relapsed: 0"
  ) %>%
  group_by(patient) %>%
  summarise(
    last_followup = max(start, na.rm = TRUE),
    .groups = "drop"
  )

events_combined2 <- events_combined %>%
  left_join(lastfu_tbl, by = "patient") %>%
  mutate(
    is_ongoing = is_interval &
      event == "Chemotherapy" &
      ( is.na(end) | end == last_followup )
  )

# Turn off is_ongoing for VA-02 since treatment stopped
events_combined2 <- events_combined2 %>%
  left_join(lastfu_tbl, by = "patient") %>%
  mutate(
      is_ongoing = if_else(patient == "VA-02", FALSE, is_ongoing)
  )

events_combined2 <- events_combined2 %>%
  mutate(
    point_type = case_when(
      event %in% c("BM sample collection", "cfDNA sample collection") ~ "Sample collection",
      event %in% c("MRD (MFC)", "MRD (clonoSEQ)", "MRD (clinical)")      ~ "Clinical MRD assay",
      TRUE                                                             ~ NA_character_
    ),
    # Events at the index baseline describe presentation, not on-study relapse.
    # Show red progression squares only strictly after baseline.
    is_progression = event == "Progression" & start_day > 0,
    is_transplant  = event == "Transplant"
    )  

events_combined2 <- events_combined2 %>%
  mutate(
    # if event is NA, replace with "Baseline"
    event = if_else(is.na(event), "Baseline", event)
  )

events_combined2 <- events_combined2 %>%
  mutate(
    across(
      starts_with("is_"),       # pick every column whose name starts with "is_"
      ~ replace_na(.x, FALSE)   # turn any NA into FALSE
    )
  )

# 2) Shapes for the two point types
shape_map_simple <- c(
  "Sample collection" = 16,  # circle
  "Clinical MRD assay"         = 17,   #  triangle,
  "Progression" = 15   # ← square
)

events_combined2 <- events_combined2 %>%
  mutate(
    chemo_group_simple = case_when(
      event != "Chemotherapy" ~ NA_character_,
      # IMiDs
      str_detect(details, regex("Lenalidomide|Rev|Pom(Dex|D)|Iberdomide", TRUE)) ~ "IMiD-based",
      # Proteasome inhibitors
      str_detect(details, regex("CY[Bb]OR|Bort|Carfil|KD|Ixa|Ixazomib", TRUE))   ~ "PI-based",
      # Antibodies
      str_detect(details, regex("Dara|Daratumumab|Elranatamab", TRUE))           ~ "Antibody-based",
      # fallback
      TRUE                                                                      ~ "Other"
    )
  )

# Previously treated 
pretx_list <- c("EK-09", "IMG-098", "SPORE_0009")
events_combined2 <- events_combined2 %>%
  mutate(
    previously_treated = patient %in% pretx_list
  )

# Limit early chemo bars
events_combined2 <- events_combined2 %>%
  mutate(
    start_day_plot = if_else(
      event == "Chemotherapy",
      pmax(start_day, 0),  # cap at 0 only for chemo
      start_day_plot               # otherwise leave as‐is
    ),
    end_day_plot = if_else(
      event == "Chemotherapy",
      pmax(end_day, 0),    # cap at 0 only for chemo
      start_day_plot                 # otherwise leave as‐is
    )
  )

## Add ongoing to patient who stopped maintenance early 
#events_combined2 <- events_combined2 %>%
#  mutate(
#    is_ongoing = if_else(
#      patient == "VA-02" & event == "Last follow-up",
#      TRUE,
#      is_ongoing
#    )
#  )

chemo_cols_simple <- c(
  "IMiD-based"     = "#009E73",  # the same green you use for MRD+/maintenance curves
  "PI-based"       = "#E69F00",  # orange used for the cohort bars
  "Antibody-based" = "#0072B2",  # blue used for the panel A ROC lines
  "Other"          = "#999999"   # neutral grey as a fallback
)

## Other 
chemo_cols_simple <- c(
  "IMiD-based"     = "#35B779FF",  # the same green you use for MRD+/maintenance curves
  "PI-based"       = "#E69F00FF",  # orange used for the cohort bars
  "Antibody-based" = "#9467bd",  # purple used for antibody-based treatment
  "Other"          = "#999999"   # neutral grey as a fallback
)

# Original
chemo_cols_simple <- c(
  "IMiD-based"         = "#1b9e77",
  "PI-based"           = "#d95f02",
  "Antibody-based"     = "#7570b3",
  "Other"              = "#999999"
)

## Much lighter
chemo_cols_simple <- c(
  "IMiD-based"     = "#A6D9C7",  # pale seafoam
  "PI-based"       = "#FDBB84",  # peach
  "Antibody-based" = "#B3B5D7",  # lavender-gray
  "Other"          = "#E0E0E0"   # soft gray
)

## Change to weeks 
events_combined2 <- events_combined2 %>%
  mutate(
    start_week_plot = start_day_plot/7,
    end_week_plot = end_day_plot/7
    )

## Months
events_combined2 <- events_combined2 %>%
  mutate(
    start_month_plot = start_day_plot/30.44,
    end_month_plot = end_day_plot/30.44
  )

write_csv(
  events_combined2 %>%
    filter(str_detect(as.character(.data$patient), "^SPORE")) %>%
    select(
      patient,
      event,
      start,
      end,
      baseline_date,
      start_day,
      end_day,
      start_month_plot,
      end_month_plot,
      is_interval,
      details
    ),
  file.path(swim_support_dir, "spore_plot_coordinate_audit.csv")
)

revision_swim_patients <- eligible_spring2026_swim_patients()
write_csv(
  events_combined2 %>%
    filter(.data$patient %in% revision_swim_patients) %>%
    select(
      patient,
      event,
      start,
      end,
      baseline_date,
      start_day,
      end_day,
      start_month_plot,
      end_month_plot,
      is_interval,
      details
    ),
  file.path(swim_support_dir, "oicr_revision_plot_coordinate_audit.csv")
)
rm(revision_swim_patients)

xmax <- max(events_combined2$end_week_plot, na.rm=TRUE) * 1.05
max_months <- ceiling(xmax * 7 / 30.44)  # if xmax is in weeks, convert to months
month_breaks <- seq(0, max_months, by = 3)
minor_breaks <- seq(0, max_months, by = 3)

## Consolidate Points
events_combined2 <- events_combined2 %>%
  mutate(
    event_type = case_when(
      is_transplant      ~ "Transplant",
      is_progression     ~ "Progression",
      is_ongoing         ~ "Ongoing",
      point_type == "Sample collection" ~ "Sample collection",
      point_type == "Clinical MRD assay" ~ "Clinical MRD assay",
      TRUE               ~ NA_character_
    )
  )



### Add the MRD status 
### Also try to show it by MRD status at one year maintenance 

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")

dat <- file 

# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ the master data
  left_join(cohort_df, by = "Patient") 

dat <- dat %>% 
  mutate(
    MRD_truth = case_when(
      !is.na(Adaptive_Binary)            ~ Adaptive_Binary,             # use clonoSEQ if available
      is.na(Adaptive_Binary)
      & !is.na(Flow_Binary)              ~ Flow_Binary,                 # else use MFC 
      TRUE                               ~ NA_real_                     # missing if neither assay run
    )
  )

## Keep only MRD timepoints to optimize timepoint 
dat <- dat %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>% 
  filter(!is.na(MRD_truth)) %>% 
  select(Patient, Timepoint, timepoint_info, Cohort, Adaptive_Frequency, Adaptive_Binary, Flow_pct_cells, Flow_Binary, MRD_truth)



# 1) Reduce to the 1-year maintenance rows and rename assays
dat_1yr <- dat %>%
  filter(timepoint_info == "1yr maintenance") %>%
  transmute(
    patient = Patient,
    MFC      = Flow_Binary,
    clonoSEQ = Adaptive_Binary
  ) %>%
  group_by(patient) %>%
  slice_head(n = 1) %>%   # just in case there are duplicates
  ungroup()

# 2) Ensure every patient in the figure is represented, even if the patient has
#    no 1-year maintenance MRD result. Earlier interactive runs sometimes left
#    `patient_order` as a character vector, but clean command-line runs use the
#    tibble created above. Normalize both forms here so the downstream join is
#    deterministic and does not depend on RStudio session state.
figure_patients <- if (is.data.frame(patient_order)) {
  unique(as.character(patient_order$patient))
} else {
  unique(as.character(patient_order))
}
all_patients <- tibble(patient = figure_patients)
dat_1yr <- all_patients %>%
  left_join(dat_1yr, by = "patient")

# 4) Palette
mrd_cols <- c("MRD+" = "#111111", "MRD-" = "#666666", "Missing" = "#F7F7F7")

## Now make new patient order 

# Build per-patient 1-yr row (you already did this)
# dat_1yr has columns: patient, MFC, clonoSEQ
# If you used my earlier pick-within-window, substitute that table here instead.



## Redo also with the cfWGS MRD 
cfWGS_MRD <- readRDS("Output_tables_2025/all_patients_with_BM_and_blood_calls_updated6_full.rds")

## Add it in to the other one
cfWGS_MRD_1yr <- cfWGS_MRD %>%
  filter(timepoint_info == "1yr maintenance") %>%
  transmute(
    patient = Patient,
    cfWGS_blood      = Blood_zscore_only_sites_call,
    cfWGS_blood_prob      = Blood_zscore_only_sites_prob,
    cfWGS_BM = BM_zscore_only_detection_rate_call,
    cfWGS_BM_prob = BM_zscore_only_detection_rate_prob
  ) %>%
  group_by(patient) %>%
  slice_head(n = 1) %>%   # just in case there are duplicates
  ungroup()

## Add that in as well 
dat_1yr <- dat_1yr %>% 
  left_join(cfWGS_MRD_1yr)

## Get continuous probability 
dat_1yr_prob <- cfWGS_MRD_1yr

# 3) Pivot to long format for two tiles (MFC, clonoSEQ)
df_mrd_long <- dat_1yr %>%
  pivot_longer(cols = c(MFC, clonoSEQ, cfWGS_BM, cfWGS_blood),
               names_to = "assay", values_to = "status") %>%
  mutate(
    assay  = factor(assay, levels = c("MFC", "clonoSEQ", "cfWGS_BM", "cfWGS_blood")),
    status = case_when(
      status == 1 ~ "MRD+",
      status == 0 ~ "MRD-",
      TRUE        ~ "Missing"
    ),
    status = factor(status, levels = c("MRD+", "MRD-", "Missing"))
  )


## Cluster the order 
assay_cols <- c("MFC", "clonoSEQ", "cfWGS_blood", "cfWGS_BM")
# 
# ord_df <- dat_1yr %>%
#   left_join(cohort_df, by = c("patient" = "Patient")) %>%
#   mutate(
#     # Conservative “worst-of-two” decision rule
#     any_pos = (MFC == 1) | (clonoSEQ == 1) | (cfWGS_blood == 1) | (cfWGS_BM == 1),
#     any_neg = (MFC == 0) | (clonoSEQ == 0) | (cfWGS_blood == 0) | (cfWGS_BM == 0),
#     ord_status = case_when(
#       any_pos                    ~ "MRD+",
#       !any_pos & any_neg         ~ "MRD-",
#       TRUE                       ~ "Missing"
#     ),
#     ord_status = factor(ord_status, levels = c("MRD+", "MRD-", "Missing")),
#     
#     # Helpful tie-breakers:
#     # 1) Prefer patients with clonoSEQ data (higher assay sensitivity) over MFC-only within the same status
#     has_clono = !is.na(clonoSEQ),
#     # 2) Cohort ordering with Train above Test (adjust to the labels)
#     Cohort    = fct_relevel(Cohort, c("Train", "Frontline", "Test", "Non-frontline"))
#   )

## For sorting on all MRD 
# ord_df <- dat_1yr %>%
#   left_join(cohort_df, by = c("patient" = "Patient")) %>%
#   mutate(Cohort = fct_relevel(Cohort, c("Train","Frontline","Test","Non-frontline"))) %>%
#   mutate(
#     pmax_val = do.call(pmax, c(across(all_of(assay_cols)), na.rm = TRUE)),
#     pmin_val = do.call(pmin, c(across(all_of(assay_cols)), na.rm = TRUE)),
#     any_pos  = if_else(is.infinite(pmax_val), FALSE, pmax_val == 1),
#     any_neg  = if_else(is.infinite(pmin_val), FALSE, pmin_val == 0),
#     ord_status = case_when(
#       any_pos            ~ "MRD+",
#       !any_pos & any_neg ~ "MRD-",
#       TRUE               ~ "Missing"
#     ),
#     ord_status  = factor(ord_status, levels = c("MRD+","MRD-","Missing")),
#     has_BM      = !is.na(cfWGS_BM),
#     across(all_of(assay_cols), \(z) case_when(z == 1 ~ 1, z == 0 ~ 0.5, TRUE ~ 0), .names = "sc_{.col}"),
#     status_score = rowSums(across(starts_with("sc_"))),
#     n_pos        = rowSums(across(all_of(assay_cols), \(z) z == 1), na.rm = TRUE),
#     n_measured   = rowSums(across(all_of(assay_cols), \(z) !is.na(z)))
#   ) %>%
#   select(-pmax_val, -pmin_val)

## For just sorting on cfWGS BM, then cfDNA
# --- build ord_df as you already do ---
ord_df <- dat_1yr %>%
  left_join(cohort_df, by = c("patient" = "Patient")) %>%
  # Exclude support-only patients (for example SPORE_0008) that have no
  # manuscript cohort assignment before constructing the Figure 1A order.
  filter(!is.na(Cohort)) %>%
  mutate(Cohort = forcats::fct_relevel(Cohort, c("Training","Frontline","Test","Non-frontline"))) %>%
  mutate(
    pmax_val = do.call(pmax, c(across(all_of(assay_cols)), na.rm = TRUE)),
    pmin_val = do.call(pmin, c(across(all_of(assay_cols)), na.rm = TRUE)),
    any_pos  = if_else(is.infinite(pmax_val), FALSE, pmax_val == 1),
    any_neg  = if_else(is.infinite(pmin_val), FALSE, pmin_val == 0),
    ord_status = case_when(
      any_pos            ~ "MRD+",
      !any_pos & any_neg ~ "MRD-",
      TRUE               ~ "Missing"
    ),
    ord_status  = factor(ord_status, levels = c("MRD+","MRD-","Missing")),
    has_BM      = !is.na(cfWGS_BM),
    across(all_of(assay_cols), \(z) case_when(z == 1 ~ 1, z == 0 ~ 0.5, TRUE ~ 0), .names = "sc_{.col}"),
    status_score = rowSums(across(starts_with("sc_"))),
    n_pos        = rowSums(across(all_of(assay_cols), \(z) z == 1), na.rm = TRUE),
    n_measured   = rowSums(across(all_of(assay_cols), \(z) !is.na(z))),
    # BM-first ranking: 1→top, 0→next, NA→last
    bm_rank = dplyr::case_when(cfWGS_BM == 1 ~ 0L,
                               cfWGS_BM == 0 ~ 1L,
                               TRUE          ~ 2L),
    cfDNA_rank = dplyr::case_when(cfWGS_blood == 1 ~ 0L,
                               cfWGS_blood == 0 ~ 1L,
                               TRUE          ~ 2L)
    
  ) %>%
  select(-pmax_val, -pmin_val)

# make sure Cohort has the right order once
cohort_levels <- c("Training","Frontline","Test","Non-frontline")
ord_df <- ord_df %>%
  mutate(Cohort = forcats::fct_relevel(Cohort, cohort_levels))

patient_order_mrd <- ord_df %>%
  arrange(
    bm_rank,               # BM+: 1 first, then BM-: 0, then NA
    cfDNA_rank,
    ord_status,            # then MRD+ → MRD− → Missing
    desc(status_score),
    desc(has_BM),
    desc(n_pos),
    desc(n_measured),
    as.integer(Cohort),
    patient
  ) %>%
  distinct(patient) %>%
  pull(patient)

# Final y-axis order (MRD+ at top, then MRD-, then Missing; within each, clonoSEQ-present first)
# patient_order_mrd <- ord_df %>%
#   arrange(ord_status, desc(has_clono), Cohort, patient) %>%
#   pull(patient) %>%
#   unique()
# 

# For old table version
# patient_order_mrd <- ord_df %>%
#   arrange(
#     ord_status,            # MRD+ → MRD− → Missing (already a factor in that order)
#     desc(status_score),
#     desc(has_BM),          
#     desc(n_pos),
#     desc(n_measured),
#     as.integer(Cohort),    # factor index = desired cohort order
#     patient
#   ) %>%
#   distinct(patient) %>%
#   pull(patient)


## Reorder 
# 1) Compute total follow‑up per patient
patient_order_tbl <- events_combined2 %>% filter(!(patient == "VA-02" & event == "Last follow-up")) %>% # remove this since not plotted
  group_by(Cohort, patient) %>%
  summarise(
    # use end_day if available, otherwise start_day
    followup = max(coalesce(end_day, start_day), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # within each cohort, sort descending so longest follow‑up is first
  arrange(Cohort, desc(followup))

# 2) Pull out the ordered patient vector
patient_order <- patient_order_tbl %>% pull(patient)

## Change back to tumor fraction 
patient_order_tbl_cVAF <- patient_order_cVAF %>%
  group_by(cohort, patient) %>%
  # within each cohort, sort descending so longest follow‑up is first
  arrange(cohort, desc(pct_for_plot)) %>% 
  pull(patient)

patient_order <- patient_order_tbl_cVAF 


## Or for MRD 
patient_order <- c(
  patient_order_mrd,
  setdiff(as.character(patient_order_tbl$patient), patient_order_mrd)
) %>%
  unique() %>%
  intersect(cohort_df$Patient)

if (!all(patient_order %in% cohort_df$Patient)) {
  stop("Figure 1A patient order contains a patient outside the manuscript cohort.", call. = FALSE)
}

# 3) Re‑factor the patient column
events_combined2 <- events_combined2 %>%
  mutate(
    patient = factor(patient, levels = patient_order)
  )  


### Edit the patient labels for de-identification in figure
df <- cohort_df
make_patient_id_map <- function(df, patient_col = "Patient") {
  p <- sym(patient_col)
  
  df %>%
    distinct(!!p, .keep_all = FALSE) %>%         # keep first occurrence, preserve row order
    mutate(
      .prefix = case_when(
        str_starts(!!p, "IMG")   ~ "IMG",
        str_starts(!!p, "SPORE") ~ "SPORE",
        TRUE                     ~ "M4"
      )
    ) %>%
    group_by(.prefix) %>%
    mutate(
      .idx   = row_number(),
      # width: 2 digits unless group size >= 100, then compute digits of n()
      .width = if_else(dplyr::n() >= 100L,
                       as.integer(floor(log10(dplyr::n())) + 1L),
                       2L)
    ) %>%
    ungroup() %>%
    mutate(New_ID = paste0(.prefix, "-", str_pad(.idx, width = .width, pad = "0"))) %>%
    select(!!p, New_ID) %>%
    rename(Patient = !!p)
}


id_map <- make_patient_id_map(df)

# sanity checks
stopifnot(anyDuplicated(id_map$Patient) == 0)
stopifnot(anyDuplicated(id_map$New_ID)  == 0)

cohort_df_anon <- cohort_df %>% left_join(id_map, by = "Patient")

## Re-export supp table with new IDs 
all_events_tmp <- all_events %>%
  left_join(id_map, by = c("patient" = "Patient")) %>%
  mutate(patient = New_ID) %>%
  select(names(all_events))  # keep original column order

write_csv(all_events_tmp %>% select(-any_of("details_2")), "Final Tables and Figures/Supp_Table_1_all_events_for_swim_plot_combined_updated_with_new_IDs.csv")

# Save as CSV
write.csv(id_map,
          file = "id_map.csv",
          row.names = FALSE)

# Save as RDS
saveRDS(id_map,
        file = "id_map.rds")

id_map <- readRDS("id_map.rds")

## Get labels vector
# build label lookup (id_map: Patient, New_ID)
lab_vec <- id_map$New_ID; names(lab_vec) <- id_map$Patient

lab_fun <- function(y) {
  yy  <- as.character(y)
  out <- unname(lab_vec[yy])             # remove names → pure character
  nas <- is.na(out); out[nas] <- yy[nas] # fallback to original
  out
}

# This is the order the plot actually uses (it is reversed later):
y_levels <- rev(levels(events_combined2$patient))

# Apply audited duplicate-event masks before marking progression layers.
events_combined2 <- events_combined2 %>%
  apply_manual_swim_plot_event_masks("spore0012_duplicate_progression_mask")


# 2) Numeric y for the progression layer that matches the plotted order (top row = 1)
prog_marks <- events_combined2 %>%
  dplyr::filter(is_progression) %>%
  dplyr::mutate(y_num = match(patient, y_levels))  # fast + robust

# 3) Plot
p_swim <- ggplot() +
  # Chemo intervals (coloured)
  geom_segment(
    data = events_combined2 %>% filter(is_interval, event == "Chemotherapy"),
    aes(x = start_month_plot, xend = end_month_plot, y = patient, yend = patient, colour = chemo_group_simple),
    linewidth = 5, lineend = "round"
  ) +
  # Other intervals (black)
  geom_segment(
    data = events_combined2 %>% filter(is_interval, event != "Chemotherapy"),
    aes(x = start_month_plot, xend = end_month_plot, y = patient, yend = patient),
    colour = "black", linewidth = 5, lineend = "round"
  ) +
  # Progression square (filled)
  geom_point(
    data = prog_marks %>% mutate(point_type = "Progression"),
    aes(x = start_month_plot, y = y_num, shape = point_type),
    size        = 2.6,
    colour      = "red",
    inherit.aes = TRUE,
    show.legend = FALSE
  ) +
  # Point events (ONLY two shapes)
  geom_point(
    data = events_combined2 %>% filter(!is_interval, !is.na(point_type)),
    aes(x = start_month_plot, y = patient, shape = point_type),
    size = 2.6, stroke = 0.35, colour = "black", fill = "white"
  ) +
  # Progression tick (vertical line across the bar)
  # geom_text(
  #   data   = events_combined2 %>% filter(is_progression),
  #   aes(x = start_month_plot, y = patient),
  #   label    = "|",          # single vertical bar
  #   fontface = "bold",       # bold weight
  #   size     = 5,            # tweak to taste
  #   colour   = "red"
  # ) +
  ### Thicker relapse segment
  # geom_segment(
  #   data = prog_marks,
  #   aes(
  #     x    = start_month_plot,
  #     xend = start_month_plot,
  #     y    = y_num - 0.35,   # height matches the tile height ~0.7
  #     yend = y_num + 0.35
  #   ),
  #   colour      = "red",
  #   linewidth   = 1.75,
  #   lineend     = "round",
  #   inherit.aes = FALSE
  # ) +
  # Last follow-up arrow
  geom_segment(
    data = events_combined2 %>% filter(is_ongoing),
    aes(
      x    = end_month_plot,
      xend = end_month_plot + 50/30.44,      # arrow extends 50 days to the right
      y    = patient,
      yend = patient
    ),
    colour = "black",
    linewidth = 0.6,
    lineend   = "butt",                # so the bar doesn’t cap‑over the arrow
    arrow     = arrow(
      length = unit(2, "mm"),          # size of the arrow head
      ends   = "last",                 # only draw the head at the end
      type   = "closed"                # filled triangle
    )
  ) +
  # Transplant "T"
  geom_text(
    data = events_combined2 %>% filter(is_transplant),
    aes(x = start_month_plot, y = patient, label = "T"),
    fontface = "bold", size = 3
  ) +
  # Previously treated 
  geom_point(
    data = filter(events_combined2, previously_treated),
    aes(x = -40/30.44, y = patient),
    shape  = 8,        # asterisk/star glyph
    size   = 2,        # adjust as needed
    colour = "black"
  ) +
  # Scales
  scale_colour_manual(values = chemo_cols_simple, name = "Chemotherapy regimen", drop = FALSE) +
  scale_shape_manual(values = shape_map_simple, name = "Point events", drop = FALSE) +
  scale_fill_manual(
    values = c(
      "Progression" = "red",
      "Clinical MRD assay" = "white",
      "Sample collection"  = "white"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    "\nMonths since baseline",
    limits = c(-1.4, max(month_breaks)),
    breaks = month_breaks,
    labels = as.character,
    minor_breaks = minor_breaks, # minor ticks every 3 months
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(labels = NULL, limits = y_levels) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),  # ← turns off vertical grid lines
    panel.grid.minor.x = element_blank(),  # ← turns off minor vertical grid lines
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(size = 12),
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.spacing.x   = unit(0.3, "cm"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(title = "Swim Plot of Treatment Timelines", y = NULL)

p_swim <- p_swim +
  scale_y_discrete(
    # by default it will use the patient factor levels
    limits = rev(levels(events_combined2$patient))
  ) 

## Change order 
p_swim <- p_swim + 
  scale_y_discrete(
    limits = rev(levels(events_combined2$patient)),
    breaks = NULL    # no ticks or labels
  )

ggsave(file.path(swim_support_dir, "Figure1A_draft_swim_plot_without_annotations.png"),
       p_swim,
       width  = 15,
       height = 10,
       dpi    = 500)

### Add back cohort and paired status 
# The source cohort table uses Frontline/Non-frontline, whereas older ordering
# tables used Training/Test. Normalize this annotation to the labels shown in
# the manuscript before applying the scale so Frontline rows cannot fall through to the
# default missing-value grey.
cohort_cols <- c(
  "Train" = "#1f77b4",
  "Test"  = "#e6550d"
)

## Re-factor to correct order and derive stable cohort display labels.
patient_order_combined <- patient_order_combined %>%
  select(-any_of("cohort")) %>%
  left_join(cohort_df %>% select(Patient, Cohort), by = c("patient" = "Patient")) %>%
  mutate(
    patient = factor(patient, levels = patient_order),
    cohort = case_when(
      Cohort %in% c("Training", "Frontline") ~ "Train",
      Cohort %in% c("Test", "Non-frontline") ~ "Test",
      TRUE ~ NA_character_
    ),
    cohort = factor(cohort, levels = c("Train", "Test"))
  )

if (anyNA(patient_order_combined$cohort)) {
  stop(
    "Figure 1A cohort annotation contains patients without a recognized cohort label: ",
    paste(patient_order_combined$patient[is.na(patient_order_combined$cohort)], collapse = ", "),
    call. = FALSE
  )
}

ann_cohort <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = cohort)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = cohort_cols, name = "Cohort") +
  scale_x_continuous(name = "Cohort", limits = c(0.5, 1.5), expand = c(0,0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_order),  # character levels, not numeric
    labels = lab_fun,
    expand = c(0,0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_text(size = 10, hjust = 1, lineheight = 0.8),
    plot.margin  = margin(0,0,0,0)
  )

# ann_cohort <- ggplot(patient_order_combined,
#                      aes(x = 1, y = patient, fill = cohort)) +
#   geom_tile(width = 0.9, height = 0.9) +
#   scale_fill_manual(values = cohort_cols, name = "Cohort") +
#   scale_x_continuous(name   = "Cohort", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
#  # scale_y_discrete(expand = c(0, 0),   limits = rev(patient_order)) +
#   scale_y_discrete(
#     limits = rev(patient_order),
#     labels = function(y) {  # anonymized labels
#       out <- lab_vec[y]
#       out[is.na(out)] <- y[is.na(out)]   # fallback: original if any ID missing
#       out
#     },
#     expand = c(0, 0)
#   ) +
#   theme_minimal(base_size = 8) +
#   theme(
#     panel.grid      = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.title.x    = element_text(size = 10, hjust = 0.5, lineheight = 0.8),
#     axis.ticks.x    = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.text.y     = element_text(size = 10, hjust = 1, lineheight = 0.8), # For patients
#     plot.margin     = margin(0, 0, 0, 0)
#   )

paired_cols <- c(
  "Paired"     = "#9467bd",  # purple
  "BM only"    = "#1f77b4",  # blue (same as Front‑line cohort)
  "Blood only" = "#CC6677",  # red (same as LEN‑based chemo)
  "No baseline BAM" = "#F7F7F7"
)


ann_paired <- ggplot(patient_order_combined,
                     aes(x = 1, y = patient, fill = factor(sample_type))) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = paired_cols, name = "Baseline samples available") +
  scale_x_continuous(name   = "Samples\navailable", limits = c(0.5, 1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(
    limits = rev(patient_order),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.title.y    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.margin     = margin(0, 0, 0, 0)
  )


# ### Make the tumor fraction annotation
# ## Continue from here
# # --- prep: order + helpers ---
# patient_order_cVAF <- patient_order_cVAF %>%
#   mutate(
#     patient = factor(patient, levels = patient_order_tbl_cVAF),
#     sign = case_when(
#       is.na(pct_for_plot)       ~ "missing",
#       pct_for_plot >= 0         ~ "increase",
#       TRUE                      ~ "decrease"
#     )
#   )
# 
# # symmetric limits around 0 (at least ±100 so the specified axis breaks fit where possible)
# rng     <- max(abs(patient_order_cVAF$pct_for_plot), na.rm = TRUE)
# max_abs <- max(100, ceiling(rng / 10) * 10)
# limits_x <- c(-max_abs, max_abs)
# 
# # specified axis breaks, clipped to the plotting range
# breaks_wanted <- c(-100, -50, 0, 50, 100)
# breaks_x <- breaks_wanted[breaks_wanted >= limits_x[1] & breaks_wanted <= limits_x[2]]
# 
# # a small x for left annotation (inside plot area)
# x_annot <- limits_x[1] + 0.02 * diff(limits_x)
# 
# ann_tf <- ggplot(patient_order_cVAF,
#                  aes(x = pct_for_plot, y = patient, group = 1)) +
#   # draw change line from 0 to each point
#   geom_segment(aes(x = 0, xend = pct_for_plot, yend = patient),
#                linewidth = 0.6, color = "grey70", na.rm = TRUE) +
#   # point color by source metric
#   geom_point(aes(color = source_metric), size = 2, na.rm = TRUE) +
#   # dotted line at zero
#   geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, color = "grey40") +
#   # x axis symmetric and labelled in %
#   scale_x_continuous(
#     name   = expression(Delta~"ctDNA cVAF"),
#     limits = c(-100, 100),        # adjust if the range is larger
#     breaks = seq(-100, 100, 50),
#     labels = scales::label_number(accuracy = 1, suffix = "%"),
#     expand = c(0, 0)
#   ) +
#   # use the patient order
#   scale_y_discrete(
#     limits = rev(patient_order_tbl_cVAF),
#     expand = c(0, 0)
#   ) +
#   # color palette for clarity
#   scale_color_manual(
#     name = "cVAF Source",
#     values = c("bm" = "black", "blood" = "grey55"),
#     labels = c("bm" = "Bone marrow", "blood" = "Blood"),
#     na.translate = FALSE
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_minimal(base_size = 9) +
#   theme(
#     panel.grid      = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.text.y     = element_text(size = 9, hjust = 1, lineheight = 0.85),
#     plot.margin     = margin(5, 15, 5, 5),
#     legend.position = "NA",
#     legend.title    = element_text(size = 9),
#     legend.text     = element_text(size = 8)
#   )
# 
# ann_tf
# 
# ## Hide y-axis 
# ann_tf <- ann_tf + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# 
# ann_tf




# 5) Build the annotation band

df_mrd_long <- df_mrd_long %>%
  mutate(patient = factor(patient, levels = patient_order))

ann_mrd <- ggplot(df_mrd_long %>% filter(assay %in% c("cfWGS_BM", "cfWGS_blood")),
                  aes(x = assay, y = patient, fill = status)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = mrd_cols, name = "MRD at 1yr-maintenance") +
  scale_x_discrete(
    name = "MRD (1yr)",
    position = "bottom",
    labels = c(
      "MFC"         = "MFC",
      "clonoSEQ"    = "clonoSEQ",
      "cfWGS_BM"    = "cfWGS (BM informed)",
      "cfWGS_blood" = "cfWGS (BM naive)"
    )
  ) +
  scale_y_discrete(
    limits = rev(patient_order_mrd),  # keep the defined order
    labels = lab_fun,
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1),  # rotated labels
    plot.margin  = margin(0, 0, 0, 0)
  )

# ### Old way
# ann_tf <- ggplot(patient_order_cVAF,
#                  aes(x = pct_for_plot, y = patient, group = 1)) +
#   geom_path(colour = "grey70", size = 0.4) +
#   geom_point(colour = "grey20", size = 2) +
#   scale_x_continuous(
#     name   = "Delta\nctDNA cVAF",
#     limits = c(min(patient_order_cVAF$pct_for_plot), max(patient_order_cVAF$pct_for_plot, na.rm = TRUE) * 1.2),
#     expand = c(0, 0)
#   ) +
#   scale_y_discrete(
#     limits = rev(patient_order_tbl_cVAF),
#     expand = c(0, 0)    # shrink vertical padding between rows
#   ) +
#   theme_minimal(base_size = 8) +
#   theme(
#     panel.grid      = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.text.y     = element_text(size = 10, hjust = 1, lineheight = 0.8),
#     plot.margin     = margin(0, 1, 0, 1)
#   )

# ### Show the actual probability in BM if they have or else blood 
# ## Continue from here
# # --- prep: order + helpers ---
# patient_order_cVAF <- patient_order_cVAF %>%
#   mutate(
#     patient = factor(patient, levels = patient_order_tbl_cVAF),
#     sign = case_when(
#       is.na(pct_for_plot)       ~ "missing",
#       pct_for_plot >= 0         ~ "increase",
#       TRUE                      ~ "decrease"
#     )
#   )
# 
# # symmetric limits around 0 (at least ±100 so the specified axis breaks fit where possible)
# rng     <- max(abs(patient_order_cVAF$pct_for_plot), na.rm = TRUE)
# max_abs <- max(100, ceiling(rng / 10) * 10)
# limits_x <- c(-max_abs, max_abs)
# 
# # specified axis breaks, clipped to the plotting range
# breaks_wanted <- c(-100, -50, 0, 50, 100)
# breaks_x <- breaks_wanted[breaks_wanted >= limits_x[1] & breaks_wanted <= limits_x[2]]
# 
# # a small x for left annotation (inside plot area)
# x_annot <- limits_x[1] + 0.02 * diff(limits_x)
# 
# ann_tf <- ggplot(patient_order_cVAF,
#                  aes(x = pct_for_plot, y = patient, group = 1)) +
#   # draw change line from 0 to each point
#   geom_segment(aes(x = 0, xend = pct_for_plot, yend = patient),
#                linewidth = 0.6, color = "grey70", na.rm = TRUE) +
#   # point color by source metric
#   geom_point(aes(color = source_metric), size = 2, na.rm = TRUE) +
#   # dotted line at zero
#   geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, color = "grey40") +
#   # x axis symmetric and labelled in %
#   scale_x_continuous(
#     name   = "MRD Probability",
#     limits = c(-100, 100),        # adjust if the range is larger
#     breaks = seq(-100, 100, 50),
#     labels = scales::label_number(accuracy = 1, suffix = "%"),
#     expand = c(0, 0)
#   ) +
#   # use the patient order
#   scale_y_discrete(
#     limits = rev(patient_order_tbl_cVAF),
#     expand = c(0, 0)
#   ) +
#   # color palette for clarity
#   scale_color_manual(
#     name = "cVAF Source",
#     values = c("bm" = "black", "blood" = "grey55"),
#     labels = c("bm" = "Bone marrow", "blood" = "Blood"),
#     na.translate = FALSE
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_minimal(base_size = 9) +
#   theme(
#     panel.grid      = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.ticks.y    = element_blank(),
#     axis.text.y     = element_text(size = 9, hjust = 1, lineheight = 0.85),
#     plot.margin     = margin(5, 15, 5, 5),
#     legend.position = "NA",
#     legend.title    = element_text(size = 9),
#     legend.text     = element_text(size = 8)
#   )
# 
# ann_tf
# 
# ## Hide y-axis
# ann_tf <- ann_tf + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
# 
# ann_tf


## Assemble final plot
final_plot <- ann_cohort + ann_mrd + p_swim +
 # plot_layout(widths = c(0.015, 0.06, 0.925),
  plot_layout(widths = c(0.015, 0.03, 0.925),
              guides = "collect") &
  theme(
    # position & spacing
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    # unify title/text/key sizes
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.8, "lines")
  ) &
  guides(
    fill   = guide_legend(
      nrow          = 1,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    colour = guide_legend(
      nrow          = 4,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8)
    ),
    shape  = guide_legend(
      nrow          = 3,
      byrow         = TRUE,
      title.position = "top",
      label.theme   = element_text(size = 7),
      title.theme   = element_text(size = 8),
      override.aes  = list(size = 2.6)  # match the point size
    )
  )

final_plot

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_updated10.png",
       final_plot,
       width  = 15,
       height = 10,
       dpi    = 500)

ggsave("Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated10A.png",
       final_plot,
       width  = 16.5,
       height = 10,
       dpi    = 500)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 1A
#
# What this is:
#   The wide treatment/sample timing swim plot with annotation tracks.
#
# Why it is here:
#   This is the script-generated component used for final Main Figure 1A. The
#   assembled Figure 1 PDF is created outside this script, but this PNG is the
#   traceable plotted source component.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Figure1A_swimplot_with_3_annotations_wide_updated10A.png",
  artifact_id = "FIG1A",
  role = "figure_panel_png",
  description = "Wide treatment and sample timing swim plot used as Main Figure 1A.",
  script_name = "2_1_Part2_Cohort_Swim_Plot.R"
)



## Add legend 

# 1) Build a dummy data.frame with one row per symbol
legend_df <- data.frame(
  y     = c(3, 2, 1),
  label = c("Transplant", "Progression", "Ongoing")
)

# 2) Make the legend plot
p_symbols <- ggplot(legend_df, aes(y = y)) +
  # Title
  ggtitle("Symbols") +
  # 2a) Transplant: a bold "T"
  geom_text(
    data = subset(legend_df, label=="Transplant"),
    aes(x = 0, label = "T"),
    fontface = "bold",
    size     = 6
  ) +
  # 2b) Progression: red vertical line
  geom_text(
    data = subset(legend_df, label=="Progression"),
    aes(x = 0, y = y),
    label    = "|",          # single vertical bar
    fontface = "bold",       # bold weight
    size     = 6,            # match the other legend sizes
    colour   = "red"
  ) +
  # 2c) Ongoing: black right arrow
  geom_segment(
    data = subset(legend_df, label=="Ongoing"),
    aes(x = 0, xend = 0.6, yend = y),
    arrow   = arrow(length = unit(4, "mm"), ends="last", type="closed"),
    colour  = "black",
    size    = 1
  ) +
  # text labels
  geom_text(
    aes(x = 1.2, label = label),
    hjust = 0,
    size  = 5
  ) +
  # clean up
  scale_y_continuous(limits = c(0.5, 3.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.2, 2), expand = c(0,0)) +
  theme_void() +
  theme(
    plot.title      = element_text(size=14, face="bold", hjust=0),
    plot.margin     = margin(5,5,5,5)
  )

# 3) (Optional) view it
print(p_symbols)











#### Optional standalone legend export
# This block is not used to make the mapped Figure 1A component. It is retained
# only as a helper for exporting standalone legends during manual assembly.
export_standalone_swimplot_legend <- FALSE

if (isTRUE(export_standalone_swimplot_legend)) {
# 1. Dummy data for chemo colours
df_chemo <- data.frame(
  chemo_group = factor(names(chemo_cols), levels = names(chemo_cols)),
  x = seq_along(names(chemo_cols)),
  y = 1
)

# 2. Dummy data for point-event shapes
df_shape <- data.frame(
  event = factor(names(shape_map), levels = names(shape_map)),
  x = seq_along(names(shape_map)),
  y = 2
)

# 3. Build dummy plot to generate full legend
legend_plot <- ggplot() +
  geom_point(
    data = df_chemo,
    aes(x = x, y = y, colour = chemo_group),
    size = 5
  ) +
  geom_point(
    data = df_shape,
    aes(x = x, y = y, shape = event),
    size = 5
  ) +
  scale_colour_manual(
    name   = "Chemotherapy regimen",
    values = chemo_cols,
    drop   = FALSE
  ) +
  scale_shape_manual(
    name   = "Point events",
    values = shape_map,
    drop   = FALSE
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.x = unit(0.3, "cm"),
    legend.spacing.y = unit(0.3, "cm")
  )

# 4. Extract/export the legend grob if cowplot is available.
if (requireNamespace("cowplot", quietly = TRUE)) {
  full_legend <- cowplot::get_legend(legend_plot)
  
  # Draw to the current device for interactive inspection.
  grid::grid.newpage()
  grid::grid.draw(full_legend)
  
  # Save standalone legend helpers.
  ggsave(file.path(swim_support_dir, "swimplot_full_legend.png"), full_legend,
         width = 14, height = 2, dpi = 500)
  ggsave(file.path(swim_support_dir, "full_swimplot_legend.pdf"), full_legend,
         width = 8, height = 2, device = cairo_pdf)
  ggsave(file.path(swim_support_dir, "full_swimplot_legend.png"), full_legend,
         width = 8, height = 2, dpi = 300)
} else {
  message("Skipping optional standalone swim-plot legend export because cowplot is not installed.")
}
} else {
  message("Skipping optional standalone swim-plot legend export; Figure 1A uses the main swim-plot output.")
}
