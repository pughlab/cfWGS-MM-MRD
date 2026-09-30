###############################################################################
# 4_2_Compare_subclonal_evolution.R
#
# Purpose:
#   Copy-number-based assessment of emergent subclonal CNA events between
#   paired baseline and relapse cfDNA WGS samples (30-40x). Identifies CNAs
#   present at relapse but absent (or sub-threshold) at baseline, generating
#   per-patient longitudinal CNA plots and an emergent-event summary table.
#
# Inputs:
#   - Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#       (integrated feature table; output of 2_0)
#   - cohort_assignment_table_updated.rds
#       (M4/SPORE/IMMAGINE cohort labels for patient filtering)
#
# Outputs:
#   - Subclonal_evolution_plots.pdf  (per-patient CNA track plots)
#   - Final Tables and Figures/Emergent_CNA_event_total.csv
#   - Final Tables and Figures/Emergent_CNA_event_gains.csv
#   - Final Tables and Figures/Emergent_CNA_loss_events.csv
#       (one row per patient, progression sample, and changed CNA arm)
#   - Final Tables and Figures/Subclonal_evolution_revision_test_cohort_*.csv
#       (cohort-expansion evaluability, calls, changes, and summary tables)
#   - Final Tables and Figures/Tumor_fraction_summary_*.csv
#       (auxiliary per-patient tumour-fraction summaries)
#
# Units of analysis:
#   - CNA event tables contain one row per patient, progression sample, and CNA
#     arm whose binary call changed relative to that patient's baseline sample.
#   - Per-patient plots contain all dated baseline-to-progression plasma samples
#     for patients with both baseline and progression cfDNA.
#   - Tumour-fraction summary tables contain one row per eligible patient.
#
# Dependencies:
#   tidyverse, lubridate
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/4_2_Compare_subclonal_evolution.R
#
# Manuscript outputs created/updated:
#   - Extended Data Figure 10A-B: emergent CNA/subclonal evolution support
#     tables and frozen final ichorCNA figure components when available.
#
# Pipeline role:
#   The current final figure panels use genome-wide ichorCNA visualizations
#   assembled outside this R script. This script regenerates the in-repository
#   event tables that identify which CNAs are newly detected at relapse, and
#   places those support tables beside the final figure artifacts.
#
# Manuscript-use boundary:
#   The emergent/gained/lost CNA tables are the manuscript-supporting outputs
#   for Extended Data Figure 10. Subclonal_evolution_plots.pdf and the
#   tumour-fraction timing tables are auxiliary
#   descriptive analyses; they are not the final ED10 panels.
#
# Author:    Dory Abelman
# Last update: August 2025
###############################################################################
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#

## ---- 0. USER SETTINGS -------------------------------------------------------
in_rds        <- "Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds"
out_plot_pdf  <- "Subclonal_evolution_plots.pdf"
out_events_csv<- "Emergent_CNA_events.csv"
cohort_df <- readRDS("cohort_assignment_table_updated.rds")
outdir <- "Final Tables and Figures/"

## ---- 1. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# Shared manuscript-output helpers.
# Extended Data Figure 10A/B uses external ichorCNA component plots assembled
# outside this R script. This script contributes the in-repo support
# CSVs and preserves the frozen final PDF when it is present locally.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

# Shared Spring 2026 revision metadata helpers. These keep the test-cohort
# sensitivity analysis aligned with the same metadata authority and exclusion
# table used by the rest of the revision-aware workflow.
.shared_helper <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.shared_helper)) {
  .shared_helper <- "helpers.R"
}
source(.shared_helper)
rm(.shared_helper)

## ---- 2. Load data -----------------------------------------------------------
all_df <- if (grepl("\\.rds$", in_rds, ignore.case = TRUE)) {
  readRDS(in_rds)
} else {
  read_csv(in_rds)          # assumes header row
}

event_cols <- c("del1p", "amp1q", "del13q", "del17p")

require_columns(
  all_df,
  c(
    "Patient", "Sample", "Sample_ID", "Sample_type", "Study",
    "timepoint_info", "Date_of_sample_collection",
    "Tumor_Fraction", "Num_days_to_closest_relapse",
    event_cols
  ),
  "Integrated feature table for subclonal CNA analysis"
)
require_columns(cohort_df, "Patient", "Cohort assignment table")

## ---- 3. Keep cfDNA & harmonise flags ----------------------------------------
cfDNA_df <- all_df %>%
  filter(str_detect(Sample_type, regex("plasma", ignore_case = TRUE))) %>%
  filter(Patient %in% cohort_df$Patient) %>%
  mutate(
    # Baseline ≈ diagnosis / T0
    is_baseline = timepoint_info %in% c("Diagnosis", "Baseline"),
    # Relapse / progression
    is_relapse  = timepoint_info %in% c("Relapse", "Progression"),
    Sample_Date = as_date(Date_of_sample_collection)
  ) 

## ---- 4. Keep patients with baseline *and* relapse cfDNA ---------------------
pts_keep <- cfDNA_df %>%
  group_by(Patient) %>%
  summarise(has_bl = any(is_baseline),
            has_rl = any(is_relapse),
            .groups = "drop") %>%
  filter(has_bl & has_rl) %>%
  pull(Patient)

cfDNA_df <- cfDNA_df %>% filter(Patient %in% pts_keep)
cfDNA_df <- cfDNA_df %>% mutate(enrollment_source = "Original cohort")

## ---- 5. Simple FGA proxy + numeric CNA flags --------------------------------
cfDNA_df <- cfDNA_df %>%
  mutate(across(all_of(event_cols), ~ as.numeric(.x))) %>%   # TRUE → 1, FALSE → 0
  rowwise() %>%
  mutate(FGA_proxy = sum(c_across(all_of(event_cols)), na.rm = TRUE) /
           length(event_cols)) %>%
  ungroup()

## ---- 6. Emergent CNA table --------------------------------------------------
# DEFINITION: An "emergent" CNA is one that is ABSENT (0) at baseline and
# PRESENT (1) at relapse.  The delta is computed per-patient as
#   relapse_value - baseline_value.
# A delta of +1 therefore captures 0->1 (newly acquired), whereas:
#   delta of  0 = no change (present at both or absent at both)
#   delta of -1 = lost at relapse (present at baseline, absent at relapse)
# Only delta == 1 events are retained as "Emergent" and written to CSV.
# NOTE: event_cols are binary (0/1) flags produced by 1_4_Process_CNA_Data.R.
emergent_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - .x[which(is_baseline)][1], .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(Event   = str_remove(Event, "^delta_"),
         Emergent = Delta == 1) %>%
  filter(Emergent) %>%
  select(Patient, Sample, Event) %>%
  arrange(Patient, Event)

edfig10_total_events_csv <- file.path(outdir, "Emergent_CNA_event_total.csv")
write_csv(emergent_tbl, edfig10_total_events_csv)


# Gains:
gains_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - first(.x), .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(
    Event    = str_remove(Event, "^delta_"),
    Emergent = Delta == 1
  ) %>%
  filter(Emergent) %>%
  select(Patient, Sample, Event)

edfig10_gain_events_csv <- file.path(outdir, "Emergent_CNA_event_gains.csv")
write_csv(gains_tbl, edfig10_gain_events_csv)

# Losses:
losses_tbl <- cfDNA_df %>%
  arrange(Patient, Sample_Date) %>%
  group_by(Patient) %>%
  mutate(across(all_of(event_cols),
                ~ .x - first(.x), .names = "delta_{col}")) %>%
  filter(is_relapse) %>%
  pivot_longer(starts_with("delta_"),
               names_to  = "Event",
               values_to = "Delta") %>%
  mutate(
    Event = str_remove(Event, "^delta_"),
    Lost  = Delta == -1
  ) %>%
  filter(Lost) %>%
  select(Patient, Sample, Event)

edfig10_loss_events_csv <- file.path(outdir, "Emergent_CNA_loss_events.csv")
write_csv(losses_tbl, edfig10_loss_events_csv)

## ---- 6B. Newly eligible test-cohort CNA analysis and QC ----------------------
# The Spring 2026 OICR revision samples expand the main subclonal-evolution
# cohort below. The revision-only files retained in this block are support
# and QC views; they are no longer a separate sensitivity denominator.
empty_revision_event_tbl <- tibble(
  Patient = character(),
  Sample_ID = character(),
  Sample = character(),
  Event = character(),
  baseline_value = numeric(),
  relapse_value = numeric(),
  Delta = numeric(),
  baseline_tumor_fraction = numeric(),
  followup_tumor_fraction = numeric(),
  Date_of_sample_collection = as.Date(character()),
  Num_days_to_closest_relapse = numeric()
)

revision_cfDNA_evaluable <- cfDNA_df[0, ]

make_revision_cna_change_table <- function(revision_cfDNA_df, event_cols) {
  if (nrow(revision_cfDNA_df) == 0) return(empty_revision_event_tbl)

  baseline_long <- revision_cfDNA_df %>%
    filter(is_baseline) %>%
    arrange(Patient, Sample_Date) %>%
    group_by(Patient) %>%
    slice(1) %>%
    ungroup() %>%
    select(Patient, baseline_tumor_fraction = Tumor_Fraction, all_of(event_cols)) %>%
    pivot_longer(
      cols = all_of(event_cols),
      names_to = "Event",
      values_to = "baseline_value"
    )

  relapse_long <- revision_cfDNA_df %>%
    filter(is_relapse) %>%
    arrange(Patient, Sample_Date) %>%
    select(
      Patient, Sample_ID, Sample, Date_of_sample_collection,
      Num_days_to_closest_relapse,
      followup_tumor_fraction = Tumor_Fraction,
      all_of(event_cols)
    ) %>%
    pivot_longer(
      cols = all_of(event_cols),
      names_to = "Event",
      values_to = "relapse_value"
    )

  relapse_long %>%
    inner_join(baseline_long, by = c("Patient", "Event")) %>%
    mutate(
      baseline_value = as.numeric(baseline_value),
      relapse_value = as.numeric(relapse_value),
      Delta = relapse_value - baseline_value
    ) %>%
    select(
      Patient, Sample_ID, Sample, Event, baseline_value, relapse_value, Delta,
      baseline_tumor_fraction, followup_tumor_fraction,
      Date_of_sample_collection, Num_days_to_closest_relapse
    ) %>%
    arrange(Patient, Sample_ID, Event)
}

revision_metadata <- load_spring2026_revision_metadata(required = FALSE)
if (!is.null(revision_metadata)) {
  require_columns(
    revision_metadata,
    c("Patient", "Sample_ID", "Study", "Sample_type", "timepoint_info"),
    "Spring 2026 revision metadata for subclonal CNA sensitivity analysis"
  )

  revision_patient_ids <- sort(unique(revision_metadata$Patient))
  revision_sample_ids <- unique(revision_metadata$Sample_ID)

  revision_all_df <- all_df %>%
    filter(
      Study == "IMMAGINE_revision_OICR",
      Patient %in% revision_patient_ids,
      Sample_ID %in% revision_sample_ids
    )

  revision_cfDNA_all <- revision_all_df %>%
    filter(str_detect(Sample_type, regex("plasma", ignore_case = TRUE))) %>%
    mutate(
      is_baseline = timepoint_info %in% c("Diagnosis", "Baseline"),
      is_relapse = timepoint_info %in% c("Relapse", "Progression"),
      Sample_Date = as_date(Date_of_sample_collection)
    )

  revision_evaluability <- tibble(Patient = revision_patient_ids) %>%
    left_join(
      revision_all_df %>%
        group_by(Patient) %>%
        summarise(
          feature_rds_rows = dplyr::n(),
          feature_rds_samples = n_distinct(Sample_ID),
          .groups = "drop"
        ),
      by = "Patient"
    ) %>%
    left_join(
      revision_cfDNA_all %>%
        group_by(Patient) %>%
        summarise(
          plasma_cfDNA_rows = dplyr::n(),
          baseline_cfDNA_rows = sum(is_baseline, na.rm = TRUE),
          relapse_cfDNA_rows = sum(is_relapse, na.rm = TRUE),
          maintenance_cfDNA_rows = sum(timepoint_info == "Maintenance", na.rm = TRUE),
          has_baseline_cfDNA = any(is_baseline),
          has_relapse_cfDNA = any(is_relapse),
          .groups = "drop"
        ),
      by = "Patient"
    ) %>%
    mutate(
      across(
        c(feature_rds_rows, feature_rds_samples, plasma_cfDNA_rows,
          baseline_cfDNA_rows, relapse_cfDNA_rows, maintenance_cfDNA_rows),
        ~ replace_na(.x, 0L)
      ),
      across(
        c(has_baseline_cfDNA, has_relapse_cfDNA),
        ~ replace_na(.x, FALSE)
      ),
      evaluable_for_revision_cna_subclonal_analysis =
        has_baseline_cfDNA & has_relapse_cfDNA,
      reason_not_evaluable = case_when(
        evaluable_for_revision_cna_subclonal_analysis ~ NA_character_,
        plasma_cfDNA_rows == 0 ~ "No plasma cfDNA rows in integrated feature table",
        !has_baseline_cfDNA & !has_relapse_cfDNA ~ "Missing both baseline and relapse/progression plasma cfDNA",
        !has_baseline_cfDNA ~ "Missing baseline plasma cfDNA",
        !has_relapse_cfDNA ~ "Missing relapse/progression plasma cfDNA",
        TRUE ~ "Not evaluable under current baseline/relapse cfDNA rules"
      )
    ) %>%
    arrange(Patient)

  revision_pts_keep <- revision_evaluability %>%
    filter(evaluable_for_revision_cna_subclonal_analysis) %>%
    pull(Patient)

  revision_cfDNA_evaluable <- revision_cfDNA_all %>%
    filter(Patient %in% revision_pts_keep) %>%
    mutate(enrollment_source = "Newly eligible test-cohort case") %>%
    mutate(across(all_of(event_cols), ~ as.numeric(.x))) %>%
    rowwise() %>%
    mutate(FGA_proxy = sum(c_across(all_of(event_cols)), na.rm = TRUE) /
             length(event_cols)) %>%
    ungroup()

  missing_revision_cna_calls <- revision_cfDNA_evaluable %>%
    filter(if_any(all_of(event_cols), is.na)) %>%
    distinct(Patient, Sample_ID)
  if (nrow(missing_revision_cna_calls) > 0L) {
    stop(
      "Revision patients classified as CNA-evaluable have missing target-arm calls: ",
      paste(
        paste0(missing_revision_cna_calls$Patient, "/",
               missing_revision_cna_calls$Sample_ID),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  revision_cna_calls <- revision_cfDNA_evaluable %>%
    select(
      Patient, Sample_ID, Sample, timepoint_info, Date_of_sample_collection,
      Num_days_to_closest_relapse, Tumor_Fraction, FGA_proxy, all_of(event_cols)
    ) %>%
    arrange(Patient, Date_of_sample_collection, Sample_ID)

  revision_cna_changes <- make_revision_cna_change_table(
    revision_cfDNA_evaluable,
    event_cols
  )

  revision_emergent_tbl <- revision_cna_changes %>%
    filter(Delta == 1) %>%
    arrange(Patient, Event, Sample_ID)

  revision_gain_events_tbl <- revision_emergent_tbl

  revision_loss_events_tbl <- revision_cna_changes %>%
    filter(Delta == -1) %>%
    arrange(Patient, Event, Sample_ID)

  revision_summary_tbl <- tibble(
    metric = c(
      "revision_metadata_rows_after_exclusions",
      "revision_metadata_patients_after_exclusions",
      "revision_rows_in_integrated_feature_rds",
      "revision_patients_in_integrated_feature_rds",
      "revision_plasma_cfDNA_rows_in_integrated_feature_rds",
      "revision_patients_evaluable_for_baseline_relapse_cfDNA_cna_change",
      "revision_relapse_or_progression_samples_evaluated",
      "revision_emergent_cna_event_rows",
      "revision_patients_with_emergent_cna_events",
      "revision_lost_cna_event_rows",
      "revision_patients_with_lost_cna_events",
      "revision_lost_call_followup_tumor_fraction_min",
      "revision_lost_call_followup_tumor_fraction_max"
    ),
    value = as.character(c(
      nrow(revision_metadata),
      length(revision_patient_ids),
      nrow(revision_all_df),
      n_distinct(revision_all_df$Patient),
      nrow(revision_cfDNA_all),
      length(revision_pts_keep),
      n_distinct(revision_cna_changes$Sample_ID),
      nrow(revision_emergent_tbl),
      n_distinct(revision_emergent_tbl$Patient),
      nrow(revision_loss_events_tbl),
      n_distinct(revision_loss_events_tbl$Patient),
      if (nrow(revision_loss_events_tbl) > 0L) {
        min(revision_loss_events_tbl$followup_tumor_fraction, na.rm = TRUE)
      } else {
        NA_real_
      },
      if (nrow(revision_loss_events_tbl) > 0L) {
        max(revision_loss_events_tbl$followup_tumor_fraction, na.rm = TRUE)
      } else {
        NA_real_
      }
    ))
  ) %>%
    bind_rows(
      tibble(
        metric = "interpretation",
        value = if (nrow(revision_emergent_tbl) == 0 &&
                    nrow(revision_loss_events_tbl) > 0) {
          paste0(
            "No emergent del1p/amp1q/del13q/del17p CNA events were detected ",
            "among revision patients evaluable under baseline-plus-relapse/",
            "progression plasma cfDNA rules. ",
            nrow(revision_loss_events_tbl),
            " apparent 1-to-0 transitions were observed across ",
            n_distinct(revision_loss_events_tbl$Patient),
            " patients at follow-up tumor fractions of ",
            sprintf(
              "%.1f%%-%.1f%%",
              100 * min(revision_loss_events_tbl$followup_tumor_fraction, na.rm = TRUE),
              100 * max(revision_loss_events_tbl$followup_tumor_fraction, na.rm = TRUE)
            ),
            "; these are reported as loss of detectability and do ",
            "not establish biological clonal loss."
          )
        } else if (nrow(revision_emergent_tbl) == 0) {
          "No emergent del1p/amp1q/del13q/del17p CNA events were detected among revision patients evaluable under baseline-plus-relapse/progression plasma cfDNA rules."
        } else {
          "Emergent del1p/amp1q/del13q/del17p CNA events were detected among revision patients evaluable under baseline-plus-relapse/progression plasma cfDNA rules."
        }
      )
    )

  revision_evaluability_csv <- file.path(
    outdir,
    "Subclonal_evolution_revision_test_cohort_evaluability.csv"
  )
  revision_cna_calls_csv <- file.path(
    outdir,
    "Subclonal_evolution_revision_test_cohort_cna_calls.csv"
  )
  revision_cna_changes_csv <- file.path(
    outdir,
    "Subclonal_evolution_revision_test_cohort_cna_changes.csv"
  )
  revision_summary_csv <- file.path(
    outdir,
    "Subclonal_evolution_revision_test_cohort_summary.csv"
  )
  revision_total_events_csv <- file.path(
    outdir,
    "Emergent_CNA_event_total_revision_test_cohort.csv"
  )
  revision_gain_events_csv <- file.path(
    outdir,
    "Emergent_CNA_event_gains_revision_test_cohort.csv"
  )
  revision_loss_events_csv <- file.path(
    outdir,
    "Emergent_CNA_loss_events_revision_test_cohort.csv"
  )

  write_csv(revision_evaluability, revision_evaluability_csv)
  write_csv(revision_cna_calls, revision_cna_calls_csv)
  write_csv(revision_cna_changes, revision_cna_changes_csv)
  write_csv(revision_summary_tbl, revision_summary_csv)
  write_csv(revision_emergent_tbl, revision_total_events_csv)
  write_csv(revision_gain_events_tbl, revision_gain_events_csv)
  write_csv(revision_loss_events_tbl, revision_loss_events_csv)

  ms_copy_artifact(
    source_path = revision_evaluability_csv,
    artifact_id = "EDFIG10A_B_A",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10 revision/test-cohort sensitivity: patient-level evaluability for baseline-plus-relapse/progression plasma cfDNA CNA-change analysis.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_cna_calls_csv,
    artifact_id = "EDFIG10A_B_A",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10 revision/test-cohort sensitivity: per-sample CNA calls for patients evaluable for baseline-plus-relapse/progression plasma cfDNA CNA-change analysis.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_cna_changes_csv,
    artifact_id = "EDFIG10A_B_A",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10 revision/test-cohort sensitivity: all baseline-to-relapse CNA deltas before filtering to emergent or lost events.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_summary_csv,
    artifact_id = "EDFIG10A_B_A",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10 revision/test-cohort sensitivity summary for CNA-change evaluability and event counts.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_total_events_csv,
    artifact_id = "EDFIG10A_B_A",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10A revision/test-cohort sensitivity: emergent CNA events among evaluable Spring 2026 OICR revision plasma pairs.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_gain_events_csv,
    artifact_id = "EDFIG10A_B_B",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10B revision/test-cohort sensitivity: gained CNA events among evaluable Spring 2026 OICR revision plasma pairs.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  ms_copy_artifact(
    source_path = revision_loss_events_csv,
    artifact_id = "EDFIG10A_B_B",
    role = "supporting_data_csv",
    description = "Extended Data Figure 10B revision/test-cohort sensitivity: lost CNA events among evaluable Spring 2026 OICR revision plasma pairs.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )

  message(
    "Revision/test-cohort CNA sensitivity analysis: ",
    length(revision_pts_keep), " evaluable patients, ",
    nrow(revision_emergent_tbl), " emergent CNA events, ",
    nrow(revision_loss_events_tbl), " lost CNA events."
  )
} else {
  warning(
    "Spring 2026 revision metadata was not available or revision integration is disabled; ",
    "revision/test-cohort subclonal CNA sensitivity outputs were not generated."
  )
}

## ---- 6C. Expand the main analysis cohort -----------------------------------
# Explicitly combine the original eligible cases with every newly eligible
# Spring 2026 test-cohort case. Enrollment source is retained for auditing,
# but all downstream plots, counts, and event tables use this unified cohort.
cfDNA_df <- bind_rows(cfDNA_df, revision_cfDNA_evaluable) %>%
  arrange(Patient, Sample_Date, Sample_ID) %>%
  distinct(Patient, Sample_ID, .keep_all = TRUE)

pts_keep <- sort(unique(cfDNA_df$Patient))

expanded_cohort_membership <- cfDNA_df %>%
  distinct(Patient, enrollment_source) %>%
  arrange(Patient)
write_csv(
  expanded_cohort_membership,
  file.path(outdir, "Subclonal_evolution_expanded_cohort_membership.csv")
)

# Recompute the main event tables after expansion. The baseline row is selected
# explicitly rather than assuming the first chronologically sorted row is a
# valid baseline specimen.
make_main_event_table <- function(dat, delta_to_keep) {
  dat %>%
    arrange(Patient, Sample_Date) %>%
    group_by(Patient) %>%
    mutate(
      across(
        all_of(event_cols),
        ~ .x - .x[which(is_baseline)[1]],
        .names = "delta_{col}"
      )
    ) %>%
    filter(is_relapse) %>%
    pivot_longer(
      starts_with("delta_"), names_to = "Event", values_to = "Delta"
    ) %>%
    mutate(Event = str_remove(Event, "^delta_")) %>%
    filter(Delta == delta_to_keep) %>%
    ungroup() %>%
    select(Patient, Sample, Event) %>%
    arrange(Patient, Event)
}

emergent_tbl <- make_main_event_table(cfDNA_df, 1)
gains_tbl <- emergent_tbl
losses_tbl <- make_main_event_table(cfDNA_df, -1)

write_csv(emergent_tbl, edfig10_total_events_csv)
write_csv(gains_tbl, edfig10_gain_events_csv)
write_csv(losses_tbl, edfig10_loss_events_csv)

message(
  "Expanded main subclonal-evolution cohort: ",
  n_distinct(cfDNA_df$Patient), " patients (",
  sum(expanded_cohort_membership$enrollment_source == "Original cohort"),
  " original + ",
  sum(expanded_cohort_membership$enrollment_source == "Newly eligible test-cohort case"),
  " newly eligible)."
)

# MANUSCRIPT OUTPUT: Extended Data Figure 10A/B support data
# The current final figure panels are externally assembled from ichorCNA
# genome-wide CNA plots. The CSVs generated here are the in-repo support
# tables that identify emergent/gained/lost CNA events used to support those
# panels. They are copied beside the final figure artifacts so readers can
# trace the figure back to the tabular calls regenerated by this script.
ms_copy_artifact(
  source_path = edfig10_total_events_csv,
  artifact_id = "EDFIG10A_B_A",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10A support data: emergent CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

ms_copy_artifact(
  source_path = edfig10_gain_events_csv,
  artifact_id = "EDFIG10A_B_B",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10B support data: gained CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

ms_copy_artifact(
  source_path = edfig10_loss_events_csv,
  artifact_id = "EDFIG10A_B_B",
  role = "supporting_data_csv",
  description = "Extended Data Figure 10B support data: lost CNA events in paired baseline/progression cfDNA samples.",
  script_name = "4_2_Compare_subclonal_evolution.R"
)

edfig10_final_pdf_candidates <- c(
  file.path("Manuscript_Exports", "02_extended_data_figures", "Extended_Data_Figure_10", "final_artifacts", "Extended_Data_Figure_10.pdf"),
  file.path("Figures_Exported", "Final_Feb2026", "Extended_Data_Figure_10.pdf"),
  file.path("reproducible_workflow", "outputs", "frozen", "extended_figures", "Extended_Data_Figure_10.pdf")
)
edfig10_final_pdf <- edfig10_final_pdf_candidates[file.exists(edfig10_final_pdf_candidates)][1]

if (!is.na(edfig10_final_pdf)) {
  ms_copy_artifact(
    source_path = edfig10_final_pdf,
    artifact_id = "EDFIG10A_B_A",
    role = "external_final_figure_pdf",
    description = "Extended Data Figure 10 final PDF: externally assembled ichorCNA component figure preserved as a frozen manuscript artifact.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )
  ms_copy_artifact(
    source_path = edfig10_final_pdf,
    artifact_id = "EDFIG10A_B_B",
    role = "external_final_figure_pdf",
    description = "Extended Data Figure 10 final PDF: externally assembled ichorCNA component figure preserved as a frozen manuscript artifact.",
    script_name = "4_2_Compare_subclonal_evolution.R"
  )
} else {
  warning(
    "No frozen Extended Data Figure 10 PDF was found. The in-repo support CSVs were generated, ",
    "but the final panel assembly depends on external ichorCNA component plots."
  )
}

## ---- 7. Per-patient plots (Tumour fraction + FGA proxy) ---------------------
plot_input <- cfDNA_df %>%
  filter(!is.na(Sample_Date), is.finite(as.numeric(Sample_Date)))

skipped_plot_rows <- cfDNA_df %>%
  filter(is.na(Sample_Date) | !is.finite(as.numeric(Sample_Date)))
if (nrow(skipped_plot_rows) > 0) {
  write_csv(skipped_plot_rows, file.path(outdir, "Subclonal_evolution_skipped_plot_rows_missing_dates.csv"))
}

plot_list <- plot_input %>%
  mutate(
    is_relapse = factor(is_relapse,
                        levels = c(FALSE, TRUE),
                        labels = c("Non-relapse", "Relapse"))
  ) %>%
  dplyr::group_nest(Patient) %>%        # <-- use dplyr::group_nest()
  mutate(
    plot = purrr::map2(
      data, Patient,
      ~ ggplot(.x, aes(x = Sample_Date)) +
        geom_line(aes(y = Tumor_Fraction),
                  colour = "steelblue", linewidth = 1) +
        geom_point(aes(y = Tumor_Fraction,
                       shape = is_relapse,
                       fill  = is_relapse),
                   colour = "black", size = 3) +
        geom_line(aes(y = FGA_proxy),
                  colour = "firebrick", linewidth = 1,
                  linetype = "dashed") +
        geom_point(aes(y = FGA_proxy),
                   colour = "firebrick", size = 3,
                   shape = 23, fill = "white") +
        scale_y_continuous(
          name     = "Tumour fraction (ichorCNA)",
          sec.axis = sec_axis(~., name = "FGA proxy (0-1)")
        ) +
        scale_shape_manual(values = c(21, 24)) +
        scale_fill_manual(values = c("grey80", "firebrick")) +
        labs(
          title = .y,                # patient_id now comes from .y
          x     = "Sample date"
        ) +
        theme_bw(base_size = 11) +
        theme(
          legend.position = "bottom",
          plot.title      = element_text(face = "bold")
        )
    )
  ) %>%
  filter(map_int(data, nrow) > 0)

pdf(out_plot_pdf, width = 7.5, height = 5.5)
if (nrow(plot_list) > 0) walk(plot_list$plot, print)
dev.off()


### Summarize sample availability and emergent-CNA timing
# ---- 8. Cohort and emergent-CNA numerical summaries -------------------------

# 1) Basic counts
df_keep        <- cfDNA_df
n_samples      <- nrow(df_keep)
patients       <- sort(unique(df_keep$Patient))
n_patients     <- length(patients)

# 2) Samples per patient
spp            <- df_keep %>%
  group_by(Patient) %>%
  summarise(n = dplyr::n(), .groups = "drop")
avg_spp        <- median(spp$n, na.rm = TRUE)
range_spp      <- range(spp$n, na.rm = TRUE)

# 3) Emergent‐CNA patients
emergent_pts   <- emergent_tbl %>% distinct(Patient) %>% pull(Patient)
n_emergent     <- length(emergent_pts)
pct_emergent   <- n_emergent / n_patients * 100

# 4) Days before progression for emergent events
#    (we look up Num_days_to_closest_relapse from cfDNA_df for those samples)
days_vec       <- cfDNA_df %>%
  filter(Sample %in% emergent_tbl$Sample) %>%
  pull(Num_days_to_closest_relapse)

mean_days      <- mean(days_vec, na.rm = TRUE)
range_days     <- range(days_vec, na.rm = TRUE)
iqr_days       <- IQR(days_vec, na.rm = TRUE)



## Additional summaries

# 1) Add days-to-progression to emergent_tbl
emergent_tbl_days <- emergent_tbl %>%
  left_join(
    cfDNA_df %>% 
      select(Patient, Sample, Num_days_to_closest_relapse),
    by = c("Patient", "Sample")
  )

# 2) Overall summary
overall_days <- emergent_tbl_days$Num_days_to_closest_relapse
overall_summary <- tibble(
  n_events   = nrow(emergent_tbl_days),
  mean_days  = mean(overall_days, na.rm = TRUE),
  median_days= median(overall_days, na.rm = TRUE),
  range_days = paste(range(overall_days, na.rm = TRUE), collapse = "–"),
  iqr_days   = IQR(overall_days, na.rm = TRUE)
)

print(overall_summary)

# 3) Per–event type summary (optional)
per_event_summary <- emergent_tbl_days %>%
  group_by(Event) %>%
  summarise(
    n_events    = dplyr::n(),
    mean_days   = mean(Num_days_to_closest_relapse, na.rm = TRUE),
    median_days = median(Num_days_to_closest_relapse, na.rm = TRUE),
    range_days  = paste(range(Num_days_to_closest_relapse, na.rm = TRUE), collapse = "–"),
    iqr_days    = IQR(Num_days_to_closest_relapse, na.rm = TRUE),
    .groups     = "drop"
  )

print(per_event_summary)


### See second closest draw time

# 1) Identify patients of interest (e.g. those with emergent events)
patients_of_interest <- emergent_tbl %>%
  distinct(Patient) %>%
  pull(Patient)

# 2) For each patient, find the second‐closest cfDNA draw to progression
second_closest <- cfDNA_df %>%
  filter(Patient %in% patients_of_interest) %>%
  # keep only non‐relapse draws (so relapse itself with zero days doesn’t count)
  filter(Num_days_to_closest_relapse > 0) %>%
  group_by(Patient) %>%
  filter(timepoint_info != "Relapse") %>%
  arrange(Num_days_to_closest_relapse) %>%
  slice(1) %>%  # closest distance
  ungroup() %>%
  select(
    Patient,
    Sample_ID,
    Sample_Date       = Date_of_sample_collection,
    Days_before_relapse = Num_days_to_closest_relapse,
    Timepoint         = timepoint_info
  )

# Compute median, range, and IQR:
second_closest %>%
  summarise(
    median_days = median(Days_before_relapse, na.rm = TRUE),
    range_days  = paste(range(Days_before_relapse, na.rm = TRUE), collapse = "–"),
    iqr_days    = IQR(Days_before_relapse, na.rm = TRUE)
  ) %>%
  print()

# If you prefer separate base‐R values:
days <- second_closest$Days_before_relapse
cat("Median:", median(days, na.rm = TRUE), "\n")
cat("Range:", paste(range(days, na.rm = TRUE), collapse = "–"), "\n")
cat("IQR:", IQR(days, na.rm = TRUE), "\n")

print(second_closest)



### Now see how early tumor fraction rose 

# 1) Patients who progressed (same as before)
pts_prog <- emergent_tbl %>% distinct(Patient) %>% pull(Patient)

# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # Nadir is restricted to valid samples collected on/before progression.
    tf_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) min(Tumor_Fraction[idx]) else NA_real_
    },
    days_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) {
        Num_days_to_closest_relapse[idx[which.min(Tumor_Fraction[idx])]]
      } else {
        NA_real_
      }
    },
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)




# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)



write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_new_CNA_patients.csv"))


## Do now for just those who didn't show high risk CNAs 
# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(!Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # Nadir is restricted to valid samples collected on/before progression.
    tf_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) min(Tumor_Fraction[idx]) else NA_real_
    },
    days_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) {
        Num_days_to_closest_relapse[idx[which.min(Tumor_Fraction[idx])]]
      } else {
        NA_real_
      }
    },
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)




# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)



write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_new_CNA_patients.csv"))




### redo for all patients 
### Now see how early tumor fraction rose 

# 1) Patients who progressed (same as before)
pts_prog <- cfDNA_df %>% distinct(Patient) %>% pull(Patient)

# 2) Per‐patient TF metrics
tf_summary <- cfDNA_df %>%
  filter(Patient %in% pts_prog) %>%
  group_by(Patient) %>%
  summarise(
    # baseline TF (first non‐relapse “Diagnosis” sample)
    tf_baseline = Tumor_Fraction[which(is_baseline & !is_relapse)][1],
    # relapse TF (first relapse sample)
    tf_relapse  = Tumor_Fraction[which(is_relapse)][1],
    # Nadir is restricted to valid samples collected on/before progression.
    tf_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) min(Tumor_Fraction[idx]) else NA_real_
    },
    days_nadir = {
      idx <- which(
        !is_relapse & Num_days_to_closest_relapse >= 0 &
          is.finite(Tumor_Fraction)
      )
      if (length(idx)) {
        Num_days_to_closest_relapse[idx[which.min(Tumor_Fraction[idx])]]
      } else {
        NA_real_
      }
    },
    # magnitude of rise
    tf_rise     = tf_relapse - tf_nadir,
    .groups     = "drop"
  )

# 3) Compute medians, IQRs, and ranges
baseline_q <- quantile(tf_summary$tf_baseline, probs = c(0.25, 0.75), na.rm = TRUE) * 100
relapse_q  <- quantile(tf_summary$tf_relapse,  probs = c(0.25, 0.75), na.rm = TRUE) * 100
rise_q     <- quantile(tf_summary$tf_rise,     probs = c(0.25, 0.75), na.rm = TRUE) * 100
days_q     <- quantile(tf_summary$days_nadir,  probs = c(0.25, 0.75), na.rm = TRUE)



# 3) Compute medians and ranges
baseline_range <- range(tf_summary$tf_baseline, na.rm = TRUE) * 100
relapse_range  <- range(tf_summary$tf_relapse,  na.rm = TRUE) * 100
rise_range     <- range(tf_summary$tf_rise,     na.rm = TRUE) * 100
days_range     <- range(tf_summary$days_nadir,  na.rm = TRUE)




write_csv(tf_summary, file.path(outdir, "Tumor_fraction_summary_all_CNA_patients.csv"))
