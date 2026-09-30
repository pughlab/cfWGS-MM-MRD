# ==============================================================================
# 2_3_Feature_Concordance_And_Mutation_Counts.R
#
# Purpose:
#   1. Select one diagnosis/baseline analysis record per patient, including the
#      documented exceptions and modality-specific source rows used below.
#   2. Compute FISH ↔ WGS concordance (overall and by ctDNA fraction) for
#      copy-number alterations and IGH translocations.
#   3. Compute matched BM ↔ cfDNA mutation-set concordance and summarise
#      baseline mutation counts by cohort.
#   4. Calculate Spearman correlations between mutation burden and
#      clinical/fragmentomic features and fit the exploratory BM count model.
#   5. Generate figures used in the manuscript:
#        • Boxplots of mutation counts and cfDNA tumor fraction by cohort
#        • Scatterplots of mutation burden vs. tumor fraction, fragment-size score, albumin
#        • Dumbbell plots of event-level concordance, sensitivity & specificity by ctDNA fraction
#
# Inputs:
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds
#   - Output_tables_2025_updated/patient_mutation_counts.csv (from script 2_2)
#   - Jan2025_exported_data/CNA_at_FISH_sites_combined.rds
#   - Jan2025_exported_data/Sample_ploidy_from_sequenza_400.rds
#   - Jan2025_exported_data/mutation_export_updated_more_info2.rds
#   - Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#   - combined_clinical_data_updated_April2025.csv
#   - Output_tables_2025_updated/merged_mut.rds, merged_trans.rds, and
#     merged_CNA.rds for the detailed BM-versus-cfDNA performance section
#   - Cohort assignments loaded through load_final_cohort_assignment()
#
# Outputs:
#   - FISH/WGS concordance tables in Output_tables_2025_updated/
#   - Concordance components for Supplementary Table 2 and the feature-
#     correlation CSV for Supplementary Table 3 in Final Tables and Figures/
#   - Figures in Final Tables and Figures/Baseline_concordance/
#   - R objects (RDS) for downstream concordance/source-data reuse
#   - Source tables supporting Extended Data Figure 2 panels A-C and E-F.
#     The numerical source for the manually assembled panel D is exported by
#     5_1_Export_Locked_Figure_Source_Data.R.
#
# Required packages:
#   tidyverse, purrr, stringr, writexl, glue, Hmisc, broom,
#   ggpubr, patchwork, viridis, scales
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_3_Feature_Concordance_And_Mutation_Counts.R
#
# Manuscript outputs created/updated:
#   - Extended Data Figure 2A-C and 2E-F: baseline concordance, mutation burden,
#     and feature-correlation panels. The script can stage the retained full
#     Extended Data Figure 2 PDF, but it does not draw panel D.
#   - Supplementary Table 2: this script supplies the BM/cfDNA performance,
#     FISH/WGS agreement, FISH-probe, and per-sample-call components. The
#     delivered six-sheet workbook also contains the feature catalogue and
#     per-variant VAF sheets assembled from upstream outputs.
#   - Supplementary Table 3: baseline mutation-count and feature-correlation
#     source table.
#
# Units of analysis:
#   - Baseline clinical/mutation summaries: one selected row per patient.
#     The documented CA-02 consolidation and IMG-142/IMG-235 close-date pairing
#     can fill complementary fields from more than one visit, so these are
#     integrated analysis records and not always single physical specimens.
#   - FISH/WGS concordance: one evaluable patient × event × WGS source row.
#   - Mutation-set concordance: one matched patient × baseline timepoint row;
#     confusion counts are over unique mutation identities within that row.
#     Its true-negative count is defined relative to the complete observed
#     mutation-identity universe in this section. A patient with no mutation row
#     in either compartment cannot enter this mutation-row-derived comparison.
#   - Feature correlations: pairwise-complete selected baseline patient rows;
#     n_pairs is therefore allowed to differ across variable pairs.
#
# Final tables used in the manuscript:
#   The cleanly numbered tables are stored in
#   Final docs/Final Tables and Supplementary Tables. Historical filenames in
#   this script are retained as provenance and intermediate compatibility
#   outputs; table status should not be inferred from those names alone.
#
# Pipeline role:
#   This script quantifies how well cfDNA and BM WGS recover clinically reported
#   FISH/cytogenetic events, then summarizes mutation burden and its association
#   with tumour fraction and fragmentomic features. Concordance is stratified by
#   ctDNA fraction because low tumour fraction is a known biological and technical
#   limit on detecting CNAs and structural variants in plasma.
#
# ==============================================================================
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#
# Analyst note:
#   BAM archive/unarchive helper tables are operational diagnostics only. They
#   are skipped by default so routine manuscript regeneration does not require
#   storage-location spreadsheets or create root-level helper files. To regenerate
#   those support files intentionally, run with:
#     CFWGS_RUN_BAM_ARCHIVE_DIAGNOSTICS=true Rscript Scripts_2025/Final_Scripts/2_3_Feature_Concordance_And_Mutation_Counts.R
#

library(tidyverse)   # dplyr, tidyr, readr, etc.
library(purrr)       # for pmap_dfr
library(stringr)     # str_detect, str_to_lower, etc.
library(writexl)     # write_xlsx()
library(glue)        # for building sentence output
library(Hmisc)    # for rcorr()
library(broom)    # for tidy()
library(purrr)    # for map_df()
library(tibble)   # for tibble()
library(readxl)   # read_excel() for BAM storage helper tables
library(ggpubr)   # stat_compare_means() for baseline concordance plots
library(scales)   # percent_format() for plot axes
library(patchwork) # plot_spacer() and plot layouts
library(viridis)  # viridis() colors for concordance plots

# Shared helper for final manuscript-organized outputs.
# The script continues to write its historical outputs. The helper additionally
# copies each final manuscript component into
# Scripts_2025/Final_Scripts/final_manuscript_objects with labels such as
# Extended_Data_Figure_2A and Supplementary_Table_3.
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

file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")



##### PART 1: See concordance to FISH 
## 1.  PARAMETERS  ---------------------------------------------------
## ------------------------------------------------------------------
# tf_cut is the fixed 5% ctDNA-fraction boundary used for the stratified
# concordance summaries. FISH/CNA sections classify exactly 0.05 as high TF
# (`>= tf_cut`). The earlier mutation-set section retains its historical
# `> tf_cut` comparison, so a sample exactly at 0.05 would be classified
# differently between those sections. This script documents but does not alter
# that established behavior.
tf_cut    <- 0.05            # FISH/CNA high-TF stratum uses >= 0.05
baseline  <- c("Diagnosis","Baseline")   # recognise baseline labels
run_bam_archive_diagnostics <- tolower(Sys.getenv("CFWGS_RUN_BAM_ARCHIVE_DIAGNOSTICS", "false")) %in%
  c("true", "t", "1", "yes", "y")

## ------------------------------------------------------------------
## 2.  STARTING DATA  ------------------------------------------------
## ------------------------------------------------------------------
dat <- file   # <- the tibble

## ------------------------------------------------------------------
## 3.  KEEP BASELINE SAMPLES AND ENSURE NO DUPLICATES ---------------
## ------------------------------------------------------------------
# 1. Subset to only Diagnosis or Baseline timepoints
dat_tb <- dat %>% 
  filter(
    str_to_lower(timepoint_info) %in% str_to_lower(baseline)
  )


# 2. Check for duplicate rows per patient in this subset
dup_patients <- dat_tb %>%
  count(Patient) %>%
  filter(n > 1) %>%
  pull(Patient)

# ensure Date is Date class
dat_tb <- dat_tb %>%
  mutate(Date = as.Date(Date))

# 1) Remove CA-03 timepoint 02
dat_tb2 <- dat_tb %>%
  filter(!(Patient == "CA-03" & timepoint_info == "02"))

# 2) Consolidate the two CA-02 rows
resp_CA02 <- dat_tb2 %>%
  filter(Patient == "CA-02") %>%
  # order so timepoint “01” comes before “02”
  arrange(factor(timepoint_info, levels = c("01","02"))) %>%
  summarise(across(everything(), ~{
    vals <- .
    # first non-NA in order
    first_val <- vals[which(!is.na(vals))[1]]
    # any other non-NA ≠ first_val
    other    <- vals[!is.na(vals) & vals != first_val][1]
    # if first_val is Unknown/Other but other is a “real” value, use other
    if (!is.na(first_val) &&
        first_val %in% c("Unknown","Other") &&
        !is.na(other) &&
        !other %in% c("Unknown","Other")) {
      other
    } else {
      first_val
    }
  }))

# 3) Drop all CA-02 originals
dat_tb3 <- dat_tb2 %>% filter(Patient != "CA-02")

# 4) Apply audited manual baseline-row selection for duplicate candidates.
dat_tb4 <- dat_tb3 %>%
  filter_manual_baseline_row_selection("spore0009_baseline_only")

# 5) Re-bind the collapsed CA-02 row
dat_tb_final <- bind_rows(dat_tb4, resp_CA02) %>%
  arrange(Patient)

# Preserve the baseline/diagnosis cfDNA feature row independently of the row
# later selected for BM-informed patient-level concordance. A patient can have
# modality-specific baseline sources at different visits. In particular,
# SPORE_0012 contributes its baseline cfDNA mutation catalogue, ichorCNA tumour
# fraction, and fragment-size score at T1, whereas T4 is the manually selected
# BM baseline. Collapsing to the T4 row before retaining these cfDNA features
# previously paired the T1 mutation count with missing T4 cfDNA covariates.
baseline_cfDNA_feature_catalogue <- dat_tb_final %>%
  filter(!is.na(.data$Blood_Mutation_Count)) %>%
  transmute(
    Patient = .data$Patient,
    cfDNA_feature_Sample_Code = as.character(.data$Sample_Code),
    cfDNA_feature_Timepoint = as.character(.data$Timepoint),
    cfDNA_feature_Date = as.Date(.data$Date),
    cfDNA_feature_mutation_count_original = as.integer(.data$Blood_Mutation_Count),
    cfDNA_feature_tumour_fraction = as.numeric(
      .data$WGS_Tumor_Fraction_Blood_plasma_cfDNA
    ),
    cfDNA_feature_FS = as.numeric(.data$FS)
  )

duplicate_cfDNA_feature_sources <- baseline_cfDNA_feature_catalogue %>%
  count(.data$Patient, name = "n_cfDNA_feature_sources") %>%
  filter(.data$n_cfDNA_feature_sources != 1L)
if (nrow(duplicate_cfDNA_feature_sources)) {
  stop(
    "Baseline cfDNA feature source is not unique by patient: ",
    paste(duplicate_cfDNA_feature_sources$Patient, collapse = ", "),
    call. = FALSE
  )
}

dat_base <- dat_tb_final 

## Remove dup
dat_base <- dat_base %>%
  filter(!(Patient == "CA-03" & Timepoint == "02"))

# Preserve all existing one-row baseline patients and add complementary BM and
# cfDNA fields from a second baseline/diagnosis row only when the two specimen
# dates are within 30 days. This resolves T0/T1 label mismatches for IMG-142
# (5 days) and IMG-235 (7 days) without removing or redefining any patient that
# was already represented by a complete single baseline row.
baseline_rows_for_date_pairing <- dat_base %>%
  mutate(
    .baseline_row_id = row_number(),
    Date = as.Date(Date)
  )

additive_baseline_date_pairs <- inner_join(
  baseline_rows_for_date_pairing %>%
    filter(!is.na(BM_Mutation_Count), !is.na(Date)) %>%
    transmute(
      Patient,
      BM_row_id = .baseline_row_id,
      BM_Sample_Code = as.character(Sample_Code),
      BM_Timepoint = as.character(Timepoint),
      BM_date = Date
    ),
  baseline_rows_for_date_pairing %>%
    filter(!is.na(Blood_Mutation_Count), !is.na(Date)) %>%
    transmute(
      Patient,
      cfDNA_row_id = .baseline_row_id,
      cfDNA_Sample_Code = as.character(Sample_Code),
      cfDNA_Timepoint = as.character(Timepoint),
      cfDNA_date = Date
    ),
  by = "Patient"
) %>%
  filter(BM_row_id != cfDNA_row_id) %>%
  mutate(abs_days = abs(as.integer(BM_date - cfDNA_date))) %>%
  filter(abs_days <= 30L) %>%
  arrange(Patient, abs_days, BM_date, cfDNA_date, BM_row_id, cfDNA_row_id) %>%
  group_by(Patient) %>%
  slice_head(n = 1L) %>%
  ungroup()

expected_date_added_patients <- c("IMG-142", "IMG-235")
missing_expected_date_additions <- setdiff(
  expected_date_added_patients,
  additive_baseline_date_pairs$Patient
)
if (length(missing_expected_date_additions)) {
  stop(
    "Expected <=30-day complementary baseline pair(s) were not recovered: ",
    paste(missing_expected_date_additions, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(additive_baseline_date_pairs)) {
  collapse_pair_fill_missing <- function(pair_row) {
    pair_row <- as.list(pair_row)
    pair_ids <- c(pair_row$BM_row_id, pair_row$cfDNA_row_id)
    pair_data <- baseline_rows_for_date_pairing %>%
      filter(.baseline_row_id %in% pair_ids) %>%
      arrange(Date, .baseline_row_id)
    primary <- pair_data[1, , drop = FALSE]
    secondary <- pair_data[2, , drop = FALSE]

    for (column_name in setdiff(names(primary), ".baseline_row_id")) {
      primary_value <- primary[[column_name]]
      secondary_value <- secondary[[column_name]]
      primary_missing <- is.na(primary_value)
      if (is.character(primary_value)) {
        primary_missing <- primary_missing | !nzchar(primary_value)
      }
      primary[[column_name]][primary_missing] <- secondary_value[primary_missing]
    }
    primary
  }

  collapsed_date_pair_rows <- purrr::map_dfr(
    seq_len(nrow(additive_baseline_date_pairs)),
    ~ collapse_pair_fill_missing(additive_baseline_date_pairs[.x, , drop = FALSE])
  )

  dat_base <- bind_rows(
    baseline_rows_for_date_pairing %>%
      filter(!Patient %in% additive_baseline_date_pairs$Patient),
    collapsed_date_pair_rows
  ) %>%
    select(-.baseline_row_id) %>%
    arrange(Patient, Date, Timepoint, Sample_Code)
} else {
  dat_base <- baseline_rows_for_date_pairing %>% select(-.baseline_row_id)
}

dir.create(file.path("Output_tables_2025", "clinical_support"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(
  additive_baseline_date_pairs %>%
    mutate(pairing_rule = "add_to_existing_patient_set_if_within_30_days"),
  file.path(
    "Output_tables_2025",
    "clinical_support",
    "feature_concordance_additive_baseline_pairs_within_30d_audit.csv"
  )
)

baseline_duplicate_audit <- dat_base %>%
  group_by(Patient) %>%
  mutate(n_patient_baseline_candidates = n()) %>%
  ungroup() %>%
  filter(n_patient_baseline_candidates > 1) %>%
  mutate(
    manually_designated_baseline = .data$Patient == "SPORE_0012" &
      as.character(.data$Sample_Code) == "SPORE_0012_T4" &
      as.character(.data$Timepoint) == "4" &
      as.character(.data$timepoint_info) == "Baseline",
    patient_baseline_timepoint_rank = case_when(
      .data$manually_designated_baseline ~ 0L,
      str_detect(as.character(Timepoint), regex("^T?0$", ignore_case = TRUE)) ~ 0L,
      str_detect(as.character(Timepoint), regex("^T?1$|^0?1$", ignore_case = TRUE)) ~ 1L,
      TRUE ~ 2L
    ),
    patient_baseline_has_feature_evidence = !is.na(BM_Mutation_Count) | !is.na(Blood_Mutation_Count)
  ) %>%
  group_by(Patient) %>%
  arrange(
    desc(manually_designated_baseline),
    is.na(Date),
    Date,
    patient_baseline_timepoint_rank,
    desc(patient_baseline_has_feature_evidence),
    Timepoint,
    .by_group = TRUE
  ) %>%
  mutate(
    selected_for_patient_baseline_concordance = row_number() == 1L,
    patient_baseline_selection_reason = case_when(
      selected_for_patient_baseline_concordance & manually_designated_baseline ~
        "selected manually designated SPORE_0012 T4 BM baseline for BM-informed baseline analyses",
      selected_for_patient_baseline_concordance ~ "selected earliest dated baseline/diagnosis candidate, preferring T0/T1 and rows with WGS mutation evidence",
      TRUE ~ "not selected for patient-level baseline concordance to enforce one baseline/diagnosis row per patient"
    )
  ) %>%
  ungroup()

if (nrow(baseline_duplicate_audit) > 0L) {
  dir.create(file.path("Output_tables_2025", "clinical_support"), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(
    baseline_duplicate_audit,
    file.path("Output_tables_2025", "clinical_support", "feature_concordance_patient_baseline_duplicate_selection_audit.csv")
  )
  dat_base <- dat_base %>%
    left_join(
      baseline_duplicate_audit %>%
        select(Patient, Sample_Code, Timepoint, selected_for_patient_baseline_concordance),
      by = c("Patient", "Sample_Code", "Timepoint")
    ) %>%
    filter(is.na(selected_for_patient_baseline_concordance) | selected_for_patient_baseline_concordance) %>%
    select(-selected_for_patient_baseline_concordance)
}

# 6) Quick duplicate check
dups <- dat_base %>%
  count(Patient) %>%
  filter(n > 1)

if (nrow(dups)) {
  warning("Still multiple Diagnosis/Baseline rows for: ",
          paste(unique(dups$Patient), collapse = ", "))
} else {
  message("All patients now have at most one Diagnosis/Baseline row.")
}



### Filter to the ones interested in 
## Pull from previous export 
cohort_df <- load_final_cohort_assignment()
keep_patients <- cohort_df$Patient

## Keep only interested patients 
dat_base <- dat_base %>% filter(Patient %in% keep_patients)


# -----------------------------------------------------------
# 2.  DEFINE COHORTS  ------------
# -----------------------------------------------------------
dat_base <- dat_base %>%
  left_join(cohort_df, by = "Patient")


dat_base <- dat_base %>%
  mutate(Cohort = case_when(
    Cohort == "Frontline"     ~ "Frontline induction-transplant",
    TRUE                      ~ Cohort
  ))

dat_base$cohort <- dat_base$Cohort ## for consistency 

# Keep mutation-burden summaries and ED2E/F consistent with the regenerated
# heatmap and ED2G catalogues. The integrated aggregate contains legacy count
# fields, including several cfDNA values counted before rsID removal. The
# heatmap export is the canonical patient-level no-rsID count table and uses the
# exact ED2G catalogue for evaluable paired baseline patients.
mutation_count_path <- file.path(
  "Output_tables_2025_updated", "patient_mutation_counts.csv"
)
if (!file.exists(mutation_count_path)) {
  stop(
    "Missing canonical baseline mutation-count table: ", mutation_count_path,
    ". Run 2_2_Baseline_demographics_by_WGS_heatmap_updated.R first.",
    call. = FALSE
  )
}
canonical_baseline_mutation_counts <- readr::read_csv(
  mutation_count_path,
  show_col_types = FALSE
) %>%
  select(
    .data$Patient,
    canonical_BM_Mutation_Count = .data$BM_Mutation_Count,
    canonical_Blood_Mutation_Count = .data$Blood_Mutation_Count,
    .data$BM_Mutation_Count_Source,
    .data$Blood_Mutation_Count_Source
  ) %>%
  distinct()
duplicate_canonical_mutation_counts <- canonical_baseline_mutation_counts %>%
  count(.data$Patient) %>%
  filter(.data$n > 1L)
if (nrow(duplicate_canonical_mutation_counts)) {
  stop(
    "Canonical baseline mutation-count table is not unique by patient: ",
    paste(duplicate_canonical_mutation_counts$Patient, collapse = ", "),
    call. = FALSE
  )
}
dat_base <- dat_base %>%
  left_join(canonical_baseline_mutation_counts, by = "Patient") %>%
  mutate(
    BM_Mutation_Count = coalesce(
      as.integer(.data$canonical_BM_Mutation_Count),
      as.integer(.data$BM_Mutation_Count)
    ),
    Blood_Mutation_Count = coalesce(
      as.integer(.data$canonical_Blood_Mutation_Count),
      as.integer(.data$Blood_Mutation_Count)
    )
  ) %>%
  select(
    -.data$canonical_BM_Mutation_Count,
    -.data$canonical_Blood_Mutation_Count
  )

# Align cfDNA-specific covariates to the same baseline/diagnosis cfDNA source
# that supplied the mutation catalogue, even when a different row is retained
# as the patient's BM baseline. Prefer the modality-specific source values;
# fall back to the selected patient-level row only when no cfDNA source row is
# available in the baseline candidate table.
dat_base <- dat_base %>%
  rename(
    selected_row_cfDNA_tumour_fraction =
      .data$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    selected_row_cfDNA_FS = .data$FS
  ) %>%
  left_join(baseline_cfDNA_feature_catalogue, by = "Patient") %>%
  mutate(
    WGS_Tumor_Fraction_Blood_plasma_cfDNA = coalesce(
      .data$cfDNA_feature_tumour_fraction,
      .data$selected_row_cfDNA_tumour_fraction
    ),
    FS = coalesce(
      .data$cfDNA_feature_FS,
      .data$selected_row_cfDNA_FS
    )
  )

cfDNA_feature_alignment_audit <- dat_base %>%
  select(
    .data$Patient,
    selected_patient_row_Sample_Code = .data$Sample_Code,
    selected_patient_row_Date = .data$Date,
    .data$cfDNA_feature_Sample_Code,
    .data$cfDNA_feature_Timepoint,
    .data$cfDNA_feature_Date,
    .data$cfDNA_feature_mutation_count_original,
    canonical_cfDNA_mutation_count = .data$Blood_Mutation_Count,
    .data$selected_row_cfDNA_tumour_fraction,
    aligned_cfDNA_tumour_fraction =
      .data$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    .data$selected_row_cfDNA_FS,
    aligned_cfDNA_FS = .data$FS
  )

readr::write_csv(
  cfDNA_feature_alignment_audit,
  file.path(
    "Output_tables_2025", "clinical_support",
    "feature_concordance_cfDNA_baseline_feature_alignment_audit.csv"
  )
)

spore0012_cfDNA_alignment <- cfDNA_feature_alignment_audit %>%
  filter(.data$Patient == "SPORE_0012")
if (
  nrow(spore0012_cfDNA_alignment) != 1L ||
    spore0012_cfDNA_alignment$cfDNA_feature_Sample_Code != "SPORE_0012_T1" ||
    spore0012_cfDNA_alignment$canonical_cfDNA_mutation_count != 3793L ||
    spore0012_cfDNA_alignment$aligned_cfDNA_tumour_fraction != 0
) {
  stop(
    "SPORE_0012 cfDNA baseline alignment failed; expected T1, ",
    "3,793 mutations, and ichorCNA tumour fraction 0.",
    call. = FALSE
  )
}

cfDNA_mutation_tf_missing <- dat_base %>%
  filter(
    !is.na(.data$Blood_Mutation_Count),
    is.na(.data$WGS_Tumor_Fraction_Blood_plasma_cfDNA)
  ) %>%
  pull(.data$Patient)
if (length(cfDNA_mutation_tf_missing)) {
  stop(
    "Baseline cfDNA mutation count lacks a matched ichorCNA tumour fraction for: ",
    paste(cfDNA_mutation_tf_missing, collapse = ", "),
    call. = FALSE
  )
}

img181_count_check <- dat_base %>%
  filter(.data$Patient == "IMG-181") %>%
  select(.data$Patient, .data$BM_Mutation_Count, .data$Blood_Mutation_Count)
if (
  nrow(img181_count_check) != 1L ||
    img181_count_check$BM_Mutation_Count != 2505L ||
    img181_count_check$Blood_Mutation_Count != 888L
) {
  stop(
    "IMG-181 canonical baseline mutation-count check failed; expected BM=2505 and cfDNA=888.",
    call. = FALSE
  )
}

# Recover baseline FISH calls that can be lost when the integrated aggregate is
# rebuilt from clinical metadata without assay-specific FISH columns.  Calls in
# the aggregate are used; the deeper IMMAGINE patient/master tables
# only fill missing values.  The audit and hard checks prevent the known
# IMG-060/IMG-098/IMG-181 positive events from silently disappearing again.
fish_fields <- c("T_4_14", "T_11_14", "T_14_16", "DEL_17P", "DEL_1P", "AMP_1Q")
fish_before <- dat_base %>%
  dplyr::select(.data$Patient, dplyr::all_of(fish_fields))
fish_catalogue <- load_revision_inclusive_baseline_fish_calls(file)
fish_call_sources <- attr(fish_catalogue, "call_sources")

dat_base <- dat_base %>%
  dplyr::select(-dplyr::all_of(fish_fields)) %>%
  dplyr::left_join(fish_catalogue, by = "Patient")

fish_recovery_audit <- fish_before %>%
  tidyr::pivot_longer(-.data$Patient, names_to = "feature", values_to = "call_before") %>%
  dplyr::left_join(
    dat_base %>%
      dplyr::select(.data$Patient, dplyr::all_of(fish_fields)) %>%
      tidyr::pivot_longer(-.data$Patient, names_to = "feature", values_to = "call_after"),
    by = c("Patient", "feature")
  ) %>%
  dplyr::left_join(fish_call_sources, by = c("Patient", "feature", "call_after" = "call")) %>%
  dplyr::filter(is.na(.data$call_before) & !is.na(.data$call_after)) %>%
  dplyr::arrange(.data$Patient, .data$feature)

dir.create("Output_tables_2025/feature_concordance_support", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  fish_recovery_audit,
  "Output_tables_2025/feature_concordance_support/baseline_fish_call_recovery_audit.csv"
)

expected_positive_calls <- tibble::tribble(
  ~Patient,  ~feature,
  "IMG-060", "DEL_1P",
  "IMG-060", "DEL_17P",
  "IMG-098", "T_4_14",
  "IMG-181", "DEL_17P"
)
missing_expected_positive_calls <- expected_positive_calls %>%
  dplyr::left_join(
    dat_base %>%
      dplyr::select(.data$Patient, dplyr::all_of(fish_fields)) %>%
      tidyr::pivot_longer(-.data$Patient, names_to = "feature", values_to = "call"),
    by = c("Patient", "feature")
  ) %>%
  dplyr::filter(.data$call != "Positive" | is.na(.data$call))
if (nrow(missing_expected_positive_calls)) {
  stop(
    "Required baseline FISH-positive calls are missing after clinical-source recovery: ",
    paste0(missing_expected_positive_calls$Patient, "/", missing_expected_positive_calls$feature, collapse = ", "),
    call. = FALSE
  )
}

## Edit low confidence call
dat_base <- dat_base %>%
  mutate(across(
    starts_with("WGS_IGH_"),
    ~ if_else(grepl("ZC-02", Patient), 0L, .)
  ))

## ------------------------------------------------------------------
## 4.  MAPPINGS  -----------------------------------------------------
## ------------------------------------------------------------------
map <- tribble(
  ~fish,      ~wgs_bm,                  ~wgs_cf,                          ~type,
  "DEL_1P",   "WGS_del1p_BM_cells",     "WGS_del1p_Blood_plasma_cfDNA",   "CNA",
  "AMP_1Q",   "WGS_amp1q_BM_cells",     "WGS_amp1q_Blood_plasma_cfDNA",   "CNA",
  "DEL_17P",  "WGS_del17p_BM_cells",    "WGS_del17p_Blood_plasma_cfDNA",  "CNA",
  "T_4_14",   "WGS_IGH_FGFR3_BM_cells", "WGS_IGH_FGFR3_Blood_plasma_cfDNA","Translocation",
  "T_11_14",  "WGS_IGH_CCND1_BM_cells", "WGS_IGH_CCND1_Blood_plasma_cfDNA","Translocation",
  "T_14_16",  "WGS_IGH_MAF_BM_cells",   "WGS_IGH_MAF_Blood_plasma_cfDNA",  "Translocation"
)

## ------------------------------------------------------------------
## 5.  TIDY TO LONG FORMAT  -----------------------------------------
## ------------------------------------------------------------------
long <- map %>% 
  pmap_dfr(function(fish, wgs_bm, wgs_cf, type){
    
    dat_base %>% 
      # grab cohort here
      select(
        Patient, 
        cohort,
        Sample_Code, 
        !!sym(fish), 
        !!sym(wgs_bm), 
        !!sym(wgs_cf),
        WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      rename(
        fish_call = !!sym(fish),
        wgs_bm    = !!sym(wgs_bm),
        wgs_cf    = !!sym(wgs_cf),
        tf        = WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      mutate(
        event      = fish,
        type       = type,
        # standardise calls…
        fish_call = case_when(
          str_detect(str_to_lower(fish_call), "pos|^1$|true") ~ 1,
          str_detect(str_to_lower(fish_call), "neg|^0$|false") ~ 0,
          TRUE ~ NA_real_
        ),
        # make WGS calls logical 0/1
        across(
          c(wgs_bm, wgs_cf),
          ~ case_when(
            str_detect(str_to_lower(.), "^1$|true")  ~ 1,
            str_detect(str_to_lower(.), "^0$|false") ~ 0,
            TRUE                                   ~ NA_real_
          )
        ),
        # tumour‐fraction bucket
        tf_group = case_when(
          is.na(tf)      ~ "tf_unknown",
          tf >= tf_cut   ~ "high_tf",
          TRUE           ~ "low_tf"
        )
      ) %>% 
      pivot_longer(
        cols      = c(wgs_bm, wgs_cf),
        names_to  = "wgs_source",
        values_to = "wgs_call"
      ) %>% 
      mutate(
        wgs_source = recode(
          wgs_source,
          wgs_bm = "BM_cells",
          wgs_cf = "cfDNA"
        )
      )
    
  })



## ------------------------------------------------------------------
## 6.  FUNCTION TO SUMMARISE ACCURACY  -------------------------------
## ------------------------------------------------------------------
summarise_concord <- function(df){
  df %>% 
    filter(!is.na(fish_call) & !is.na(wgs_call)) %>%       # both performed
    summarise(
      n          = dplyr::n(),
      tp         = sum(fish_call == 1 & wgs_call == 1),
      tn         = sum(fish_call == 0 & wgs_call == 0),
      fp         = sum(fish_call == 0 & wgs_call == 1),
      fn         = sum(fish_call == 1 & wgs_call == 0),
      concord    = (tp + tn) / n,
      sens       = tp / (tp + fn),
      spec       = tn / (tn + fp)
    )
}

## ------------------------------------------------------------------
## 7.  OVERALL, & BY TF GROUP  --------------------------------------
## ------------------------------------------------------------------

# A)  per‑COHORT × TF‑GROUP  (all combinations)
conc_cohort_tf <- long %>%                     # <‑‑ the “long” object you built
  group_by(event, type, wgs_source,
           cohort,                       # Frontline induction-transplant / pre‑treated
           tf_group)   %>% summarise_concord() %>%                   # high_tf / low_tf / tf_unknow %>% summarise_concord() %>% 
  ungroup()

# B)  per‑COHORT (ignore TF bucket)  → store tf_group == "all"
conc_cohort_overall <- long %>%
  group_by(event, type, wgs_source, cohort) %>% 
  summarise_concord() %>% 
  mutate(tf_group = "all") %>%            # sentinel level
  ungroup()

# C)  per‑TF‑GROUP *irrespective* of cohort (optional – drop if you don’t need it)
conc_tf_overall <- long %>%
  group_by(event, type, wgs_source, tf_group) %>% 
  summarise_concord() %>% 
  mutate(cohort = "All evaluable") %>%    # pooled across cohorts, stratified by TF
  ungroup()

# D) pooled across both cohort and tumour-fraction strata.  This explicit row
# prevents the all-evaluable overall estimate from being confused with the
# all-evaluable `tf_unknown` subset in the reported tables.
conc_all_overall <- long %>%
  group_by(event, type, wgs_source) %>%
  summarise_concord() %>%
  mutate(
    cohort = "All evaluable",
    tf_group = "all"
  ) %>%
  ungroup()

# E)  bind them all and order nicely
concordance_tbl <- bind_rows(
  conc_cohort_overall,
  conc_cohort_tf,
  conc_tf_overall,
  conc_all_overall
) %>% 
  arrange(type, event, wgs_source, cohort, tf_group)

# have a quick look
print(concordance_tbl, n = 20)

### Get overall concordance for translocations and CNAs 
# Translocation concordance on bone‐marrow WGS vs FISH - used in manuscript
transloc_BM <- long %>% 
  filter(type == "Translocation", wgs_source == "BM_cells") %>% 
  summarise_concord()

# Translocation concordance on cfDNA WGS vs FISH
transloc_cf <- long %>% 
  filter(type == "Translocation", wgs_source == "cfDNA") %>% 
  summarise_concord()

# CNA concordance on bone‐marrow WGS vs FISH
cna_BM <- long %>% 
  filter(type == "CNA", wgs_source == "BM_cells") %>% 
  summarise_concord()

# CNA concordance on cfDNA WGS vs FISH
cna_cf <- long %>% 
  filter(type == "CNA", wgs_source == "cfDNA") %>% 
  summarise_concord()



## ------------------------------------------------------------------
## 8. SAVE RESULTS  -----------------------------------------
## ------------------------------------------------------------------
writexl::write_xlsx(concordance_tbl, "Output_tables_2025_updated/FISH_WGS_concordance_with_cohort_updated5.xlsx")



#### Now redo the above, but looking at the specific FISH probe locations instead 
# Load FISH CNA data 
FISH_CNA_combined <- readRDS("Jan2025_exported_data/CNA_at_FISH_sites_combined.rds")

# Keep only Bone Marrow cells and plasma cfDNA samples
FISH_CNA_combined <- FISH_CNA_combined %>%
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"))

patient_key_cols <- intersect(c("Patient", "Patient.x", "Patient.y", "patient"), names(FISH_CNA_combined))
if (length(patient_key_cols) == 0) {
  stop("CNA_at_FISH_sites_combined.rds does not contain a patient identifier column.")
}
FISH_CNA_combined <- FISH_CNA_combined %>%
  mutate(Patient = coalesce(!!!rlang::syms(patient_key_cols)))

## Add ploidy to this
Ploidy_in_BM <- readRDS("Jan2025_exported_data/Sample_ploidy_from_sequenza_400.rds")

colnames(Ploidy_in_BM) <- c(
  "Sample_ID",             # Unique sequencing sample identifier
  "Ploidy_Estimate",       # Estimated ploidy from Sequenza confints file
  "Baseline_Intercept",    # Baseline intercept (if from regression or model fit)
  "Baseline_Source",       # Indicates source of baseline (e.g., BM, cfDNA, etc.)
  "Bam_clean_tmp",              # Cleaned BAM filename or path
  "Patient",            # Patient identifier linking multiple samples
  "Timepoint"        # Sample collection timepoint (Diagnosis, Relapse, etc.)
)


## Pivot wider 
FISH_CNA_wide <- FISH_CNA_combined %>%
  mutate(
    source = case_when(
      Sample_type == "Blood_plasma_cfDNA" ~ "blood",
      Sample_type == "BM_cells"           ~ "BM",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(source)) %>%
  group_by(Patient, Timepoint, source) %>%
  arrange(desc(coalesce(Tumor_Fraction, -Inf)), Date_of_sample_collection) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    Patient, Timepoint, source,
    probe_call_amp1q, probe_call_del17p, probe_call_del1p,
    is_altered_at_probe_amp1q, is_altered_at_probe_del17p, is_altered_at_probe_del1p
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(
      probe_call_amp1q, probe_call_del17p, probe_call_del1p,
      is_altered_at_probe_amp1q, is_altered_at_probe_del17p, is_altered_at_probe_del1p
    ),
    names_glue = "{.value}_{source}"
  ) %>%
  relocate(
    Patient, Timepoint,
    probe_call_amp1q_blood,  probe_call_amp1q_BM,
    probe_call_del17p_blood, probe_call_del17p_BM,
    probe_call_del1p_blood,  probe_call_del1p_BM,
    is_altered_at_probe_amp1q_blood,  is_altered_at_probe_amp1q_BM,
    is_altered_at_probe_del17p_blood, is_altered_at_probe_del17p_BM,
    is_altered_at_probe_del1p_blood,  is_altered_at_probe_del1p_BM
  )

# Check for duplicates
dupes <- FISH_CNA_wide %>%
  count(Patient, Timepoint) %>%
  filter(n > 1)

# Print duplicates
if (nrow(dupes) > 0) {
  message("⚠️ Found duplicate Patient-Timepoint combinations:")
  print(dupes)
} else {
  message("✅ No duplicates found.")
}


## Add it to the other dataframe 
dat_base_with_FISH <- dat_base %>% 
  left_join(FISH_CNA_wide) %>% 
  left_join(Ploidy_in_BM)


## ------------------------------------------------------------------
## 4.  MAPPINGS  -----------------------------------------------------
## ------------------------------------------------------------------
map <- tribble(
  ~fish,      ~wgs_bm,                  ~wgs_cf,                          ~type,
  "DEL_1P",   "is_altered_at_probe_del1p_BM",     "is_altered_at_probe_del1p_blood",   "CNA",
  "AMP_1Q",   "is_altered_at_probe_amp1q_BM",     "is_altered_at_probe_amp1q_blood",   "CNA",
  "DEL_17P",  "is_altered_at_probe_del17p_BM",    "is_altered_at_probe_del17p_blood",  "CNA",
)

## ------------------------------------------------------------------
## 5.  TIDY TO LONG FORMAT  -----------------------------------------
## ------------------------------------------------------------------
long <- map %>% 
  pmap_dfr(function(fish, wgs_bm, wgs_cf, type){
    
    dat_base_with_FISH %>% 
      # grab cohort here
      select(
        Patient, 
        cohort,
        Sample_Code, 
        !!sym(fish), 
        !!sym(wgs_bm), 
        !!sym(wgs_cf),
        WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      rename(
        fish_call = !!sym(fish),
        wgs_bm    = !!sym(wgs_bm),
        wgs_cf    = !!sym(wgs_cf),
        tf        = WGS_Tumor_Fraction_Blood_plasma_cfDNA
      ) %>% 
      mutate(
        event      = fish,
        type       = type,
        # standardise calls…
        fish_call = case_when(
          str_detect(str_to_lower(fish_call), "pos|^1$|true") ~ 1,
          str_detect(str_to_lower(fish_call), "neg|^0$|false") ~ 0,
          TRUE ~ NA_real_
        ),
        # make WGS calls logical 0/1
        across(
          c(wgs_bm, wgs_cf),
          ~ case_when(
            str_detect(str_to_lower(.), "^1$|true")  ~ 1,
            str_detect(str_to_lower(.), "^0$|false") ~ 0,
            TRUE                                   ~ NA_real_
          )
        ),
        # tumour‐fraction bucket
        tf_group = case_when(
          is.na(tf)      ~ "tf_unknown",
          tf >= tf_cut   ~ "high_tf",
          TRUE           ~ "low_tf"
        )
      ) %>% 
      pivot_longer(
        cols      = c(wgs_bm, wgs_cf),
        names_to  = "wgs_source",
        values_to = "wgs_call"
      ) %>% 
      mutate(
        wgs_source = recode(
          wgs_source,
          wgs_bm = "BM_cells",
          wgs_cf = "cfDNA"
        )
      )
    
  })




## ------------------------------------------------------------------
## 7.  OVERALL, & BY TF GROUP  --------------------------------------
## ------------------------------------------------------------------

# A)  per‑COHORT × TF‑GROUP  (all combinations)
conc_cohort_tf2 <- long %>%                     # <‑‑ the “long” object you built
  group_by(event, type, wgs_source,
           cohort,                       # Frontline induction-transplant / pre‑treated
           tf_group)   %>% summarise_concord() %>%                   # high_tf / low_tf / tf_unknow %>% summarise_concord() %>% 
  ungroup()

# B)  per‑COHORT (ignore TF bucket)  → store tf_group == "all"
conc_cohort_overall2 <- long %>%
  group_by(event, type, wgs_source, cohort) %>% 
  summarise_concord() %>% 
  mutate(tf_group = "all") %>%            # sentinel level
  ungroup()

# C)  per‑TF‑GROUP *irrespective* of cohort (optional – drop if you don’t need it)
conc_tf_overall2 <- long %>%
  group_by(event, type, wgs_source, tf_group) %>% 
  summarise_concord() %>% 
  mutate(cohort = "All evaluable") %>%    # pooled across cohorts, stratified by TF
  ungroup()

# D) pooled across both cohort and tumour-fraction strata, matching the
# arm-level table above and preserving the true all-evaluable overall result.
conc_all_overall2 <- long %>%
  group_by(event, type, wgs_source) %>%
  summarise_concord() %>%
  mutate(
    cohort = "All evaluable",
    tf_group = "all"
  ) %>%
  ungroup()

# E)  bind them all and order nicely
concordance_tbl_at_FISH_probe <- bind_rows(
  conc_cohort_overall2,
  conc_cohort_tf2,
  conc_tf_overall2,
  conc_all_overall2
) %>% 
  arrange(type, event, wgs_source, cohort, tf_group)

# have a quick look
print(concordance_tbl_at_FISH_probe, n = 20)

## ------------------------------------------------------------------
## 8. SAVE RESULTS  -----------------------------------------
## ------------------------------------------------------------------
writexl::write_xlsx(concordance_tbl_at_FISH_probe, "Output_tables_2025_updated/FISH_WGS_concordance_at_FISH_probes_CNA.xlsx")










## Put everything together for manuscript 

# ------------------------------------------------------------------
# helper : grab one concordance number
# ------------------------------------------------------------------
get_conc <- function(df,
                     aberration,          # "Translocation" / "CNA"
                     source,              # "BM_cells" / "cfDNA"
                     tf_grp    = "all",   # "all" / "high_tf" / "low_tf" / "tf_unknown"
                     cohort_val = NA      # "Frontline induction-transplant" / "Non-frontline"
){
  out <- df %>% 
    filter(
      type       == aberration,
      wgs_source == source,
      ( is.na(cohort_val) & is.na(cohort) ) | (cohort == cohort_val),
      tf_group   == tf_grp
    ) %>% 
    pull(concord)
  if(length(out)==0) NA_real_ else out[1]
}

# -----------------------------------------------------------------------------
# A) build the full grid of parameters we want
# -----------------------------------------------------------------------------
cohort_vals  <- sort(unique(na.omit(concordance_tbl$cohort)))
sources      <- unique(concordance_tbl$wgs_source)
aberrations  <- c("Translocation","CNA")
tf_groups    <- c("all","high_tf","low_tf")

param_grid <- expand_grid(
  cohort     = cohort_vals,
  wgs_source = sources,
  aberration = aberrations,
  tf_group   = tf_groups
)

# -----------------------------------------------------------------------------
# B) compute concordance for each combination
# -----------------------------------------------------------------------------
results2 <- param_grid %>%
  rowwise() %>%
  mutate(
    concordance = get_conc(
      concordance_tbl,
      aberration, wgs_source,
      tf_grp    = tf_group,
      cohort_val= cohort
    )
  ) %>%
  ungroup()

print(results2)

# -----------------------------------------------------------------------------
# 4) write everything out
# -----------------------------------------------------------------------------
if(!dir.exists("Output_tables_2025")) {
  dir.create("Output_tables_2025")
}

### This is the overall concordance to FISH
write.csv(
  results2,
  "Output_tables_2025_updated/concordance_by_cohort_and_source_and_tf_to_FISH_updated5.csv",
  row.names = FALSE
)






#### Now see the proportion of cases with evidence of disease, calculated in earlier script
## For reporting in manuscript
evidence_summary <- dat_base %>%
  summarise(
    BM_non_na      = sum(!is.na(WGS_Evidence_of_Disease_BM_cells)),
    BM_positive    = sum(WGS_Evidence_of_Disease_BM_cells == 1, na.rm = TRUE),
    BM_pct         = BM_positive / BM_non_na * 100,
    
    Blood_non_na   = sum(!is.na(WGS_Evidence_of_Disease_Blood_plasma_cfDNA)),
    Blood_positive = sum(WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 1, na.rm = TRUE),
    Blood_pct      = Blood_positive / Blood_non_na * 100
  )







###### PART 2: See mutation overlap based on the specific base change 
mutation_data_total <- readRDS("Jan2025_exported_data/mutation_export_updated_more_info2.rds")
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds")
.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

combined_clinical_data_updated <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)


# 1) Annotate the mutation table with clinical metadata
mut_feat <- mutation_data_total %>%
  left_join(
    All_feature_data %>%
      select(Sample, Patient, Sample_type, Tumor_Fraction, Timepoint),
    by = "Sample"
  ) %>%
  mutate(
    Patient      = coalesce(Patient.x, Patient.y),
    Timepoint    = coalesce(Timepoint.x, Timepoint.y),
    Sample_type  = coalesce(Sample_type.x, Sample_type.y)
  ) %>%
  select(
    -Patient.x, -Patient.y,
    -Timepoint.x, -Timepoint.y,
    -Sample_type.x, -Sample_type.y
  ) %>%
  # restrict to only the two sample types of interest
  filter(Sample_type %in% c("BM_cells","Blood_plasma_cfDNA")) 

# 2) Identify only patient×timepoints that have both BM & cfDNA  
matched_pts <- combined_clinical_data_updated %>%
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  distinct(Patient, Timepoint, Sample_type) %>%
  group_by(Patient, Timepoint) %>%
  filter(n_distinct(Sample_type) == 2) %>%
  ungroup() %>%
  select(Patient, Timepoint) %>% 
  unique()

# 3) Keep only those matched cases but add doulbe negatives
mut_matched <- mut_feat %>%
  semi_join(matched_pts, by = c("Patient","Timepoint")) %>%
  inner_join(matched_pts,  by = c("Patient","Timepoint"))

mut_matched$Date_of_sample_collection <-
  as.Date(mut_matched$Date_of_sample_collection)

combined_clinical_data_updated$Date_of_sample_collection <-
  as.Date(combined_clinical_data_updated$Date_of_sample_collection)

# Filter to the manuscript cohort and add only the clinical field used below.
# Explicit keys are essential here: the expanded QC mutation export retains
# caller-dependent annotation columns (including date fields), and an implicit
# join would incorrectly treat every same-named annotation column as a key.
mut_matched <- mut_matched %>%
  left_join(cohort_df, by = "Patient") %>%
  left_join(
    combined_clinical_data_updated %>%
      select(Patient, Timepoint, Sample_type, timepoint_info) %>%
      distinct(),
    by = c("Patient", "Timepoint", "Sample_type"),
    suffix = c("", ".clinical")
  ) %>%
  mutate(
    timepoint_info = coalesce(
      as.character(timepoint_info),
      as.character(timepoint_info.clinical)
    )
  ) %>%
  select(
    -timepoint_info.clinical
  )
mut_matched <- mut_matched %>% filter(timepoint_info %in% c("Baseline", "Diagnosis")) # get baseline
mut_matched <- mut_matched %>% filter(!is.na(Cohort))

# 4) Build per‐patient×timepoint sets of genes
mut_sets <- mut_matched %>%
  group_by(Patient, Timepoint, Cohort, Sample_type) %>%
  summarise(
    muts = list(unique(Mutation_cDNA)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Sample_type,
    values_from = muts,
    values_fill = list(muts = list(character(0)))  # in case one arm has zero
  )

# Command-line runs may produce a matched baseline subset where one sample type
# has no mutation rows after upstream filtering. In an interactive
# RStudio session this was easy to miss because objects from earlier runs could
# remain in memory. Define the expected list columns explicitly so downstream
# concordance statistics always compare BM-vs-cfDNA mutation sets with a clear
# empty-set convention.
expected_mutation_sample_types <- c("BM_cells", "Blood_plasma_cfDNA")
for (sample_type_col in expected_mutation_sample_types) {
  if (!sample_type_col %in% names(mut_sets)) {
    mut_sets[[sample_type_col]] <- rep(list(character(0)), nrow(mut_sets))
  }
}

normalise_mutation_set <- function(x) {
  x <- unlist(x, use.names = FALSE)
  x <- as.character(x[!is.na(x)])
  unique(x)
}

mut_sets <- mut_sets %>%
  mutate(
    BM_cells = map(BM_cells, normalise_mutation_set),
    Blood_plasma_cfDNA = map(Blood_plasma_cfDNA, normalise_mutation_set)
  )

tf_cut <- 0.05   # 5 %

# 1) grab only the cfDNA rows
cfDNA_tf_all <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  select(Patient, Timepoint, Tumor_Fraction) %>%
  distinct()   # in case you have duplicate

# grab cfDNA TF for every matched Pt×TP
cfDNA_tf_all <- cfDNA_tf_all %>%
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction)          ~ "tf_unknown",
      Tumor_Fraction >  tf_cut       ~ "high_tf",
      TRUE                     ~ "low_tf"
    )
  )

# add to the mutation sets
mut_sets <- mut_sets %>%
  left_join(cfDNA_tf_all, by = c("Patient","Timepoint"))

# Mutation-set specificity requires an explicit mutation universe. Interactive
# sessions sometimes had N_muts in memory; clean command-line runs do not. Use
# the observed baseline matched mutation universe from this section, which
# preserves the intended comparison without relying on hidden workspace state.
N_muts <- length(unique(mut_matched$Mutation_cDNA[!is.na(mut_matched$Mutation_cDNA)]))
if (!is.finite(N_muts) || N_muts == 0) {
  stop("Cannot compute mutation-set specificity: no baseline matched mutations define the mutation universe.", call. = FALSE)
}

## 5)  per-row concordance 
concordance_row <- mut_sets %>%
  rowwise() %>%
  mutate(
    tp  = length(intersect(BM_cells, Blood_plasma_cfDNA)),
    fn  = length(setdiff(BM_cells, Blood_plasma_cfDNA)),
    fp  = length(setdiff(Blood_plasma_cfDNA, BM_cells)),
    tn  = N_muts - length(union(BM_cells, Blood_plasma_cfDNA)),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / length(union(BM_cells, Blood_plasma_cfDNA))
  ) %>%
  ungroup()

## 6)  global concordance **within each TF bucket** and **overall**
concordance_tf <- concordance_row %>%
  group_by(Cohort) %>%
  mutate(tf_group = replace_na(tf_group, "tf_unknown")) %>%
  group_by(tf_group, Cohort) %>%
  summarise(
    tp  = sum(tp),
    fn  = sum(fn),
    fp  = sum(fp),
    tn  = sum(tn),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / (tp + fn + fp),   # global Jaccard
    .groups = "drop"
  )

# add an “overall” row
concordance_overall <- concordance_row %>%
  group_by(Cohort) %>%
  summarise(
    tp = sum(tp),
    fn = sum(fn),
    fp = sum(fp),
    tn = sum(tn),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    jaccard     = tp / (tp + fn + fp)
  ) %>%
  mutate(tf_group = "all")

mean_jaccard <- concordance_row %>%
  group_by(Cohort) %>%
  summarise(
    avg_jaccard = mean(jaccard, na.rm = TRUE),
    sd_jaccard  = sd(jaccard, na.rm = TRUE),
    n_samples   = dplyr::n()
  )

concordance_global <- bind_rows(concordance_tf, concordance_overall) %>%
  mutate(
    concordance = (tp + tn) / (tp + tn + fp + fn)
  ) %>%
  select(tf_group, tp, fn, fp, tn,
         sensitivity, specificity,
         jaccard, concordance, Cohort)

cat("\nGlobal concordance by tumour-fraction bucket:\n")
print(concordance_global)


## Export
# 1) CSV
outdir <- "Final Tables and Figures/Baseline_concordance/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  concordance_global,
  file = file.path(outdir, "mutation_concordance_global.csv"),
  row.names = FALSE
)


### This above does not include double negatives 


## Optional operational diagnostic: BAM archive/unarchive helper tables.
##
## This block does not contribute to Extended Data Figure 2 or Supplementary
## Tables 2/3. It was used during data-management review to identify BAM files
## that might need to be retrieved from storage. It is now opt-in so the routine
## manuscript pipeline does not depend on `All_bam_storage_locations.xlsx` and
## does not create support files in the project root.
if (isTRUE(run_bam_archive_diagnostics)) {
  bam_storage_path <- "All_bam_storage_locations.xlsx"
  if (!file.exists(bam_storage_path)) {
    stop(
      "BAM archive diagnostics requested, but `All_bam_storage_locations.xlsx` was not found.",
      call. = FALSE
    )
  }

  bam_diag_dir <- file.path("Output_tables_2025_updated", "support_only_bam_archive_diagnostics")
  dir.create(bam_diag_dir, recursive = TRUE, showWarnings = FALSE)

  mutation_igv_review <- mutation_data_total %>%
    left_join(cohort_df, by = "Patient") %>%
    filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
    filter(!is.na(Cohort))

  write_csv(mutation_igv_review, file.path(bam_diag_dir, "Mutation_iGV_verification.csv"))

  all_bam_storage <- read_excel(bam_storage_path)

  mutation_igv_review <- mutation_igv_review %>%
    mutate(Bam = paste0(Sample, ".bam")) %>%
    left_join(all_bam_storage, by = "Bam")

  missing_path <- mutation_igv_review %>%
    filter(is.na(Path) & is.na(Location))

  all_bam_storage <- all_bam_storage %>%
    mutate(Bam_nosuffix = str_remove_all(Bam, "_WG|_PG"))

  joined <- missing_path %>%
    mutate(Bam_nosuffix = str_remove_all(Bam, "_WG|_PG")) %>%
    left_join(
      all_bam_storage %>%
        select(Bam_storage = Bam, Path, Location, Bam_nosuffix, Cohort),
      by = "Bam_nosuffix"
    )

  write_csv(joined, file.path(bam_diag_dir, "Bam_to_unarchive1.csv"))
  write_csv(mutation_igv_review, file.path(bam_diag_dir, "Bam_to_unarchive2.csv"))

  unmatched <- joined %>% filter(is.na(Bam_storage))

  if (nrow(unmatched) > 0 &&
      requireNamespace("stringdist", quietly = TRUE) &&
      "Bam_noWG" %in% names(all_bam_storage) &&
      any(!is.na(all_bam_storage$Bam_noWG))) {
    similar_matches <- unmatched %>%
      rowwise() %>%
      mutate(
        best_match = {
          dists <- stringdist::stringdist(Bam, all_bam_storage$Bam_noWG)
          i <- which.min(dists)
          all_bam_storage$Bam[i]
        },
        min_dist = {
          dists <- stringdist::stringdist(Bam, all_bam_storage$Bam_noWG)
          min(dists)
        },
        best_path = {
          dists <- stringdist::stringdist(Bam, all_bam_storage$Bam_noWG)
          i <- which.min(dists)
          all_bam_storage$Path[i]
        },
        best_location = {
          dists <- stringdist::stringdist(Bam, all_bam_storage$Bam_noWG)
          i <- which.min(dists)
          all_bam_storage$Location[i]
        }
      ) %>%
      ungroup()

    write_csv(similar_matches, file.path(bam_diag_dir, "Bam_to_unarchive_fuzzy_matches.csv"))
  } else {
    message("Skipping optional BAM fuzzy-match diagnostic: no unmatched BAMs, no stringdist package, or no Bam_noWG column.")
  }

  if (exists("filtered") && all(c("Sample", "Tumor_Sample_Barcode") %in% names(filtered))) {
    export_df <- filtered[, c("Sample", "Tumor_Sample_Barcode")]
    readr::write_csv(as.data.frame(export_df), file.path(bam_diag_dir, "bam_list_and_barcodes.csv"))
  } else {
    message("Skipping optional bam_list_and_barcodes.csv diagnostic: object `filtered` is not defined in this command-line run.")
  }
} else {
  message("Skipping support-only BAM archive diagnostics. Set CFWGS_RUN_BAM_ARCHIVE_DIAGNOSTICS=true to regenerate them.")
}






#### PART 3: Mutation counts and other info 

##### Now get the mutation counts 
# 1) Summary of mutations detected at baseline
baseline_summary <- dat_base %>%
  group_by(cohort) %>%
  summarise(
    n               = dplyr::n(),
    mean_BM         = mean(BM_Mutation_Count,   na.rm = TRUE),
    median_BM       = median(BM_Mutation_Count, na.rm = TRUE),
    sd_BM           = sd(BM_Mutation_Count,     na.rm = TRUE),
    range_BM        = paste0(min(BM_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(BM_Mutation_Count, na.rm = TRUE)),
    mean_Blood      = mean(Blood_Mutation_Count,   na.rm = TRUE),
    median_Blood    = median(Blood_Mutation_Count, na.rm = TRUE),
    sd_Blood        = sd(Blood_Mutation_Count,     na.rm = TRUE),
    range_Blood     = paste0(min(Blood_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(Blood_Mutation_Count, na.rm = TRUE))
  )

print(baseline_summary)

### Add sentences and stats 
# 2) Descriptive sentences––––––––––––––––––––––––––––––––––––––––
sentences <- baseline_summary %>%
  transmute(
    sentence = glue(
      "In the {cohort} cohort (n = {n}), the mean bone marrow mutation count was ",
      "{round(mean_BM,1)} (SD = {round(sd_BM,1)}, range {range_BM}), and the mean ",
      "blood mutation count was {round(mean_Blood,1)} (SD = {round(sd_Blood,1)}, ",
      "range {range_Blood})."
    )
  ) %>%
  pull(sentence)

# Print them to console (you can copy-paste this block into Word)
cat(sentences, sep = "\n\n")

# 3) Overall summary across all cohorts ––––––––––––––––––––––––––––––––––––
overall_summary <- dat_base %>%
  summarise(
    n               = dplyr::n(),
    mean_BM         = mean(BM_Mutation_Count,   na.rm = TRUE),
    median_BM       = median(BM_Mutation_Count, na.rm = TRUE),
    sd_BM           = sd(BM_Mutation_Count,     na.rm = TRUE),
    range_BM        = paste0(min(BM_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(BM_Mutation_Count, na.rm = TRUE)),
    mean_Blood      = mean(Blood_Mutation_Count,   na.rm = TRUE),
    median_Blood    = median(Blood_Mutation_Count, na.rm = TRUE),
    sd_Blood        = sd(Blood_Mutation_Count,     na.rm = TRUE),
    range_Blood     = paste0(min(Blood_Mutation_Count, na.rm = TRUE),
                             "–",
                             max(Blood_Mutation_Count, na.rm = TRUE))
  )

print(overall_summary)

# 4) Descriptive sentence for overall data –––––––––––––––––––––––––––––––––
overall_sentence <- glue(
  "Across all cohorts (n = {overall_summary$n}), the mean bone marrow mutation count was ",
  "{round(overall_summary$mean_BM,1)} (SD = {round(overall_summary$sd_BM,1)}, range {overall_summary$range_BM}), ",
  "and the mean blood mutation count was {round(overall_summary$mean_Blood,1)} ",
  "(SD = {round(overall_summary$sd_Blood,1)}, range {overall_summary$range_Blood})."
)

cat(overall_sentence)


#––– 3) Statistical testing––––––––––––––––––––––––––––––––––––––––––––

run_two_group_test <- function(data, value_col, group_col, method = c("t.test", "wilcox.test")) {
  method <- match.arg(method)
  test_df <- data %>%
    select(value = all_of(value_col), group = all_of(group_col)) %>%
    filter(!is.na(.data$value), !is.na(.data$group))
  group_counts <- test_df %>%
    count(.data$group, name = "n") %>%
    arrange(.data$group)

  has_two_groups <- nrow(group_counts) == 2
  enough_data <- has_two_groups && if (method == "t.test") {
    all(group_counts$n >= 2)
  } else {
    all(group_counts$n >= 1)
  }

  if (!enough_data) {
    return(list(
      result = NULL,
      p_value = NA_real_,
      note = paste0(
        "Skipped ", method, " for ", value_col,
        ": expected two cohorts with ",
        ifelse(method == "t.test", "at least two", "at least one"),
        " non-missing observations each; observed ",
        paste(paste0(group_counts$group, " n=", group_counts$n), collapse = "; ")
      )
    ))
  }

  formula <- stats::as.formula("value ~ group")
  result <- if (method == "t.test") {
    stats::t.test(formula, data = test_df)
  } else {
    stats::wilcox.test(formula, data = test_df)
  }
  list(result = result, p_value = result$p.value, note = NA_character_)
}

# BM mutation counts: Frontline vs Non-frontline
bm_ttest <- run_two_group_test(dat_base, "BM_Mutation_Count", "cohort", "t.test")
bm_wilcox <- run_two_group_test(dat_base, "BM_Mutation_Count", "cohort", "wilcox.test")

# Blood mutation counts: Frontline vs Non-frontline
blood_ttest <- run_two_group_test(dat_base, "Blood_Mutation_Count", "cohort", "t.test")
blood_wilcox <- run_two_group_test(dat_base, "Blood_Mutation_Count", "cohort", "wilcox.test")

# Show results when the cohort sizes support the requested test
bm_ttest$result
bm_wilcox$result
blood_ttest$result
blood_wilcox$result

#––– 4) Optional: extract p-values–––––––––––––––––––––––––––––––––––
pvals <- tibble(
  assay = c("BM", "BM (Wilcox)", "Blood", "Blood (Wilcox)"),
  p_value = c(bm_ttest$p_value,
              bm_wilcox$p_value,
              blood_ttest$p_value,
              blood_wilcox$p_value),
  note = c(bm_ttest$note,
           bm_wilcox$note,
           blood_ttest$note,
           blood_wilcox$note)
)
pvals


#### Assess other correlations with number of mutations detected at baseline

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 1) Correlation matrix (Spearman) + p‐values across all numeric vars  
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# Select only the numeric columns included in the final table. The selected_row_* fields are
# temporary patient-row helpers used before modality-specific cfDNA alignment;
# retaining them would duplicate the cfDNA TF/FS variables used in the final table and
# expose the obsolete n = 61 association in Supplementary Table 3.
num_df <- dat_base %>%
  select(where(is.numeric)) %>%
  select(-any_of(c(
    "selected_row_cfDNA_tumour_fraction",
    "selected_row_cfDNA_FS"
  )))

# rcorr returns:
#  • r : correlation matrix
#  • P : p-value matrix
#  • n : matrix of counts used in each pairwise test
rc  <- rcorr(as.matrix(num_df), type = "spearman")

r_mat <- rc$r
p_mat <- rc$P
n_mat <- rc$n

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# Flatten and include the sample count for each pair
flatten_corr <- function(r_mat, p_mat, n_mat) {
  idx <- which(lower.tri(p_mat), arr.ind = TRUE)
  tibble(
    var1     = rownames(p_mat)[idx[,1]],
    var2     = colnames(p_mat)[idx[,2]],
    rho      = r_mat[idx],
    p_val    = p_mat[idx],
    n_pairs  = n_mat[idx]
  )
}

all_corrs <- flatten_corr(r_mat, p_mat, n_mat)

# Apply Benjamini-Hochberg adjustment across the exploratory correlation matrix.
all_corrs <- all_corrs %>%
  mutate(
    p_adj = p.adjust(p_val, method = "BH")
  )

## Export this
# Keep the historical misspelled "Suplementary" filename because the manuscript
# source map and retained final table use that provenance path. The updated2 file
# is the confirmed source for renamed Supplementary Table 3; the updated file is
# retained as a compatibility export for earlier script runs.
write.csv(all_corrs %>% filter(!is.na(p_adj)), file = "Final Tables and Figures/Suplementary_Table_2_All_Feature_Correlations_updated2.csv")
write.csv(all_corrs %>% filter(!is.na(p_adj)), file = "Final Tables and Figures/Suplementary_Table_2_All_Feature_Correlations_updated.csv")

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 3
#
# What this is:
#   Feature-correlation table comparing baseline mutation burden and related
#   cfWGS/clinical variables, with Benjamini-Hochberg adjusted p-values.
#
# Why it is here:
#   The historical filename says "Suplementary_Table_2", but the audited
#   manuscript source map identifies this export as final Supplementary Table 3.
#
# Current provenance note:
#   The revision-inclusive table regenerated from the current combined cohort is
#   used for the final table. Packaging scripts must use this regenerated file rather
#   than the older retained manuscript CSV.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Suplementary_Table_2_All_Feature_Correlations_updated2.csv",
  artifact_id = "STABLE3",
  role = "regenerated_table_csv",
  description = paste(
    "Regenerated feature-correlation table for Supplementary Table 3;",
    "historical filename retained for provenance."
  ),
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

# now significant by raw p-value
sig_raw <- all_corrs %>% filter(p_val < 0.05)

# significant by FDR
sig_fdr <- all_corrs %>% filter(p_adj < 0.05)

# view both
print(sig_raw)   # exploratory list
print(sig_fdr)   # more stringent list




#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 3) Specific Spearman tests of mutation count vs. tumor fraction + FS  
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
specific_tests <- list(
  BM_vs_BM_TF = cor.test(
    dat_base$BM_Mutation_Count,
    dat_base$WGS_Tumor_Fraction_BM_cells,
    method = "spearman",
    use = "complete.obs"
  ),
  Blood_vs_Blood_TF = cor.test(
    dat_base$Blood_Mutation_Count,
    dat_base$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    method = "spearman",
    use = "complete.obs"
  ),
  BM_vs_FS = cor.test(
    dat_base$BM_Mutation_Count,
    dat_base$FS,
    method = "spearman",
    use = "complete.obs"
  ),
  Blood_vs_FS = cor.test(
    dat_base$Blood_Mutation_Count,
    dat_base$FS,
    method = "spearman",
    use = "complete.obs"
  )
)

# gather into one tidy table
tidy_specific <- map_df(specific_tests, tidy, .id = "comparison")
print(tidy_specific)


### Other tests 
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 1) ISS stage (ordinal) → BM mutation count
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# overall difference
kruskal.test(BM_Mutation_Count ~ ISS_STAGE, data = dat_base)
kruskal.test(Blood_Mutation_Count ~ ISS_STAGE, data = dat_base)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 2) Cytogenetic risk (binary: high vs standard)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
cyto_levels <- unique(na.omit(dat_base$Cytogenetic_Risk))
if (length(cyto_levels) == 2) {
  wilcox.test(BM_Mutation_Count ~ Cytogenetic_Risk,
              data = dat_base,
              na.action = na.exclude)
} else {
  message(
    "Skipping exploratory Cytogenetic_Risk Wilcoxon test: expected 2 observed levels, found ",
    length(cyto_levels), "."
  )
}


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 3) Ig subtype (multilevel categorical)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
kruskal.test(BM_Mutation_Count ~ Subtype,
             data = dat_base)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 4) Spearman correlations for continuous predictors
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
library(Hmisc)
cont_vars <- c("WGS_Tumor_Fraction_BM_cells",
               "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
               "FS",
               "AGE",
               "Plasma_pct",
               "dFLC",
               "LDH")
m <- rcorr(
  as.matrix(dat_base %>% select(BM_Mutation_Count, all_of(cont_vars))),
  type = "spearman"
)
# view rho matrix and p-values
m$r
m$P


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# 5) Simple multivariable model
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––  
# log-transform counts + adjust for stage + TF + FS + age
mod <- lm(log1p(BM_Mutation_Count) ~ 
            factor(ISS_STAGE) +
            WGS_Tumor_Fraction_BM_cells +
            FS +
            AGE,
          data = dat_base)
summary(mod)

### Not significant



### Add some plots

# Create a directory for figures
if (!dir.exists("Final Tables and Figures/Baseline_concordance")) dir.create("Final Tables and Figures/Baseline_concordance")

### Updated style 
# palette
cohort_cols <- c(
  `Frontline induction-transplant` = "#3182bd",
  `Non-frontline`    = "#e6550d"
)

format_p <- function(p) {
  if (p < 0.01) {
    "<0.01"
  } else {
    sprintf("%.2f", p)
  }
}

# Extended Data Figure 2 support - boxplots of baseline BM vs cfDNA mutation counts.
# Older output filenames in this block may still contain "Figure2" labels.
plot_df <- build_baseline_mutation_count_plot_data(
  dat_base,
  baseline_candidates = dat_tb_final
)

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  plot_df,
  file.path(
    "Output_tables_2025", "clinical_support",
    "extended_data_figure_2e_mutation_count_plot_audit.csv"
  )
)
edfig2e_generated_source_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "generated", "figure_components", "Extended_Data_Figure_2", "panel_E"
)
dir.create(edfig2e_generated_source_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  plot_df %>% dplyr::select("Patient", "cohort", "Assay", "MutCount"),
  file.path(
    edfig2e_generated_source_dir,
    "Extended_Data_Figure_2E_mutation_counts_by_cohort_source_data.csv"
  )
)

p1 <- ggplot(plot_df, aes(cohort, MutCount, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  facet_wrap(~Assay) +
  stat_compare_means(method = "wilcox.test",
                     label    = "p.format",
                     label.y  = max(plot_df$MutCount) * 1.05) +
  scale_fill_manual(values = cohort_cols) +
  labs(
    title    = "Baseline mutation counts by cohort",
   # subtitle = "Primary vs Test, in BM and cfDNA",
    x        = "Cohort",
    y        = "Number of mutations"
  ) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training Cohort",
      "Non-frontline"                  = "Test Cohort"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 12,  hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    strip.text    = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2F_boxplot.png", p1, width = 5, height = 4, dpi = 500)

library(ggpubr)  # for stat_compare_means

p1 <- ggplot(plot_df, aes(cohort, MutCount, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  facet_wrap(~Assay) +
  # add the Wilcoxon bracket + star
  stat_compare_means(
    comparisons = list(c("Frontline induction-transplant", "Non-frontline")),
    method      = "wilcox.test",
    label       = "p.signif",     # will show * / ** / ***  
    tip.length  = 0.02,           # how far past the box ends the little ticks go
    bracket.size = 0.4            # thickness of the bracket line
  ) +
  
  scale_fill_manual(values = cohort_cols) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training Cohort",
      "Non-frontline"                  = "Test Cohort"
    )
  ) +
  labs(
    title = "Baseline mutation counts by cohort",
    x     = "Cohort",
    y     = "Number of mutations"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 14,  hjust = 0.5),
    strip.text     = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2F_boxplot_with_bracket.png", p1, width = 5, height = 4.25, dpi = 600)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 2E
#
# What this is:
#   Baseline mutation-count boxplot by cohort, with the original Wilcoxon
#   comparison bracket and training/test cohort labels.
#
# Why it is here:
#   The final assembled Extended Data Figure 2 uses this PNG as panel E.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Baseline_concordance/Figure2F_boxplot_with_bracket.png",
  artifact_id = "EDFIG2E",
  role = "figure_panel_png",
  description = "Baseline mutation-count boxplot used as Extended Data Figure 2E.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)



#### Now redo but percent high tumor fraction instead 
### This is 2B
plot_df <- dat_base %>%
  select(Patient, cohort, WGS_Tumor_Fraction_Blood_plasma_cfDNA) 

edfig2a_source_path <- file.path(
  outdir,
  "Extended_Data_Figure_2A_tumor_fraction_by_cohort_source_data.csv"
)
readr::write_csv(
  plot_df %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)),
  edfig2a_source_path
)

edfig2a_component_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "generated", "figure_components", "Extended_Data_Figure_2", "panel_A"
)
dir.create(edfig2a_component_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  plot_df %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)),
  file.path(
    edfig2a_component_dir,
    "Extended_Data_Figure_2A_tumor_fraction_by_cohort_source_data.csv"
  )
)

p2 <- ggplot(plot_df, aes(cohort, WGS_Tumor_Fraction_Blood_plasma_cfDNA, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, colour = "black", size = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, colour = "black") +
  # add the Wilcoxon bracket + star
  stat_compare_means(
    comparisons = list(c("Frontline induction-transplant", "Non-frontline")),
    method      = "wilcox.test",
    label       = "p.signif",     # will show * / ** / ***  
    tip.length  = 0.02,           # how far past the box ends the little ticks go
    bracket.size = 0.4            # thickness of the bracket line
  ) +
   scale_y_continuous(
    labels = percent_format(accuracy = 1) #,
  #  limits = c(0, 1)       # if you want the axis to go from 0% to 100%
  ) +
  scale_fill_manual(values = cohort_cols) +
  scale_x_discrete(
    labels = c(
      "Frontline induction-transplant" = "Training Cohort",
      "Non-frontline"                  = "Test Cohort"
    )
  ) +
  labs(
    title = "cfDNA tumor fraction by cohort",
    x     = "Cohort",
    y     = "cfDNA tumor fraction (ichorCNA)"
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black")+
  theme_classic(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", size = 12,  hjust = 0.5),
    strip.text     = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2B_boxplot_with_bracket_tumor_fraction.png", p2, width = 4, height = 4, dpi = 600)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 2A
#
# What this is:
#   cfDNA tumor-fraction boxplot by cohort, including the 5% tumor-fraction
#   reference line used throughout this concordance analysis.
#
# Why it is here:
#   The final assembled Extended Data Figure 2 uses this PNG as panel A.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Baseline_concordance/Figure2B_boxplot_with_bracket_tumor_fraction.png",
  artifact_id = "EDFIG2A",
  role = "figure_panel_png",
  description = "cfDNA tumor-fraction by cohort boxplot used as Extended Data Figure 2A.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = edfig2a_source_path,
  artifact_id = "EDFIG2A",
  role = "source_data_csv",
  description = "Patient-level cfDNA tumor fractions plotted in Extended Data Figure 2A.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)


# Extended Data Figure 2 support - BM vs cfDNA mutation-count scatter.
rho_test <- cor.test(dat_base$BM_Mutation_Count,
                     dat_base$Blood_Mutation_Count,
                     method = "spearman")
rho  <- round(rho_test$estimate, 2)
pval_raw <- rho_test$p.value
pval     <- format_p(pval_raw)

p2 <- ggplot(dat_base, aes(BM_Mutation_Count, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho, "\np = ", pval),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation burden: BM vs cfDNA",
    x     = "BM mutation count",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(outdir, "Figure2B_scatter_BM_vs_cfDNA.png"), p2, width = 4, height = 4, dpi = 500)


# Extended Data Figure 2 support - cfDNA mutation count vs ichorCNA tumour fraction.
tf_test <- cor.test(dat_base$Blood_Mutation_Count,
                    dat_base$WGS_Tumor_Fraction_Blood_plasma_cfDNA,
                    method = "spearman")
rho_tf <- round(tf_test$estimate, 2)
p_tf   <- signif(tf_test$p.value, 2)
p_tf     <- format_p(p_tf)

p3 <- ggplot(dat_base,
             aes(WGS_Tumor_Fraction_Blood_plasma_cfDNA, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_tf, "\np = ", p_tf),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs tumour fraction",
    x     = "cfDNA tumour fraction (ichorCNA)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(outdir, "Figure2C_scatter_tf.png"), p3, width = 4, height = 4, dpi = 500)


# Extended Data Figure 2 support - cfDNA mutation count vs fragment-size score (FS).
fs_test <- cor.test(dat_base$Blood_Mutation_Count,
                    dat_base$FS, method = "spearman")
rho_fs <- round(fs_test$estimate, 2)
p_fs   <- signif(fs_test$p.value, 2)
p_fs   <- format_p(p_fs)

p4 <- ggplot(dat_base,
             aes(FS, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_fs, "\np = ", p_fs),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs fragment-size score",
    x     = "Fragment-size score (FS)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(outdir, "Figure2D_scatter_FS.png"), p4, width = 4, height = 4, dpi = 500)


# Extended Data Figure 2 support - cfDNA mutation count vs serum albumin.
alb_test <- cor.test(dat_base$Blood_Mutation_Count,
                     dat_base$Albumin, method = "spearman")
rho_alb <- round(alb_test$estimate, 2)
p_alb   <- signif(alb_test$p.value, 2)
p_alb     <- format_p(p_alb)

p5 <- ggplot(dat_base,
             aes(Albumin, Blood_Mutation_Count, color = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("ρ = ", rho_alb, "\np = ", p_alb),
           hjust = 1.1, vjust = 1.1, size = 4) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  labs(
    title = "Mutation count vs serum albumin",
    x     = "Serum albumin (g/L)",
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(outdir, "Figure2E_scatter_albumin.png"), p5, width = 4, height = 4, dpi = 500)


## Combine  
combined <- p1 + plot_spacer() + p2 + p3 + p4 + p5 +
  plot_layout(widths  = c(1.7, .2, 1, 1, 1, 1), 
              nrow = 1, # 
              heights = c(1)) +   # bottom row slightly shorter
  plot_annotation(
    title = "Baseline Concordance and Clinical Correlates",
    theme = theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave("Final Tables and Figures/Baseline_concordance/Figure2D_combined.png",
       combined,
       width  = 16,
       height = 4,
       dpi    = 600)



### As facet 
# 1) define which x‐vars go in which panel
panels <- tibble(
  var       = c(
    "BM_Mutation_Count",
    "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
    "FS",
    "Albumin"
  ),
  panel_lab = c(
    "BM mutation count",
    "cfDNA tumour fraction\n(ichorCNA)",
    "Fragment-size score (FS)",
    "Serum albumin (g/L)"
  )
)

# 2) pivot dat_base into long form
df_long <- dat_base %>%
  pivot_longer(
    cols      = panels$var,
    names_to  = "var",
    values_to = "x"
  ) %>%
  left_join(panels, by="var")

# 3) compute Spearman ρ and p for each panel
corr_df <- df_long %>%
  group_by(panel_lab) %>%
  summarise(
    n_complete = sum(complete.cases(x, Blood_Mutation_Count)),
    rho = cor(
      x, Blood_Mutation_Count,
      method = "spearman",
      use    = "complete.obs"   # <- NEW
    ),
    p   = cor.test(
      x, Blood_Mutation_Count,
      method = "spearman",
      use    = "complete.obs"   # <- cor.test() understands it too
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_text = if_else(p < 0.01, "p < 0.01", sprintf("p = %.2f", p)),
    label  = sprintf("ρ = %.2f\n%s", rho, p_text)
  )

# 4) now make the faceted scatter
p_combined <- ggplot(df_long, aes(x = x, y = Blood_Mutation_Count, colour = cohort)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed") +
  facet_wrap(~ panel_lab, scales = "free_x", nrow = 1) +
  # add the per‐panel ρ/p annotation
  geom_text(
    data = corr_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = cohort_cols,   # the existing colours
    labels = c(
      "Frontline induction-transplant" = "Training Cohort",
      "Non-frontline"                  = "Test Cohort"
    ),
    name = "Cohort"
  ) +
  labs(
    title = "Clinical correlates of cfDNA mutation burden",
    x     = NULL,
    y     = "cfDNA mutation count"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 11),
    axis.text       = element_text(size = 9),
    legend.position = "top"
  )

# 5) save
ggsave("Final Tables and Figures/Baseline_concordance/Figure2D_facetted_scatter.png",
       p_combined,
       width  = 10,  # for a 4‐panel row
       height = 4,
       dpi    = 600)

# Keep the plotted rows and displayed Spearman statistics synchronized with
# the ED2F image. These sidecars are regenerated from the same `df_long` and
# `corr_df` objects used above, so future cohort additions cannot leave stale
# source tables behind a newly rendered panel.
edfig2f_source_path <- file.path(
  outdir,
  "Extended_Data_Figure_2F_clinical_correlates_source_data.csv"
)
edfig2f_summary_path <- file.path(
  outdir,
  "Extended_Data_Figure_2F_clinical_correlates_spearman_summary.csv"
)

readr::write_csv(
  df_long %>%
    select(Patient, cohort, Blood_Mutation_Count, var, x, panel_lab),
  edfig2f_source_path
)
readr::write_csv(
  corr_df %>%
    select(panel_lab, n_complete, rho, p, p_text, label),
  edfig2f_summary_path
)

# Refresh the established generated-component filenames as well as the compact
# manuscript-object copies written below. Some source-data assembly workflows
# still discover ED2F through these longer legacy filenames.
edfig2f_component_dir <- file.path(
  "Scripts_2025", "Final_Scripts", "final_manuscript_objects",
  "generated", "figure_components", "Extended_Data_Figure_2", "panel_F"
)
dir.create(edfig2f_component_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  df_long %>%
    select(Patient, cohort, Blood_Mutation_Count, var, x, panel_lab),
  file.path(
    edfig2f_component_dir,
    "Extended_Data_Figure_2F_clinical_correlates_source_data.csv"
  )
)
readr::write_csv(
  corr_df %>%
    select(panel_lab, n_complete, rho, p, p_text, label),
  file.path(
    edfig2f_component_dir,
    "Extended_Data_Figure_2F_clinical_correlates_spearman_summary.csv"
  )
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 2F
#
# What this is:
#   Faceted scatterplot of cfDNA mutation burden against baseline clinical and
#   fragmentomic correlates, with Spearman labels and cohort colors.
#
# Why it is here:
#   The final assembled Extended Data Figure 2 uses this PNG as panel F.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Baseline_concordance/Figure2D_facetted_scatter.png",
  artifact_id = "EDFIG2F",
  role = "figure_panel_png",
  description = "Faceted mutation-burden correlation scatterplot used as Extended Data Figure 2F.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = edfig2f_source_path,
  artifact_id = "EDFIG2F",
  role = "source_data_csv",
  description = "Source data plotted in Extended Data Figure 2F.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = edfig2f_summary_path,
  artifact_id = "EDFIG2F",
  role = "summary_csv",
  description = "Complete-case Spearman statistics displayed in Extended Data Figure 2F.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)






#### Now make a dumbbell plot to show concordance between BM and cfDNA
# 1) reshape, keep only frontline + TF strata
event_tf_conc <- concordance_tbl %>%
  filter(cohort == "Frontline induction-transplant",
         tf_group %in% c("high_tf","low_tf")) %>%
  transmute(
    event    = event,
    TF_group = factor(tf_group,
                      levels = c("high_tf","low_tf"),
                      labels = c("High TF","Low TF")),
    sample   = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    conc     = concord * 100      # percent
  )

# 2) extract BM concordance at High TF, to define ordering
bm_high <- event_tf_conc %>%
  filter(TF_group=="High TF", sample=="BM") %>%
  arrange(conc) %>%
  pull(event)

# 3) apply that ordering to the factor
event_tf_conc <- event_tf_conc %>%
  mutate(
    event = factor(event, levels = bm_high)
  )

# Get overall
overall_conc <- concordance_tbl %>%
  filter(
    cohort   == "Frontline induction-transplant",
    tf_group == "all"
  ) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    sample = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    conc   = concord * 100                       # percent
  )

# Overall sensitivity across all tumour-fraction strata.
sens_overall <- long %>%
  filter(cohort   == "Frontline induction-transplant") %>%
  group_by(event, wgs_source) %>%
  summarise_concord() %>%         # gives sens, etc
  ungroup() %>%
  transmute(
    event      = factor(event, levels = bm_high),
    sample     = if_else(wgs_source=="BM_cells","BM","cfDNA"),
    sens_pct   = sens * 100       # percent
  )


# 4) plot
p_tf <- ggplot(event_tf_conc,
               aes(x = conc, y = event, group = event)) +
  # grey horizontal connector
  geom_line(color="grey80", size=0.6) +
  # two TF‐group points
  geom_point(aes(colour = TF_group), size = 3) +
  # one facet per sample
  facet_wrap(~sample, nrow = 1) +
  scale_colour_viridis_d(
    option = "D", end = 0.8,
    name = "Tumour fraction"
  ) +
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0,100, by=25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "SV/CNA concordance with FISH by tumour fraction",
    x     = "Concordance with FISH",
    y     = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    panel.spacing   = unit(1, "lines")
  )


# 5) save
ggsave("Final Tables and Figures/Baseline_concordance/Fig2B_event_concordance_by_TF_updated4.png", p_tf,
       width = 5, height = 4, dpi = 600)


# Add overall sensitivity as a star marker on the tumour-fraction concordance
# plot. This lets the panel show both concordance by tumour fraction and the
# aggregate detection sensitivity for each event.
tf_plot_df <- bind_rows(
  # high/low TF
  event_tf_conc %>%
    rename(value = conc) %>%
    transmute(event, sample, Measure = TF_group, value),
  # overall sensitivity
  sens_overall %>%
    transmute(event, sample,
              Measure = factor("Overall sensitivity",
                               levels=c("High TF","Low TF","Overall sensitivity")),
              value = sens_pct)
) %>%
  mutate(
    Measure = factor(Measure, 
                     levels=c("High TF","Low TF","Overall sensitivity"))
  )

# 2) the plot
pretty_events <- c(
  T_4_14   = "t(4;14)",
  T_11_14  = "t(11;14)",
  T_14_16  = "t(14;16)",
  AMP_1Q   = "amp(1q)",
  DEL_17P  = "del(17p)",
  DEL_1P   = "del(1p)"
)

p_tf_sens2 <- ggplot(tf_plot_df, aes(x = value, y = event, group = event)) +
  
  # grey connector only for the TF strata
  # geom_line(
  #   data = filter(tf_plot_df, Measure %in% c("High TF","Low TF")),
  #   aes(x = value, y = event, group = event),
  #   colour = "grey80", size = 0.6
  # ) +
  # 
  # all three measures as points
  geom_point(aes(colour = Measure, shape = Measure),
             size = 3, stroke = 1) +
  
  # facet by BM vs cfDNA
  facet_wrap(~ sample, nrow = 1) +
  
  # Nice labels
  scale_y_discrete(
    labels = pretty_events
  ) +
  # single legend with both colour + shape
  scale_colour_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = viridis(2, end = 0.8)[1],
      "Low TF"              = viridis(2, end = 0.8)[2]
   #   "Overall sensitivity" = "black"
    )
  ) +
  scale_shape_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = 16,
      "Low TF"              = 16
   #   "Overall sensitivity" = 8
  # "Overall sensitivity" = NA
    )
  ) +
  
  # nice breathing room at 0% and 100%
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0,100,25), "%"),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  labs(
  #  title = "SV/CNA concordance with FISH (●) and overall sensitivity",
   title = "Structural and copy-number variant concordance with FISH",
    x     = "Percent",
    y     = NULL
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

# 3) save
ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2B_event_concordance_with_sensitivity_updated5.png",
  p_tf_sens2, width = 5.5, height = 4, dpi = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 2C
#
# What this is:
#   Structural-variant and copy-number FISH concordance plot stratified by high
#   versus low tumor fraction, faceted by BM and cfDNA source.
#
# Why it is here:
#   The audited manuscript map assigns this PNG to final Extended Data Figure
#   2C, even though the historical filename still contains "Fig2B".
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Baseline_concordance/Fig2B_event_concordance_with_sensitivity_updated5.png",
  artifact_id = "EDFIG2B_D_C",
  role = "training_cohort_figure_panel_png",
  description = "Training-cohort-only FISH concordance panel retained as a secondary analysis.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

# -------------------------------------------------------------------------
# Alternate ED2C-style panel: all evaluable diagnosis/baseline patients
#
# The mapped manuscript ED2C panel above intentionally preserves the original
# training-cohort-only analysis. This alternate version pools all patients that
# reached the baseline concordance table after the same diagnosis/baseline,
# duplicate-resolution, and cohort-assignment filters used upstream.
# -------------------------------------------------------------------------
event_tf_conc_all_evaluable <- concordance_tbl %>%
  filter(cohort == "All evaluable", tf_group %in% c("high_tf", "low_tf")) %>%
  transmute(
    event = event,
    TF_group = factor(
      tf_group,
      levels = c("high_tf", "low_tf"),
      labels = c("High TF", "Low TF")
    ),
    sample = if_else(wgs_source == "BM_cells", "BM", "cfDNA"),
    conc = concord * 100,
    n = n,
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn
  )

bm_high_all_evaluable <- event_tf_conc_all_evaluable %>%
  filter(TF_group == "High TF", sample == "BM") %>%
  arrange(conc) %>%
  pull(event)

event_tf_conc_all_evaluable <- event_tf_conc_all_evaluable %>%
  mutate(event = factor(event, levels = bm_high_all_evaluable))

sens_overall_all_evaluable <- long %>%
  group_by(event, wgs_source) %>%
  summarise_concord() %>%
  ungroup() %>%
  transmute(
    event = factor(event, levels = bm_high_all_evaluable),
    sample = if_else(wgs_source == "BM_cells", "BM", "cfDNA"),
    Measure = factor(
      "Overall sensitivity",
      levels = c("High TF", "Low TF", "Overall sensitivity")
    ),
    value = sens * 100,
    n = n,
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn
  )

tf_plot_df_all_evaluable <- bind_rows(
  event_tf_conc_all_evaluable %>%
    transmute(event, sample, Measure = TF_group, value = conc, n, tp, tn, fp, fn),
  sens_overall_all_evaluable
) %>%
  mutate(
    Measure = factor(
      Measure,
      levels = c("High TF", "Low TF", "Overall sensitivity")
    )
  )

readr::write_csv(
  tf_plot_df_all_evaluable,
  file.path(
    outdir,
    "Extended_Data_Figure_2C_all_evaluable_baseline_FISH_concordance_source_data.csv"
  )
)

# Keep the aggregate sensitivity rows in the source-data export for auditability,
# but do not display them in the concordance panel.
tf_plot_df_all_evaluable_display <- tf_plot_df_all_evaluable %>%
  filter(Measure %in% c("High TF", "Low TF")) %>%
  mutate(Measure = droplevels(Measure))

p_tf_sens2_all_evaluable <- ggplot(
  tf_plot_df_all_evaluable_display,
  aes(x = value, y = event, group = event)
) +
  geom_point(aes(colour = Measure, shape = Measure), size = 3, stroke = 1) +
  facet_wrap(~ sample, nrow = 1) +
  scale_y_discrete(labels = pretty_events) +
  scale_colour_manual(
    name = "Concordance",
    values = c(
      "High TF" = viridis(2, end = 0.8)[1],
      "Low TF" = viridis(2, end = 0.8)[2]
    )
  ) +
  scale_shape_manual(
    name = "Concordance",
    values = c(
      "High TF" = 16,
      "Low TF" = 16
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0, 100, 25), "%"),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  labs(
    title = "Structural and copy-number variant concordance with FISH",
    x = "Percent",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.box = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

all_evaluable_ed2c_path <- file.path(
  outdir,
  "Fig2B_event_concordance_with_sensitivity_all_evaluable_baseline.png"
)
ggsave(all_evaluable_ed2c_path, p_tf_sens2_all_evaluable, width = 5.5, height = 4, dpi = 600)

# Preserve the narrower training-cohort-only panel under an explicit filename,
# then make the long-standing manuscript filename point to the revision-
# inclusive all-evaluable panel. This prevents downstream assemblers or users
# from silently selecting a scientifically different denominator based only on
# the historical filename.
training_only_ed2c_path <- file.path(
  outdir,
  "Fig2B_event_concordance_with_sensitivity_training_cohort.png"
)
file.copy(
  file.path(outdir, "Fig2B_event_concordance_with_sensitivity_updated5.png"),
  training_only_ed2c_path,
  overwrite = TRUE
)
file.copy(
  all_evaluable_ed2c_path,
  file.path(outdir, "Fig2B_event_concordance_with_sensitivity_updated5.png"),
  overwrite = TRUE
)
tf_plot_df <- tf_plot_df_all_evaluable

ms_copy_artifact(
  source_path = all_evaluable_ed2c_path,
  artifact_id = "EDFIG2B_D_C",
  role = "all_evaluable_figure_panel_png",
  description = "All-evaluable training-plus-test baseline/diagnosis FISH concordance panel for Extended Data Figure 2C.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = file.path(
    outdir,
    "Extended_Data_Figure_2C_all_evaluable_baseline_FISH_concordance_source_data.csv"
  ),
  artifact_id = "EDFIG2B_D_C",
  role = "all_evaluable_source_data_csv",
  description = "Source data for the alternate all-evaluable baseline/diagnosis FISH concordance panel.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

# The revision-inclusive, training-plus-test panel is used as
# Extended Data Figure 2C.  Re-export it under the primary roles after retaining
# the historical training-only companion above so downstream assemblers cannot
# silently select the narrower panel.
ms_copy_artifact(
  source_path = all_evaluable_ed2c_path,
  artifact_id = "EDFIG2B_D_C",
  role = "figure_panel_png",
  description = "Primary Extended Data Figure 2C: all-evaluable training-plus-test FISH concordance by tumor fraction.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = file.path(
    outdir,
    "Extended_Data_Figure_2C_all_evaluable_baseline_FISH_concordance_source_data.csv"
  ),
  artifact_id = "EDFIG2B_D_C",
  role = "source_data_csv",
  description = "Primary Extended Data Figure 2C source data across all evaluable training and test baseline specimens.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)


# Calculate sensitivity for CNAs and translocations separately, by WGS source
# Add type_group column and filter for events we care about
long_metrics <- long %>%
  mutate(type_group = case_when(
    grepl("^IGH_", event) ~ "Translocation",
    type == "CNA" ~ "CNA",
    TRUE ~ type
  )) %>%
  filter(!is.na(fish_call) & !is.na(wgs_call)) %>%
  filter(type_group %in% c("CNA", "Translocation"))

# Function to compute metrics for any grouping
compute_metrics <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars)), .add = FALSE) %>%
    summarise(
      tp = sum(wgs_call == 1 & fish_call == 1),
      tn = sum(wgs_call == 0 & fish_call == 0),
      fp = sum(wgs_call == 1 & fish_call == 0),
      fn = sum(wgs_call == 0 & fish_call == 1),
      sensitivity  = tp / (tp + fn),
      specificity  = tn / (tn + fp),
      concordance  = (tp + tn) / (tp + tn + fp + fn),
      ppv          = tp / (tp + fp),
      npv          = tn / (tn + fn),
      .groups = "drop"
    )
}

# Metrics per variant category (CNA vs Translocation)
metrics_by_category <- compute_metrics(long_metrics, c("type_group", "wgs_source"))

# Metrics per feature (event)
metrics_by_feature <- compute_metrics(long_metrics, c("type_group", "event", "wgs_source"))


# Export to CSVs
write_csv(metrics_by_category, "Exported_data_tables_clinical/Supp_table_wgs_fish_metrics_by_category.csv")
write_csv(metrics_by_feature, "Exported_data_tables_clinical/Supp_table_wgs_fish_metrics_by_feature.csv")


## Put together 
metrics_by_category <- metrics_by_category %>%
  mutate(
    event = dplyr::case_when(
      type_group == "CNA" ~ "All_CNAs",
      type_group == "Translocation" ~ "All_Translocations",
      TRUE ~ "All_Other"
    )
  ) %>%
  relocate(event, .after = type_group)

metrics_combined <- bind_rows(
  metrics_by_category,
  metrics_by_feature
) %>%
  arrange(type_group, factor(event, levels = c("All_CNAs","All_Translocations")), wgs_source)

## Round 
metrics_combined <- metrics_combined %>%
  mutate(
    dplyr::across(
      c(sensitivity, specificity, concordance, ppv, npv),
      ~ round(.x, 3)
    )
  )

write_csv(metrics_combined, "Exported_data_tables_clinical/Supp_table_2_WGS_vs_FISH_metrics_combined3.csv")



### Now go more into depth for the cfDNA-BM concordance and plot
dir <- "Output_tables_2025_updated/"
merged_mut   <- read_rds(file.path(dir, "merged_mut.rds"))
merged_trans <- read_rds(file.path(dir, "merged_trans.rds"))
merged_CNA   <- read_rds(file.path(dir, "merged_CNA.rds"))

tf_cutoff <- 0.05   # 5% ctDNA threshold

# 1) Build a per‐event × TF‐group performance table
## First CNA
perf_by_tf <- merged_CNA %>%
  # Revision-inclusive analysis: retain all cohort-assigned patients with an
  # evaluable matched baseline BM/cfDNA pair. Historically this block filtered
  # to the frontline cohort, which silently excluded eligible test-cohort pairs
  # from Supplementary Table 2 and Extended Data Figure 2B.
  filter(
    timepoint_info_blood  == "Baseline"
  ) %>%
  # assign each sample to Low / High TF
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction_blood)    ~ NA_character_,
      Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
      TRUE                            ~ "Low TF"
    )
  ) %>%
  
  # pivot the five CNAs into long form
  pivot_longer(
    cols = matches("^(del1p|amp1q|del13q|del17p|hyperdiploid)_(BM|blood)$"),
    names_to      = c("event","source"),
    names_pattern = "(.*)_(BM|blood)$",
    values_to     = "call"
  ) %>%
  # normalize calls to logical
  mutate(
    call   = call == "Yes",
    source = if_else(source=="BM",   "BM", "cfDNA")
  ) %>%
  # spread BM vs blood side by side
  pivot_wider(
    id_cols    = c(Patient, event, tf_group),
    names_from = source,
    values_from= call
  ) %>%
  
  # now summarise per‐event × TF‐group
  group_by(event, tf_group) %>%
  filter(!is.na(cfDNA), !is.na(BM)) %>% # remove when either is NA
  summarise(
    n            = dplyr::n(),                                    # samples
    tp           = sum(cfDNA & BM,   na.rm=TRUE),
    tn           = sum(!cfDNA & !BM, na.rm=TRUE),
    fp           = sum(cfDNA & !BM,  na.rm=TRUE),
    fn           = sum(!cfDNA & BM,  na.rm=TRUE),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  )

# 2) Add the “All”‐TF row for each event
perf_all <- perf_by_tf %>%
  filter(!is.na(tf_group)) %>%         # drop any NA‐TF rows
  group_by(event) %>%
  summarise(
    n            = sum(n),
    tp           = sum(tp),
    tn           = sum(tn),
    fp           = sum(fp),
    fn           = sum(fn),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  ) %>%
  mutate(tf_group = "All")

perf_tf_complete <- bind_rows(perf_by_tf, perf_all) %>%
  mutate(
    tf_group = factor(tf_group, levels = c("Low TF","High TF","All"))
  )


### Now redo for translocations 
merged_trans <- merged_trans %>%
  mutate(
    Patient = str_remove(Patient_Timepoint, "_Baseline$")
  )

# Patient/event-level audit for the revision-inclusive matched-pair analysis.
# This makes the denominator and the contribution of each cohort directly
# traceable rather than recoverable only from aggregate confusion matrices.
translocation_pair_audit <- merged_trans %>%
  filter(timepoint_info_blood == "Baseline") %>%
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction_blood) ~ NA_character_,
      Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
      TRUE ~ "Low TF"
    )
  ) %>%
  pivot_longer(
    cols = matches("^(IGH_MAF|IGH_MYC|IGH_CCND1|IGH_FGFR3)_(BM|blood)$"),
    names_to = c("event", "source"),
    names_pattern = "(.*)_(BM|blood)$",
    values_to = "call"
  ) %>%
  pivot_wider(
    id_cols = c(Patient, cohort, Tumor_Fraction_blood, tf_group, event),
    names_from = source,
    values_from = call
  ) %>%
  filter(!is.na(BM), !is.na(blood)) %>%
  mutate(
    pair_class = case_when(
      BM == "Yes" & blood == "Yes" ~ "TP",
      BM == "No" & blood == "No" ~ "TN",
      BM == "No" & blood == "Yes" ~ "FP",
      BM == "Yes" & blood == "No" ~ "FN",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(cohort, Patient, event)

if (!any(translocation_pair_audit$cohort == "Non-frontline")) {
  stop("Revision-inclusive translocation audit contains no test-cohort pairs.")
}
readr::write_csv(
  translocation_pair_audit,
  "Final Tables and Figures/Baseline_concordance/revision_inclusive_baseline_translocation_pair_audit.csv"
)

perf_by_tf <- merged_trans %>%
  # Revision-inclusive analysis: retain eligible training- and test-cohort
  # matched baseline pairs (see the corresponding CNA block above).
  filter(
    timepoint_info_blood  == "Baseline"
  ) %>%
  # assign each sample to Low / High TF
  mutate(
    tf_group = case_when(
      is.na(Tumor_Fraction_blood)    ~ NA_character_,
      Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
      TRUE                            ~ "Low TF"
    )
  ) %>%
  
  # pivot the five CNAs into long form
  pivot_longer(
    cols = matches("^(IGH_MAF|IGH_MYC|IGH_CCND1|IGH_FGFR3)_(BM|blood)$"),
    names_to      = c("event","source"),
    names_pattern = "(.*)_(BM|blood)$",
    values_to     = "call"
  ) %>%
  # normalize calls to logical
  mutate(
    call   = call == "Yes",
    source = if_else(source=="BM",   "BM", "cfDNA")
  ) %>%
  # spread BM vs blood side by side
  pivot_wider(
    id_cols    = c(Patient, event, tf_group),
    names_from = source,
    values_from= call
  ) %>%
  
  # now summarise per‐event × TF‐group
  group_by(event, tf_group) %>%
  filter(!is.na(cfDNA), !is.na(BM)) %>% # remove when either is NA
  summarise(
    n            = dplyr::n(),                                    # samples
    tp           = sum(cfDNA & BM,   na.rm=TRUE),
    tn           = sum(!cfDNA & !BM, na.rm=TRUE),
    fp           = sum(cfDNA & !BM,  na.rm=TRUE),
    fn           = sum(!cfDNA & BM,  na.rm=TRUE),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  )

# 2) Add the “All”‐TF row for each event
perf_all <- perf_by_tf %>%
  filter(!is.na(tf_group)) %>%         # drop any NA‐TF rows
  group_by(event) %>%
  summarise(
    n            = sum(n),
    tp           = sum(tp),
    tn           = sum(tn),
    fp           = sum(fp),
    fn           = sum(fn),
    sensitivity  = tp/(tp + fn),
    specificity  = tn/(tn + fp),
    concordance  = (tp + tn)/n,
    .groups      = "drop"
  ) %>%
  mutate(tf_group = "All")

perf_tf_complete_trans <- bind_rows(perf_by_tf, perf_all) %>%
  mutate(
    tf_group = factor(tf_group, levels = c("Low TF","High TF","All"))
  )

## Now bind the two rows together 
perf_tf_complete <- bind_rows(perf_tf_complete, perf_tf_complete_trans)

## Add the mutations 
concordance_global$event <- "Mutations"
tmp <- concordance_global %>%
  # Pool the frozen training and expanded test cohorts by summing their
  # confusion-matrix counts, then recompute all derived metrics. Averaging the
  # cohort-specific rates would weight cohorts equally rather than weighting
  # the underlying evaluable mutation events.
  group_by(tf_group, event) %>%
  summarise(
    tp = sum(tp, na.rm = TRUE),
    fn = sum(fn, na.rm = TRUE),
    fp = sum(fp, na.rm = TRUE),
    tn = sum(tn, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sensitivity = if_else((tp + fn) > 0, tp / (tp + fn), NA_real_),
    specificity = if_else((tn + fp) > 0, tn / (tn + fp), NA_real_),
    jaccard = if_else((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_),
    concordance = if_else(
      (tp + tn + fp + fn) > 0,
      (tp + tn) / (tp + tn + fp + fn),
      NA_real_
    )
  )
tmp <- tmp %>% 
  mutate(
    tf_group = case_when(
      tf_group == "high_tf" ~ "High TF",
      tf_group == "low_tf"  ~ "Low TF",
      tf_group == "all"     ~ "All",
      TRUE             ~ tf_group
    )
  )

perf_tf_complete <- bind_rows(perf_tf_complete, tmp)


## Now make plot
# 1) reshape, keep only frontline + TF strata
event_conc <- perf_tf_complete %>%
  filter(tf_group %in% c("High TF","Low TF")) %>%
  transmute(
    event    = event,
    TF_group = factor(tf_group,
                      levels = c("High TF","Low TF"),
                      labels = c("High TF","Low TF")),
    conc     = concordance * 100      # percent
  )

# 2) extract BM concordance at High TF, to define ordering
bm_high <- event_conc %>%
  filter(TF_group=="High TF") %>%
  arrange(conc) %>%
  pull(event)

# 3) apply that ordering to the factor
event_conc <- event_conc %>%
  mutate(
    event = factor(event, levels = bm_high)
  )

# Get overall
overall_conc <- perf_tf_complete %>%
  filter(tf_group %in% c("All")) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    conc     = concordance * 100      # percent
  )

# Overall sensitivity across all tumour-fraction strata.
sens_overall <- perf_tf_complete %>%
  filter(tf_group %in% c("All")) %>%
  transmute(
    event  = factor(event, levels = bm_high),    # same ordering
    sens_pct   = sensitivity* 100      # percent
  )
  
# Add overall sensitivity as a star marker on the BM-vs-cfDNA concordance plot.
tf_plot_df_BM <- bind_rows(
  # high/low TF
  event_conc %>%
    rename(value = conc) %>%
    transmute(event, Measure = TF_group, value),
  # overall sensitivity
  sens_overall %>%
    transmute(event,
              Measure = factor("Overall sensitivity",
                               levels=c("High TF","Low TF","Overall sensitivity")),
              value = sens_pct)
) %>%
  mutate(
    Measure = factor(Measure, 
                     levels=c("High TF","Low TF","Overall sensitivity"))
  )

# 2) the plot
pretty_events <- c(
  amp1q        = "amp(1q)",
  del13q       = "del(13q)",
  del17p       = "del(17p)",
  del1p        = "del(1p)",
  hyperdiploid = "hyperdiploid",
  IGH_CCND1    = "t(11;14) IGH-CCND1",
  IGH_FGFR3    = "t(4;14) IGH-FGFR3",
  IGH_MAF      = "t(14;16) IGH-MAF",
  IGH_MYC      = "t(8;14) IGH-MYC",
  Mutations      = "Mutations"
  
)

p_tf_sens <- ggplot(tf_plot_df_BM, aes(x = value, y = event, group = event)) +
  
  # grey connector only for the TF strata
  geom_line(
    data = filter(tf_plot_df_BM, Measure %in% c("High TF","Low TF")),
    aes(x = value, y = event, group = event),
    colour = "grey80", size = 0.6
  ) +
  
  # all three measures as points
  geom_point(aes(colour = Measure, shape = Measure),
             size = 3, stroke = 1) +

  # Nice labels
  scale_y_discrete(
    labels = pretty_events
  ) +
  # single legend with both colour + shape
  scale_colour_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = viridis(2, end = 0.8)[1],
      "Low TF"              = viridis(2, end = 0.8)[2],
      "Overall sensitivity" = "black"
    )
  ) +
  scale_shape_manual(
    name   = "Concordance",
    values = c(
      "High TF"             = 16,
      "Low TF"              = 16,
      "Overall sensitivity" = 8
    )
  ) +
  
  # nice breathing room at 0% and 100%
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0,100,25), "%"),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  labs(
    title = "BM vs cfDNA SV/CNA concordance (●) and sensitivity",
    x     = "Percent",
    y     = NULL
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

# 3) save
ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2C_event_concordance_between_BM_and_cfDNA_with_sensitivity4_updated2.png",
  p_tf_sens, width = 5, height = 4, dpi = 600
)

bad <- tf_plot_df_BM %>%
  dplyr::filter(is.na(value) | value < 0 | value > 100)
bad




# Final three-panel event-level summary for Extended Data Figure 2:
# concordance, sensitivity, and specificity by tumour-fraction stratum.
# 1.  Reshape:  Concordance, Sensitivity, Specificity  ×  TF-group
perf_long <- perf_tf_complete %>%
  filter(tf_group %in% c("High TF", "Low TF")) %>%         # keep only strata
  pivot_longer(
    cols      = c(concordance, sensitivity, specificity),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Percent  = Value * 100,
    TF_group = factor(tf_group, levels = c("High TF", "Low TF"))
  )


# 2.  Event ordering (by High-TF Concordance)
event_order <- perf_long %>%
  filter(Metric == "concordance", TF_group == "High TF") %>%
  arrange(Percent) %>%
  pull(event)

perf_long <- perf_long %>%
  mutate(event = factor(event, levels = event_order))

# Source data used for Extended Data Figure 2B. The
# primary performance object is revision-inclusive (all evaluable baseline
# training and test pairs), so export that exact object under the historical
# source-data filename consumed by the locked-source workbook builder.
edfig2b_primary_source_path <- file.path(
  outdir,
  "Extended_Data_Figure_2B_BM_cfDNA_performance_byTF_plot_source_data.csv"
)
readr::write_csv(perf_long, edfig2b_primary_source_path)

# 4.  Plot
p_3panel <- ggplot(perf_long,
                   aes(x = Percent, y = event, group = event)) +
  # grey connector between High / Low TF
 # geom_line(aes(group = interaction(event, Metric)),
 #           colour = "grey80", linewidth = 0.6) +
  # points for the two strata
  geom_point(aes(colour = TF_group, shape = TF_group),
             size = 3, stroke = 0.8) +
  facet_wrap(~ Metric, nrow = 1,
             labeller = labeller(Metric = c(
               concordance = "Concordance",
               sensitivity = "Sensitivity",
               specificity = "Specificity"
             ))) +
  scale_y_discrete(labels = pretty_events) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_colour_manual(
    values = c("High TF" = viridis(2, end = 0.8)[1],
               "Low TF"  = viridis(2, end = 0.8)[2]),
    name   = "Tumour fraction"
  ) +
  scale_shape_manual(
    values = c("High TF" = 16, "Low TF" = 16),
    name   = "Tumour fraction"
  ) +
  labs(
    title = "BM vs cfDNA performance by ctDNA fraction",
    x     = "Percent",
    y     = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 8),
    legend.position = "top",
    legend.box      = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

ggsave(
  "Final Tables and Figures/Baseline_concordance/Fig2C_BM_cfDNA_conc_sens_spec_byTF_updated6.png",
  p_3panel, width = 5.5, height = 4, dpi = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 2B
#
# What this is:
#   Three-panel BM-vs-cfDNA concordance/sensitivity/specificity summary by high
#   versus low tumor fraction.
#
# Why it is here:
#   The audited manuscript map assigns this PNG to final Extended Data Figure
#   2B, even though the historical filename still contains "Fig2C".
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Baseline_concordance/Fig2C_BM_cfDNA_conc_sens_spec_byTF_updated6.png",
  artifact_id = "EDFIG2B_D_B",
  role = "figure_panel_png",
  description = "BM-vs-cfDNA concordance, sensitivity, and specificity panel used as Extended Data Figure 2B.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

# -------------------------------------------------------------------------
# Independent ED2B reconstruction: all evaluable diagnosis/baseline patients
#
# The primary performance object above is already revision-inclusive. This
# second construction independently recomputes the same all-evaluable scope
# from the matched CNA, translocation, and mutation objects and writes an
# explicitly named audit companion; it is not a separate training-only result.
# -------------------------------------------------------------------------
summarise_bm_cfdna_event_performance <- function(data, event_regex) {
  data %>%
    filter(timepoint_info_blood == "Baseline") %>%
    mutate(
      tf_group = case_when(
        is.na(Tumor_Fraction_blood) ~ NA_character_,
        Tumor_Fraction_blood >= tf_cutoff ~ "High TF",
        TRUE ~ "Low TF"
      )
    ) %>%
    pivot_longer(
      cols = matches(event_regex),
      names_to = c("event", "source"),
      names_pattern = "(.*)_(BM|blood)$",
      values_to = "call"
    ) %>%
    mutate(
      call = call == "Yes",
      source = if_else(source == "BM", "BM", "cfDNA")
    ) %>%
    pivot_wider(
      id_cols = c(Patient, event, tf_group),
      names_from = source,
      values_from = call
    ) %>%
    group_by(event, tf_group) %>%
    filter(!is.na(cfDNA), !is.na(BM)) %>%
    summarise(
      n = dplyr::n(),
      tp = sum(cfDNA & BM, na.rm = TRUE),
      tn = sum(!cfDNA & !BM, na.rm = TRUE),
      fp = sum(cfDNA & !BM, na.rm = TRUE),
      fn = sum(!cfDNA & BM, na.rm = TRUE),
      sensitivity = tp / (tp + fn),
      specificity = tn / (tn + fp),
      concordance = (tp + tn) / n,
      .groups = "drop"
    )
}

add_all_tf_summary <- function(perf_by_tf) {
  bind_rows(
    perf_by_tf,
    perf_by_tf %>%
      filter(!is.na(tf_group)) %>%
      group_by(event) %>%
      summarise(
        n = sum(n),
        tp = sum(tp),
        tn = sum(tn),
        fp = sum(fp),
        fn = sum(fn),
        sensitivity = tp / (tp + fn),
        specificity = tn / (tn + fp),
        concordance = (tp + tn) / n,
        .groups = "drop"
      ) %>%
      mutate(tf_group = "All")
  ) %>%
    mutate(tf_group = factor(tf_group, levels = c("Low TF", "High TF", "All")))
}

perf_tf_complete_all_evaluable <- bind_rows(
  add_all_tf_summary(
    summarise_bm_cfdna_event_performance(
      merged_CNA,
      "^(del1p|amp1q|del13q|del17p|hyperdiploid)_(BM|blood)$"
    )
  ),
  add_all_tf_summary(
    summarise_bm_cfdna_event_performance(
      merged_trans,
      "^(IGH_MAF|IGH_MYC|IGH_CCND1|IGH_FGFR3)_(BM|blood)$"
    )
  )
)

mutation_perf_by_tf_all_evaluable <- concordance_row %>%
  mutate(
    tf_group = case_when(
      replace_na(tf_group, "tf_unknown") == "high_tf" ~ "High TF",
      replace_na(tf_group, "tf_unknown") == "low_tf" ~ "Low TF",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(event = "Mutations", tf_group) %>%
  summarise(
    n = dplyr::n(),
    tp = sum(tp),
    tn = sum(tn),
    fp = sum(fp),
    fn = sum(fn),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    concordance = (tp + tn) / (tp + tn + fp + fn),
    .groups = "drop"
  )

perf_tf_complete_all_evaluable <- bind_rows(
  perf_tf_complete_all_evaluable,
  add_all_tf_summary(mutation_perf_by_tf_all_evaluable)
)

perf_long_all_evaluable <- perf_tf_complete_all_evaluable %>%
  filter(tf_group %in% c("High TF", "Low TF")) %>%
  pivot_longer(
    cols = c(concordance, sensitivity, specificity),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Percent = Value * 100,
    TF_group = factor(tf_group, levels = c("High TF", "Low TF"))
  )

event_order_all_evaluable <- perf_long_all_evaluable %>%
  filter(Metric == "concordance", TF_group == "High TF") %>%
  arrange(Percent) %>%
  pull(event)

perf_long_all_evaluable <- perf_long_all_evaluable %>%
  mutate(event = factor(event, levels = event_order_all_evaluable))

readr::write_csv(
  perf_long_all_evaluable,
  file.path(
    outdir,
    "Extended_Data_Figure_2B_all_evaluable_baseline_BM_cfDNA_performance_source_data.csv"
  )
)

p_3panel_all_evaluable <- ggplot(
  perf_long_all_evaluable,
  aes(x = Percent, y = event, group = event)
) +
  geom_point(aes(colour = TF_group, shape = TF_group), size = 3, stroke = 0.8) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    labeller = labeller(
      Metric = c(
        concordance = "Concordance",
        sensitivity = "Sensitivity",
        specificity = "Specificity"
      )
    )
  ) +
  scale_y_discrete(labels = pretty_events) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_colour_manual(
    values = c(
      "High TF" = viridis(2, end = 0.8)[1],
      "Low TF" = viridis(2, end = 0.8)[2]
    ),
    name = "Tumour fraction"
  ) +
  scale_shape_manual(
    values = c("High TF" = 16, "Low TF" = 16),
    name = "Tumour fraction"
  ) +
  labs(
    title = "BM vs cfDNA performance by ctDNA fraction",
    x = "Percent",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.box = "horizontal",
    panel.spacing.x = unit(1.2, "lines")
  )

all_evaluable_ed2b_path <- file.path(
  outdir,
  "Fig2C_BM_cfDNA_conc_sens_spec_byTF_all_evaluable_baseline.png"
)
ggsave(all_evaluable_ed2b_path, p_3panel_all_evaluable, width = 5.5, height = 4, dpi = 600)

ms_copy_artifact(
  source_path = all_evaluable_ed2b_path,
  artifact_id = "EDFIG2B_D_B",
  role = "all_evaluable_figure_panel_png",
  description = "Alternate all-evaluable baseline/diagnosis BM-vs-cfDNA concordance, sensitivity, and specificity panel for Extended Data Figure 2B.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)
ms_copy_artifact(
  source_path = file.path(
    outdir,
    "Extended_Data_Figure_2B_all_evaluable_baseline_BM_cfDNA_performance_source_data.csv"
  ),
  artifact_id = "EDFIG2B_D_B",
  role = "all_evaluable_source_data_csv",
  description = "Source data for the alternate all-evaluable baseline/diagnosis BM-vs-cfDNA performance panel.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

# Manuscript note: Extended Data Figure 2D is the VA-09 BM WGS chr1
# depth-ratio/copy-number line plot with the 1q FISH-probe interval marked.
# The panel is still preserved as part of the manually assembled final PDF,
# but its numerical source is no longer unresolved: the frozen Sequenza chr1
# segments and FISH-probe interval are exported by
# 5_1_Export_Locked_Figure_Source_Data.R.
edfig2_final_pdf_candidates <- c(
  file.path("Manuscript_Exports", "02_extended_data_figures", "Extended_Data_Figure_2", "final_artifacts", "Extended_Data_Figure_2.pdf"),
  file.path("Figures_Exported", "Final_Feb2026", "Extended_Data_Figure_2.pdf"),
  file.path("reproducible_workflow", "outputs", "frozen", "extended_figures", "Extended_Data_Figure_2.pdf")
)
edfig2_final_pdf <- edfig2_final_pdf_candidates[file.exists(edfig2_final_pdf_candidates)][1]

if (!is.na(edfig2_final_pdf)) {
  ms_copy_artifact(
    source_path = edfig2_final_pdf,
    artifact_id = "EDFIG2B_D_D",
    role = "manual_final_figure_pdf",
    description = "Extended Data Figure 2D: VA-09 BM WGS chr1 depth-ratio/copy-number line plot with the 1q FISH-probe interval, preserved as part of the frozen Extended Data Figure 2 PDF.",
    script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
  )
} else {
  warning(
    "Extended Data Figure 2D is a manually assembled panel and no frozen Extended Data Figure 2 PDF was found locally."
  )
}



# Get overall summary
# --- 0) Define event groups (adjust if the names differ) ---------------------
cna_events           <- c("amp1q","del13q","del17p","del1p","hyperdiploid")
translocation_events <- c("IGH_CCND1","IGH_FGFR3","IGH_MAF","IGH_MYC")
mutation_events      <- c("Mutations")

# --- 1) Per-event metrics  -----------------------------------
perf_event_metrics <- perf_tf_complete %>%
  mutate(
    # ensure counts are integers and define total_n robustly
    total_n = tp + tn + fp + fn,
    prevalence = ifelse(total_n > 0, (tp + fn) / total_n, NA_real_),
    ppv       = ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_),  # precision
    npv       = ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_),
    accuracy  = ifelse(total_n > 0, (tp + tn) / total_n, NA_real_),
    f1_score  = ifelse((ppv + sensitivity) > 0, 2 * (ppv * sensitivity) / (ppv + sensitivity), NA_real_),
    fdr       = ifelse((tp + fp) > 0, fp / (tp + fp), NA_real_),
    forate    = ifelse((tn + fn) > 0, fn / (tn + fn), NA_real_),  # false omission rate
    mcc = ifelse(
      (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0,
      (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),
      NA_real_
    ),
    # Jaccard on the confusion matrix
    jaccard = ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_)
  ) %>%
  transmute(
    category = case_when(
      event %in% cna_events ~ "CNA",
      event %in% translocation_events ~ "Translocation",
      event %in% mutation_events ~ "Mutation",
      TRUE ~ "Other"
    ),
    event,
    tf_group,
    n  = total_n,
    tp, tn, fp, fn,
    sensitivity, specificity, accuracy, ppv, npv, f1_score,
    fdr, forate, mcc, jaccard
  )

# --- 2) Category-level metrics (label as All_*) -------------------------------
perf_summary <- perf_tf_complete %>%
  mutate(
    category = case_when(
      event %in% cna_events ~ "CNA",
      event %in% translocation_events ~ "Translocation",
      event %in% mutation_events ~ "Mutation",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  group_by(category, tf_group) %>%
  summarise(
    tp = sum(tp, na.rm = TRUE),
    tn = sum(tn, na.rm = TRUE),
    fp = sum(fp, na.rm = TRUE),
    fn = sum(fn, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  mutate(
    n           = tp + tn + fp + fn,
    sensitivity = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
    specificity = ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_),
    ppv         = ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_),
    npv         = ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_),
    accuracy    = ifelse(n > 0, (tp + tn) / n, NA_real_),
    f1_score    = ifelse((ppv + sensitivity) > 0, 2 * (ppv * sensitivity) / (ppv + sensitivity), NA_real_),
    fdr         = ifelse((tp + fp) > 0, fp / (tp + fp), NA_real_),
    forate      = ifelse((tn + fn) > 0, fn / (tn + fn), NA_real_),
    mcc = ifelse(
      (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0,
      (tp*tn - fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),
      NA_real_
    ),
    jaccard     = ifelse((tp + fp + fn) > 0, tp / (tp + fp + fn), NA_real_),
    event = dplyr::case_when(
      category == "CNA" ~ "All_CNAs",
      category == "Translocation" ~ "All_Translocations",
      category == "Mutation" ~ "All_Mutations",
      TRUE ~ "All_Other"
    )
  ) %>%
  ungroup() %>%
  select(category, event, tf_group, n, tp, tn, fp, fn,
         sensitivity, specificity, accuracy, ppv, npv, f1_score,
         fdr, forate, mcc, jaccard)

# --- 3) Bind category-level + per-event into one table ------------------------
perf_combined <- bind_rows(perf_summary, perf_event_metrics) %>%
  arrange(match(category, c("Translocation","CNA","Mutation","Other")),
          match(tf_group, c("High TF","Low TF","high_tf","low_tf","BM","Blood","All","NA","NA ")),
          event)

# Round
perf_combined <- perf_combined %>%
  mutate(across(c(sensitivity, specificity, accuracy, ppv, npv, f1_score,
                  fdr, forate, mcc, jaccard),
                ~ round(.x, 3)))
# --- 4) Export ----------------------------------------------------------------
write_csv(perf_combined,
          "Final Tables and Figures/Supplentary_table_BM_vs_cfDNA_performance_by_event_and_category3.csv") ## this makes second part of supp table 2
saveRDS(perf_combined,
        "Final Tables and Figures/Supplentary_table_BM_vs_cfDNA_performance_by_event_and_category3.rds")


## above combined to supplementary table 2

## Reassemble here 
# Packages
library(openxlsx)

# Helper to turn literal "#NUM!" into NA for every column
fix_num <- function(df) {
  df %>%
    mutate(across(everything(),
                  ~ na_if(as.character(.x), "#NUM!")))
}

# ---- Panel A ----
panel_a <- read_csv(
  "Final Tables and Figures/Supp_table_concordance_summary_BM_cfDNA_updated.csv",
  na = c("", "NA", "#NUM!"),        # catch #NUM! while reading
  show_col_types = FALSE
) %>% fix_num()

# ---- Panels B & C ----
panel_b <- fix_num(as.data.frame(perf_combined))
panel_c <- fix_num(as.data.frame(metrics_combined))

# ---- Build workbook ----
wb <- createWorkbook()
addWorksheet(wb, "BM_vs_cfDNA_by_cohort")
addWorksheet(wb, "BM_vs_cfDNA_performance_detaile")
addWorksheet(wb, "FISH_vs_BM_and_cfDNA")

writeData(wb, "BM_vs_cfDNA_by_cohort",           panel_a)
writeData(wb, "BM_vs_cfDNA_performance_detaile", panel_b)
writeData(wb, "FISH_vs_BM_and_cfDNA",            panel_c)

# Nice-to-haves: auto column widths + freeze header rows
for (s in sheets <- c("BM_vs_cfDNA_by_cohort",
                      "BM_vs_cfDNA_performance_detaile",
                      "FISH_vs_BM_and_cfDNA")) {
  setColWidths(wb, sheet = s, cols = 1:200, widths = "auto")
  freezePane(wb, sheet = s, firstRow = TRUE)
}

# ---- Save workbook ----
out_xlsx <- "Final Tables and Figures/Supplementary_Table_2_BM_cfDNA_Concordance_updated.xlsx"
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Wrote: ", normalizePath(out_xlsx))

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 2
#
# What this is:
#   Multi-sheet workbook summarizing BM-vs-cfDNA concordance and performance,
#   plus FISH-vs-WGS metrics used for the Extended Data Figure 2 analysis.
#
# Why it is here:
#   This is the script-generated workbook counterpart for final Supplementary
#   Table 2. Historical filename drift is preserved in the source-map notes.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = out_xlsx,
  artifact_id = "STABLE2",
  role = "workbook_xlsx",
  description = "BM/cfDNA concordance workbook used as Supplementary Table 2.",
  script_name = "2_3_Feature_Concordance_And_Mutation_Counts.R"
)

### Export additional important things 
# =====================================================================
# FINAL EXPORTS – put this block right before the script ends
# =====================================================================
outdir <- "Final Tables and Figures/Baseline_concordance/"

## 1. R-objects (RDS) --------------------------------------------------
saveRDS(long,                 file.path(outdir, "FISH_WGS_long_call_table.rds"))
saveRDS(perf_tf_complete,     file.path(outdir, "BM_cfDNA_performance_byTF.rds"))
saveRDS(perf_long,            file.path(outdir, "BM_cfDNA_performance_long_forPlot.rds"))
saveRDS(perf_summary,            file.path(outdir, "BM_cfDNA_performance_summary.rds"))
saveRDS(tf_plot_df,            file.path(outdir, "FISH_BM_cfDNA_performance.rds"))

## 2. Simple CSV / XLSX tables ----------------------------------------
readr::write_csv(baseline_summary,
                 file.path(outdir, "baseline_mutation_summary.csv"))

readr::write_csv(pvals,
                 file.path(outdir, "baseline_mutation_count_pvals.csv"))
readr::write_csv(mean_jaccard,
                 file.path(outdir, "mean_jaccard_baseline.csv"))


## Now get the actual per-sample calls 
dat_small <- dat_base_with_FISH %>%
  select(
    # --- Clinical context ---
    Patient,
    Timepoint,
    timepoint_info,
    Cohort,
    WGS_Tumor_Fraction_Blood_plasma_cfDNA,
    Ploidy_Estimate,

    # --- FISH Calls ---
    DEL_1P, AMP_1Q, DEL_17P, T_4_14, T_11_14, T_14_16,
    
    # --- Probe-level and alteration calls ---
    probe_call_amp1q_blood,
    probe_call_amp1q_BM,
    probe_call_del17p_blood,
    probe_call_del17p_BM,
    probe_call_del1p_blood,
    probe_call_del1p_BM,
    is_altered_at_probe_amp1q_blood,
    is_altered_at_probe_amp1q_BM,
    is_altered_at_probe_del17p_blood,
    is_altered_at_probe_del17p_BM,
    is_altered_at_probe_del1p_blood,
    is_altered_at_probe_del1p_BM,
    
    # --- WGS arm-level CNAs ---
    WGS_del1p_BM_cells,
    WGS_del1p_Blood_plasma_cfDNA,
    WGS_amp1q_BM_cells,
    WGS_amp1q_Blood_plasma_cfDNA,
    WGS_del17p_BM_cells,
    WGS_del17p_Blood_plasma_cfDNA,
    
    # --- WGS translocations ---
    WGS_IGH_FGFR3_BM_cells,
    WGS_IGH_FGFR3_Blood_plasma_cfDNA,
    WGS_IGH_CCND1_BM_cells,
    WGS_IGH_CCND1_Blood_plasma_cfDNA,
    WGS_IGH_MAF_BM_cells,
    WGS_IGH_MAF_Blood_plasma_cfDNA
  ) %>%
  rename(
    # --- Core clinical ---
    Ploidy_estimate_BM = Ploidy_Estimate,
    
    # --- FISH Calls ---
    FISH_Call_DEL_1P   = DEL_1P,
    FISH_Call_AMP_1Q   = AMP_1Q,
    FISH_Call_DEL_17P  = DEL_17P,
    FISH_Call_T_4_14   = T_4_14,
    FISH_Call_T_11_14  = T_11_14,
    FISH_Call_T_14_16  = T_14_16,
    
    # --- Probe-level (add WGS_ prefix) ---
    WGS_probe_call_amp1q_blood  = probe_call_amp1q_blood,
    WGS_probe_call_amp1q_BM     = probe_call_amp1q_BM,
    WGS_probe_call_del17p_blood = probe_call_del17p_blood,
    WGS_probe_call_del17p_BM    = probe_call_del17p_BM,
    WGS_probe_call_del1p_blood  = probe_call_del1p_blood,
    WGS_probe_call_del1p_BM     = probe_call_del1p_BM,
    
    WGS_is_altered_at_probe_amp1q_blood  = is_altered_at_probe_amp1q_blood,
    WGS_is_altered_at_probe_amp1q_BM     = is_altered_at_probe_amp1q_BM,
    WGS_is_altered_at_probe_del17p_blood = is_altered_at_probe_del17p_blood,
    WGS_is_altered_at_probe_del17p_BM    = is_altered_at_probe_del17p_BM,
    WGS_is_altered_at_probe_del1p_blood  = is_altered_at_probe_del1p_blood,
    WGS_is_altered_at_probe_del1p_BM     = is_altered_at_probe_del1p_BM,
    
    # --- Arm-level CNAs ---
    WGS_arm_del1p_BM_cells              = WGS_del1p_BM_cells,
    WGS_arm_del1p_Blood_plasma_cfDNA    = WGS_del1p_Blood_plasma_cfDNA,
    WGS_arm_amp1q_BM_cells              = WGS_amp1q_BM_cells,
    WGS_arm_amp1q_Blood_plasma_cfDNA    = WGS_amp1q_Blood_plasma_cfDNA,
    WGS_arm_del17p_BM_cells             = WGS_del17p_BM_cells,
    WGS_arm_del17p_Blood_plasma_cfDNA   = WGS_del17p_Blood_plasma_cfDNA
  )

# Fix cohort 
dat_small <- dat_small %>%
  mutate(
    Cohort = case_when(
      Cohort == "Frontline induction-transplant" ~ "Training Cohort",
      Cohort == "Non-frontline" ~ "Test Cohort",
      TRUE ~ Cohort  # keep other labels unchanged if any exist
    )
  )

## Rename patient IDs 
# --- 0) Load ID map
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

# Make a joinable version that uses the same Patient key as the output (New_ID if available)
dat_small <- dat_small %>%
  left_join(id_map, by = c("Patient" = "Patient")) %>%
  mutate(Patient = coalesce(New_ID, Patient)) 

# Export arm-level CNA calls as explicit binary integers for Supplementary
# Table 2F. Preserve unevaluable/missing compartments as NA rather than
# converting them to 0, because 0 means an evaluated negative call.
supp_table_2f_binary_arm_columns <- c(
  "WGS_arm_del1p_BM_cells",
  "WGS_arm_del1p_Blood_plasma_cfDNA",
  "WGS_arm_amp1q_BM_cells",
  "WGS_arm_amp1q_Blood_plasma_cfDNA",
  "WGS_arm_del17p_BM_cells",
  "WGS_arm_del17p_Blood_plasma_cfDNA"
)
missing_supp_table_2f_binary_columns <- setdiff(
  supp_table_2f_binary_arm_columns,
  names(dat_small)
)
if (length(missing_supp_table_2f_binary_columns)) {
  stop(
    "Supplementary Table 2F is missing expected arm-level WGS columns: ",
    paste(missing_supp_table_2f_binary_columns, collapse = ", "),
    call. = FALSE
  )
}
dat_small <- dat_small %>%
  mutate(
    across(
      all_of(supp_table_2f_binary_arm_columns),
      ~ if_else(is.na(.x), NA_integer_, as.integer(.x))
    )
  )


## 3. One XLSX workbook with the key performance tables ----------------
## Keep the analysis output and the source consumed by the submission-package
## assembler synchronized. Previously the assembler read an older workbook in
## Output_tables_2025, allowing the frontline-only table to persist even after
## the revision-inclusive analysis had been generated elsewhere.
fish_concordance_for_publication <- concordance_tbl %>%
  mutate(cohort = dplyr::coalesce(cohort, "All evaluable"))
fish_probe_concordance_for_publication <- concordance_tbl_at_FISH_probe %>%
  mutate(cohort = dplyr::coalesce(cohort, "All evaluable"))

supplementary_table_2_sheets <- list(
  A_BM_cfDNA_perf_byTF = perf_combined,
  B_FISH_vs_Sample_Concordance = fish_concordance_for_publication,
  C_FISH_vs_Sample_At_Probe = fish_probe_concordance_for_publication,
  D_Individual_Calls = dat_small %>% select(-New_ID)
)
supplementary_table_2_sheets <- relabel_publication_workbook_tables(
  supplementary_table_2_sheets,
  "Supplementary Table 2"
)

supplementary_table_2_analysis_path <- file.path(
  outdir,
  "Supplementary_Table_2_SV_CNA_performance_summary_updated2.xlsx"
)
supplementary_table_2_package_source_path <- file.path(
  "Output_tables_2025",
  "Supplementary_Table_2_SV_CNA_performance_summary_updated2.xlsx"
)

writexl::write_xlsx(
  supplementary_table_2_sheets,
  path = supplementary_table_2_analysis_path
)
writexl::write_xlsx(
  supplementary_table_2_sheets,
  path = supplementary_table_2_package_source_path
)

message("✓ Additional outputs written to ", outdir)
