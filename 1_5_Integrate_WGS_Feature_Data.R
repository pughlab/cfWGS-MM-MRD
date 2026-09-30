# =============================================================================
# Script: 1_5_Integrate_WGS_Feature_Data.R
#
# Description:
#   Integrates copy-number (CNA), translocation, tumor-fraction, and mutation
#   features into a unified `All_feature_data` table for downstream analysis.
#   Steps:
#     1. Load and harmonize clinical metadata (sample IDs).
#     2. Read pre-exported CNA and translocation tables.
#     3. Read ichorCNA tumor-fraction data.
#     4. Load mutation calls from combined MAFs (BM and blood) via `read.maf()`.
#     5. Merge CNA + translocation with metadata and tumor fraction.
#     6. Subset MAF objects to `myeloma_genes`, classify VAF and mutation types,
#        and summarize per sample.
#     7. Join mutation summary into the feature table, fill missing values.
#     8. Apply the final rule-based `Evidence_of_Disease` definition, followed
#        by the retained manual IGV-review override.
#     9. Save the mutation helpers and integrated feature tables as RDS and TSV.
#
# Inputs:
#   • combined_clinical_data_updated_April2025.csv
#   • Jan2025_exported_data/cna_data_ichorCNA.rds
#   • Jan2025_exported_data/cna_data_from_sequenza_400_updated.rds
#   • Jan2025_exported_data/translocation_data_cytoband_updated.rds
#   • Oct 2024 data/tumor_fraction_cfWGS.txt
#   • combined_maf_temp_blood_Jan2025.maf
#   • combined_maf_temp_bm_Jan2025.maf (preferred) or
#     combined_maf_temp_bm_May2025.maf (legacy fallback)
#   • Jan2025_exported_data/FISH_data_from_sequenza_400_updated.rds
#   • Jan2025_exported_data/FISH_probe_calls_bin_cytoband_ichorCNA.rds
#   • Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm
#   • Optional historical mutation recovery source:
#     Jan2025_exported_data/mutation_export_updated_more_info.rds
#
# Analysis unit:
#   The final table contains one retained row per biological sample, keyed by
#   Patient + Sample_ID + Sample_type + Timepoint + timepoint_info. When multiple
#   BAM rows map to that key, the deterministic ranking near the end keeps the
#   row with the strongest disease signal and exports the removed rows for audit.
#
# Active downstream outputs:
#   • Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#   • Jan2025_exported_data/All_feature_data_Sep2025_updated2.txt
#       - primary integrated WGS feature table consumed by 2_0, 2_2, 2_3,
#         2_4, 3_1, and downstream manuscript analyses.
#   • Jan2025_exported_data/CNA_translocation_Sep2025_updated2.rds
#   • Jan2025_exported_data/CNA_translocation_Sep2025_updated2.txt
#       - CNA/translocation-only helper consumed by baseline heatmap scripts.
#   • Jan2025_exported_data/CNA_at_FISH_sites_combined.rds
#   • Jan2025_exported_data/CNA_at_FISH_sites_combined.txt
#       - FISH-probe CNA helper consumed by 2_3 concordance analyses.
#   • Jan2025_exported_data/mutation_export_updated2.rds
#   • Jan2025_exported_data/mutation_export_updated.txt
#       - compact myeloma-panel mutation helper table.
#   • Jan2025_exported_data/mutation_export_updated_more_info2.rds
#   • Jan2025_exported_data/mutation_export_updated_more_info.txt
#       - expanded mutation/QC helper consumed by 2_3.
#
# Support/QC outputs:
#   • Output_tables_2025/feature_integration_support/samples_missing_metadata_after_cna_translocation_join.csv
#   • Output_tables_2025/feature_integration_support/excluded_bm_ichor_cna_rows.csv
#   • Output_tables_2025/feature_integration_support/excluded_bm_ichor_fish_cna_rows.csv
#   • Output_tables_2025/feature_integration_support/cna_translocation_identity_column_reconciliation_audit.csv
#   • Output_tables_2025/feature_integration_support/pre_dedup_exact_cna_translocation_rows.csv
#   • Output_tables_2025/feature_integration_support/active_cna_translocation_duplicate_sample_rows.csv
#   • Output_tables_2025/feature_integration_support/verified_baseline_mutation_rows_preserved_from_previous_helper.csv
#   • Output_tables_2025/feature_integration_support/mutation_samples_missing_cna_translocation_rows.csv
#   • Output_tables_2025/feature_integration_support/all_feature_data_biological_replicates_removed.csv
#   • Output_tables_2025/feature_integration_support/evidence_rule_del13q_sensitivity_check.csv
#
# Dependencies:
#   library(dplyr)
#   library(tidyr)
#   library(readr)
#   library(stringr)
#   library(readxl)
#   library(maftools)
#   library(purrr)
#
# Usage:
#   source("1_5_Integrate_WGS_Feature_Data.R")
#
# How to run from the repository root:
#   Rscript 1_5_Integrate_WGS_Feature_Data.R
#
# Failure behavior:
#   The script stops before writing the final integrated table when a retained
#   mutation-positive sample has no matching CNA/translocation identity row.
#   Resolve the upstream sample key and compare the rebuilt table with the
#   retained output rather than converting an unknown CNA/translocation row to
#   a negative call.
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script merges WGS-derived mutation, CNA,
#     translocation, and tumor-fraction features used by baseline heatmaps,
#     concordance tables, model training, and longitudinal analyses.
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#
# Pipeline role and reproducibility note:
#   This is the main WGS feature-integration step. It should be rerun whenever
#   upstream mutation, CNA, translocation, tumor-fraction, or reviewed IGV
#   translocation inputs change. It performs deterministic joins and rule-based
#   feature calls; no model fitting, threshold selection, or stochastic sampling
#   is performed here.
#

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(readxl)
library(maftools)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

# Define export directories.
export_dir <- "Jan2025_exported_data"
support_table_dir <- file.path("Output_tables_2025", "feature_integration_support")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(support_table_dir, recursive = TRUE, showWarnings = FALSE)

support_file <- function(filename) {
  file.path(support_table_dir, filename)
}

reconcile_joined_identity_columns <- function(df, columns, audit_filename) {
  audit_rows <- purrr::map_dfr(columns, function(base_col) {
    candidates <- intersect(c(base_col, paste0(base_col, ".x"), paste0(base_col, ".y")), names(df))
    if (length(candidates) <= 1L) {
      return(tibble())
    }

    values <- df %>%
      mutate(.row_id = row_number()) %>%
      select(any_of(c("Sample", "Bam", "Patient", "Sample_ID", "Timepoint", "timepoint_info")),
             .row_id, all_of(candidates)) %>%
      pivot_longer(
        cols = all_of(candidates),
        names_to = "source_column",
        values_to = "source_value",
        values_transform = list(source_value = as.character)
      ) %>%
      filter(!is.na(source_value), nzchar(source_value))

    values %>%
      group_by(.row_id) %>%
      filter(n_distinct(source_value) > 1L) %>%
      ungroup() %>%
      mutate(canonical_column = base_col) %>%
      select(canonical_column, everything())
  })

  readr::write_csv(audit_rows, support_file(audit_filename))

  for (base_col in columns) {
    candidates <- intersect(c(base_col, paste0(base_col, ".x"), paste0(base_col, ".y")), names(df))
    if (length(candidates) > 1L) {
      df[[base_col]] <- dplyr::coalesce(!!!df[candidates])
      df <- df %>% select(-any_of(setdiff(candidates, base_col)))
    }
  }

  df
}

require_files <- function(paths, description) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing required ", description, ":\n  ",
      paste(missing_paths, collapse = "\n  "),
      call. = FALSE
    )
  }
  invisible(paths)
}

require_files(
  c(
    "combined_clinical_data_updated_April2025.csv",
    file.path(export_dir, "cna_data_ichorCNA.rds"),
    file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"),
    file.path(export_dir, "translocation_data_cytoband_updated.rds"),
    "Oct 2024 data/tumor_fraction_cfWGS.txt",
    file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"),
    file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"),
    file.path(export_dir, "Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm")
  ),
  "feature-integration inputs"
)

# Load clinical metadata.
# This table provides the patient, timepoint, sample type, and BAM identifiers
# needed to connect genomic caller outputs back to the manuscript cohorts.
metada_df_mutation_comparison <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
) %>%
  mutate(
    # Tumor_Sample_Barcode is the caller-neutral sample key used to merge
    # metadata with mutation, CNA, and translocation outputs. It removes assay
    # tokens and processing suffixes while retaining the biological sample ID.
    Tumor_Sample_Barcode = Bam %>%
      str_remove_all("_PG|_WG") %>%
      str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""),
    Bam_clean_tmp = str_remove(Bam, "\\.bam$")
  ) %>%
  mutate(
    Sample_ID = if_else(
      (is.na(Sample_ID) | !nzchar(as.character(Sample_ID))) &
        !is.na(Patient) & nzchar(as.character(Patient)) &
        !is.na(Timepoint) & nzchar(as.character(Timepoint)) &
        !is.na(Sample_type) & nzchar(as.character(Sample_type)),
      paste(Patient, paste0("T", Timepoint), Sample_type, sep = "_"),
      as.character(Sample_ID)
    )
  )

# Load CNA, translocation, tumor-fraction, and mutation inputs.
# cna_data_ichorCNA.rds  - arm-level binary calls from 1_4_Process_CNA_Data.R
# cna_data_from_sequenza - arm-level calls from 1_4A using the Sequenza CNV caller;
#   Sequenza provides more accurate purity/ploidy estimates
#   and takes precedence over ichorCNA when both callers processed the same BAM.
# translocation_data     - Ig translocation binary flags from 1_3
# tumor_fraction_cfWGS   - ichorCNA-estimated plasma tumor fraction per BAM;
#   used in the Evidence_of_Disease classifier and as a continuous feature
#   in downstream models (2_0, 3_1, 4_1)
cna_data           <- readRDS(file.path(export_dir, "cna_data_ichorCNA.rds"))
cna_data_sequenza <- readRDS(file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"))
translocation_data <- readRDS(file.path(export_dir, "translocation_data_cytoband_updated.rds"))
tumor_fraction <- read_tsv("Oct 2024 data/tumor_fraction_cfWGS.txt")
spring2026_ichor_params <- spring2026_revision_files(
  # Spring 2026 ichorCNA params are optional revision inputs. When present, they
  # add tumor-fraction/ploidy estimates for newly integrated patient samples.
  # M4CHIP dilution rows are excluded here because dilution tumor fractions are
  # handled by the dilution-series workflow, not the main patient feature table.
  "",
  "^Ichor_CNA_combined_params_summary_tumor_fraction_sex[.]tsv$"
)[1]
if (!is.na(spring2026_ichor_params) && file.exists(spring2026_ichor_params)) {
  spring2026_tumor_fraction <- read_tsv(spring2026_ichor_params, show_col_types = FALSE) %>%
    filter(!str_detect(file, "^M4CHIP_")) %>%
    transmute(
      # file contains the ichorCNA params filename; removing .params.txt returns
      # the BAM-like key used by the historical tumor_fraction table.
      Bam = str_remove(file, "[.]params[.]txt$"),
      Tumor_fraction = suppressWarnings(as.numeric(tumor_fraction)),
      Ploidy = suppressWarnings(as.numeric(ploidy))
    )
  tumor_fraction <- bind_rows(tumor_fraction, spring2026_tumor_fraction) %>%
    # If a BAM already exists in the historical table, keep that first entry.
    # The revision params extend missing rows rather than overriding historical
    # manually reviewed tumor-fraction calls.
    distinct(Bam, .keep_all = TRUE)
}

# Load mutation MAF objects.
# The upstream mutation-processing script (1_2_Process_Mutation_Data.R) writes
# project-local combined MAF files. The command-line manuscript workflow expects
# those local files to be present after Stage 1_2 so this script does not depend
# on a user-specific external drive or workstation path.
resolve_first_existing <- function(label, candidates) {
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "Could not find ", label, ". Checked:\n  ",
      paste(candidates, collapse = "\n  "),
      "\nRun 1_2_Process_Mutation_Data.R first or provide the expected MAF input.",
      call. = FALSE
    )
  }
  hit
}

blood_maf_path <- resolve_first_existing(
  "combined cfDNA/blood MAF",
  c(
    "combined_maf_temp_blood_Jan2025.maf"
  )
)
bm_maf_path <- resolve_first_existing(
  "combined baseline BM MAF",
  c(
    "combined_maf_temp_bm_Jan2025.maf",
    "combined_maf_temp_bm_May2025.maf"
  )
)
maf_object_blood <- read.maf(blood_maf_path)
maf_object_bm    <- read.maf(bm_maf_path)


### Integrate arm-level CNA and IgH translocation features
# This section produces `CNA_translocation`, the genomic feature table before
# mutation summaries are added. It is exported at the end because some
# downstream checks use the CNA/translocation-only intermediate directly.
#saveRDS(CNA_translocation, "CNA_translocation_original_Feb2025.rds")

# Define caller-source rules by sample matrix before merging CNA calls.
# Analyses reported in the manuscript use Sequenza for BM-cell CNA calls when available
# and ichorCNA for cfDNA/non-BM CNA calls. BM rows found only in ichorCNA are
# retained as fallback CNA evidence rather than being forced missing. Only
# ichorCNA rows for BAMs that also have Sequenza calls are dropped, so Sequenza
# takes precedence without losing BM samples lacking Sequenza output.
metadata_sample_types <- metada_df_mutation_comparison %>%
  transmute(
    Sample = as.character(Bam_clean_tmp),
    Sample_type = as.character(Sample_type),
    Patient = as.character(Patient),
    Sample_ID = as.character(Sample_ID),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info)
  ) %>%
  filter(!is.na(Sample), nzchar(Sample)) %>%
  distinct(Sample, .keep_all = TRUE)

bm_sample_keys <- metadata_sample_types %>%
  filter(Sample_type == "BM_cells") %>%
  distinct(Sample)

# 1) Rename Sequenza ID to Sample so downstream stays consistent
# Sequenza (1_4A) stores the sample identifier in Bam_clean_tmp rather
# than Sample because it derives from a different pipeline input file.
# Renaming here ensures it joins correctly on the 'Sample' key used
# throughout all downstream scripts (1_5 onward).
cna_seq_renamed <- cna_data_sequenza %>% select(-any_of("Sample")) %>%
  rename(Sample = Bam_clean_tmp) %>%
  # Align columns to the ichorCNA-derived table before row-binding.
  select(any_of(names(cna_data)))

# 2) Drop ichorCNA rows only when they are superseded by Sequenza.
# Where both callers processed the same BAM, Sequenza calls take precedence and
# the ichorCNA row is dropped to prevent duplicate Sample keys. BM-cell rows
# without a Sequenza replacement are retained so samples with valid CNA evidence
# are not misrepresented as CNA-missing downstream.
overlap_samples <- intersect(cna_data$Sample, cna_seq_renamed$Sample)

cna_ichor_bm_rows_excluded <- cna_data %>%
  semi_join(bm_sample_keys, by = "Sample") %>%
  filter(Sample %in% overlap_samples) %>%
  left_join(metadata_sample_types, by = "Sample") %>%
  arrange(Patient, Sample_ID, Sample)

readr::write_csv(
  cna_ichor_bm_rows_excluded,
  support_file("excluded_bm_ichor_cna_rows.csv")
)

cna_ichor_filtered <- cna_data %>%
  filter(!Sample %in% overlap_samples)

# 3) Combine CNA calls (Sequenza takes precedence)
cna_combined <- bind_rows(cna_ichor_filtered, cna_seq_renamed)

# (Optional) In case of any lingering duplicates, keep first occurrence
# cna_combined <- cna_combined %>% distinct(Sample, .keep_all = TRUE)

# 4) Join to translocation data (expects translocation_data$Sample to match Bam_clean_tmp)
CNA_translocation <- full_join(cna_combined, translocation_data, by = "Sample")

# QC
message("Combined CNA rows: ", nrow(cna_combined),
        " | After translocation join: ", nrow(CNA_translocation))

## Attach clinical/sample metadata to the CNA/translocation rows.
# First, remove the '.bam' suffix from the 'Bam' column in
# 'metada_df_mutation_comparison'.
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix

# Perform the left join
CNA_translocation <- left_join(CNA_translocation, 
                               metada_df_mutation_comparison, 
                               by = c("Sample" = "Bam_clean_tmp"))
CNA_translocation <- reconcile_joined_identity_columns(
  CNA_translocation,
  columns = c("Date_of_sample_collection", "Study", "Sample_ID"),
  audit_filename = "cna_translocation_identity_column_reconciliation_audit.csv"
)

## Add the tumor fraction info 
# Tumor_Fraction is the ichorCNA estimate of the fraction of cfDNA
# molecules derived from tumor. A single BAM may appear multiple times
# in tumor_fraction if it was run with different ichorCNA parameter sets
# (e.g. centromere-masked vs full-genome). Taking the maximum recovers
# the largest reported estimate under the historical project rule. This choice
# favors sensitivity and can be optimistic when rows represent alternative
# parameter fits rather than true replicates; the source rows should be retained
# for audit.
# Tumor_Fraction drives Evidence_of_Disease Tier 3 and is used as a
# continuous predictor in 4_1_Survival_Analysis.R.
# Keep only the max Tumor_Fraction for each Bam in tumor_fraction
tumor_fraction_max <- tumor_fraction %>%
  group_by(Bam) %>%
  summarise(Tumor_Fraction = max(Tumor_fraction, na.rm = TRUE))  # Ensure to handle NA values


CNA_translocation <- left_join(CNA_translocation, 
                               tumor_fraction_max, 
                               by = c("Sample" = "Bam"))



### Build the FISH-probe/cytoband CNA helper table
# This repeats the Sequenza-over-ichor merge for calls at clinically relevant
# probe/cytoband sites. The result is not the main All_feature_data table, but
# it supports FISH concordance and copy-number helper analyses.

# 1) Rename Sequenza ID to Sample so downstream stays consistent
FISH_sequenza <- readRDS(file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"))
FISH_ichor <- readRDS(file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"))

FISH_sequenza <- FISH_sequenza %>% select(-any_of("Sample")) %>%
  rename(Sample = Bam_clean_tmp)

# 2) Drop ichorCNA probe rows that are superseded by Sequenza
overlap_samples <- intersect(FISH_ichor$Sample, FISH_sequenza$Sample)

FISH_ichor_bm_rows_excluded <- FISH_ichor %>%
  semi_join(bm_sample_keys, by = "Sample") %>%
  filter(Sample %in% overlap_samples) %>%
  left_join(metadata_sample_types, by = "Sample") %>%
  arrange(Patient, Sample_ID, Sample)

readr::write_csv(
  FISH_ichor_bm_rows_excluded,
  support_file("excluded_bm_ichor_fish_cna_rows.csv")
)

FISH_ichor_filtered <- FISH_ichor %>%
  filter(!Sample %in% overlap_samples)

# 3) Combine CNA calls (Sequenza takes precedence)
FISH_CNA_combined <- bind_rows(FISH_ichor_filtered, FISH_sequenza)

## Attach clinical/sample metadata to the probe-level CNA rows.
# First, remove the '.bam' suffix from the 'Bam' column in
# 'metada_df_mutation_comparison'.
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix

# Perform the left join
FISH_CNA_combined <- left_join(FISH_CNA_combined, 
                               metada_df_mutation_comparison, 
                               by = c("Sample" = "Bam_clean_tmp"))

## Add the tumor fraction info 
# Keep only the max Tumor_Fraction for each Bam in tumor_fraction
tumor_fraction_max <- tumor_fraction %>%
  group_by(Bam) %>%
  summarise(Tumor_Fraction = max(Tumor_fraction, na.rm = TRUE))  # Ensure to handle NA values


FISH_CNA_combined <- left_join(FISH_CNA_combined, 
                               tumor_fraction_max, 
                               by = c("Sample" = "Bam"))

# Recalculate binary probe-level WGS alteration labels.
# Extended loss_labels includes LOSS and CNLOH (vs. only HETD/HOMD in
# 1_4_Process_CNA_Data.R) under the historical project rule. CNLOH is a
# copy-neutral allelic imbalance, not a deletion, and may be negative by FISH;
# downstream concordance must therefore interpret this flag as a WGS abnormality
# at the locus rather than literal agreement with a deletion probe.
gain_labels <- c("GAIN","AMP","HLAMP")
loss_labels <- c("HOMD","HETD","LOSS","CNLOH")

FISH_CNA_combined <- FISH_CNA_combined %>%
  mutate(
    # normalize case
    across(c(probe_call_amp1q, probe_call_del17p, probe_call_del1p),
           ~ toupper(as.character(.))),
    
    # 1 if gain label for amp1q; else 0 (including NA/NEUT)
    is_altered_at_probe_amp1q  = as.integer(!is.na(probe_call_amp1q)  &
                                              probe_call_amp1q  %in% gain_labels),
    
    # 1 if loss label for 17p and 1p; else 0
    is_altered_at_probe_del17p = as.integer(!is.na(probe_call_del17p) &
                                              probe_call_del17p %in% loss_labels),
    is_altered_at_probe_del1p  = as.integer(!is.na(probe_call_del1p)  &
                                              probe_call_del1p  %in% loss_labels)
  )


## Export the FISH-probe/cytoband CNA helper table.
saveRDS(FISH_CNA_combined, file = file.path(export_dir, "CNA_at_FISH_sites_combined.rds"))
write.table(FISH_CNA_combined,
            file = file.path(export_dir, "CNA_at_FISH_sites_combined.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Active FISH-probe CNA helper written: ",
        file.path(export_dir, "CNA_at_FISH_sites_combined.rds"))



# Filter rows where Tumor_Sample_Barcode is NA
## Can get to these later if decide to include**
samples_with_na_barcode <- CNA_translocation %>%
  filter(is.na(Tumor_Sample_Barcode)) %>% 
  unique()

readr::write_csv(
  samples_with_na_barcode,
  support_file("samples_missing_metadata_after_cna_translocation_join.csv")
)
message("Support QC table written: ",
        support_file("samples_missing_metadata_after_cna_translocation_join.csv"))

## Keep unmatched genomic rows out of the active feature table.
#
# The full-join above intentionally captures CNA/translocation calls even when
# metadata are missing, but rows without Patient/Sample_ID/timepoint identity
# cannot be joined correctly to clinical MRD, longitudinal samples, or
# manuscript cohorts. Leaving them in `All_feature_data` creates silent junk
# rows that downstream scripts may ignore inconsistently. The unmatched rows
# remain fully exported in the support QC table above for investigation.
CNA_translocation <- CNA_translocation %>%
  filter(
    !is.na(Patient),
    !is.na(Sample_ID),
    !is.na(Timepoint),
    !is.na(timepoint_info)
  )

pre_dedup_exact_cna_translocation_rows <- CNA_translocation %>%
  add_count(across(everything()), name = "exact_duplicate_rows") %>%
  filter(exact_duplicate_rows > 1L) %>%
  distinct() %>%
  arrange(Sample, Patient, Sample_ID)

readr::write_csv(
  pre_dedup_exact_cna_translocation_rows,
  support_file("pre_dedup_exact_cna_translocation_rows.csv")
)

## Remove exact duplicated BAM-level rows before auditing true active duplicate
# sample keys. Non-identical duplicate Sample rows remain visible below and in
# the downstream dataflow audit.
CNA_translocation <- CNA_translocation %>%
  distinct()

active_duplicate_sample_rows <- CNA_translocation %>%
  group_by(Sample) %>%
  filter(!is.na(Sample), dplyr::n() > 1L) %>%
  ungroup() %>%
  arrange(Sample, Patient, Sample_ID)

readr::write_csv(
  active_duplicate_sample_rows,
  support_file("active_cna_translocation_duplicate_sample_rows.csv")
)

##### Extract mutation data for specified genes
# These 30 recurrently mutated genes constitute the MM-specific panel
# used for somatic mutation detection. Selection criteria:
#   (a) established high-risk markers: TP53, RB1, ATM;
#   (b) MAPK pathway drivers (most common class in relapsed MM): KRAS, NRAS, BRAF;
#   (c) Ig translocation partner genes: CCND1, FGFR3, MMSET, MYC, BCL2;
#   (d) RNA stability / NF-kB regulators frequently acquired at relapse.
# Limiting subsetMaf to this panel speeds processing and ensures only
# clinically interpretable mutations inform the Evidence_of_Disease flag.
myeloma_genes <- c(
  "TP53",    # ~10-15%; high-risk MM
  "KRAS",    # ~20-25%; MAPK/ERK pathway
  "NRAS",    # ~20-25%; MAPK/ERK pathway
  "BRAF",    # ~5-10%; MAPK/ERK pathway
  "FAM46C",  # ~10-15%; RNA stability
  "DIS3",    # ~10-15%; RNA degradation
  "CYLD",    # ~5-10%; NF-κB regulator
  "ATM",     # ~5%; DNA damage repair
  "CCND1",   # ~15-20%; t(11;14), cyclin D1
  "MYC",     # ~15-20%; MYC translocations
  "RB1",     # ~5-10%; cell cycle control
  "TRAF3",   # ~5%; NF-κB regulator
  "IRF4",    # ~5%; plasma cell differentiation
  "FGFR3",   # ~10-15%; t(4;14), receptor tyrosine kinase
  "MMSET",   # ~10-15%; t(4;14), epigenetics
  "BCL2",    # ~15-20%; t(11;14), venetoclax target
  "IKZF1",   # ~5%; transcription regulation
  "IKZF3",   # ~5%; transcription regulation
  "CDKN2C",  # ~5-10%; cell cycle regulation
  "KDM6A",   # ~5%; epigenetics
  "SETD2",   # ~5%; histone modification
  "PTEN",    # ~5%; tumor suppressor
  "XBP1",    # ~5%; plasma cell differentiation
  "MAX",     # ~5%; MYC regulatory partner
  "SP140",   # ~5%; immune dysregulation
  "NFKBIA",  # ~5%; NF-κB inhibitor
  "NFKB2",   # ~5%; NF-κB activator
  "PRDM1",   # ~5%; plasma cell differentiation
  "EGR1",    # ~5%; early growth response
  "LTB"      # <5%; rare but part of NF-κB signaling
)


## Combine BM and blood mutation data and export the primary mutation table.
maf_subset <- subsetMaf(maf = maf_object_bm, genes = myeloma_genes, includeSyn = FALSE)
maf_subset_blood <- subsetMaf(maf = maf_object_blood, genes = myeloma_genes, includeSyn = FALSE)


# Bone-marrow
# filter(t_depth > 10): retains variant calls with at least 11 total reads at
# the locus; calls at exactly depth 10 are excluded. This threshold is an
# operational support filter; it does not establish analytical sensitivity and
# should be interpreted with the sample's assay depth and caller QC.
temp_bm <- maf_subset@data %>%
  filter(t_depth > 10) %>%                           # only well-supported calls
  mutate(
    Sample          = sub("\\.bam$", "", Bam),       # drop .bam
    Mutation_cDNA   = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic= paste(Chromosome,
                            Start_Position,
                            Reference_Allele,
                            Tumor_Seq_Allele2,
                            sep = "_"),
    # Mutation_Type: simplified 4-class encoding used in the Evidence_of_Disease
    # logic and in the baseline WGS heatmap (Extended Data Figure 1).
    #   Truncating  - nonsense or frameshift (almost always loss-of-function)
    #   Missense    - single amino-acid substitution or in-frame indel
    #   Splice_Site - affects canonical donor/acceptor; typically loss-of-function
    #   Other       - UTR, intronic, silent; excluded from most downstream analyses
    Mutation_Type   = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")       ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")           ~ "Missense",
      Variant_Classification == "Splice_Site"                          ~ "Splice_Site",
      TRUE                                                             ~ "Other"
    )
  ) %>%
  select(
    Tumor_Sample_Barcode, Sample, Hugo_Symbol,
    Mutation_cDNA, Mutation_Genomic,
    Mutation_Type, t_depth, VAF
  ) %>%
  distinct()

# Blood
temp_blood <- maf_subset_blood@data %>%
  filter(t_depth > 10) %>%
  mutate(
    Sample          = sub("\\.bam$", "", Bam),
    Mutation_cDNA   = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic= paste(Chromosome,
                            Start_Position,
                            Reference_Allele,
                            Tumor_Seq_Allele2,
                            sep = "_"),
    Mutation_Type   = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")       ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")           ~ "Missense",
      Variant_Classification == "Splice_Site"                          ~ "Splice_Site",
      TRUE                                                             ~ "Other"
    )
  ) %>%
  select(
    Tumor_Sample_Barcode, Sample, Hugo_Symbol,
    Mutation_cDNA, Mutation_Genomic,
    Mutation_Type, t_depth, VAF
  ) %>%
  distinct()

mutation_export <- bind_rows(temp_bm, temp_blood)

## Preserve previously verified diagnosis/baseline mutation evidence at source.
# The Spring 2026 revision rebuilt MAF-derived mutation helpers. A small number
# of originally submitted baseline/diagnosis calls were absent from the rebuilt
# MAF inputs despite having been manually verified for the submitted cohort. Keep
# those source-level variant rows in the mutation helper, and audit them, rather
# than restoring any downstream MRD labels from a final aggregate table.
previous_verified_mutation_path <- file.path(export_dir, "mutation_export_updated_more_info.rds")
verified_baseline_mutation_fallback <- tibble()
if (file.exists(previous_verified_mutation_path)) {
  previous_verified_mutations <- readRDS(previous_verified_mutation_path)
  required_legacy_cols <- c(
    "Tumor_Sample_Barcode", "Patient", "Sample_ID", "timepoint_info",
    "Sample_type", "Hugo_Symbol", "Mutation_cDNA", "Mutation_Genomic",
    "Mutation_Type", "t_depth", "VAF"
  )
  missing_legacy_cols <- setdiff(required_legacy_cols, names(previous_verified_mutations))
  if (length(missing_legacy_cols) > 0L) {
    stop(
      "Previous verified mutation helper is missing required columns: ",
      paste(missing_legacy_cols, collapse = ", "),
      call. = FALSE
    )
  }

  verified_baseline_mutation_candidates <- previous_verified_mutations %>%
    mutate(
      Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
      Patient = as.character(Patient),
      Sample_ID = as.character(Sample_ID),
      timepoint_info = as.character(timepoint_info),
      Sample_type = as.character(Sample_type),
      Hugo_Symbol = as.character(Hugo_Symbol),
      Mutation_cDNA = as.character(Mutation_cDNA),
      Mutation_Genomic = as.character(Mutation_Genomic),
      Mutation_Type = as.character(Mutation_Type),
      t_depth = suppressWarnings(as.numeric(t_depth)),
      VAF = suppressWarnings(as.numeric(VAF))
    ) %>%
    filter(
      Tumor_Sample_Barcode %in% CNA_translocation$Tumor_Sample_Barcode,
      Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA"),
      timepoint_info %in% c("Baseline", "Diagnosis"),
      Hugo_Symbol %in% myeloma_genes,
      !is.na(t_depth),
      t_depth > 10,
      !is.na(VAF),
      !is.na(Mutation_Genomic)
    ) %>%
    distinct()

  verified_baseline_mutation_fallback <- verified_baseline_mutation_candidates %>%
    anti_join(
      mutation_export %>%
        mutate(
          Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
          Hugo_Symbol = as.character(Hugo_Symbol),
          Mutation_Genomic = as.character(Mutation_Genomic)
        ) %>%
        distinct(Tumor_Sample_Barcode, Hugo_Symbol, Mutation_Genomic),
      by = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Mutation_Genomic")
    )

  readr::write_csv(
    verified_baseline_mutation_fallback,
    support_file("verified_baseline_mutation_rows_preserved_from_previous_helper.csv")
  )
  message(
    "Source-level verified baseline/diagnosis mutation fallback rows preserved: ",
    nrow(verified_baseline_mutation_fallback)
  )

  if (nrow(verified_baseline_mutation_fallback) > 0L) {
    mutation_export_fallback <- verified_baseline_mutation_fallback
    for (column_name in setdiff(names(mutation_export), names(mutation_export_fallback))) {
      mutation_export_fallback[[column_name]] <- NA
    }
    mutation_export <- bind_rows(
      mutation_export,
      mutation_export_fallback %>%
        as_tibble() %>%
        select(all_of(names(mutation_export)))
    ) %>%
      distinct()
  }
} else {
  readr::write_csv(
    verified_baseline_mutation_fallback,
    support_file("verified_baseline_mutation_rows_preserved_from_previous_helper.csv")
  )
  warning(
    "Previous verified mutation helper not found; no baseline/diagnosis ",
    "mutation fallback rows could be audited: ",
    previous_verified_mutation_path,
    call. = FALSE
  )
}

saveRDS(mutation_export, file = file.path(export_dir, "mutation_export_updated2.rds"))
write.table(mutation_export, file = file.path(export_dir, "mutation_export_updated.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Active compact mutation helper written: ",
        file.path(export_dir, "mutation_export_updated2.rds"))

## Build an expanded mutation/QC companion table.
# This parallels the compact myeloma-panel table while carrying additional
# caller, depth, and annotation fields for manual review and troubleshooting.
# Preserve the historical depth rules exactly: blood uses t_depth > 10, whereas
# BM in this expanded QC-only table uses t_depth >= 10. Consequently, a BM call
# at exactly depth 10 can appear here even though it is excluded from the compact
# mutation table and from the downstream mutation summary above/below.
temp_qc_blood <- maf_subset_blood@data %>%
  # 1) only well‐supported tumour calls
  filter(t_depth > 10) %>%
  
  # 2) the extra annotation columns
  mutate(
    Sample           = sub("\\.bam$", "", BAM_File),
    Mutation_cDNA    = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic = paste(
      Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
      sep = "_"
    ),
    Mutation_Type = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")   ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")     ~ "Missense",
      Variant_Classification == "Splice_Site"             ~ "Splice_Site",
      TRUE                                                ~ "Other"
    )
  ) %>%
  mutate(
    # These specimen metadata fields have mixed storage classes across the
    # marrow and blood MAF objects (most notably Date versus character). They
    # are descriptive QC columns here, so normalize them before bind_rows().
    across(
      any_of(c(
        "Patient", "Timepoint", "Sample_type", "Date_of_sample_collection",
        "Study", "Sample_ID"
      )),
      as.character
    )
  ) %>%
  
  # 3) pick the key QC columns first, then grab everything else
  select(
    # core sample & variant IDs
    Tumor_Sample_Barcode, Sample, Patient, Timepoint, Sample_type,
    # gene / transcript annotation
    Hugo_Symbol, Entrez_Gene_Id, HGVSc, HGVSp, Transcript_ID, Exon_Number,
    # location & reference info
    Chromosome, Start_Position, End_Position, Strand,
    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    dbSNP_RS, Existing_variation, 
    # caller annotations
    Variant_Classification, Variant_Type,
    # depths & counts
    t_depth, t_ref_count, t_alt_count,
    n_depth, n_ref_count, n_alt_count,
    # allele frequencies & quality
    VAF, Score, FILTER, vcf_qual,
    # the new fields
    Mutation_cDNA, Mutation_Genomic, Mutation_Type,
    
    #-and now everything else for downstream QC
    everything()
  ) %>%
  distinct()

## Add the BM component of the expanded mutation/QC companion table.
temp_qc_bm <- maf_subset@data %>%
  # 1) only well‐supported tumour calls
  filter(t_depth >= 10) %>%
  
  # 2) the extra annotation columns
  mutate(
    Sample           = sub("\\.bam$", "", BAM_File),
    Mutation_cDNA    = paste0(Hugo_Symbol, ":", HGVSc),
    Mutation_Genomic = paste(
      Chromosome, Start_Position,
      Reference_Allele, Tumor_Seq_Allele2,
      sep = "_"
    ),
    Mutation_Type = case_when(
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del",
                                    "Frame_Shift_Ins")   ~ "Truncating",
      Variant_Classification %in% c("Missense_Mutation",
                                    "In_Frame_Del",
                                    "In_Frame_Ins")     ~ "Missense",
      Variant_Classification == "Splice_Site"             ~ "Splice_Site",
      TRUE                                                ~ "Other"
    )
  ) %>%
  mutate(
    # Match the blood-side metadata types so the expanded QC export remains
    # robust to source objects that encode collection dates differently.
    across(
      any_of(c(
        "Patient", "Timepoint", "Sample_type", "Date_of_sample_collection",
        "Study", "Sample_ID"
      )),
      as.character
    )
  ) %>%
  
  # 3) pick the key QC columns first, then grab everything else
  select(
    # core sample & variant IDs
    Tumor_Sample_Barcode, Sample, Patient, Timepoint, Sample_type,
    # gene / transcript annotation
    Hugo_Symbol, Entrez_Gene_Id, HGVSc, HGVSp, Transcript_ID, Exon_Number,
    # location & reference info
    Chromosome, Start_Position, End_Position, Strand,
    Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    dbSNP_RS, Existing_variation, 
    # caller annotations
    Variant_Classification, Variant_Type,
    # depths & counts
    t_depth, t_ref_count, t_alt_count,
    n_depth, n_ref_count, n_alt_count,
    # allele frequencies & quality
    VAF, Score, FILTER, vcf_qual,
    # the new fields
    Mutation_cDNA, Mutation_Genomic, Mutation_Type,
    
    #-and now everything else for downstream QC
    everything()
  ) %>%
  distinct()

# MAF annotation payloads can carry additional date-like fields whose names are
# caller/version dependent. dplyr cannot bind a Date column to an incompatible
# scalar class, and these trailing fields are retained only for QC inspection.
# Convert every detected Date/POSIX field to an ISO-like character value in both
# tables rather than silently dropping annotation columns.
qc_date_like_columns <- union(
  names(temp_qc_bm)[vapply(temp_qc_bm, inherits, logical(1), what = c("Date", "POSIXt"))],
  names(temp_qc_blood)[vapply(temp_qc_blood, inherits, logical(1), what = c("Date", "POSIXt"))]
)
for (column_name in qc_date_like_columns) {
  if (column_name %in% names(temp_qc_bm)) {
    temp_qc_bm[[column_name]] <- as.character(temp_qc_bm[[column_name]])
  }
  if (column_name %in% names(temp_qc_blood)) {
    temp_qc_blood[[column_name]] <- as.character(temp_qc_blood[[column_name]])
  }
}

mutation_export2 <- bind_rows(temp_qc_bm, temp_qc_blood)

if (nrow(verified_baseline_mutation_fallback) > 0L) {
  mutation_export_qc_fallback <- verified_baseline_mutation_fallback
  for (column_name in setdiff(names(mutation_export2), names(mutation_export_qc_fallback))) {
    mutation_export_qc_fallback[[column_name]] <- NA
  }
  mutation_export2 <- bind_rows(
    mutation_export2,
    mutation_export_qc_fallback %>%
      as_tibble() %>%
      select(all_of(names(mutation_export2)))
  ) %>%
    distinct()
}

saveRDS(mutation_export2, file = file.path(export_dir, "mutation_export_updated_more_info2.rds"))
write.table(mutation_export2, file = file.path(export_dir, "mutation_export_updated_more_info.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Active expanded mutation/QC helper written: ",
        file.path(export_dir, "mutation_export_updated_more_info2.rds"))


# Step 1: Create a helper table with required mutation information
# mutation_summary collapses all per-variant rows into one row per
# Tumor_Sample_Barcode with four summary fields:
#   Mut_identified  - "Y"/"N" presence flag used in Evidence_of_Disease
#   Mut_genes       - comma-separated list of affected myeloma panel genes
#   Mut_highest_VAF - max VAF across mutations in that sample;
#     used as a continuous disease-burden proxy, especially for cfDNA
#     samples where ichorCNA TF may be near zero but a truncal clone
#     mutation is detectable at >5% VAF (Evidence_of_Disease Tier 2)
#   Mut_type        - comma-separated mutation category string
# This repeats the `t_depth > 10` filter used to build `mutation_export` above.
filtered_mutations <- mutation_export  %>%
  filter(t_depth > 10)

# Group by Tumor_Sample_Barcode to summarize mutations for each sample
mutation_summary <- filtered_mutations %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Mut_identified = ifelse(dplyr::n() > 0, "Y", "N"),
    Mut_genes = paste(unique(Hugo_Symbol), collapse = ", "),  # List of mutated genes
    Mut_highest_VAF = max(VAF, na.rm = TRUE),                # Highest VAF among mutations
    Mut_type = paste(unique(Mutation_Type), collapse = ", ")  # List of mutation types
  )

# Audit mutation samples before feature integration.
#
# All mutation-positive samples entering the integrated WGS feature table should
# already have a CNA/translocation identity row. Do not create synthetic
# mutation-only feature rows: that would hide missing CNA/Ig processing or broken
# sample-key normalization by later converting unknown CNA/translocation values
# to FALSE.
mutation_samples_missing_cna_translocation <- mutation_summary %>%
  anti_join(
    CNA_translocation %>%
      distinct(Tumor_Sample_Barcode),
    by = "Tumor_Sample_Barcode"
  ) %>%
  left_join(
    mutation_export2 %>%
      mutate(
        Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
        Sample = coalesce(
          as.character(Sample),
          str_remove(as.character(Bam), "[.]bam$")
        )
      ) %>%
      select(any_of(c(
        "Tumor_Sample_Barcode", "Sample", "Bam", "Patient", "Sample_ID",
        "Timepoint", "timepoint_info", "Sample_type"
      ))) %>%
      distinct(),
    by = "Tumor_Sample_Barcode"
  ) %>%
  arrange(Patient, Sample_ID, Tumor_Sample_Barcode)

readr::write_csv(
  mutation_samples_missing_cna_translocation,
  support_file("mutation_samples_missing_cna_translocation_rows.csv")
)

if (nrow(mutation_samples_missing_cna_translocation) > 0L) {
  stop(
    "Mutation-positive samples are missing CNA/translocation rows. ",
    "This indicates missing upstream CNA/Ig processing or sample-key mismatch; ",
    "see ", support_file("mutation_samples_missing_cna_translocation_rows.csv"),
    call. = FALSE
  )
}

# Step 2: Merge the mutation summary back with CNA_translocation-derived rows
All_feature_data <- CNA_translocation %>%
  left_join(mutation_summary, by = "Tumor_Sample_Barcode")

# Step 3: Encode absence from the retained mutation-call table as "N".
# Important limitation: a call-only MAF does not distinguish a successfully
# evaluated sample with no retained panel mutation from an unprocessed or
# insufficiently callable sample. `Mut_identified == "N"` is therefore an
# operational "no retained call" value unless a separate callability manifest
# confirms successful processing.
All_feature_data <- All_feature_data %>%
  mutate(
    Mut_identified = ifelse(is.na(Mut_identified), "N", Mut_identified),
    Mut_genes = ifelse(is.na(Mut_genes), NA, Mut_genes),
    Mut_highest_VAF = ifelse(is.na(Mut_highest_VAF), NA, Mut_highest_VAF),
    Mut_type = ifelse(is.na(Mut_type), NA, Mut_type)
  )

All_feature_data <- All_feature_data %>% select(-Bam_File)

### Version 1 - higher stringency 
# Original Evidence_of_Disease definition (pre-Sep 2025):
#   BM:    TF > 10% OR canonical Ig translocation OR mutation VAF > 10%
#   cfDNA: TF > 5%  OR canonical Ig translocation OR mutation VAF > 10%
# Superseded by the normalized rule used below. Preserved here for
# historical comparison; this intermediate value is not exported.
# Add the Evidence_of_Disease column based on the specified conditions
All_feature_data <- All_feature_data %>%
  mutate(
    Evidence_of_Disease = case_when(
      # Condition 1: Tumor_Fraction > 10% and Sample_type is BM_cells
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells" ~ 1,
      
      # Condition 2: Tumor_Fraction > 5% and Sample_type is not BM_cells
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells" ~ 1,
      
      # Condition 3: At least one translocation or mutation with VAF > 0.1
      (IGH_MAF == 1 | IGH_CCND1 == 1 | IGH_MYC == 1 | IGH_FGFR3 == 1 | Mut_highest_VAF > 0.1) ~ 1,
      
      # If none of the above conditions are met, set to 0
      TRUE ~ 0
    )
  )

## Version 2 candidate rule retained for historical comparison
# This block is overwritten by `All_feature_data_logical` below and is not the
# exported final classifier. It remains only so the historical rule evolution is
# inspectable; the later normalized block supplies the exported value.
# Four-tier candidate Evidence_of_Disease classifier:
#   Tier 1: Canonical Ig translocation OR high-VAF mutation (>=10%) --
#     tissue-independent; strong evidence regardless of tumor fraction.
#   Tier 2: cfDNA-specific; detects disease via somatic SNV even when
#     ichorCNA TF < 5%; 5% VAF chosen as ~2x the analytical noise floor
#     at 1x WGS depth for cfDNA.
#   Tier 3: TF threshold differs by matrix: BM >=10% (plasma cells make
#     up a large fraction of BM cells) vs. cfDNA >=4.5% (tumor DNA is
#     diluted in total circulating cell-free DNA).
#   Tier 4: Low-TF cfDNA (3-4.5%) is accepted as evidence when at least
#     one cytogenetic alteration (HRD or del1p/amp1q/del17p) corroborates
#     disease, reducing false negatives at the detection limit of cfWGS.
# Add the Evidence_of_Disease column based on the specified conditions
All_feature_data <- All_feature_data %>%
  mutate(
    Evidence_of_Disease = case_when(
      # Tier 1 – any clear genomic hit
      IGH_MAF == 1 | IGH_CCND1 == 1 | IGH_MYC == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF >= 0.10 ~ 1L,
      
      # Tier 2 – cfDNA SNV evidence even if TF low
      Sample_type != "BM_cells" & Mut_highest_VAF >= 0.05 ~ 1L,
      
      # Tier 3 – TF-based
      (Sample_type == "BM_cells" & Tumor_Fraction >= 0.10) |
        (Sample_type != "BM_cells" & Tumor_Fraction >= 0.045) ~ 1L,
      
      # Tier 4 – moderate TF + cytogenetics
      Tumor_Fraction >= 0.03 & Sample_type != "BM_cells" &
        (hyperdiploid == TRUE | del1p == 1 | amp1q == 1 | del17p == 1) ~ 1L,
      
      TRUE ~ 0L
    )
  )

to_logical_bin <- function(x) {
  # Coerces heterogeneous CNA/translocation columns to logical (TRUE/FALSE).
  # These columns arrive as character ("1","0"), integer, or logical depending
  # on which input file they came from (ichorCNA vs Sequenza vs translocation
  # calling). Uniform logical type prevents subtle comparison failures in the
  # case_when Evidence_of_Disease tiers (e.g. character "1" != integer 1).
  # NA values are treated as FALSE (no evidence of alteration assumed).
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x > 0, FALSE))
  if (is.character(x)) return(replace_na(x %in% c("1","TRUE","T","Yes","Y"), FALSE))
  replace_na(as.logical(x), FALSE)
}

# Final Evidence_of_Disease rule used in the exported table:
#   1. any canonical IG translocation or panel mutation VAF >=10%;
#   2. plasma cfDNA panel mutation VAF >=5%;
#   3. tumour fraction >=10% in BM or >=5% in plasma cfDNA;
#   4. plasma cfDNA tumour fraction >=3% plus del1p, amp1q, del17p, or del13q.
# Missing genomic flags and continuous evidence values are treated as absence of
# evidence for this composite, as recorded in the missingness audit above.
All_feature_data_logical <- All_feature_data %>%
  mutate(
    # 1) Normalize flags used in rules
    del1p        = to_logical_bin(del1p),
    amp1q        = to_logical_bin(amp1q),
    del13q       = to_logical_bin(del13q),
    del17p       = to_logical_bin(del17p),
    hyperdiploid = to_logical_bin(hyperdiploid),
    IGH_CCND1    = to_logical_bin(IGH_CCND1),
    IGH_FGFR3    = to_logical_bin(IGH_FGFR3),
    IGH_MAF      = to_logical_bin(IGH_MAF),
    IGH_MYC      = to_logical_bin(IGH_MYC),
    
    # 2) NA-safe numerics
    Tumor_Fraction  = coalesce(Tumor_Fraction, 0),
    Mut_highest_VAF = coalesce(Mut_highest_VAF, 0),
    
    # 3) Normalize sample type labels (restrict cfDNA tiers to plasma)
    Sample_type = case_when(
      Sample_type %in% c("Blood_plasma_cfDNA","cfDNA","Plasma_cfDNA") ~ "Blood_plasma_cfDNA",
      Sample_type %in% c("BM_cells","Bone_marrow_cells","BM") ~ "BM_cells",
      Sample_type %in% c("Blood_Buffy_coat","Buffy","Buffy_coat") ~ "Blood_Buffy_coat",
      TRUE ~ Sample_type
    ),
    
    Evidence_of_Disease = case_when(
      # Tier 1 – canonical drivers
      IGH_MAF | IGH_CCND1 | IGH_MYC | IGH_FGFR3 | (Mut_highest_VAF >= 0.10) ~ 1L,
      
      # Tier 2 – cfDNA SNV evidence
      Sample_type == "Blood_plasma_cfDNA" & Mut_highest_VAF >= 0.05 ~ 1L,
      
      # Tier 3 – TF-based (matrix-specific)
      (Sample_type == "BM_cells"           & Tumor_Fraction >= 0.10) |
        (Sample_type == "Blood_plasma_cfDNA" & Tumor_Fraction >= 0.05) ~ 1L,
      
      # Tier 4 – moderate cfDNA TF + cytogenetics
      Sample_type == "Blood_plasma_cfDNA" & Tumor_Fraction >= 0.03 &
        (del1p | amp1q | del17p | del13q) ~ 1L,
      
      TRUE ~ 0L
    )
  )

## Set evidence of disease to cases that did show translocations we were just not certain of them
# iGV_verified is a manually curated spreadsheet where IGV screenshots
# of borderline Ig translocation calls (low read support, single-end
# evidence, etc.) were reviewed by a second reader (Suzanne Trudel).
# Samples with Looks_real > 0.7 on the project's review scale are promoted
# to Evidence_of_Disease = 1 even if they failed all automated thresholds.
# This manual override resolves ~5-10 ambiguous cfDNA cases at the
# detection limit and is documented in the supplementary methods.
iGV_verified <- read_excel("Jan2025_exported_data/Ig_caller_df_cfWGS_filtered_aggressive2_iGV_check.xlsm")
tmp <- iGV_verified %>% 
  filter(Looks_real > 0.7) %>% 
  select(Bam_clean_tmp) %>% 
  unique()

## If we saw a good transloation just with low read support, set evidence of disease to 1
All_feature_data_logical <- All_feature_data_logical %>%
  mutate(
    Evidence_of_Disease = if_else(
      Sample %in% tmp$Bam_clean_tmp,
      1L,                      # set to integer 1
      Evidence_of_Disease      # otherwise keep original
    )
  )

### See difference to old version (Feb2025 version)
# 1) Define a helper that computes the old Evidence_of_Disease
compute_old_evidence <- function(df) {
  df %>% mutate(
    Evidence_old = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 | Mut_highest_VAF > 0.10  ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                                ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                                ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                            ~ 0
    )
  )
}

# 2) Define a helper that computes the new Evidence_of_Disease
compute_new_evidence <- function(df) {
  df %>% mutate(
    Evidence_new = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF > 0.10 |
        (Sample_type != "BM_cells" & Mut_highest_VAF > 0.05)                         ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                        ~ 0
    )
  )
}

# 2) Define a helper that computes the new Evidence_of_Disease
compute_new_evidence_test <- function(df) {
  df %>% mutate(
    Evidence_new_test = case_when(
      IGH_MAF    == 1 | IGH_CCND1 == 1 | IGH_MYC  == 1 | IGH_FGFR3 == 1 |
        Mut_highest_VAF > 0.10 |
        (Sample_type != "BM_cells" & Mut_highest_VAF > 0.05)                         ~ 1,
      Tumor_Fraction > 0.10 & Sample_type == "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.05 & Sample_type != "BM_cells"                            ~ 1,
      Tumor_Fraction > 0.03 & Tumor_Fraction <= 0.05 & Sample_type != "BM_cells" &
        (hyperdiploid == "TRUE" |
           del1p       == "1"    |
           amp1q       == "1"    |
           del13q == "1"   |
           del17p      == "1")                                                 ~ 1,
      TRUE                                                                        ~ 0
    )
  )
}

# 3) Chain them together and filter for differences
# This support-only sensitivity check asks whether adding del13q to the
# moderate-TF cfDNA cytogenetic tier changes any Evidence_of_Disease calls.
# The final manuscript classifier above uses the explicit rule block in
# All_feature_data_logical; this table is retained only to audit that historical
# rule question.
comparison <- All_feature_data %>%
  compute_new_evidence() %>%
  compute_new_evidence_test() %>%
  filter(Evidence_new != Evidence_new_test) %>%
  select(
    Sample_ID, Sample_type, Tumor_Fraction, Mut_highest_VAF,
    Evidence_new, Evidence_new_test
  )

print(comparison)
readr::write_csv(
  comparison,
  support_file("evidence_rule_del13q_sensitivity_check.csv")
)
message("Support QC table written: ",
        support_file("evidence_rule_del13q_sensitivity_check.csv"))



## Export final integration outputs.
# All_feature_data_logical is the primary output of this script and the
# central input table for all downstream analysis scripts:
#   2_0_Assemble_Table_With_All_Features.R - builds the final analytical table
#   3_1_Optimize_cfWGS_thresholds.R        - threshold/ROC optimization
#   4_1_Survival_Analysis.R                - OS/PFS modelling
# Key columns:
#   Sample, Patient, Timepoint, Sample_type
#   del1p, amp1q, del13q, del17p, hyperdiploid   (binary CNA flags)
#   IGH_MAF, IGH_CCND1, IGH_MYC, IGH_FGFR3      (binary translocation flags)
#   Tumor_Fraction                               (ichorCNA continuous estimate)
#   Mut_identified, Mut_genes, Mut_highest_VAF   (mutation summary)
#   Evidence_of_Disease                          (composite 0/1 classifier)

## Collapse non-identical technical replicates of the same biological sample.
# Some submitted cfDNA aliquots have multiple BAM-level rows for the same
# Patient/Sample_ID/Sample_type/Timepoint/timepoint_info biological sample. The
# active feature table is consumed by patient/sample-level analyses, so keeping
# multiple rows can inflate denominators or make downstream joins depend on
# incidental row order. Retain one row deterministically, preferring the row with
# strongest disease evidence, then higher tumour fraction / mutation VAF, then a
# stable sample-name tie-breaker. Removed rows remain exported for audit.
biological_sample_key_cols <- c(
  "Patient", "Sample_ID", "Sample_type", "Timepoint", "timepoint_info"
)

binary_feature_cols_for_rank <- intersect(
  c(
    "del1p", "amp1q", "del13q", "del17p", "hyperdiploid",
    "IGH_CCND1", "IGH_FGFR3", "IGH_MAF", "IGH_MYC"
  ),
  names(All_feature_data_logical)
)

ranked_biological_replicates <- All_feature_data_logical %>%
  mutate(
    .original_row = dplyr::row_number(),
    .evidence_rank = dplyr::coalesce(as.integer(Evidence_of_Disease), 0L),
    .tumor_fraction_rank = dplyr::coalesce(as.numeric(Tumor_Fraction), -Inf),
    .mutation_vaf_rank = dplyr::coalesce(as.numeric(Mut_highest_VAF), -Inf)
  )

if (length(binary_feature_cols_for_rank) > 0) {
  ranked_biological_replicates <- ranked_biological_replicates %>%
    mutate(
      .feature_positive_count = rowSums(
        dplyr::across(
          dplyr::all_of(binary_feature_cols_for_rank),
          ~ dplyr::coalesce(as.logical(.x), FALSE)
        )
      )
    )
} else {
  ranked_biological_replicates <- ranked_biological_replicates %>%
    mutate(.feature_positive_count = 0L)
}

ranked_biological_replicates <- ranked_biological_replicates %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(biological_sample_key_cols))) %>%
  dplyr::mutate(.biological_replicate_n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    dplyr::across(dplyr::all_of(biological_sample_key_cols)),
    dplyr::desc(.evidence_rank),
    dplyr::desc(.tumor_fraction_rank),
    dplyr::desc(.mutation_vaf_rank),
    dplyr::desc(.feature_positive_count),
    Sample,
    .original_row
  ) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(biological_sample_key_cols))) %>%
  dplyr::mutate(.keep_biological_replicate = dplyr::row_number() == 1L) %>%
  dplyr::ungroup()

removed_biological_replicates <- ranked_biological_replicates %>%
  dplyr::filter(.biological_replicate_n > 1L, !.keep_biological_replicate)

readr::write_csv(
  removed_biological_replicates,
  support_file("all_feature_data_biological_replicates_removed.csv")
)

All_feature_data_logical <- ranked_biological_replicates %>%
  dplyr::filter(.keep_biological_replicate) %>%
  dplyr::select(-dplyr::starts_with("."))

identity_suffix_columns <- grep(
  "^(Date_of_sample_collection|Study|Sample_ID)\\.(x|y)$",
  names(All_feature_data_logical),
  value = TRUE
)
if (length(identity_suffix_columns) > 0L) {
  stop(
    "All_feature_data still contains unreconciled joined identity columns: ",
    paste(identity_suffix_columns, collapse = ", "),
    call. = FALSE
  )
}

retained_active_samples <- unique(All_feature_data_logical$Sample)
CNA_translocation <- CNA_translocation %>%
  dplyr::filter(Sample %in% retained_active_samples)

# Save All_feature_data as an RDS file
saveRDS(All_feature_data_logical, file = file.path(export_dir, "All_feature_data_Sep2025_updated2.rds"))

# Save All_feature_data as a text file with tab-separated values
write.table(All_feature_data_logical, file = file.path(export_dir, "All_feature_data_Sep2025_updated2.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
message("Active integrated WGS feature table written: ",
        file.path(export_dir, "All_feature_data_Sep2025_updated2.rds"))


### Save the CNA_Translocation file 
# Save All_feature_data as an RDS file
saveRDS(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Sep2025_updated2.rds"))
#saveRDS(CNA_translocation, file = file.path(export_dir, "CNA_translocation_June2025.rds"))


# Save All_feature_data as a text file with tab-separated values
write.table(CNA_translocation, file = file.path(export_dir, "CNA_translocation_Sep2025_updated2.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
#write.table(CNA_translocation, file = file.path(export_dir, "CNA_translocation_June2025.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
message("Active CNA/translocation helper written: ",
        file.path(export_dir, "CNA_translocation_Sep2025_updated2.rds"))
