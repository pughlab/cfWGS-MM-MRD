# =============================================================================
# 1_8_Process_Cumulative_VAFs_MRDetect.R
# Project:  cfWGS MRDetect (Winter 2025)
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_8_Process_Cumulative_VAFs_MRDetect.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes MRDetect outputs for
#     patient samples that feed model training, clinical concordance, and
#     survival analyses.
#
# Downstream use:
#   - `2_0_Assemble_Table_With_All_Features.R` reads the BM-informed and
#     blood-informed longitudinal CSVs and adds their MRDetect features to the
#     integrated analysis table.
#   - `1_8B_Export_MRDetect_patient_feature_table.R` and
#     `2_4_Longitudinal_features_analysis.R` read the complete z-scored RDS.
#   - `2_1_Part2_Cohort_Swim_Plot.R` reads the all-patient blood-informed CSV.
#   - `1_10_Estimate_MRDetect_LOD_proxy.R` and
#     `3_1_part2_Apply_cfWGS_thresholds_to_dilution_series.R` read the retained
#     historical CHARM control export.
#   The integrated outputs subsequently feed the model, concordance, survival,
#   and longitudinal manuscript figures and tables.
#
# Author:   Dory Abelman
# Date:     January 2025
# Last Updated: September 2026
#
# Purpose:
#   1. Combine historical and Spring 2026 MRDetect result files.
#   2. Correct the preserved Winter 2025 parser column order and apply audited
#      rerun replacements.
#   3. Attach mutation-panel and queried-sample metadata, then calculate
#      mutation-panel-specific healthy-control z-scores.
#   4. Export row-level MRDetect results and longitudinal BM- and blood-informed
#      summaries used by downstream integration and figure scripts.
#
# Unit of analysis:
#   One MRDetect result row is one queried BAM evaluated against one personalized
#   mutation-list VCF under one mutation-source/filter-source combination.
#
# Dependencies:
#   • readr, data.table, tidyverse (dplyr, tidyr, stringr), openxlsx
#   • ggplot2, ggbreak, patchwork, scales, conflicted
#
# Principal input files:
#   • MRDetect_output_winter_2025/MRDetect_outputs/*.csv
#   • Spring 2026 revision MRDetect exports returned by
#     spring2026_revision_files()
#   • combined_clinical_data_updated_April2025.csv plus the revision metadata
#     incorporated by read_combined_clinical_metadata_with_revision()
#   • Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#   • cohort_assignment_table_updated.rds
#   • Scripts_2025/Final_Scripts/input_overrides/
#     VA15_mrdetect_rerun_2026-07-10.csv, when present
#   • helpers.R, which resolves Spring 2026 inputs, sample identities, cohort
#     metadata, strict MRDetect parsing, and control deduplication
#
# Output Directory (created if necessary):
#   • MRDetect_output_winter_2025/Processed_R_outputs/
#   • Writes:
#       - cfWGS_Winter2025All_MRDetect_Sep2025.txt
#       - cfWGS_Winter2025All_MRDetect_Sep2025.rds
#       - cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.txt
#       - cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds
#       - All_detection_rates_baseline_and_controls_Feb2026.csv/.rds
#         (legacy filename; contains the historical CHARM control reference rows)
#       - BM_muts_plots_baseline/cfWGS_MRDetect_BM_data_updated_Sep.csv
#       - Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated Sep.csv
#       - Blood_muts_plots_baseline/
#         cfWGS MRDetect Blood data updated Sep with all patients.csv
#       - clinical_support/baseline_diagnosis_mutation_source_by_patient.csv
#       - audit CSVs describing excluded mutation panels, control coverage,
#         rerun overrides, source selection, aliases, and duplicate sources
#
# Important interpretation details:
#   • Personalized mutation-list VCFs are restricted to Baseline or Diagnosis
#     sources before the final BM- and blood-informed exports are written.
#   • Historical queried samples are normalized against historical CHARM
#     controls; Spring 2026 XPlus samples request the XPlus CHARM reference.
#   • `Mrd_by_WGS` and `Cumulative_VAF` use the retained exploratory 4.5-z-score
#     screen. Final cfWGS binary calls are defined downstream by script 3_1.
#   • The “All_detection_rates_baseline_and_controls” filename is historical;
#     the current export contains historical CHARM control rows only.
#
# Usage:
#   Run from the project root with the command shown above.
# =============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#

# ──────────────────────────────────────────────────────────────────────────────
# 1) Load libraries
# ──────────────────────────────────────────────────────────────────────────────
library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggbreak)
library(patchwork)
library(scales)
#library(conflicted)
library(lubridate)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

restrict_to_single_baseline_mutation_source <- function(dat, audit_dir, audit_prefix) {
  # ## Keep one baseline/diagnosis personalized mutation source per patient
  # Spring 2026 metadata can contain both Diagnosis and Baseline rows for the same
  # patient and mutation-source compartment. MRDetect rows are patient-sample x
  # mutation-list comparisons; if two baseline-like mutation lists are carried
  # forward for the same patient, downstream longitudinal summaries double-count
  # that patient. This helper deterministically chooses one source and writes the
  # decision table before filtering.
  baseline_labels <- c("Baseline", "Diagnosis")

  baseline_sources <- dat %>%
    filter(timepoint_info %in% baseline_labels) %>%
    mutate(
      # Prefer the canonical biological source identity when metadata exists.
      # Historical and revision MRDetect files can name the same baseline sample
      # with different VCF/BAM strings; those aliases should not be treated as
      # separate baseline mutation sources. Fall back to VCF strings only when
      # source metadata is unavailable.
      source_date = suppressWarnings(as.Date(Date_of_sample_collection_Sample_ID)),
      biological_source_key = if_else(
        !is.na(Patient) & !is.na(Mut_source) & !is.na(Sample_ID) & !is.na(source_date),
        paste(Patient, Mut_source, Sample_ID, source_date, sep = "|"),
        NA_character_
      ),
      source_key = coalesce(
        biological_source_key,
        na_if(as.character(VCF_factor), ""),
        na_if(as.character(VCF_clean), ""),
        na_if(as.character(VCF), ""),
        na_if(as.character(Sample_ID), ""),
        paste(as.character(Patient), as.character(Mut_source), as.character(Timepoint), as.character(timepoint_info), sep = "|")
      ),
      timepoint_rank = case_when(
        # Earlier baseline-like timepoints are preferred when dates cannot break
        # the tie. T0/TP0/D0/01 all represent the earliest baseline family.
        str_to_upper(Timepoint) %in% c("T0", "TP0", "D0", "0", "01", "1") ~ 1,
        str_to_upper(Timepoint) %in% c("T1", "TP1") ~ 2,
        str_detect(Timepoint, "^[Tt]?[0-9]+$") ~ as.numeric(str_remove(str_to_upper(Timepoint), "^T")) + 10,
        TRUE ~ 999
      ),
      label_rank = case_when(
        # Diagnosis is preferred over Baseline only after date/timepoint ranking.
        # This preserves explicit diagnosis mutation lists when both labels exist
        # for the same dated/timepoint context.
        timepoint_info == "Diagnosis" ~ 1,
        timepoint_info == "Baseline" ~ 2,
        TRUE ~ 9
      )
    ) %>%
    distinct(
      Patient, Mut_source, source_key, Sample_ID, Timepoint, timepoint_info,
      source_date, timepoint_rank, label_rank
    )

  duplicate_groups <- baseline_sources %>%
    count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
    filter(n_baseline_diagnosis_mutation_sources > 1)

  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(duplicate_groups) == 0L) {
    write_csv(
      baseline_sources %>%
        mutate(
          selected_source_key = NA_character_,
          selected_source = NA_character_,
          retained = NA,
          mutation_source = paste(Sample_ID, Timepoint, timepoint_info, source_date, sep = "|")
        ) %>%
        filter(FALSE),
      file.path(audit_dir, paste0(audit_prefix, "_baseline_mutation_source_deduplication_audit.csv"))
    )
    return(dat)
  }

  selected_sources <- baseline_sources %>%
    semi_join(duplicate_groups, by = c("Patient", "Mut_source")) %>%
    arrange(
      # Deterministic tie break: earliest collection date, earliest timepoint,
      # preferred label, then lexical source key. This avoids run-to-run changes
      # if input row order changes.
      Patient, Mut_source,
      coalesce(source_date, as.Date("9999-12-31")),
      timepoint_rank,
      label_rank,
      source_key
    ) %>%
    group_by(Patient, Mut_source) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      Patient, Mut_source, selected_source_key = source_key,
      selected_source = paste(Sample_ID, Timepoint, timepoint_info, source_date, sep = "|")
    )

  baseline_decisions <- baseline_sources %>%
    semi_join(duplicate_groups, by = c("Patient", "Mut_source")) %>%
    left_join(selected_sources, by = c("Patient", "Mut_source")) %>%
    mutate(
      retained = source_key == selected_source_key,
      mutation_source = paste(Sample_ID, Timepoint, timepoint_info, source_date, sep = "|")
    ) %>%
    arrange(Patient, Mut_source, desc(retained), timepoint_rank, source_key)

  write_csv(
    baseline_decisions,
    file.path(audit_dir, paste0(audit_prefix, "_baseline_mutation_source_deduplication_audit.csv"))
  )

  dat %>%
    mutate(
      source_date = suppressWarnings(as.Date(Date_of_sample_collection_Sample_ID)),
      biological_source_key = if_else(
        !is.na(Patient) & !is.na(Mut_source) & !is.na(Sample_ID) & !is.na(source_date),
        paste(Patient, Mut_source, Sample_ID, source_date, sep = "|"),
        NA_character_
      ),
      source_key = coalesce(
        biological_source_key,
        na_if(as.character(VCF_factor), ""),
        na_if(as.character(VCF_clean), ""),
        na_if(as.character(VCF), ""),
        na_if(as.character(Sample_ID), ""),
        paste(as.character(Patient), as.character(Mut_source), as.character(Timepoint), as.character(timepoint_info), sep = "|")
      )
    ) %>%
    left_join(selected_sources, by = c("Patient", "Mut_source")) %>%
    # Patients without duplicate baseline-like mutation lists pass through.
    # Patients with duplicates keep only the selected source_key.
    filter(is.na(selected_source_key) | source_key == selected_source_key) %>%
    select(-source_date, -biological_source_key, -source_key, -selected_source_key, -selected_source)
}

deduplicate_baseline_source_alias_query_rows <- function(dat, audit_dir, audit_prefix) {
  # ## Collapse VCF-name aliases only when they represent the same source/query
  # Once historical and revision IDs are mapped to one biological source, the same
  # queried BAM can occasionally have both a canonical VCF row and a historical
  # alias VCF row. Keep the canonical VCF for that exact source/query pair, but do
  # not remove alias-only longitudinal rows such as old MRDetect calls that were
  # never rerun under the revision VCF name.
  required_cols <- c(
    "Patient", "Mut_source", "Sample_ID", "Date_of_sample_collection_Sample_ID",
    "Sample_ID_Bam", "Filter_source", "VCF_clean"
  )
  if (!all(required_cols %in% names(dat))) return(dat)

  dat_keyed <- dat %>%
    mutate(
      source_date = suppressWarnings(as.Date(Date_of_sample_collection_Sample_ID)),
      biological_source_key = if_else(
        !is.na(Patient) & !is.na(Mut_source) & !is.na(Sample_ID) & !is.na(source_date),
        paste(Patient, Mut_source, Sample_ID, source_date, sep = "|"),
        paste(Patient, Mut_source, Sample_ID, VCF_clean, sep = "|")
      ),
      source_query_alias_key = paste(
        biological_source_key,
        Sample_ID_Bam,
        Filter_source,
        sep = "|"
      ),
      source_alias_priority = case_when(
        str_detect(VCF_clean, "__.+__vs__") ~ 1L,
        TRUE ~ 2L
      ),
      source_alias_row_id = row_number()
    )

  duplicate_alias_groups <- dat_keyed %>%
    count(source_query_alias_key, name = "n_source_query_alias_rows") %>%
    filter(n_source_query_alias_rows > 1)

  selected_alias_rows <- dat_keyed %>%
    group_by(source_query_alias_key) %>%
    arrange(source_alias_priority, desc(total_reads), VCF_clean, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(source_query_alias_key, selected_source_alias_row_id = source_alias_row_id)

  alias_decisions <- dat_keyed %>%
    semi_join(duplicate_alias_groups, by = "source_query_alias_key") %>%
    left_join(selected_alias_rows, by = "source_query_alias_key") %>%
    mutate(retained = source_alias_row_id == selected_source_alias_row_id) %>%
    select(
      Patient, Mut_source, Sample_ID, Sample_ID_Bam, Filter_source,
      VCF_clean, source_date, biological_source_key,
      source_alias_priority, total_reads, retained
    ) %>%
    arrange(Patient, Mut_source, Sample_ID, Sample_ID_Bam, desc(retained), source_alias_priority)

  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    alias_decisions,
    file.path(audit_dir, paste0(audit_prefix, "_baseline_source_alias_query_deduplication_audit.csv"))
  )

  dat_keyed %>%
    left_join(selected_alias_rows, by = "source_query_alias_key") %>%
    filter(is.na(selected_source_alias_row_id) | source_alias_row_id == selected_source_alias_row_id) %>%
    select(
      -source_date, -biological_source_key, -source_query_alias_key,
      -source_alias_priority, -source_alias_row_id, -selected_source_alias_row_id
    )
}

filter_to_baseline_mutation_sources <- function(dat, audit_dir, audit_prefix) {
  baseline_labels <- c("Baseline", "Diagnosis")

  excluded_sources <- dat %>%
    filter(is.na(timepoint_info) | !timepoint_info %in% baseline_labels) %>%
    distinct(
      Patient, Mut_source, Sample_ID, VCF_clean, VCF_factor, Timepoint,
      timepoint_info
    ) %>%
    arrange(Patient, Mut_source, Timepoint, VCF_clean)

  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    excluded_sources,
    file.path(audit_dir, paste0(audit_prefix, "_excluded_nonbaseline_mutation_sources.csv"))
  )

  dat %>%
    filter(timepoint_info %in% baseline_labels)
}

# resolve common conflicts
#conflicted::conflicts_prefer("dplyr::mutate")
#conflicted::conflicts_prefer("dplyr::filter")
#conflicted::conflicts_prefer("dplyr::select")
#conflicted::conflicts_prefer("dplyr::summarize")


# ──────────────────────────────────────────────────────────────────────────────
# 2) Set up paths and create output directory
# ──────────────────────────────────────────────────────────────────────────────
input_root <- "MRDetect_output_winter_2025/MRDetect_outputs/"
outdir     <- "MRDetect_output_winter_2025/Processed_R_outputs/"
project    <- "cfWGS_Winter2025"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


# ──────────────────────────────────────────────────────────────────────────────
# 3) Read and combine all CSVs, label with ‘source_file’, ‘Mut_source’ & ‘Filter_source’
# ──────────────────────────────────────────────────────────────────────────────
csv_files <- list.files(input_root, pattern = "\\.csv$", full.names = TRUE)
spring2026_patient_mrdetect_files <- spring2026_revision_files(
  # The combined patient export contains the Spring 2026 query BAMs. M4CHIP
  # query rows are removed below because they belong to the dilution workflow.
  "MRDetect_outputs",
  "^MRDetect_all_RESULTS_combined_with_source[.]csv$"
)
spring2026_final_xplus_control_files <- spring2026_revision_files(
  "MRDetect_outputs",
  "^MRDetect_all_RESULTS_combined_with_source_final_healthy_control_Xplus[.]csv$"
)
spring2026_incremental_xplus_control_files <- spring2026_revision_files(
  "MRDetect_outputs",
  "^MRDetect_all_RESULTS_combined_with_source_healthy_control_rerun[0-9]+[.]csv$"
)
# Use the complete final export when it is present. Incremental reruns
# remain a supported fallback so the script is still reproducible before final
# aggregation is available.
spring2026_xplus_control_files <- if (length(spring2026_final_xplus_control_files)) {
  spring2026_final_xplus_control_files
} else {
  spring2026_incremental_xplus_control_files
}
spring2026_mrdetect_files <- unique(c(
  spring2026_patient_mrdetect_files,
  spring2026_xplus_control_files
))
csv_files <- unique(c(csv_files, spring2026_mrdetect_files))
if (!length(csv_files)) {
  stop("No MRDetect CSV files found in historical or Spring 2026 revision inputs.", call. = FALSE)
}

read_and_label <- function(file) {
  # Historical MRDetect CSVs and Spring 2026 combined exports do not always carry
  # the same file provenance columns. Standardize to filename/source_file and add
  # input_source_file so later code can infer mutation source and audit origin.
  df <- read_mrdetect_csv_strict(file)
  if (basename(file) %in% basename(spring2026_xplus_control_files)) {
    unexpected_bams <- unique(df$BAM[!is_xplus_charm_healthy_bam(df$BAM)])
    if (length(unexpected_bams)) {
      stop(
        "Incremental XPlus healthy-control rerun contains non-allowlisted BAMs: ",
        paste(utils::head(unexpected_bams, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    df <- df %>%
      filter(str_detect(.data$source_file, "_VS_IMG-"))
  }
  if (!"filename" %in% names(df)) {
    if ("source_file" %in% names(df)) {
      df$filename <- df$source_file
    } else {
      df$filename <- basename(file)
    }
  }
  if (!"source_file" %in% names(df)) {
    df$source_file <- basename(file)
  }
  df$input_source_file <- basename(file)
  return(df)
}
all_files <- bind_rows(lapply(csv_files, read_and_label))

if (length(spring2026_final_xplus_control_files)) {
  # Do not mix partial legacy control rows with the complete all-by-all
  # final export. Patient/query rows from the other files remain unchanged.
  all_files <- all_files %>%
    filter(
      !is_xplus_charm_healthy_bam(.data$BAM) |
        .data$input_source_file %in% basename(spring2026_final_xplus_control_files)
    )
}

all_files <- all_files %>%
  filter(
    # M4CHIP/PWGVAL rows are dilution-series query BAMs. They are excluded from
    # the main patient-sample MRDetect output so they do not enter clinical
    # patient scoring, survival analyses, or main-cohort model inputs.
    !str_detect(source_file, "^M4CHIP_"),
    !str_detect(filename, "^M4CHIP_"),
    !str_detect(input_source_file, "dilution_series")
  )

# Preserve legacy non-control row handling exactly. Only XPlus healthy-control
# rows can be repeated across incremental rerun exports and therefore require
# deduplication before estimating the reference mean and SD.
xplus_control_rows <- all_files %>%
  filter(is_xplus_charm_healthy_bam(.data$BAM)) %>%
  deduplicate_mrdetect_results()
all_files <- bind_rows(
  all_files %>% filter(!is_xplus_charm_healthy_bam(.data$BAM)),
  xplus_control_rows
)
rm(xplus_control_rows)

all_files <- all_files %>%
  mutate(
    Mut_source = case_when(
      str_detect(input_source_file, "BM_muts") | str_detect(source_file, "BM_muts") ~ "BM_cells",
      str_detect(input_source_file, "Blood_muts") | str_detect(source_file, "Blood_muts") ~ "Blood",
      TRUE                                   ~ NA_character_
    ),
    Filter_source = case_when(
      str_detect(source_file, "encode_only") ~ "Encode_only",
      str_detect(source_file, "STR_encode") |
        str_detect(source_file, "encode_STR_removed") ~ "STR_encode",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  # drop 'VCF' if it already exists
  select(-any_of("VCF")) %>%
  mutate(
    # 1) remove parser artifacts before using the complete result filename as
    #    provenance/key. Some historical rows contain an extra empty CSV field,
    #    which readr preserves as a leading comma in `filename`.
    filename = filename %>%
      str_remove("^,+") %>%
      str_remove("^\\./"),
    # 2) rebuild VCF from whatever follows "_VS_" in the (cleaned) filename
    VCF = sub(".*_VS_", "\\1", filename),
    # 3) strip “.mutect2…” and “.fil…” from the VCF string itself
    VCF_clean = VCF %>%
      str_remove("\\.mutect2.*") %>%
      str_remove("\\.fil.*") %>%
      str_remove("\\.somatic.*")
  )

# fix one truncated VCF_clean
all_files <- all_files %>%
  mutate(
    VCF_clean = ifelse(
      VCF_clean == "TFRIM4_0189_Bm_P_WG_ZC-02",
      "TFRIM4_0189_Bm_P_WG_ZC-02-01-O-DNA",
      VCF_clean
    )
  )


# ──────────────────────────────────────────────────────────────────────────────
# 4) Load clinical metadata and unify ‘VCF_clean’ naming
# ──────────────────────────────────────────────────────────────────────────────
cfWGS_metadata <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv",
  include_revision_extra = TRUE
) %>%
  mutate(
    VCF_clean_merge = str_remove(Bam, "\\.filter.*"),
    # correct internal PG→WG for those five patient‐specific IDs
    VCF_clean_merge = if_else(
      VCF_clean_merge %in% c(
        "TFRIM4_0031_Bm_P_PG_M4-CA-02-01-O-DNA",
        "TFRIM4_0032_Bm_P_PG_M4-HP-01-01-O-DNA",
        "TFRIM4_0033_Bm_P_PG_M4-MJ-06-01-O-DNA",
        "TFRIM4_0034_Bm_P_PG_M4-VA-02-01-O-DNA",
        "TFRIM4_0035_Bm_P_PG_M4-VA-06-01-O-DNA"
      ),
      str_replace(VCF_clean_merge, "PG", "WG"),
      VCF_clean_merge
    )
  )

cfWGS_identity_map <- build_cfwgs_sample_identity_map(cfWGS_metadata)

spring2026_panel_metadata <- cfWGS_metadata %>%
  # Spring 2026 mutation-list names are stable mutect2_pair_id values rather than
  # always being reconstructable from BAM basenames. Build a small lookup so
  # MRDetect VCF_clean can be assigned to the correct patient/sample compartment.
  filter(!is.na(mutect2_pair_id), nzchar(mutect2_pair_id)) %>%
  transmute(
    VCF_clean = mutect2_pair_id,
    VCF_panel_patient = Patient,
    VCF_panel_sample_id = Sample_ID,
    VCF_panel_sample_type = Sample_type,
    VCF_panel_timepoint_info = timepoint_info
  ) %>%
  distinct(VCF_clean, .keep_all = TRUE)

if (nrow(spring2026_panel_metadata) > 0) {
  all_files <- all_files %>%
    left_join(spring2026_panel_metadata, by = "VCF_clean") %>%
    mutate(
      # If the VCF panel metadata says the personalized mutation list came from
      # BM or blood, trust that metadata over filename heuristics.
      Mut_source = case_when(
        VCF_panel_sample_type == "BM_cells" ~ "BM_cells",
        VCF_panel_sample_type == "Blood_plasma_cfDNA" ~ "Blood",
        TRUE ~ Mut_source
      ),
      VCF_panel_is_allowed_panel_timepoint = case_when(
        # Main MRDetect scoring only uses baseline/diagnosis personalized
        # mutation lists. Later/progression mutation lists are exported to an
        # audit file and excluded to avoid circularly defining disease using a
        # future relapse/progression tumor sample.
        is.na(VCF_panel_timepoint_info) ~ TRUE,
        VCF_panel_timepoint_info %in% c("Baseline", "Diagnosis") ~ TRUE,
        TRUE ~ FALSE
      )
    )

  excluded_disallowed_panel <- all_files %>%
    filter(
      !is.na(VCF_panel_timepoint_info),
      !VCF_panel_is_allowed_panel_timepoint
    )
  if (nrow(excluded_disallowed_panel) > 0) {
    readr::write_csv(
      excluded_disallowed_panel,
      file.path(outdir, "spring2026_mrdetect_excluded_disallowed_vcf_panels.csv")
    )
  }

  all_files <- all_files %>%
    filter(VCF_panel_is_allowed_panel_timepoint) %>%
    select(-VCF_panel_is_allowed_panel_timepoint)
}

# ──────────────────────────────────────────────────────────────────────────────
# 5) Standardize column names: replace spaces & dots with underscores
# ──────────────────────────────────────────────────────────────────────────────
# 1) remove trailing spaces
colnames(all_files) <- gsub("\\s+$", "", colnames(all_files))
# 2) replace any runs of spaces with "_"
colnames(all_files) <- gsub("\\s+", "_", colnames(all_files))
# 3) replace "." with "_"
colnames(all_files) <- gsub("\\.", "_", colnames(all_files))


# -----------------------------------------------------------------------------
# 6) MRDetect parser compatibility correction
# -----------------------------------------------------------------------------
#
# The Winter 2025 MRDetect CSV files used for the manuscript were generated by a
# parser version that assigned three count columns in the wrong order. The
# manuscript analysis has always corrected that ordering before calculating
# detection-rate features. Keep this TRUE for the preserved manuscript input set.
# If future MRDetect CSVs are produced by a corrected parser, set this to FALSE
# only after confirming the input column definitions against parser documentation.
#
# Key definitions after the correction:
#   sites_detected  = number of patient-specific mutant sites with at least one read
#   reads_detected  = number of individual reads carrying a somatic mutation
#   total_reads     = total reads inspected at those sites (denominator for detection_rate)
apply_mrdetect_parser_column_correction <- TRUE
required_mrdetect_count_cols <- c("sites_detected", "reads_detected", "total_reads")
missing_mrdetect_count_cols <- setdiff(required_mrdetect_count_cols, names(all_files))
if (length(missing_mrdetect_count_cols) > 0) {
  stop(
    "MRDetect count columns are missing after column-name standardization: ",
    paste(missing_mrdetect_count_cols, collapse = ", "),
    call. = FALSE
  )
}

if (isTRUE(apply_mrdetect_parser_column_correction)) {
  temp_sites <- all_files$sites_detected
  temp_reads <- all_files$reads_detected
  temp_total <- all_files$total_reads

  all_files <- all_files %>%
    mutate(
      sites_detected = temp_reads,
      reads_detected = temp_total,
      total_reads    = temp_sites
    )

  rm(temp_sites, temp_reads, temp_total)
}

# -----------------------------------------------------------------------------
# 6b) Apply audited MRDetect rerun replacements
# -----------------------------------------------------------------------------
# Reruns are stored as a small, version-controlled table keyed by the complete
# historical result filename. This avoids editing the large raw batch exports
# and makes every replacement reviewable. Counts in this table are already in
# their corrected biological meaning (i.e., after the parser rotation above).
mrdetect_override_file <- file.path(
  "Scripts_2025", "Final_Scripts", "input_overrides",
  "VA15_mrdetect_rerun_2026-07-10.csv"
)

if (file.exists(mrdetect_override_file)) {
  mrdetect_overrides <- readr::read_csv(
    mrdetect_override_file,
    show_col_types = FALSE
  )

  override_value_cols <- c(
    "sites_checked", "reads_checked", "sites_detected", "reads_detected",
    "total_reads", "detection_rate",
    "detection_rate_as_reads_detected_over_reads_checked",
    "detection_rate_as_reads_detected_over_total_reads",
    "sites_detection_rate"
  )
  required_override_cols <- c(
    "filename", "BAM", "VCF_panel_sample_id", override_value_cols,
    "eligible_for_main_baseline_analysis", "eligibility_reason",
    "source_directory"
  )
  missing_override_cols <- setdiff(required_override_cols, names(mrdetect_overrides))
  if (length(missing_override_cols) > 0L) {
    stop(
      "MRDetect override table is missing required columns: ",
      paste(missing_override_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(mrdetect_overrides$filename)) {
    stop("MRDetect override filenames must be unique.", call. = FALSE)
  }

  expected_rates <- mrdetect_overrides %>%
    transmute(
      filename,
      detection_rate_expected = sites_detected / reads_checked,
      reads_checked_rate_expected = reads_detected / reads_checked,
      total_reads_rate_expected = reads_detected / total_reads,
      sites_rate_expected = sites_detected / sites_checked
    )
  rate_check <- mrdetect_overrides %>%
    left_join(expected_rates, by = "filename") %>%
    filter(
      abs(detection_rate - detection_rate_expected) > 1e-10 |
        abs(detection_rate_as_reads_detected_over_reads_checked - reads_checked_rate_expected) > 1e-10 |
        abs(detection_rate_as_reads_detected_over_total_reads - total_reads_rate_expected) > 1e-10 |
        abs(sites_detection_rate - sites_rate_expected) > 1e-10
    )
  if (nrow(rate_check) > 0L) {
    stop("MRDetect override rate columns do not agree with their count columns.", call. = FALSE)
  }
  if (any(
    mrdetect_overrides$sites_detected > mrdetect_overrides$sites_checked |
      mrdetect_overrides$reads_detected > mrdetect_overrides$reads_checked |
      mrdetect_overrides$reads_checked > mrdetect_overrides$total_reads
  )) {
    stop("MRDetect override counts violate denominator constraints.", call. = FALSE)
  }

  override_match_counts <- all_files %>%
    count(filename, name = "n_input_matches") %>%
    right_join(mrdetect_overrides %>% select(filename), by = "filename")
  if (any(is.na(override_match_counts$n_input_matches)) ||
      any(override_match_counts$n_input_matches < 1L)) {
    stop(
      "Each MRDetect override must match at least one historical input row. ",
      "Check the complete filename keys. Unmatched: ",
      paste(
        override_match_counts$filename[is.na(override_match_counts$n_input_matches)],
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  override_values <- mrdetect_overrides %>%
    select(filename, all_of(override_value_cols)) %>%
    rename_with(~ paste0(.x, "__override"), all_of(override_value_cols))

  all_files <- all_files %>%
    left_join(override_values, by = "filename")
  for (value_col in override_value_cols) {
    override_col <- paste0(value_col, "__override")
    all_files[[value_col]] <- dplyr::coalesce(
      all_files[[override_col]], all_files[[value_col]]
    )
  }
  all_files <- all_files %>%
    select(-ends_with("__override"))

  override_audit <- mrdetect_overrides %>%
    left_join(override_match_counts, by = "filename") %>%
    mutate(override_file = mrdetect_override_file)
  readr::write_csv(
    override_audit,
    file.path(outdir, "mrdetect_rerun_override_audit.csv")
  )
  message("Applied ", nrow(mrdetect_overrides), " audited MRDetect rerun replacements.")
}

# Validate every active row, not only the optional rerun replacements. These
# checks document the corrected column meanings and stop a future parser/input
# change from silently altering the MRDetect rates. They do not modify values.
mrdetect_rate_cols <- c(
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads", "sites_detection_rate"
)
if (any(!vapply(all_files[mrdetect_rate_cols], is.numeric, logical(1)))) {
  stop("MRDetect count and rate columns must all be numeric.", call. = FALSE)
}
if (any(!is.finite(as.matrix(all_files[mrdetect_rate_cols]))) ||
    any(as.matrix(all_files[mrdetect_rate_cols]) < 0)) {
  stop("MRDetect count and rate columns must be finite and non-negative.", call. = FALSE)
}
if (any(all_files$sites_checked == 0 | all_files$reads_checked == 0 |
        all_files$total_reads == 0)) {
  stop("MRDetect rate denominators must be greater than zero.", call. = FALSE)
}
if (any(
  all_files$sites_detected > all_files$sites_checked |
    all_files$reads_detected > all_files$reads_checked |
    all_files$reads_checked > all_files$total_reads
)) {
  stop("MRDetect counts violate their denominator constraints.", call. = FALSE)
}
mrdetect_rate_tolerance <- 1e-10
if (any(
  abs(all_files$detection_rate - all_files$sites_detected / all_files$reads_checked) > mrdetect_rate_tolerance |
    abs(all_files$detection_rate_as_reads_detected_over_reads_checked -
          all_files$reads_detected / all_files$reads_checked) > mrdetect_rate_tolerance |
    abs(all_files$detection_rate_as_reads_detected_over_total_reads -
          all_files$reads_detected / all_files$total_reads) > mrdetect_rate_tolerance |
    abs(all_files$sites_detection_rate - all_files$sites_detected / all_files$sites_checked) > mrdetect_rate_tolerance
)) {
  stop("MRDetect stored rates do not agree with their corrected count columns.", call. = FALSE)
}
rm(mrdetect_rate_cols, mrdetect_rate_tolerance)



# ──────────────────────────────────────────────────────────────────────────────
# 7) Annotate each row with ‘Study’ based on BAM prefix
# ──────────────────────────────────────────────────────────────────────────────# BAM file naming conventions used to identify cohort:
#   TFRI...  = M4 cohor)
#   MY...    = MyP cohort 
#   EGA...   = Landau et al. external validation cohort (from EGA repository)
#   HCCCFD.. = HCC patients used as non-MM cancer controls
#   HCC-...  = HCC healthy donors used as cancer-free controls
#   TGL...   = CHARM healthy donors (TGL49 panel; PRIMARY HEALTHY CONTROL POPULATION
#              used for z-score normalization throughout this analysis)
#   SPORE... = SPORE cohort (secondary MM cohort)
all_files <- all_files %>%
  mutate(
    Study = case_when(
      is_xplus_charm_healthy_bam(BAM) ~ "XPLUS_CHARM_healthy",
      str_detect(BAM, "^TFRI")       ~ "M4",
      str_detect(BAM, "^MY")         ~ "MyP",
      str_detect(BAM, "^EGA")        ~ "Landau",
      str_detect(BAM, "^HCCCFD")     ~ "HCC",
      str_detect(BAM, "^HCC\\-")      ~ "HCC_healthy",
      str_detect(BAM, "^TGL")        ~ "CHARM_healthy",
      str_detect(BAM, "^SPORE")      ~ "SPORE",
      TRUE                            ~ NA_character_
    )
  )


# ──────────────────────────────────────────────────────────────────────────────
# 8) Merge MRDetect outputs with clinical metadata
# ──────────────────────────────────────────────────────────────────────────────
tmp_meta <- cfWGS_identity_map %>%
  # Metadata is first joined by VCF_clean_merge, the historical BAM-derived key.
  # Revision rows that use mutect2_pair_id are filled in the fallback join below.
  filter(!is.na(metadata_join_vcf_clean), nzchar(metadata_join_vcf_clean)) %>%
  arrange(metadata_join_priority) %>%
  distinct(metadata_join_vcf_clean, .keep_all = TRUE) %>%
  select(-Bam, -any_of("VCF_clean_merge")) %>%
  rename(VCF_clean_merge = metadata_join_vcf_clean) %>%
  rename(Study_VCF = Study)

Merged_MRDetect <- all_files %>%
  left_join(tmp_meta, by = c("VCF_clean" = "VCF_clean_merge"))

if ("mutect2_pair_id" %in% names(tmp_meta)) {
  # Spring 2026 MRDetect VCF_clean often equals mutect2_pair_id. The fallback
  # join fills metadata columns that were missing after the historical BAM-key
  # join without overwriting successfully joined historical rows.
  tmp_meta_pair <- tmp_meta %>%
    filter(!is.na(mutect2_pair_id), nzchar(mutect2_pair_id)) %>%
    distinct(mutect2_pair_id, .keep_all = TRUE)

  pair_joined <- all_files %>%
    left_join(tmp_meta_pair, by = c("VCF_clean" = "mutect2_pair_id"))

  fill_cols <- intersect(names(Merged_MRDetect), names(pair_joined))
  metadata_cols <- setdiff(fill_cols, names(all_files))
  for (col in metadata_cols) {
    Merged_MRDetect[[col]] <- dplyr::coalesce(Merged_MRDetect[[col]], pair_joined[[col]])
  }
}

## Fix issue with ones that don't copy over

# ──────────────────────────────────────────────────────────────────────────────
# 9) Add BAM‐level sample info (Sample_ID, Patient, Sample_type, timepoint_info)
# ──────────────────────────────────────────────────────────────────────────────
# Make a helper to normalise BAM strings
normalize_bam <- function(x) {
  # BAM labels differ across callers by optional PG/WG tokens and by whether the
  # ".filter." token is present. Normalize only for joining; keep the original
  # BAM fields in outputs so provenance remains inspectable.
  x %>%
    # remove internal _PG_ or _WG_ tokens
    str_replace_all("_[PW]G_", "_") %>%
    str_replace("^([PW]G)_", "") %>%
    str_replace("_([PW]G)$", "") %>%
    # remove the .filter. segment if present
    str_replace("\\.filter\\.", ".")
}

# Build bam_info with normalised key
bam_info <- cfWGS_identity_map %>%
  filter(!is.na(metadata_join_bam), nzchar(metadata_join_bam)) %>%
  mutate(Bam_norm = normalize_bam(metadata_join_bam)) %>%
  arrange(metadata_join_priority) %>%
  distinct(Bam_norm, .keep_all = TRUE) %>%
  transmute(
    Bam = metadata_join_bam,
    Bam_norm,
    Sample_ID,
    Patient,
    Study,
    # Platform follows the physical sequencing run, not whether the biological
    # specimen overlaps an existing cohort row. Keep the legacy column name for
    # downstream compatibility; it now marks every revision-run XPlus BAM,
    # including the IMG-159 technical repeat.
    revision_new_xplus_bam = Study == "IMMAGINE_revision_OICR",
    Sample_type,
    timepoint_info
  ) %>%
  rename_with(~ paste0(.x, "_Bam"),
              .cols = c(Sample_ID, Patient, Study, revision_new_xplus_bam, Sample_type, timepoint_info))

# Add the same normalized key to Merged_MRDetect
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(BAM_norm = str_replace_all(BAM, "_[PW]G_", "_") |> 
           str_replace("_([PW]G)$", "") |> 
           str_replace("\\.filter\\.", ".") |> 
           str_replace("^([PW]G)_", ""))

# Join using the normalized key
Merged_MRDetect <- Merged_MRDetect %>%
  left_join(bam_info, by = c("BAM_norm" = "Bam_norm"))

# optional: drop the helper columns if you don’t need them in the final table
Merged_MRDetect <- Merged_MRDetect %>% select(-BAM_norm)

## Old way
# bam_info <- cfWGS_metadata %>%
#   select(Bam, Sample_ID, Patient, Sample_type, timepoint_info) %>%
#   rename_with(~ paste0(.x, "_Bam"), .cols = -Bam)
# 
# Merged_MRDetect <- Merged_MRDetect %>%
#   left_join(bam_info, by = c("BAM" = "Bam"))
# 

Merged_MRDetect %>%
  filter(is.na(Patient_Bam)) %>%
  pull(BAM) %>% 
  unique() # diagnostic only: non-empty output means BAM metadata did not join.
  
# ──────────────────────────────────────────────────────────────────────────────
# 10) Classify matched vs. unmatched plasma; force CHARM_healthy→cfDNA
# ──────────────────────────────────────────────────────────────────────────────
# MRDetect runs each cfDNA sample (BAM) against the VCF from one patient's
# bone-marrow mutations. That VCF is the "patient" the sample is being tested against.
#   Matched_plasma   = BAM and VCF come from the SAME patient (the true signal)
#   Unmatched_plasma = BAM and VCF come from DIFFERENT patients
#                      (used as cross-patient negative controls)
# All CHARM_healthy donors are always run against a patient's VCF, so their
# Sample_type_Bam field is set to Blood_plasma_cfDNA here.
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(
    plotting_type       = if_else(Patient_Bam == Patient, "Matched_plasma", "Unmatched_plasma"),
    Sample_type_Bam     = if_else(
      Study %in% c("CHARM_healthy", "XPLUS_CHARM_healthy"),
      "Blood_plasma_cfDNA",
      Sample_type_Bam
    )
  )


# ──────────────────────────────────────────────────────────────────────────────
# 11) Subset to cfDNA timepoints only.
# The queried BAM may be a longitudinal cfDNA sample, but the tumor-mutation VCF
# panel above is restricted to Diagnosis/Baseline labels. Relapse-context samples
# that are the baseline for a follow-up series should be labelled Baseline in the
# metadata before entering this analysis.
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect_cfDNA <- Merged_MRDetect %>%
  filter(Sample_type_Bam == "Blood_plasma_cfDNA",
         timepoint_info_Bam %in% c("Baseline","Diagnosis","Progression"))
# This subset is retained as an interactive diagnostic object. The saved raw,
# z-scored, BM-informed, and blood-informed outputs below are built from
# `Merged_MRDetect`/`Merged_MRDetect_zscore`, not from this object.


# ──────────────────────────────────────────────────────────────────────────────
# 12) Compute selected healthy-control means & SDs → join back for z-scores
# ──────────────────────────────────────────────────────────────────────────────
Merged_MRDetect$VCF_factor <- factor(Merged_MRDetect$VCF, levels = unique(Merged_MRDetect$VCF))

# Z-SCORE NORMALIZATION STRATEGY:
# For each patient VCF (VCF_factor), we compute a healthy-control reference
# distribution using healthy donors run against that same VCF.
# This is very important: each patient has a different mutation set (different VCF),
# so the background noise level differs per patient. By grouping on VCF_factor,
# we get a VCF-specific mean and SD for healthy controls, ensuring that the
# z-score for a patient sample is relative to healthy donors tested against
# that exact same set of patient mutations.
# Historical queried samples use the original CHARM background. Spring 2026
# revision queried samples request the XPlus CHARM background only when they are
# new revision samples rather than original-repo BAMs carried in the revision
# metadata for harmonized naming.
#   z = (patient_rate - mean_HC_rate) / sd_HC_rate
# High z-score (>>2) → detection rate significantly above healthy-control noise,
# suggesting true circulating tumor DNA.
Merged_MRDetect <- Merged_MRDetect %>%
  mutate(
    healthy_reference_requested = case_when(
      # Spring 2026 revision queried BAMs were sequenced on the NovaSeq X Plus
      # workflow, so their MRDetect z-scores should use XPlus CHARM controls.
      # Biological overlap does not change the sequencing platform.
      revision_new_xplus_bam_Bam %in% TRUE ~ "XPLUS_CHARM_healthy",
      Study == "XPLUS_CHARM_healthy" ~ "XPLUS_CHARM_healthy",
      TRUE ~ "CHARM_healthy"
    )
  )

healthy_reference_candidates <- Merged_MRDetect %>%
  filter(Study %in% c("CHARM_healthy", "XPLUS_CHARM_healthy", "HCC_healthy")) %>%
  filter(
    .data$Study != "XPLUS_CHARM_healthy" |
      is_xplus_mrdetect_reference_bam(.data$BAM)
  ) %>%
  mutate(
    healthy_reference_tier = case_when(
      # Keep XPlus CHARM controls separate from legacy CHARM controls. Historical
      # rows request CHARM_healthy below; Spring 2026 revision rows request
      # XPLUS_CHARM_healthy and are not silently normalized to regular controls.
      Study == "CHARM_healthy" ~ "CHARM_healthy",
      Study == "XPLUS_CHARM_healthy" ~ "XPLUS_CHARM_healthy",
      TRUE ~ "non_mm_healthy_fallback"
    )
  )

healthy_reference_counts <- healthy_reference_candidates %>%
  count(
    VCF_factor, Mut_source, Filter_source, healthy_reference_tier,
    name = "n_healthy_reference_rows"
  )

xplus_control_coverage_audit <- healthy_reference_candidates %>%
  filter(.data$healthy_reference_tier == "XPLUS_CHARM_healthy") %>%
  group_by(.data$VCF_factor, .data$Mut_source, .data$Filter_source) %>%
  summarise(
    n_xplus_control_bams = n_distinct(.data$BAM),
    xplus_control_bams = paste(sort(unique(basename(.data$BAM))), collapse = ";"),
    n_expected_allowlist_bams = length(xplus_mrdetect_reference_bams()),
    missing_allowlist_bams = paste(
      setdiff(xplus_mrdetect_reference_bams(), unique(basename(.data$BAM))),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(.data$VCF_factor, .data$Mut_source, .data$Filter_source)
write_csv(
  xplus_control_coverage_audit,
  file.path(outdir, "spring2026_mrdetect_xplus_control_coverage_audit.csv")
)

healthy_reference_requests <- Merged_MRDetect %>%
  distinct(VCF_factor, Mut_source, Filter_source, healthy_reference_requested)

healthy_reference_choice <- healthy_reference_requests %>%
  tidyr::crossing(
    healthy_reference_tier = c(
      "CHARM_healthy",
      "XPLUS_CHARM_healthy",
      "non_mm_healthy_fallback"
    )
  ) %>%
  left_join(
    healthy_reference_counts,
    by = c("VCF_factor", "Mut_source", "Filter_source", "healthy_reference_tier")
  ) %>%
  mutate(
    n_healthy_reference_rows = coalesce(n_healthy_reference_rows, 0L),
    healthy_reference_priority = case_when(
      healthy_reference_requested == "XPLUS_CHARM_healthy" &
        healthy_reference_tier == "XPLUS_CHARM_healthy" ~ 1L,
      healthy_reference_requested == "CHARM_healthy" &
        healthy_reference_tier == "CHARM_healthy" ~ 1L,
      healthy_reference_requested == "CHARM_healthy" &
        healthy_reference_tier == "non_mm_healthy_fallback" ~ 2L,
      TRUE ~ 99L
    )
  ) %>%
  filter(healthy_reference_priority < 99L) %>%
  group_by(VCF_factor, Mut_source, Filter_source, healthy_reference_requested) %>%
  arrange(healthy_reference_priority, desc(n_healthy_reference_rows), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    healthy_reference_status = if_else(
      n_healthy_reference_rows > 0,
      "available",
      "missing"
    ),
    healthy_reference_tier = if_else(
      n_healthy_reference_rows > 0,
      healthy_reference_tier,
      NA_character_
    )
  )

write_csv(
  healthy_reference_choice,
  file.path(outdir, "spring2026_mrdetect_healthy_reference_choice_audit.csv")
)

zscore_lookup <- healthy_reference_candidates %>%
  # Compute VCF-specific healthy-control means and SDs after the tier choice.
  # This keeps each patient mutation-list background estimate separate.
  select(-healthy_reference_requested) %>%
  inner_join(
    healthy_reference_choice %>%
      filter(healthy_reference_status == "available") %>%
      select(
        VCF_factor, Mut_source, Filter_source,
        healthy_reference_requested, healthy_reference_tier
      ),
    by = c("VCF_factor", "Mut_source", "Filter_source", "healthy_reference_tier")
  ) %>%
  select(VCF_factor, Mut_source, Filter_source, healthy_reference_requested,
         detection_rate,
         detection_rate_as_reads_detected_over_reads_checked,
         detection_rate_as_reads_detected_over_total_reads,
         sites_detection_rate) %>%
  group_by(VCF_factor, Mut_source, Filter_source, healthy_reference_requested) %>%
  summarize(
    mean_det_charm                     = mean(detection_rate, na.rm = TRUE),
    sd_det_charm                       = sd(detection_rate, na.rm = TRUE),
    mean_det_checked_charm             = mean(detection_rate_as_reads_detected_over_reads_checked, na.rm = TRUE),
    sd_det_checked_charm               = sd(detection_rate_as_reads_detected_over_reads_checked, na.rm = TRUE),
    mean_det_total_charm               = mean(detection_rate_as_reads_detected_over_total_reads, na.rm = TRUE),
    sd_det_total_charm                 = sd(detection_rate_as_reads_detected_over_total_reads, na.rm = TRUE),
    mean_sites_charm                   = mean(sites_detection_rate, na.rm = TRUE),
    sd_sites_charm                     = sd(sites_detection_rate, na.rm = TRUE)
  ) %>%
  ungroup()

Merged_MRDetect_zscore <- Merged_MRDetect %>%
  left_join(
    healthy_reference_choice %>%
      select(
        VCF_factor, Mut_source, Filter_source,
        healthy_reference_requested, healthy_reference_tier,
        n_healthy_reference_rows, healthy_reference_status
      ),
    by = c("VCF_factor", "Mut_source", "Filter_source", "healthy_reference_requested")
  ) %>%
  left_join(
    zscore_lookup,
    by = c("VCF_factor", "Mut_source", "Filter_source", "healthy_reference_requested")
  ) %>%
  mutate(
    # These z-scores are descriptive MRDetect features. If the selected healthy
    # reference has zero/NA SD, the resulting z-score remains NA rather than
    # inventing a finite value.
    detection_rate_zscore_charm                   = (detection_rate - mean_det_charm) / sd_det_charm,
    detection_rate_zscore_reads_checked_charm     = (detection_rate_as_reads_detected_over_reads_checked - mean_det_checked_charm) / sd_det_checked_charm,
    detection_rate_zscore_total_reads_charm       = (detection_rate_as_reads_detected_over_total_reads - mean_det_total_charm) / sd_det_total_charm,
    sites_rate_zscore_charm                       = (sites_detection_rate - mean_sites_charm) / sd_sites_charm,
    # Difference vs healthy controls (absolute values)
    detection_rate_diff_vs_hc                     = detection_rate - mean_det_charm,
    detection_rate_checked_diff_vs_hc             = detection_rate_as_reads_detected_over_reads_checked - mean_det_checked_charm,
    detection_rate_total_diff_vs_hc               = detection_rate_as_reads_detected_over_total_reads - mean_det_total_charm,
    sites_rate_diff_vs_hc                         = sites_detection_rate - mean_sites_charm,
    # Difference vs healthy controls (as percentages)
    detection_rate_checked_diff_vs_hc_pct         = (detection_rate_checked_diff_vs_hc) * 100
  )


# ──────────────────────────────────────────────────────────────────────────────
# 13) Identify and drop any duplicate runs (keep max total_reads)
# ──────────────────────────────────────────────────────────────────────────────
dup_keys <- c("BAM","Mut_source","Filter_source","VCF","VCF_clean","Study",
              "Sample_ID_Bam","Patient_Bam","VCF_factor")
Merged_MRDetect_zscore <- Merged_MRDetect_zscore %>%
  # If multiple MRDetect rows represent the same BAM/VCF/source/filter run, keep
  # the row with the most total reads because it has the largest denominator for
  # estimating detection rate.
  group_by(across(all_of(dup_keys))) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup()


# ──────────────────────────────────────────────────────────────────────────────
# 14) Save both “raw” and “z-scored” MRDetect tables
# ──────────────────────────────────────────────────────────────────────────────
base_name_raw    <- paste0(project, "All_MRDetect_Sep2025")
base_name_zscore <- paste0(project, "All_MRDetect_with_Zscore_Sep2025")

write.table(Merged_MRDetect,
            file = file.path(outdir, paste0(base_name_raw, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect,
        file = file.path(outdir, paste0(base_name_raw, ".rds")))

write.table(Merged_MRDetect_zscore,
            file = file.path(outdir, paste0(base_name_zscore, ".txt")),
            sep = "\t", row.names = FALSE)

saveRDS(Merged_MRDetect_zscore,
        file = file.path(outdir, paste0(base_name_zscore, ".rds")))

message("→ MRDetect tables (raw & z-scored) written to: ", outdir)


# ──────────────────────────────────────────────────────────────────────────────
# 14a) Export the historical CHARM healthy-control reference rows
# ──────────────────────────────────────────────────────────────────────────────
# Despite its retained legacy filename, this table contains CHARM_healthy rows
# only. It is the historical control input used by scripts 1_10 and 3_1_part2;
# changing its membership here would change downstream calibration analyses.
detection_rates_all_samples <- Merged_MRDetect_zscore %>%
  filter(
    # Historical CHARM controls use the baseline-family labels 01/T0.
    (Timepoint %in% c("01", "T0") & Study == "CHARM_healthy")
  ) %>%
  select(
    Study, Patient, Patient_Bam, Sample_ID_Bam, BAM, VCF, VCF_clean,
    Sample_type_Bam, timepoint_info_Bam, Mut_source, Filter_source,
    sites_checked, reads_checked, sites_detected, reads_detected, total_reads,
    detection_rate, 
    detection_rate_as_reads_detected_over_reads_checked,
    detection_rate_as_reads_detected_over_total_reads,
    sites_detection_rate,
    detection_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm,
    detection_rate_zscore_total_reads_charm,
    sites_rate_zscore_charm
  ) %>%
  arrange(Study, Patient_Bam, Mut_source, Filter_source)

# Export as CSV and RDS
write_csv(detection_rates_all_samples,
          file = file.path(outdir, "All_detection_rates_baseline_and_controls_Feb2026.csv"))
saveRDS(detection_rates_all_samples,
        file = file.path(outdir, "All_detection_rates_baseline_and_controls_Feb2026.rds"))

message("→ Historical CHARM healthy-control detection rates written to: ", outdir)
message("  Total samples exported: ", nrow(detection_rates_all_samples))


## See what is NA
# show BAM file paths where Patient_Bam is NA
Merged_MRDetect_zscore %>%
  filter(is.na(Patient_Bam)) %>%
  pull(BAM) %>% unique()

# ──────────────────────────────────────────────────────────────────────────────
# 15) Additional Processing - BM muts
# ──────────────────────────────────────────────────────────────────────────────
### Now Process Further
### First BM Muts

# 1) Start from the z-scored MRDetect table
df <- Merged_MRDetect_zscore %>%
  filter(plotting_type == "Matched_plasma",
         Mut_source      == "BM_cells",
         Filter_source   == "STR_encode")

# Rename df's columns to the names num_cols vector expects
rename_map <- c(
  mean_detection_rate_charm                       = "mean_det_charm",
  sd_detection_rate_charm                         = "sd_det_charm",
  mean_detection_rate_reads_checked_charm         = "mean_det_checked_charm",
  sd_detection_rate_reads_checked_charm           = "sd_det_checked_charm",
  mean_detection_rate_total_reads_charm           = "mean_det_total_charm",
  sd_detection_rate_total_reads_charm             = "sd_det_total_charm",
  mean_sites_rate_charm                           = "mean_sites_charm",
  sd_sites_rate_charm                             = "sd_sites_charm"
)

df <- df %>% rename(!!!rename_map)

# 2) Collapse duplicate BAM runs by averaging all numeric QC/MRD metrics
num_cols <- c(
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "mean_detection_rate_charm", "sd_detection_rate_charm",
  "mean_detection_rate_reads_checked_charm", "sd_detection_rate_reads_checked_charm",
  "mean_detection_rate_total_reads_charm", "sd_detection_rate_total_reads_charm",
  "mean_sites_rate_charm", "sd_sites_rate_charm",
  "detection_rate_zscore_charm","sites_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm"
)

df <- df %>%
  group_by(
    Sample_ID, Sample_ID_Bam, VCF, VCF_clean, Study, Patient,
    Date_of_sample_collection, Sample_type, Timepoint, Study_VCF,
    timepoint_info, Sample_type_Bam, timepoint_info_Bam,
    Mut_source, Filter_source, plotting_type, VCF_factor
  ) %>%
  summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3) Recompute MRD status & cumulative VAF
# Mrd_by_WGS threshold of 4.5 z-scores: sites_rate_zscore_charm measures how many
# standard deviations above the CHARM healthy-control mean this sample's
# site-detection rate falls. A z-score > 4.5 was chosen as a conservative positive
# threshold to maximize specificity at this screening step. The final model-based
# probability thresholds in script 3_1 are more precisely optimized to the
# frontline training data and should be used for clinical calls.
#
# Cumulative_VAF: the mutation detection rate (reads_detected / reads_checked)
# for MRD-positive samples only; set to 0 for negatives. This is the primary
# quantitative cfDNA feature passed to downstream analysis scripts.
df <- df %>%
  mutate(
    Mrd_by_WGS    = if_else(sites_rate_zscore_charm > 4.5, "Positive", "Negative"),
    Cumulative_VAF = if_else(Mrd_by_WGS == "Positive",
                             detection_rate_as_reads_detected_over_reads_checked,
                             0)
  )

## Check for dups 
# Diagnostics: where are the dupes?
cfWGS_metadata %>%
  count(Sample_ID) %>% filter(n > 1) %>% arrange(desc(n)) -> meta_dupes

df %>%
  count(Sample_ID, Sample_ID_Bam) %>% filter(n > 1) %>% arrange(desc(n)) -> df_pair_dupes

message("# meta dupes (by Sample_ID): ", nrow(meta_dupes))
message("# df pair dupes (by Sample_ID, Sample_ID_Bam): ", nrow(df_pair_dupes))

# Make metadata unique per Sample_ID (break ties deterministically)
cc <- cfWGS_metadata %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection)) %>%
  arrange(Sample_ID, desc(Date_of_sample_collection)) %>%
  group_by(Sample_ID) %>%
  # keep only ONE row per Sample_ID (latest date; if ties, take the first after arrange)
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Sample_ID, Date_of_sample_collection, Timepoint)

# sanity check
stopifnot(nrow(cc) == dplyr::n_distinct(cc$Sample_ID))

# Also ensure df has unique pairs before joining dates
df_pairs <- df %>%
  select(Sample_ID, Sample_ID_Bam, Patient) %>%
  distinct()

# 4) Join dates for both IDs
start_dates <- df_pairs %>%
  left_join(cc, by = "Sample_ID") %>%
  rename(Date_of_sample_collection_Sample_ID = Date_of_sample_collection,
         Timepoint_Sample_ID                  = Timepoint) %>%
  left_join(cc, by = c("Sample_ID_Bam" = "Sample_ID")) %>%
  rename(Date_of_sample_collection_Sample_ID_Bam = Date_of_sample_collection,
         Timepoint_Sample_ID_Bam                  = Timepoint) %>%
  mutate(
    num_days  = as.numeric(
      difftime(Date_of_sample_collection_Sample_ID_Bam,
               Date_of_sample_collection_Sample_ID,
               units = "days")
    ),
    num_weeks = num_days / 7
  )

# Final guard: confirm uniqueness using the same key used for the downstream join.
start_date_dupes <- start_dates %>%
  count(Sample_ID, Sample_ID_Bam, Patient) %>%
  filter(n > 1)
if (nrow(start_date_dupes) > 0) {
  readr::write_csv(
    start_date_dupes,
    file.path(outdir, "start_date_duplicate_key_audit.csv")
  )
  stop(
    "Duplicate start-date rows remain for the Sample_ID/Sample_ID_Bam/Patient join key. ",
    "Audit written to: ", file.path(outdir, "start_date_duplicate_key_audit.csv"),
    call. = FALSE
  )
}

# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_marrow” from the active integrated WGS feature table
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds")
cohort_df <- readRDS("cohort_assignment_table_updated.rds")
good_pts <- All_feature_data %>%
  filter(Sample_type == "BM_cells",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

spring2026_good_bm_pts <- cfWGS_metadata %>%
  filter(
    Study == "IMMAGINE_revision_OICR",
    Sample_type == "BM_cells",
    timepoint_info %in% c("Diagnosis", "Baseline"),
    !is.na(mutect2_pair_id),
    nzchar(mutect2_pair_id)
  ) %>%
  pull(Patient) %>%
  unique()

good_pts <- union(good_pts, spring2026_good_bm_pts)

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_marrow = if_else(Patient %in% good_pts, "Yes", "No")
  )

# 7) Calculate START_DATE and percent_change from baseline
# Percent change is undefined when the reference rate is zero. R therefore
# records 0/0 as NA/NaN; retain that missing value rather than calling it 0%.
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(Good_baseline_marrow == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 

## Get the change since the first treatment timepoint
# The same zero-reference rule applies to treatment-relative percent change.
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- filter_to_baseline_mutation_sources(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "bm"
)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "bm"
)

## SPORE_0012 uses its T4 BM as the designated baseline mutation source. The
## MRDetect input matrix also contains older T1/T2 plasma queries for this BM
## panel; those samples predate the mutation source and therefore cannot be
## scored prospectively. Keep the T4->T5 comparison and audit/remove only the
## temporally impossible pre-T4 pairs.
spore0012_prebaseline_query_audit <- combined_data_plot %>%
  filter(
    Patient == "SPORE_0012",
    Sample_ID == "SPORE_0012_T4_BM_cells",
    !is.na(START_DATE),
    START_DATE < 0
  ) %>%
  mutate(exclusion_reason = "tested plasma predates designated T4 BM baseline")

readr::write_csv(
  spore0012_prebaseline_query_audit,
  file.path(export_dir, "spore0012_prebaseline_bm_query_exclusion_audit.csv")
)

combined_data_plot <- combined_data_plot %>%
  filter(
    !(
      Patient == "SPORE_0012" &
        Sample_ID == "SPORE_0012_T4_BM_cells" &
        !is.na(START_DATE) &
        START_DATE < 0
    )
  )

combined_data_plot <- deduplicate_baseline_source_alias_query_rows(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "bm"
)

bm_source_duplicate_audit <- combined_data_plot %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
  count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
  filter(n_baseline_diagnosis_mutation_sources > 1) %>%
  left_join(
    combined_data_plot %>%
      filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
      distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
      group_by(Patient, Mut_source) %>%
      summarise(
        mutation_sources = paste(
          sort(unique(paste(Sample_ID, Timepoint, timepoint_info, sep = "|"))),
          collapse = "; "
        ),
        .groups = "drop"
      ),
    by = c("Patient", "Mut_source")
  )

write_csv(
  bm_source_duplicate_audit,
  file.path(export_dir, "spring2026_bm_mrdetect_duplicate_baseline_source_audit.csv")
)

if (nrow(bm_source_duplicate_audit) > 0L) {
  warning(
    "Patients with more than one baseline/diagnosis BM-mutation source were found. ",
    "Audit written to: ",
    file.path(export_dir, "spring2026_bm_mrdetect_duplicate_baseline_source_audit.csv")
  )
}

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS_MRDetect_BM_data_updated_Sep.csv")
)



# ──────────────────────────────────────────────────────────────────────────────
# 15) Additional Processing - Blood muts
# ──────────────────────────────────────────────────────────────────────────────
### Next for blood-derived muts
#Merged_MRDetect_zscore <- readRDS(file = "MRDetect_output_winter_2025/Processed_R_outputs/cfWGS_Winter2025All_MRDetect_with_Zscore_Sep2025.rds")
# 1) Start from the z-scored MRDetect table
df <- Merged_MRDetect_zscore %>%
  filter(plotting_type == "Matched_plasma",
         Mut_source      == "Blood",
         Filter_source   == "STR_encode")

# Rename df's columns to the names num_cols vector expects
rename_map <- c(
  mean_detection_rate_charm                       = "mean_det_charm",
  sd_detection_rate_charm                         = "sd_det_charm",
  mean_detection_rate_reads_checked_charm         = "mean_det_checked_charm",
  sd_detection_rate_reads_checked_charm           = "sd_det_checked_charm",
  mean_detection_rate_total_reads_charm           = "mean_det_total_charm",
  sd_detection_rate_total_reads_charm             = "sd_det_total_charm",
  mean_sites_rate_charm                           = "mean_sites_charm",
  sd_sites_rate_charm                             = "sd_sites_charm"
)

df <- df %>% rename(!!!rename_map)

# 2) Collapse duplicate BAM runs by averaging all numeric QC/MRD metrics
num_cols <- c(
  "sites_checked", "reads_checked", "sites_detected", "reads_detected",
  "total_reads", "detection_rate",
  "detection_rate_as_reads_detected_over_reads_checked",
  "detection_rate_as_reads_detected_over_total_reads",
  "sites_detection_rate", "mean_detection_rate_charm", "sd_detection_rate_charm",
  "mean_detection_rate_reads_checked_charm", "sd_detection_rate_reads_checked_charm",
  "mean_detection_rate_total_reads_charm", "sd_detection_rate_total_reads_charm",
  "mean_sites_rate_charm", "sd_sites_rate_charm",
  "detection_rate_zscore_charm","sites_rate_zscore_charm",
  "detection_rate_zscore_reads_checked_charm",
  "detection_rate_zscore_total_reads_charm"
)

df <- df %>%
  group_by(
    Sample_ID, Sample_ID_Bam, VCF, VCF_clean, Study, Patient,
    Date_of_sample_collection, Sample_type, Timepoint, Study_VCF,
    timepoint_info, Sample_type_Bam, timepoint_info_Bam,
    Mut_source, Filter_source, plotting_type, VCF_factor
  ) %>%
  summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3) Recompute MRD status & cumulative VAF
df <- df %>%
  mutate(
    Mrd_by_WGS    = if_else(sites_rate_zscore_charm > 4.5, "Positive", "Negative"),
    Cumulative_VAF = if_else(Mrd_by_WGS == "Positive",
                             detection_rate_as_reads_detected_over_reads_checked,
                             0)
  ) ## Old version

# 4) Pull in clinical dates & timepoints for both Sample_ID and Sample_ID_Bam
cc <- cfWGS_metadata %>%
  group_by(Sample_ID) %>%
  slice_max(Date_of_sample_collection, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample_ID, Date_of_sample_collection, Timepoint)

# Build start_dates _including_ Patient
start_dates <- df %>%
  select(Sample_ID, Sample_ID_Bam, Patient) %>% 
  left_join(cc,                    by = "Sample_ID")           %>% 
  rename(Date_of_sample_collection_Sample_ID = Date_of_sample_collection,
         Timepoint_Sample_ID            = Timepoint)           %>% 
  left_join(cc, by = c("Sample_ID_Bam" = "Sample_ID"))         %>% 
  rename(Date_of_sample_collection_Sample_ID_Bam = Date_of_sample_collection,
         Timepoint_Sample_ID_Bam            = Timepoint)       %>% 
  mutate(
    num_days  = as.numeric(difftime(
      Date_of_sample_collection_Sample_ID_Bam,
      Date_of_sample_collection_Sample_ID,
      units = "days")),
    num_weeks = num_days / 7
  ) %>% 
  distinct()


# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_sample” from All_feature_data
good_pts <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_sample = if_else(Patient %in% good_pts, "Yes", "No")
  )

filtered_good_pts <- tibble(Patient = good_pts) %>%
  inner_join(cohort_df, by = "Patient")

filtered_good_pts %>%
  pull(Patient) %>% unique()

# 7) Calculate START_DATE and percent_change from baseline
# A zero reference rate yields an undefined percent change and remains missing.
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(Good_baseline_sample == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 


## Now edit this to show the time difference since the second timepoint 
## First calculate the percent change in detection rate and put dates

## Get the change since the first treatment timepoint
# A zero first-treatment rate likewise yields a missing relative change.
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment

## Get rid of marrows that are not at baseline 
combined_data_plot <- combined_data_plot %>%
  filter(Timepoint %in% c("01", "T0", "T1", "1"))
combined_data_plot <- combined_data_plot %>% 
  filter(Sample_type == "Blood_plasma_cfDNA")


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- filter_to_baseline_mutation_sources(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_good_baseline"
)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_good_baseline"
)

combined_data_plot <- deduplicate_baseline_source_alias_query_rows(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_good_baseline"
)

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS MRDetect Blood data updated Sep.csv")
)


## Redo for all bloods available
# 5) Attach start_dates back onto df
combined_data_plot <- df %>%
  left_join(start_dates, by = c("Sample_ID", "Sample_ID_Bam", "Patient"))

# 6) Flag “Good_baseline_sample” from All_feature_data
good_pts <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>% unique()

combined_data_plot <- combined_data_plot %>%
  mutate(
    Good_baseline_sample = if_else(Patient %in% good_pts, "Yes", "No")
  )

# 7) Calculate START_DATE and percent_change from baseline
# A zero reference rate yields an undefined percent change and remains missing.
combined_data_plot <- combined_data_plot %>%
  rename(START_DATE = num_days) %>%
  filter(#Good_baseline_sample == "Yes",
         Sample_type_Bam == "Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  mutate(
    baseline_date = min(START_DATE[timepoint_info %in% c("Baseline","Diagnosis")], na.rm = TRUE),
    baseline_rate = detection_rate_as_reads_detected_over_reads_checked[START_DATE == baseline_date][1],
    absolute_change = detection_rate_as_reads_detected_over_reads_checked - baseline_rate,
    percent_change  = if_else(!is.na(baseline_rate),
                              (absolute_change / baseline_rate) * 100,
                              NA_real_),
    Weeks_since_baseline = pmax(0, START_DATE) / 7
  ) %>%
  ungroup() 


## Now edit this to show the time difference since the second timepoint 
## First calculate the percent change in detection rate and put dates

## Get the change since the first treatment timepoint
# A zero first-treatment rate likewise yields a missing relative change.
# Define the treatment sample categories
treatment_samples <- c("Post_induction", "Post_transplant", "Maintenance", "1.5yr maintenance")

combined_data_plot <- combined_data_plot %>%
  dplyr::group_by(Patient) %>%
  dplyr::arrange(START_DATE) %>%  # Ensure entries are sorted by days since baseline within each patient
  dplyr::mutate(
    # Find the first treatment sample date if available
    first_treatment_date = ifelse(any(timepoint_info_Bam %in% treatment_samples),
                                  START_DATE[timepoint_info_Bam %in% treatment_samples][1],
                                  NA_real_),
    # Calculate the detection rate for the first treatment timepoint
    treatment_detection_rate = ifelse(!is.na(first_treatment_date),
                                      detection_rate_as_reads_detected_over_reads_checked[START_DATE == first_treatment_date],
                                      NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    # Calculate the absolute and percent change in detection rate since the first treatment timepoint
    absolute_change_detection_rate_treatment = detection_rate_as_reads_detected_over_reads_checked - treatment_detection_rate,
    percent_change_detection_rate_treatment = (absolute_change_detection_rate_treatment / treatment_detection_rate) * 100,
    # Calculate weeks since the first treatment date
    Weeks_since_first_treatment = (START_DATE - first_treatment_date) / 7
  )

# For consistency 
combined_data_plot$Weeks_since_second_start <- combined_data_plot$Weeks_since_first_treatment
combined_data_plot$percent_change_detection_rate_second_timepoint <- combined_data_plot$percent_change_detection_rate_treatment

## Get rid of marrows that are not at baseline 
combined_data_plot <- combined_data_plot %>%
  filter(Timepoint %in% c("01", "T0", "T1", "1"))
combined_data_plot <- combined_data_plot %>% 
  filter(Sample_type == "Blood_plasma_cfDNA")


# 8) Write out 
export_dir <- "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

combined_data_plot <- filter_to_baseline_mutation_sources(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_all_patients"
)

combined_data_plot <- restrict_to_single_baseline_mutation_source(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_all_patients"
)

combined_data_plot <- deduplicate_baseline_source_alias_query_rows(
  combined_data_plot,
  audit_dir = export_dir,
  audit_prefix = "blood_all_patients"
)

blood_source_duplicate_audit <- combined_data_plot %>%
  filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
  count(Patient, Mut_source, name = "n_baseline_diagnosis_mutation_sources") %>%
  filter(n_baseline_diagnosis_mutation_sources > 1) %>%
  left_join(
    combined_data_plot %>%
      filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
      distinct(Patient, Mut_source, Sample_ID, Timepoint, timepoint_info) %>%
      group_by(Patient, Mut_source) %>%
      summarise(
        mutation_sources = paste(
          sort(unique(paste(Sample_ID, Timepoint, timepoint_info, sep = "|"))),
          collapse = "; "
        ),
        .groups = "drop"
      ),
    by = c("Patient", "Mut_source")
  )

write_csv(
  blood_source_duplicate_audit,
  file.path(export_dir, "spring2026_blood_mrdetect_duplicate_baseline_source_audit.csv")
)

if (nrow(blood_source_duplicate_audit) > 0L) {
  warning(
    "Patients with more than one baseline/diagnosis blood-mutation source were found. ",
    "Audit written to: ",
    file.path(export_dir, "spring2026_blood_mrdetect_duplicate_baseline_source_audit.csv")
  )
}

readr::write_csv(
  combined_data_plot,
  file = file.path(export_dir, "cfWGS MRDetect Blood data updated Sep with all patients.csv")
)

baseline_source_inventory_from_export <- function(path, source_label) {
  if (!file.exists(path)) {
    return(tibble())
  }

  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      mutation_source_type = source_label,
      source_key = coalesce(
        na_if(as.character(VCF_factor), ""),
        na_if(as.character(VCF_clean), ""),
        na_if(as.character(VCF), ""),
        na_if(as.character(Sample_ID), "")
      )
    ) %>%
    group_by(
      Patient, mutation_source_type, source_key, Sample_ID, VCF_clean,
      Timepoint, timepoint_info, Date_of_sample_collection_Sample_ID
    ) %>%
    summarise(
      n_queried_blood_bams = n_distinct(Sample_ID_Bam, na.rm = TRUE),
      queried_timepoints = paste(
        sort(unique(na.omit(as.character(Timepoint_Sample_ID_Bam)))),
        collapse = "; "
      ),
      queried_timepoint_labels = paste(
        sort(unique(na.omit(as.character(timepoint_info_Bam)))),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    arrange(Patient, mutation_source_type, Timepoint, source_key)
}

baseline_source_inventory <- bind_rows(
  baseline_source_inventory_from_export(
    "MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/cfWGS_MRDetect_BM_data_updated_Sep.csv",
    "BM_cells"
  ),
  baseline_source_inventory_from_export(
    "MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated Sep with all patients.csv",
    "Blood_plasma_cfDNA"
  )
)

dir.create("Output_tables_2025/clinical_support", recursive = TRUE, showWarnings = FALSE)

baseline_source_inventory_ranked <- baseline_source_inventory %>%
  mutate(
    inventory_source_date = suppressWarnings(as.Date(Date_of_sample_collection_Sample_ID)),
    inventory_timepoint_rank = case_when(
      str_to_upper(Timepoint) %in% c("T0", "TP0", "D0", "0", "01", "1") ~ 1,
      str_to_upper(Timepoint) %in% c("T1", "TP1") ~ 2,
      str_detect(Timepoint, "^[Tt]?[0-9]+$") ~ as.numeric(str_remove(str_to_upper(Timepoint), "^T")) + 10,
      TRUE ~ 999
    ),
    inventory_label_rank = case_when(
      timepoint_info == "Diagnosis" ~ 1,
      timepoint_info == "Baseline" ~ 2,
      TRUE ~ 9
    )
  ) %>%
  arrange(
    Patient, mutation_source_type,
    coalesce(inventory_source_date, as.Date("9999-12-31")),
    inventory_timepoint_rank,
    inventory_label_rank,
    source_key
  ) %>%
  group_by(Patient, mutation_source_type) %>%
  mutate(retained_inventory_source = row_number() == 1L) %>%
  ungroup()

baseline_source_inventory_duplicate_audit <- baseline_source_inventory_ranked %>%
  add_count(Patient, mutation_source_type, name = "n_patient_source_inventory_rows") %>%
  filter(n_patient_source_inventory_rows > 1L) %>%
  arrange(Patient, mutation_source_type, desc(retained_inventory_source), source_key)

write_csv(
  baseline_source_inventory_duplicate_audit,
  "Output_tables_2025/clinical_support/baseline_diagnosis_mutation_source_inventory_deduplication_audit.csv"
)

baseline_source_inventory <- baseline_source_inventory_ranked %>%
  filter(retained_inventory_source) %>%
  select(
    -inventory_source_date,
    -inventory_timepoint_rank,
    -inventory_label_rank,
    -retained_inventory_source
  )

write_csv(
  baseline_source_inventory,
  "Output_tables_2025/clinical_support/baseline_diagnosis_mutation_source_by_patient.csv"
)


## Script complete.
