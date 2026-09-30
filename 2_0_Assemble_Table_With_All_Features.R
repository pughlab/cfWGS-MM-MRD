# ==============================================================================
# 2_0_Assemble_Table_With_All_Features.R
#
# Description:
#   Loads MRD assay results (flow-cytometry, Adaptive MRD counts, Rapid Novor
#   proteomic MRD, PET), clinical metadata (M4, SPORE, IMMAGINE cohorts),
#   relapse dates, fragmentomics metrics, WGS-derived tumour fractions &
#   structural features, and mutation count tables. Cleans and harmonises all
#   inputs by Patient + Timepoint (deduplicating, coalescing NAs, fixing sample
#   codes), computes derived MRD flags, and merges everything into one
#   comprehensive table. Writes CSV and RDS versions for downstream analyses.
#
# Unit of analysis:
#   - Final aggregate (`updated9`): one assembled sample record, identified by
#     Patient + Timepoint + Sample_Code + Date after source-specific cleanup.
#   - Baseline high-quality helper: one patient; BM/cfDNA availability flags
#     are collapsed across eligible baseline records.
#
# Inputs:
#   - Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#   - M4_MRD_filtered.txt                           (lab MRD data)
#   - MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/
#       cfWGS_MRDetect_BM_data_updated_Sep.csv
#   - MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/
#       "cfWGS MRDetect Blood data updated Sep with all patients.csv"
#   - Relapse_dates_M4_clean.csv
#   - combined_clinical_data_updated_April2025.csv
#   - Clinical data/SPORE/SPORE_pct_flow_extracted.xlsx
#   - Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx  (sheet 3)
#   - New OICR Submissions/derived_metadata/oicr_submission_clinical_comprehensive_rows.csv
#   - Clinical data/M4/Updated Rapid Novor Data.xlsx and PET_dataframe.xlsx
#   - Results_Fragmentomics/Key_fragmentomics_data_updated2.csv
#   - Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv
#   - Mutation_counts/mutation_counts_table_{BM,Blood}_no_RSID.txt
#   - Exported_data_tables_clinical/ M4 date/lab and relapse-date tables
#   - Optional compatibility fallbacks used only when present:
#       Final_aggregate_table_*_updated6.csv (sample identities/missing rows)
#       Final_aggregate_table_*_updated8.rds (missing Mean.Coverage values)
#   Additional revision clinical-MRD inputs are resolved conditionally in their
#   labelled sections below.
#
# Outputs:
#   - cfWGS clinical MRD values with timepoint and dates updated August 2025.csv
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated.csv/.rds
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated2.csv/.rds
#       intermediate checkpoint after timepoint/date cleanup
#   - Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.csv/.rds
#       final current revision-inclusive aggregate used downstream
#   - baseline_high_quality_patients_updated.csv/.rds
#
# Dependencies:
#   tidyverse, readxl, gridExtra, fuzzyjoin, and helpers.R
#   Run from the project root; source and output paths below are project-relative.
#
# How to run:
#   Rscript Scripts_2025/Final_Scripts/2_0_Assemble_Table_With_All_Features.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script assembles the final aggregate
#     feature/clinical table used by Table 1, longitudinal analyses, model
#     training/application, concordance scripts, and survival analyses.
#
# Author:    Dory Abelman
# Last update: September 2025
# ==============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#

# Workflow outline:
#   1. Harmonize clinical MRD measurements and sample dates across cohorts.
#   2. Attach BM and plasma MRDetect results and fragmentomics measurements.
#   3. Join clinical/FISH variables, WGS features, and mutation-count summaries.
#   4. Resolve duplicate source records while preserving sample identity.
#   5. Recompute relapse fields, create scoring/eligibility audits, and export
#      the final sample-level aggregate plus the patient-level baseline helper.



### MRD comparison figures updated Winter 2025
### Author: Dory Abelman 

library(tidyverse)
library(readxl)
library(gridExtra)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)


### Load source tables
# These inputs are produced by earlier scripts or exported from clinical/MRD
# assay pipelines. This script expects to be run from the project root so that
# all relative paths resolve correctly.
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds") 
M4_MRD_filtered <- read_tsv("M4_MRD_filtered.txt") # The MRD data from labs

path <- file.path("MRDetect_output_winter_2025", "Processed_R_outputs", "BM_muts_data")
#MRD_cfWGS_backup <- MRD_cfWGS
MRD_cfWGS_BM <- read.csv(file.path("MRDetect_output_winter_2025/Processed_R_outputs/BM_muts_plots_baseline/cfWGS_MRDetect_BM_data_updated_Sep.csv")) # MRD data from MRDetect
MRD_cfWGS_blood <- read.csv(file.path("MRDetect_output_winter_2025/Processed_R_outputs/Blood_muts_plots_baseline/cfWGS MRDetect Blood data updated Sep with all patients.csv")) # MRD data from MRDetect


Relapse_dates_M4_clean <- read_csv("Relapse_dates_M4_clean.csv") # Relapse dates
combined_clinical_data_updated <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
) ## Aggregated clinical info
SPORE_MRD <- read_excel("Clinical data/SPORE/SPORE_pct_flow_extracted.xlsx") %>% filter(!is.na(Patient))
IMMAGINE_MRD <- read_excel("Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx", sheet = 3)


### Build one clinical MRD table across M4, SPORE, and IMMAGINE
# M4 clinical MRD is already timepoint-coded. SPORE and IMMAGINE clinical MRD
# arrive by date, so they are first standardized to Patient/Date/Flow values and
# later mapped back to timepoint using the clinical sample metadata.
# Step 1: Clean and transform `% Cells` in table2
table2 <- IMMAGINE_MRD %>%
  mutate(
    # Convert `% Cells` to numeric where "Not detected" is 0
    `% Cells` = ifelse(`% Cells` == "Not detected", 0, as.numeric(`% Cells`))
  )

table2 <- table2 %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

# Step 2: Rename columns for consistency
table1 <- SPORE_MRD %>%
  rename(ID = Patient, Date = `Timepoint of interest`, Flow_Pct_Cells = `Edited_Pct_Flow`)

table2 <- table2 %>%
  rename(Date = Date, Flow_Pct_Cells = `% Cells`)

# Step 3: Merge the tables
# Join by 'ID' and 'Date'
MRD_SPORE_IMMAGINE <- bind_rows(table1, table2) %>% select(ID, Date, Flow_Pct_Cells, `MRD Status`)
MRD_SPORE_IMMAGINE <- MRD_SPORE_IMMAGINE %>%
  rename(MRD_Status = `MRD Status`) %>%
  mutate(
    MRD_Status = ifelse(is.na(Flow_Pct_Cells) | Flow_Pct_Cells > 0, "Positive", "Negative")
  )

M4_MRD_clinical <- M4_MRD_filtered %>%
  rename(MRD_Results_FLOW = MRD_Results)

## Now add to the other table 
MRD_SPORE_IMMAGINE <- MRD_SPORE_IMMAGINE %>% 
  rename(Patient = ID, Flow_pct_cells = Flow_Pct_Cells, MRD_Results_FLOW = MRD_Status)

## Multiple by 100 to be consistent with other table for FLOW 
MRD_SPORE_IMMAGINE$Flow_pct_cells <- MRD_SPORE_IMMAGINE$Flow_pct_cells*100

oicr_revision_clinical_mrd_path <- file.path(
  # Structured OICR clinical rows can include clinical MRD results for revision
  # patients. These are paired to revision plasma/cfDNA sample dates below before
  # joining into the common clinical MRD table.
  "New OICR Submissions",
  "derived_metadata",
  "oicr_submission_clinical_comprehensive_rows.csv"
)
oicr_revision_request_mrd_path <- "FINAL_sample_request_by_sample_send.csv"

oicr_revision_clinical_mrd <- tibble()
if (file.exists(oicr_revision_clinical_mrd_path)) {
  revision_metadata_for_mrd <- load_spring2026_revision_metadata(required = FALSE)

  if (!is.null(revision_metadata_for_mrd)) {
    revision_blood_sample_lookup <- revision_metadata_for_mrd %>%
      # Pair clinical MRD results to revision plasma samples by patient numeric ID
      # and MRD_Blood_Date. Blood/plasma rows only are eligible because this table
      # represents clinical MRD blood draws, not BM/buffy metadata rows.
      filter(
        Sample_type == "Blood_plasma_cfDNA",
        !is.na(patient_numeric_id),
        !is.na(Date_of_sample_collection)
      ) %>%
      transmute(
        patient_numeric_id,
        Patient,
        revision_blood_date = as.Date(Date_of_sample_collection),
        Visit_Number = as.character(Timepoint),
        Sample_Code = str_remove(
          as.character(Sample_ID),
          "(-P|-B|-O|-OZ-DNA|-OZ|_BM_cells|_Blood_plasma_cfDNA|_Blood_Buffy_coat)$"
        ),
        timepoint_info,
        revision_mrd_sample_id = Sample_ID
      ) %>%
      distinct()

    oicr_revision_clinical_mrd_audit <- read_csv(
      oicr_revision_clinical_mrd_path,
      show_col_types = FALSE
    ) %>%
      mutate(
        MRD_Blood_Date = as.Date(MRD_Blood_Date),
        MRD_Test_Date = as.Date(MRD_Test_Date)
      ) %>%
      left_join(
        revision_blood_sample_lookup,
        by = c("Patient_ID" = "patient_numeric_id", "MRD_Blood_Date" = "revision_blood_date")
      ) %>%
      mutate(
        # Keep every clinical MRD row in the audit, but integrate only rows with
        # a known MRD result and an exact revision plasma sample-date match.
        revision_mrd_pairing_status = case_when(
          is.na(MRD_Result) | MRD_Result == "Unknown" ~ "not_used_unknown_mrd_result",
          is.na(MRD_Blood_Date) ~ "not_used_missing_mrd_blood_date",
          is.na(Patient) ~ "not_used_no_matching_revision_plasma_sample",
          TRUE ~ "used_selected_revision_mrd_pair"
        )
      )

    dir.create("Output_tables_2025", showWarnings = FALSE, recursive = TRUE)
    write_csv(
      oicr_revision_clinical_mrd_audit,
      file.path("Output_tables_2025", "spring2026_oicr_clinical_mrd_pairing_audit.csv")
    )

    unmatched_oicr_revision_mrd <- oicr_revision_clinical_mrd_audit %>%
      filter(revision_mrd_pairing_status != "used_selected_revision_mrd_pair")
    if (nrow(unmatched_oicr_revision_mrd) > 0L) {
      warning(
        "Some Spring 2026 OICR clinical MRD rows were not integrated: ",
        paste(
          unique(unmatched_oicr_revision_mrd$revision_mrd_pairing_status),
          collapse = ", "
        )
      )
    }

    oicr_revision_clinical_mrd <- oicr_revision_clinical_mrd_audit %>%
      filter(revision_mrd_pairing_status == "used_selected_revision_mrd_pair") %>%
      transmute(
        # Convert OICR clinical MRD labels into the same columns used by legacy
        # M4/SPORE/IMMAGINE clinical MRD tables.
        Patient,
        Date = MRD_Blood_Date,
        Visit_Number,
        Sample_Code,
        timepoint_info,
        Flow_pct_cells = case_when(
          MRD_Result == "Absent" ~ 0,
          MRD_Result == "Present" ~ as.numeric(MRD_Residual_Pct),
          TRUE ~ NA_real_
        ),
        MRD_Results_FLOW = case_when(
          MRD_Result == "Absent" ~ "Negative",
          MRD_Result == "Present" ~ "Positive",
          TRUE ~ NA_character_
        ),
        revision_mrd_test_date = MRD_Test_Date,
        revision_mrd_blood_date = MRD_Blood_Date,
        revision_mrd_days_from_blood = MRD_Blood_Days_from_MRD,
        revision_mrd_source = "Spring 2026 OICR selected clinical MRD pairing",
        revision_mrd_decision_rationale = Decision_Rationale,
        revision_mrd_sample_id
      ) %>%
      distinct()

    if (file.exists(oicr_revision_request_mrd_path)) {
      oicr_revision_request_mrd_audit <- read_csv(
        oicr_revision_request_mrd_path,
        show_col_types = FALSE
      ) %>%
        mutate(
          Sample_Date = as.Date(Sample_Date),
          request_row_id = row_number()
        ) %>%
        filter(Request_Type %in% c("MRD Blood", "MRD-timed Blood")) %>%
        left_join(
          revision_blood_sample_lookup,
          by = c("Patient_ID" = "patient_numeric_id", "Sample_Date" = "revision_blood_date")
        ) %>%
        mutate(
          revision_request_mrd_pairing_status = case_when(
            is.na(MRD_Result) | MRD_Result == "Unknown" ~ "not_used_unknown_mrd_result",
            is.na(Sample_Date) ~ "not_used_missing_request_sample_date",
            is.na(Patient) ~ "not_used_no_exact_revision_plasma_sample_date_match",
            TRUE ~ "used_exact_request_mrd_pair"
          )
        )

      write_csv(
        oicr_revision_request_mrd_audit,
        file.path("Output_tables_2025", "spring2026_oicr_request_mrd_pairing_audit.csv")
      )

      oicr_revision_request_mrd <- oicr_revision_request_mrd_audit %>%
        filter(revision_request_mrd_pairing_status == "used_exact_request_mrd_pair") %>%
        transmute(
          Patient,
          Date = Sample_Date,
          Visit_Number,
          Sample_Code,
          timepoint_info,
          Flow_pct_cells = case_when(
            MRD_Result == "Absent" ~ 0,
            MRD_Result == "Present" ~ suppressWarnings(as.numeric(str_extract(Clinical_Context, "[0-9.]+(?=%\\))"))),
            TRUE ~ NA_real_
          ),
          MRD_Results_FLOW = case_when(
            MRD_Result == "Absent" ~ "Negative",
            MRD_Result == "Present" ~ "Positive",
            TRUE ~ NA_character_
          ),
          revision_mrd_test_date = as.Date(Sample_Date - Days_from_MRD),
          revision_mrd_blood_date = Sample_Date,
          revision_mrd_days_from_blood = Days_from_MRD,
          revision_mrd_source = "Spring 2026 OICR per-sample MRD request table",
          revision_mrd_decision_rationale = Clinical_Context,
          revision_mrd_sample_id
        ) %>%
        anti_join(
          oicr_revision_clinical_mrd %>% distinct(Patient, Sample_Code),
          by = c("Patient", "Sample_Code")
        ) %>%
        distinct()

      oicr_revision_clinical_mrd <- bind_rows(
        oicr_revision_clinical_mrd,
        oicr_revision_request_mrd
      )
    }
  }
}

M4_MRD_clinical <- M4_MRD_clinical %>% 
  rename(Patient = M4_id)

cfWGS_Clinical_MRD <- bind_rows(M4_MRD_clinical, MRD_SPORE_IMMAGINE, oicr_revision_clinical_mrd)

## Calculate preliminary composite clinical MRD fields
# These preliminary values are recalculated after the remaining Rapid Novor and
# PET records are added. The later calculation is the one retained downstream.
# These are legacy multi-assay composite helpers, not the prespecified
# isotype-specific EasyM call used in the manuscript. Script 3_1_A constructs
# that EasyM comparator using 1% for IgG and 0.05% for IgA/light-chain disease.
# The historical rule labels every row not meeting a positive condition as
# "Negative", including rows for which all component assays are missing.
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    MRD_by_clinical_testing = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 1
      Rapid_Novor > 1 ~ "Positive",
      
      # Historical fallback: all other rows, including all-missing rows, are Negative
      TRUE ~ "Negative"
    )
  )

# Create MRD_by_clinical_testing_stringent column based on the specified conditions
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    MRD_by_clinical_testing_stringent = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Adaptive MRD is positive at either recorded sensitivity
      Mrd1E6 == "Positive" ~ "Positive",
      
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 0.1
      Rapid_Novor > 0.1 ~ "Positive",
      
      # Historical fallback: all other rows, including all-missing rows, are Negative
      TRUE ~ "Negative"
    )
  )


### Add dates to M4 for use in the patient tracker 
M4_dates <- read_csv("M4_dates_Oct2024.csv") 
tmp <- M4_dates %>% 
  select(Patient, Timepoint_Code, Date) %>% 
  unique() %>% 
  rename(Visit_Number = Timepoint_Code) %>% 
  mutate(Visit_Number = gsub("R-", "R", Visit_Number))

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>% 
  left_join(tmp %>% select(Patient, Visit_Number, Date), by = c("Patient", "Visit_Number")) %>%
  mutate(Date = coalesce(Date.x, Date.y)) %>%
  select(-Date.x, -Date.y)


## Add the IDs if there is a matched sample 
tmp2 <- combined_clinical_data_updated %>% 
  #  filter(Study != "M4") %>%
  select(Patient, Date_of_sample_collection, Timepoint, Sample_ID, timepoint_info) %>% 
  unique() %>%
  mutate(
    Sample_Code = gsub("(-P|-B|-O|-OZ-DNA|-OZ|_BM_cells|_Blood_plasma_cfDNA|_Blood_Buffy_coat)$", "", Sample_ID),  # Remove specified suffixes
    Sample_Code = gsub("MyP|MyC|MyR", "IMG", Sample_Code)  # Remove specified suffixes
  ) %>%
  select(-Sample_ID) %>%  # Remove the original Sample_ID column
  unique()  # Ensure rows are unique

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>% select(-Patient_Code, -Site_Code) %>% 
  rename(Timepoint = Visit_Number)

## Merge 
# Ensure dates are in Date format
cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(Date = as.Date(Date))

tmp2 <- tmp2 %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

## Add the current central IMMAGINE/LIBERATE clinical marrow-MRD export.
##
## The historical IMMAGINE workbook above is intentionally retained because it
## contains curated assay calls used by the original analysis.  The December
## 2025 and February 2026 central exports contain additional marrow assessments
## that were not present on that legacy sheet.  We map these records by patient
## and collection date using the project-wide <=42-day clinical MRD pairing
## window. Existing integrated flow results always
## take precedence, so this block only fills genuinely missing MFC values.
central_img_mrd_2025_path <- file.path(
  "Clinical data", "IMMAGINE",
  "IMMAGINE  LIBERATE_MRD_withoutMRN_or_DOB_18Dec2025.xlsx"
)
central_img_mrd_2026_path <- file.path(
  "Clinical data", "Additional sample sheets from Esther",
  "MRD_LIBERATE_ID_19Feb2026.xlsx"
)

central_img_mrd_matched <- tibble()
if (file.exists(central_img_mrd_2025_path)) {
  central_img_patient_map <- read_excel(
    central_img_mrd_2025_path,
    sheet = "IMMAGINE & LIBERATE SAMPLES"
  ) %>%
    transmute(
      PATIENT_ID = as.character(PATIENT_ID),
      central_mrd_test_date = as.Date(COLLECTION_DATE),
      PROCEDURE_TYPE_TEXT = as.character(PROCEDURE_TYPE_TEXT),
      Patient = if_else(
        str_detect(as.character(STUDY_NAME), "^IMMAGINE\\s+[0-9]+$"),
        sprintf(
          "IMG-%03d",
          suppressWarnings(as.integer(str_extract(as.character(STUDY_NAME), "[0-9]+$")))
        ),
        NA_character_
      )
    ) %>%
    filter(
      !is.na(PATIENT_ID),
      !is.na(Patient),
      PROCEDURE_TYPE_TEXT == "Bone marrow sample"
    ) %>%
    select(-PROCEDURE_TYPE_TEXT) %>%
    distinct(PATIENT_ID, central_mrd_test_date, .keep_all = TRUE)

  central_img_mrd_2025 <- read_excel(
    central_img_mrd_2025_path,
    sheet = "MRD"
  ) %>%
    transmute(
      Patient = if_else(
        str_detect(as.character(STUDY_NAME), "^IMMAGINE\\s+[0-9]+$"),
        sprintf(
          "IMG-%03d",
          suppressWarnings(as.integer(str_extract(as.character(STUDY_NAME), "[0-9]+$")))
        ),
        NA_character_
      ),
      central_mrd_test_date = as.Date(BM_COLLCTN_DATE),
      central_mrd_status = case_when(
        MALIGNANT_PC == "Absent" ~ "Negative",
        MALIGNANT_PC == "Present" ~ "Positive",
        TRUE ~ NA_character_
      ),
      central_flow_pct_cells = case_when(
        MALIGNANT_PC == "Absent" ~ 0,
        MALIGNANT_PC == "Present" ~ suppressWarnings(as.numeric(`RESIDUAL_PC%`)),
        TRUE ~ NA_real_
      ),
      central_mrd_source = "Central IMMAGINE/LIBERATE clinical MRD export (2025-12-18)",
      central_mrd_source_priority = 1L
    )

  central_img_mrd_2026 <- tibble()
  if (file.exists(central_img_mrd_2026_path)) {
    central_img_mrd_2026 <- read_excel(central_img_mrd_2026_path) %>%
      mutate(
        PATIENT_ID = as.character(PATIENT_ID),
        central_mrd_test_date = as.Date(BONE_MARRW_COLLCTN_DATE)
      ) %>%
      left_join(
        central_img_patient_map,
        by = c("PATIENT_ID", "central_mrd_test_date"),
        relationship = "many-to-one"
      ) %>%
      transmute(
        Patient,
        central_mrd_test_date,
        central_mrd_status = case_when(
          MALIGNANT_PLASMA_CELLS == "Absent" ~ "Negative",
          MALIGNANT_PLASMA_CELLS == "Present" ~ "Positive",
          TRUE ~ NA_character_
        ),
        central_flow_pct_cells = case_when(
          MALIGNANT_PLASMA_CELLS == "Absent" ~ 0,
          MALIGNANT_PLASMA_CELLS == "Present" ~
            suppressWarnings(as.numeric(RESID_MYELOMA_CELL_PERCENT)),
          TRUE ~ NA_real_
        ),
        central_mrd_source = "Central IMMAGINE/LIBERATE clinical MRD export (2026-02-19)",
        central_mrd_source_priority = 2L
      )
  }

  central_img_mrd_raw <- bind_rows(central_img_mrd_2025, central_img_mrd_2026) %>%
    filter(
      !is.na(Patient),
      !is.na(central_mrd_test_date),
      !is.na(central_mrd_status)
    )

  central_img_internal_conflicts <- central_img_mrd_raw %>%
    group_by(Patient, central_mrd_test_date) %>%
    summarise(n_distinct_results = n_distinct(central_mrd_status), .groups = "drop") %>%
    filter(n_distinct_results > 1L)
  if (nrow(central_img_internal_conflicts) > 0L) {
    stop(
      "Central IMMAGINE/LIBERATE MRD exports contain conflicting results for ",
      nrow(central_img_internal_conflicts),
      " patient-date combinations; review before integration.",
      call. = FALSE
    )
  }

  # When both releases report the same patient/date/status, retain the newer
  # quantitative value while preserving both release labels in provenance.
  central_img_mrd <- central_img_mrd_raw %>%
    arrange(Patient, central_mrd_test_date, desc(central_mrd_source_priority)) %>%
    group_by(Patient, central_mrd_test_date, central_mrd_status) %>%
    summarise(
      central_flow_pct_cells = {
        observed_values <- central_flow_pct_cells[!is.na(central_flow_pct_cells)]
        if (length(observed_values) > 0L) observed_values[[1L]] else NA_real_
      },
      central_mrd_source = paste(unique(central_mrd_source), collapse = "; "),
      .groups = "drop"
    )

  central_img_mrd_candidates <- tmp2 %>%
    select(
      Patient,
      Date_of_sample_collection,
      Timepoint,
      Sample_Code,
      timepoint_info
    ) %>%
    distinct() %>%
    inner_join(central_img_mrd, by = "Patient") %>%
    mutate(
      central_mrd_days_from_sample = as.integer(abs(
        central_mrd_test_date - Date_of_sample_collection
      ))
    ) %>%
    filter(central_mrd_days_from_sample <= 42L)

  central_img_match_conflicts <- central_img_mrd_candidates %>%
    group_by(Sample_Code) %>%
    summarise(
      n_statuses = n_distinct(central_mrd_status),
      n_test_dates = n_distinct(central_mrd_test_date),
      .groups = "drop"
    ) %>%
    filter(n_statuses > 1L)

  if (nrow(central_img_match_conflicts) > 0L) {
    dir.create("Output_tables_2025", showWarnings = FALSE, recursive = TRUE)
    central_img_mrd_candidates %>%
      semi_join(central_img_match_conflicts, by = "Sample_Code") %>%
      arrange(Sample_Code, central_mrd_test_date) %>%
      write_csv(file.path(
        "Output_tables_2025",
        "central_imagine_liberate_mrd_42d_conflicts_for_review.csv"
      ))
  }

  central_img_mrd_matched <- central_img_mrd_candidates %>%
    anti_join(central_img_match_conflicts, by = "Sample_Code") %>%
    group_by(Sample_Code) %>%
    arrange(central_mrd_days_from_sample, desc(central_mrd_test_date)) %>%
    slice(1L) %>%
    ungroup() %>%
    transmute(
      Sample_Code,
      Patient.central = Patient,
      Date.central = Date_of_sample_collection,
      Timepoint.central = as.character(Timepoint),
      timepoint_info.central = timepoint_info,
      Flow_pct_cells.central = central_flow_pct_cells,
      MRD_Results_FLOW.central = central_mrd_status,
      revision_mrd_test_date.central = central_mrd_test_date,
      revision_mrd_blood_date.central = Date_of_sample_collection,
      revision_mrd_days_from_blood.central = central_mrd_days_from_sample,
      revision_mrd_source.central = central_mrd_source,
      revision_mrd_decision_rationale.central = paste0(
        "Nearest non-conflicting central clinical MRD match within 42-day window (",
        central_mrd_days_from_sample,
        " days from cfWGS sample)."
      ),
      revision_mrd_sample_id.central = Sample_Code
    )
}

# Join cfWGS_Clinical_MRD with tmp2 by Patient using the project-wide 42-day window
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD %>% filter(is.na(Timepoint)) %>%
  left_join(tmp2, by = "Patient") %>%
  # Calculate the absolute difference in days
  mutate(date_diff = abs(difftime(Date, Date_of_sample_collection, units = "days"))) %>%
  # Filter to keep matches within the prespecified 42-day clinical MRD window
  filter(is.na(Timepoint.x) & date_diff <= 42 | !is.na(Timepoint.x)) %>%
  # For each row where Timepoint.x or Sample_Code.x is NA, keep the closest date match
  group_by(Patient, Date) %>%
  mutate(min_diff = min(date_diff, na.rm = TRUE)) %>%
  filter(date_diff == min_diff | !is.na(Timepoint.x)) %>%
  ungroup() %>%
  # Fill in missing Sample_Code and Timepoint from tmp2 where needed
  mutate(
    Sample_Code = ifelse(is.na(Sample_Code.x), Sample_Code.y, Sample_Code.x),
    Timepoint = ifelse(is.na(Timepoint.x), Timepoint.y, Timepoint.x)
  ) %>%
  # Select final columns and remove temporary columns
  select(-Sample_Code.x, -Sample_Code.y, -Timepoint.x, -Timepoint.y, -Date_of_sample_collection, -date_diff, -min_diff) %>% 
  unique()

## Bind rows to the other 

# Make sure both data frames have Sample_Code (and Timepoint) as character
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    Sample_Code = as.character(Sample_Code),
    Timepoint   = as.character(Timepoint)
  )

cfWGS_Clinical_MRD <- cfWGS_Clinical_MRD %>%
  mutate(
    Sample_Code = as.character(Sample_Code),
    Timepoint   = as.character(Timepoint)
  )
cfWGS_Clinical_MRD_filled <- bind_rows(cfWGS_Clinical_MRD_filled, cfWGS_Clinical_MRD %>% filter(!is.na(Timepoint)))

## Add timepoint info if NA
# 1) make a clean lookup table
tmp_info <- tmp2 %>%
  ungroup() %>% 
  select(Sample_Code, timepoint_info) %>%
  distinct()

# 2) join + coalesce
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(tmp_info, by = "Sample_Code", suffix = c("", ".new")) %>%
  mutate(
    timepoint_info = coalesce(timepoint_info, timepoint_info.new)
  ) %>%
  select(-timepoint_info.new)

#cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% select(-Sample_Code)

## Add timepoint values 
# Create the mapping of timepoints to descriptions
description_mapping <- c(
  "01" = "Diagnosis",
  "03" = "Post_induction",
  "05" = "Post_transplant",
  "07" = "1yr maintenance",
  "08" = "1.5yr maintenance",
  "09" = "2yr maintenance",
  "10" = "2.5yr maintenance",
  "11" = "3yr maintenance",
  "12" = "3.5yr maintenance",
  "13" = "4yr maintenance",
  "14" = "4.5yr maintenance",
  "15" = "5yr maintenance",
  "R"  = "Relapse"
)

# Add timepoint_info based on the mapping
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    timepoint_info = if_else(
      !grepl("^(SPORE|IMG)", Patient),  # Check if Patient does not start with SPORE or IMG
      description_mapping[Timepoint],  # Apply the mapping
      timepoint_info  # Keep the existing value otherwise
    )
  )

## Fill missing flow-MRD values from the central export without overwriting any
## result already integrated from M4, SPORE, legacy IMMAGINE, or the selected
## Spring 2026 OICR pairings.
if (nrow(central_img_mrd_matched) > 0L) {
  cross_source_flow_candidates <- bind_rows(
    cfWGS_Clinical_MRD_filled %>%
      filter(!is.na(Sample_Code), !is.na(MRD_Results_FLOW)) %>%
      transmute(
        Sample_Code,
        flow_status = MRD_Results_FLOW,
        test_date = as.Date(Date),
        flow_source = coalesce(
          revision_mrd_source,
          "Historical integrated clinical MRD workbook"
        )
      ),
    central_img_mrd_candidates %>%
      transmute(
        Sample_Code,
        flow_status = central_mrd_status,
        test_date = central_mrd_test_date,
        flow_source = central_mrd_source
      )
  ) %>%
    distinct()

  cross_source_flow_conflicts <- cross_source_flow_candidates %>%
    group_by(Sample_Code) %>%
    summarise(n_statuses = n_distinct(flow_status), .groups = "drop") %>%
    filter(n_statuses > 1L)

  if (nrow(cross_source_flow_conflicts) > 0L) {
    cross_source_flow_candidates %>%
      semi_join(cross_source_flow_conflicts, by = "Sample_Code") %>%
      arrange(Sample_Code, test_date, flow_status) %>%
      write_csv(file.path(
        "Output_tables_2025",
        "clinical_mrd_42d_cross_source_conflicts_for_review.csv"
      ))
  }

  existing_flow_keys <- cfWGS_Clinical_MRD_filled %>%
    filter(!is.na(MRD_Results_FLOW) | !is.na(Flow_pct_cells)) %>%
    distinct(Sample_Code) %>%
    mutate(existing_flow_result = TRUE)

  central_img_mrd_integration_audit <- central_img_mrd_matched %>%
    left_join(existing_flow_keys, by = "Sample_Code") %>%
    mutate(
      integration_status = if_else(
        coalesce(existing_flow_result, FALSE),
        "not_used_existing_integrated_flow_result",
        "used_to_fill_missing_flow_result"
      )
    ) %>%
    select(-existing_flow_result)

  dir.create("Output_tables_2025", showWarnings = FALSE, recursive = TRUE)
  write_csv(
    central_img_mrd_integration_audit,
    file.path(
      "Output_tables_2025",
      "central_imagine_liberate_mrd_42d_integration_audit.csv"
    )
  )

  cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
    full_join(central_img_mrd_matched, by = "Sample_Code") %>%
    mutate(
      central_flow_used = is.na(MRD_Results_FLOW) &
        !is.na(MRD_Results_FLOW.central),
      Patient = coalesce(Patient, Patient.central),
      Date = coalesce(as.Date(Date), Date.central),
      Timepoint = coalesce(as.character(Timepoint), Timepoint.central),
      timepoint_info = coalesce(timepoint_info, timepoint_info.central),
      Flow_pct_cells = coalesce(Flow_pct_cells, Flow_pct_cells.central),
      MRD_Results_FLOW = coalesce(
        MRD_Results_FLOW,
        MRD_Results_FLOW.central
      ),
      revision_mrd_test_date = if_else(
        central_flow_used,
        revision_mrd_test_date.central,
        revision_mrd_test_date
      ),
      revision_mrd_blood_date = if_else(
        central_flow_used,
        revision_mrd_blood_date.central,
        revision_mrd_blood_date
      ),
      revision_mrd_days_from_blood = if_else(
        central_flow_used,
        revision_mrd_days_from_blood.central,
        revision_mrd_days_from_blood
      ),
      revision_mrd_source = if_else(
        central_flow_used,
        revision_mrd_source.central,
        revision_mrd_source
      ),
      revision_mrd_decision_rationale = if_else(
        central_flow_used,
        revision_mrd_decision_rationale.central,
        revision_mrd_decision_rationale
      ),
      revision_mrd_sample_id = if_else(
        central_flow_used,
        revision_mrd_sample_id.central,
        revision_mrd_sample_id
      ),
      unresolved_42d_flow_conflict = Sample_Code %in%
        cross_source_flow_conflicts$Sample_Code &
        !str_detect(
          coalesce(revision_mrd_source, ""),
          fixed("Spring 2026 OICR per-sample MRD request table")
        ),
      Flow_pct_cells = if_else(
        unresolved_42d_flow_conflict,
        NA_real_,
        Flow_pct_cells
      ),
      MRD_Results_FLOW = if_else(
        unresolved_42d_flow_conflict,
        NA_character_,
        MRD_Results_FLOW
      ),
      revision_mrd_decision_rationale = if_else(
        unresolved_42d_flow_conflict,
        "Not integrated: discordant clinical MFC results fall within the 42-day pairing window.",
        revision_mrd_decision_rationale
      )
    ) %>%
    select(
      -central_flow_used,
      -unresolved_42d_flow_conflict,
      -ends_with(".central")
    )
}


### Add updated Rapid Novor proteomic MRD values
# These values supersede the earlier Rapid_Novor column when both are present.
# Negative, Not Quantifiable, and NQ are treated as zero before reshaping.
Rapid_Novor_update <- read_excel("Clinical data/M4/Updated Rapid Novor Data.xlsx")

# Restructure the data
# Step 1: Standardize and Clean Non-Numeric Values
Rapid_Novor_cleaned <- Rapid_Novor_update %>%
  mutate(across(
    -`UHN Serial`, 
    ~ case_when(
      . %in% c("Negative", "Not Quantifiable", "NQ") ~ 0,  # Set these to 0
      TRUE ~ suppressWarnings(as.numeric(.))  # Convert others to numeric
    )
  )) %>%
  rename(M4_id = `UHN Serial`) # Rename UHN Serial to M4_id

# Step 2: Restructure the Data
Rapid_Novor_restructured <- Rapid_Novor_cleaned %>%
  pivot_longer(
    cols = -M4_id, # All columns except M4_id
    names_to = "Timepoint", # New column for Timepoint
    values_to = "Value" # New column for the values
  ) 

# Step 1: Remove rows where Value is NA
Rapid_Novor_restructured <- Rapid_Novor_restructured %>%
  filter(!is.na(Value))

# Step 2: Consolidate 9 and 9P into 9, and R and RP into R
Rapid_Novor_consolidated <- Rapid_Novor_restructured %>%
  mutate(Timepoint = case_when(
    Timepoint %in% c("9", "9P") ~ "9",  # Consolidate 9 and 9P into 9
    Timepoint %in% c("R", "RP") ~ "R",  # Consolidate R and RP into R
    TRUE ~ Timepoint  # Keep other timepoints as they are
  )) %>%
  group_by(M4_id, Timepoint) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")  # Average if duplicates exist

## Prep for join 
Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Patient = M4_id,  # Rename M4_id to Patient
  ) %>%
  select(-M4_id)  # Optionally remove the original M4_id column if not needed


Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Patient = case_when(
      Patient %in% c("VA-018", "VA-017", "VA-016") ~ sub("^VA-0", "VA-", Patient),  # Remove leading 0
      TRUE ~ Patient  # Keep all other values unchanged
    )
  )

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    Timepoint = if_else(Timepoint == "R", Timepoint, str_pad(Timepoint, width = 2, pad = "0")),  # Add leading 0 only if not "R"
    Sample_Code = paste0(Patient, "-", Timepoint)  # Create Sample_Code as Patient-Timepoint
  )

## Correct errors
# Create the mapping of timepoints to descriptions
description_mapping <- c(
  "01" = "Diagnosis",
  "03" = "Post_induction",
  "05" = "Post_transplant",
  "07" = "1yr maintenance",
  "08" = "1.5yr maintenance",
  "09" = "2yr maintenance",
  "10" = "2.5yr maintenance",
  "11" = "3yr maintenance",
  "12" = "3.5yr maintenance",
  "13" = "4yr maintenance",
  "14" = "4.5yr maintenance",
  "15" = "5yr maintenance",
  "R"  = "Relapse"
)

# Add timepoint_info based on the mapping
Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  mutate(
    timepoint_info = description_mapping[Timepoint]  # Map Timepoint to description
  ) 


# Check for unmapped timepoints
unmapped_timepoints <- Rapid_Novor_consolidated %>%
  filter(is.na(timepoint_info)) %>%
  pull(Timepoint) %>%
  unique()

if (length(unmapped_timepoints) > 0) {
  warning("Unmapped timepoints: ", paste(unmapped_timepoints, collapse = ", "))
}

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>%
  rename(Rapid_Novor = Value)

Rapid_Novor_consolidated <- Rapid_Novor_consolidated %>% unique()

# Perform an inner join to compare Rapid_Novor values for matching Patient and Timepoint
# Set a tolerance threshold (e.g., 0.01 for small differences)
tolerance <- 1

# Perform an inner join to compare Rapid_Novor values for matching Patient and Timepoint
# Some samples are really different
discrepancies <- Rapid_Novor_consolidated %>%
  inner_join(
    cfWGS_Clinical_MRD_filled %>% select(Patient, Timepoint, Rapid_Novor),
    by = c("Patient", "Timepoint"),
    suffix = c("_consolidated", "_clinical")
  ) %>%
  filter(!is.na(Rapid_Novor_consolidated), !is.na(Rapid_Novor_clinical)) %>% # Ensure both values are non-missing
  mutate(
    difference = abs(Rapid_Novor_consolidated - Rapid_Novor_clinical), # Calculate absolute difference
    within_tolerance = difference <= tolerance # Check if difference is within tolerance
  )


## Replace earlier Rapid Novor values with the curated updated table.
# Step 1: Remove the Rapid_Novor and Rapid_Novor_Binary columns from cfWGS_Clinical_MRD_filled
tmp <- cfWGS_Clinical_MRD_filled %>%
  select(-Rapid_Novor)

# Step 2: Check for duplicates in tmp
dups_by_pt <- tmp %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()


# Step 2: Full join with Rapid_Novor_consolidated to ensure all rows are included
# A full_join (rather than left_join) is used here because Rapid Novor proteomic
# MRD may have been run on some patients/timepoints that are not yet in the main
# cfWGS_Clinical_MRD_filled table (e.g. later-enrolled IMMAGINE patients). A
# left_join would silently drop those rows. Subsequent unique() call removes any
# row-level duplicates introduced by the join.
cfWGS_Clinical_MRD_filled <- tmp %>%
  full_join(
    Rapid_Novor_consolidated %>% 
      select(Patient, Timepoint, Rapid_Novor, Sample_Code, timepoint_info)
  ) %>% 
  unique()

# Identify duplicate rows based on Sample_Code
duplicate_rows <- cfWGS_Clinical_MRD_filled %>%
  group_by(Sample_Code) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



### Add PET MRD status
# PET is included as an additional binary clinical disease-status assay.
PET_M4 <- read_excel("Clinical data/M4/PET_dataframe.xlsx")

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% 
  full_join(PET_M4 %>% select(PET, Patient_Timepoint), by = c("Sample_Code" = "Patient_Timepoint"))


## Clean up 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_Results_FLOW = case_when(
      Flow_pct_cells == 0 ~ "Negative",
      Flow_pct_cells > 0 & MRD_Results_FLOW == "-" ~ "Positive",
      TRUE ~ MRD_Results_FLOW  # Keep the original value if neither condition applies
    )
  )

## Recalculate final composite clinical MRD fields
# These values replace the preliminary fields above and are retained in the
# final aggregate as legacy multi-assay helper fields. They are not used as the
# manuscript EasyM comparator; the isotype-specific EasyM call is generated in
# 3_1_A_Process_and_optimize_EasyM.R. The historical fallback again labels all
# rows not satisfying a positive condition as "Negative", including
# all-component-missing rows.
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_by_clinical_testing = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Mrd1E5 is Positive
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 1
      Rapid_Novor > 1 ~ "Positive",
      
      # Historical fallback: all other rows, including all-missing rows, are Negative
      TRUE ~ "Negative"
    )
  )

# The retained legacy composite uses Rapid Novor >0.05 here. This differs from
# the >0.1 preliminary block above and does not apply the manuscript's
# isotype-specific EasyM rule. It remains unchanged for backward compatibility.
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    MRD_by_clinical_testing_stringent = case_when(
      # Condition 1: MRD_Results_FLOW is Positive
      MRD_Results_FLOW == "Positive" ~ "Positive",
      
      # Condition 2: Adaptive MRD is positive at either recorded sensitivity
      Mrd1E6 == "Positive" ~ "Positive",
      
      Mrd1E5 == "Positive" ~ "Positive",
      
      # Condition 3: Rapid_Novor is greater than 0.05 (implemented legacy rule)
      Rapid_Novor > 0.05 ~ "Positive",
      
      # Historical fallback: all other rows, including all-missing rows, are Negative
      TRUE ~ "Negative"
    )
  )

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()
### Now assess correlates

# Prepare binary columns for concordance analysis with is.na checks
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    # Binary indicator for Flow: 1 if Positive, 0 if Negative, NA if missing
    Flow_Binary = case_when(
      !is.na(MRD_Results_FLOW) & MRD_Results_FLOW == "Positive" ~ 1,
      !is.na(MRD_Results_FLOW) & MRD_Results_FLOW == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for Adaptive: 1 if Mrd1E6 or Mrd1E5 is Positive, 0 if Negative, NA if missing
    # Adaptive/clonoSEQ reports sensitivity at two thresholds:
    #   Mrd1E6 = 10^-6 sensitivity (deepest; preferred for MRD-negativity calls)
    #   Mrd1E5 = 10^-5 sensitivity (less deep; used as fallback when 1E6 is NA)
    # Positive at EITHER level = MRD+.  "Indeterminate" at 1E6 + Negative at 1E5
    # is treated as MRD- (per IMWG guidance on indeterminate clonoSEQ reads).
    Adaptive_Binary = case_when(
      (!is.na(Mrd1E6) & Mrd1E6 == "Positive") | (!is.na(Mrd1E5) & Mrd1E5 == "Positive") ~ 1,
      (!is.na(Mrd1E6) & Mrd1E6 == "Indeterminate") & (!is.na(Mrd1E5) & Mrd1E5 == "Negative") ~ 0,
      (!is.na(Mrd1E6) & Mrd1E6 == "Negative") | (!is.na(Mrd1E5) & Mrd1E5 == "Negative") ~ 0,
      (!is.na(Mrd1E6) & Mrd1E6 == "Positive") ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for Rapid Novor if positive (> 0), 0 if not, NA if missing
    Rapid_Novor_Binary = case_when(
      !is.na(Rapid_Novor) & Rapid_Novor > 0 ~ 1,
      !is.na(Rapid_Novor) & Rapid_Novor <= 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for PET: 1 if Positive, 0 if Negative, NA if missing
    PET_Binary = case_when(
      !is.na(PET) & PET == "Positive" ~ 1,
      !is.na(PET) & PET == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for MRD by clinical testing (standard): 1 if Positive, 0 if Negative, NA if missing
    MRD_Clinical_Binary = case_when(
      !is.na(MRD_by_clinical_testing) & MRD_by_clinical_testing == "Positive" ~ 1,
      !is.na(MRD_by_clinical_testing) & MRD_by_clinical_testing == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Binary indicator for MRD by clinical testing (stringent): 1 if Positive, 0 if Negative, NA if missing
    MRD_Clinical_Stringent_Binary = case_when(
      !is.na(MRD_by_clinical_testing_stringent) & MRD_by_clinical_testing_stringent == "Positive" ~ 1,
      !is.na(MRD_by_clinical_testing_stringent) & MRD_by_clinical_testing_stringent == "Negative" ~ 0,
      TRUE ~ NA_real_
    )
  )

### Add cfWGS MRDetect MRD calls
# The BM and blood MRDetect files are kept separate so downstream scripts can
# compare matrix-specific cfWGS MRD performance.
# --- 1. Make a “BM” lookup table -----------------------------
tmp_bm <- MRD_cfWGS_BM %>% filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>% filter(Sample_ID != "SPORE_0009_T3_BM_cells") %>% #have better baseline
  # the columns needed for this analysis)
  select(
    Patient,
    Date  = Date_of_sample_collection_Sample_ID_Bam,
    Timepoint = Timepoint_Sample_ID_Bam,
    Mrd_by_WGS_BM    = Mrd_by_WGS,
    Cumulative_VAF_BM        = Cumulative_VAF,
    detect_rate_BM   = detection_rate_as_reads_detected_over_reads_checked,
    sites_rate_BM = sites_detection_rate,
    zscore_BM        = sites_rate_zscore_charm,
    z_score_detection_rate_BM = detection_rate_zscore_reads_checked_charm, 
    PercentChangeFromBaseline_BM       = percent_change, 
    PercentChangeAtSecondTimepoint_BM       = percent_change_detection_rate_second_timepoint
  ) %>%
  distinct()                # just in case you have duplicates

# --- 2. Make a “Blood” lookup table -------------------------
tmp_blood <- MRD_cfWGS_blood %>% filter(timepoint_info %in% c("Baseline", "Diagnosis")) %>%
  # keep only the same columns you want, but rename for clarity
  select(
    Patient,
    Date  = Date_of_sample_collection_Sample_ID_Bam,
    Timepoint = Timepoint_Sample_ID_Bam,
    Mrd_by_WGS_blood = Mrd_by_WGS,
    Cumulative_VAF_blood     = Cumulative_VAF,
    detect_rate_blood= detection_rate_as_reads_detected_over_reads_checked,
    sites_rate_blood = sites_detection_rate,
    zscore_blood     = sites_rate_zscore_charm,
    z_score_detection_rate_blood = detection_rate_zscore_reads_checked_charm, 
    PercentChangeFromBaseline_blood    = percent_change,
    PercentChangeAtSecondTimepoint_blood    = percent_change_detection_rate_second_timepoint
  ) %>%
  distinct()

#  Parse the Date columns in the look-up tables
tmp_bm <- tmp_bm %>%
  mutate(Date = as.Date(Date)) 

tmp_blood <- tmp_blood %>%
  mutate(Date = as.Date(Date))

### check for duplicates 
# 1. Find any Patient/Date/Timepoint combos that appear >1× in tmp_bm or tmp_blood
bm_dups <- tmp_bm %>% 
  count(Patient, Date, Timepoint) %>% 
  filter(n > 1)

blood_dups <- tmp_blood %>% 
  count(Patient, Date, Timepoint) %>% 
  filter(n > 1)

bm_dups
blood_dups


# --- 3. Join both into the clinical table ------------------
# Ensure all Date columns are of type Date
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date))

tmp_bm <- tmp_bm %>%
  mutate(Date = as.Date(Date))

tmp_blood <- tmp_blood %>%
  mutate(Date = as.Date(Date))
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # first bring in the BM calls
  full_join(tmp_bm,    by = c("Patient","Date","Timepoint")) %>%
  # then bring in the blood calls
  full_join(tmp_blood, by = c("Patient","Date","Timepoint"))

# --- 4. Create two binary MRD flags -------------------------

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    cfWGS_BM_Binary = case_when(
      Mrd_by_WGS_BM    == "Positive" ~ 1,
      Mrd_by_WGS_BM    == "Negative" ~ 0,
      TRUE                          ~ NA_real_
    ),
    cfWGS_blood_Binary = case_when(
      Mrd_by_WGS_blood == "Positive" ~ 1,
      Mrd_by_WGS_blood == "Negative" ~ 0,
      TRUE                          ~ NA_real_
    )
  )

# --- 5. Quick sanity-check -------------------------

# How many rows got both BM & blood info?
cfWGS_Clinical_MRD_filled %>% 
  summarize(
    total = dplyr::n(),
    with_BM    = sum(!is.na(Mrd_by_WGS_BM)),
    with_blood = sum(!is.na(Mrd_by_WGS_blood))
  )


#### Reorganize
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  select(
    # 1) Core identifiers
    Patient,
    Date,
    timepoint_info,
    Sample_Code,
    Timepoint,
    
    # 2) Clinical MRD assays
    Flow_pct_cells,
    Flow_Binary,
    MRD_Results_FLOW,
    
    PerMillionCount_Adaptive,
    Adaptive_Frequency,
    Adaptive_Binary,
    
    Rapid_Novor,
    Rapid_Novor_Binary,
    
    PET,
    PET_Binary,
    
    MRD_by_clinical_testing,
    MRD_by_clinical_testing_stringent,
    MRD_Clinical_Binary,
    MRD_Clinical_Stringent_Binary,
    
    # 3) WGS‐based MRD, BM mutations
    Mrd_by_WGS_BM,
    cfWGS_BM_Binary,
    Cumulative_VAF_BM,
    detect_rate_BM,
    zscore_BM,
    z_score_detection_rate_BM,
    PercentChangeFromBaseline_BM,
    PercentChangeAtSecondTimepoint_BM,
    
    # 4) WGS‐based MRD, blood mutations
    Mrd_by_WGS_blood,
    cfWGS_blood_Binary,
    Cumulative_VAF_blood,
    detect_rate_blood,
    zscore_blood,
    z_score_detection_rate_blood,
    PercentChangeFromBaseline_blood,
    PercentChangeAtSecondTimepoint_blood,
    
    # 5) Anything else you might have missed
    everything()
  )


cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()


## Get duplicates
# 1) Pull out all the rows where (Patient, Date) is duplicated
dups <- cfWGS_Clinical_MRD_filled %>%
  group_by(Patient, Date, Sample_Code) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

# Inspect
dups %>% 
  select(Patient, Date, Sample_Code) %>%
  arrange(Patient, Date)

## Remove and consolidate 
#cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
#  filter(
#    is.na(Sample_ID) | !str_detect(Sample_ID, "(_Blood_Buffy_coat|_BM_cells|-B$|-O$)")
#  )

### Fill missing collection dates from the M4 processing log
# The processing log can contain duplicate dates for a sample. Where possible,
# the date closest to the matching lab date is retained.
M4_processing_log_dates <- read.csv("Exported_data_tables_clinical/M4_processing_log_dates.csv")

M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(Sample_Code = paste0(Patient, "-", Timepoint))

# Compare dates for matching Sample_Codes and calculate days difference
date_comparison <- cfWGS_Clinical_MRD_filled %>%
  # make sure Sample_Code is character, not logical
  mutate(Sample_Code = as.character(Sample_Code)) %>%
  select(Sample_Code, Date) %>%
  inner_join(
    M4_processing_log_dates %>%
      # likewise ensure it’s character
      mutate(Sample_Code = as.character(Sample_Code)) %>%
      select(Sample_Code, Cleaned_Date),
    by = "Sample_Code"
  ) %>%
  mutate(
    # coerce the cleaned date to Date
    Cleaned_Date   = as.Date(Cleaned_Date),
    # calculate days difference
    Days_Difference = as.numeric(Date - Cleaned_Date)
  ) %>%
  # drop any zero or NA differences
  filter(!is.na(Days_Difference) & Days_Difference != 0)

# Find duplicate date entries
duplicate_dates <- M4_processing_log_dates %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once

# Keep only those closest to lab 
# Merge the two datasets to compare dates
M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(Cleaned_Date = as.Date(Cleaned_Date))

## From previous script
M4_Labs <- read.csv(file = "M4_labs_cleaned.csv") %>% # generated in script 1_0
  mutate(LAB_DATE = as.Date(LAB_DATE))

# Merge the two datasets with a left join to retain all processing log dates
merged_data <- M4_processing_log_dates %>%
  left_join(
    M4_Labs %>%
      select(Patient, TIMEPOINT, LAB_DATE), 
    by = c("Patient" = "Patient", "Timepoint" = "TIMEPOINT")
  )

# Calculate the absolute difference in days (assign NA for rows without LAB_DATE)
merged_data <- merged_data %>%
  mutate(Date_Difference = ifelse(is.na(LAB_DATE), NA, abs(as.numeric(Cleaned_Date - LAB_DATE))))

# For each Sample_Code, keep the row with the smallest Date_Difference (or retain rows without LAB_DATE)
closest_date_entries <- merged_data %>%
  group_by(Sample_Code) %>%
  filter(is.na(Date_Difference) | Date_Difference == min(Date_Difference, na.rm = TRUE)) %>%
  ungroup() %>%
  unique()

# Keep only relevant columns
final_processing_log <- closest_date_entries %>%
  select(Patient, Timepoint, Sample_Code, Cleaned_Date) %>%
  unique()

# Find duplicate date entries
duplicate_dates <- final_processing_log %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once


# Define the rows to exclude
duplicates_to_exclude <- data.frame(
  Sample_Code = c("MJ-03-07", "VA-05-05", "ZC-06-01"),
  Cleaned_Date = as.Date(c("2020-11-04", "2019-12-05", "2021-11-12"))
)

# Filter out the duplicates from the main dataframe
M4_processing_log_dates <- M4_processing_log_dates %>%
  anti_join(duplicates_to_exclude, by = c("Sample_Code", "Cleaned_Date"))

## Recode VA-07 since 05 was actually relapse timepoint 
M4_processing_log_dates <- M4_processing_log_dates %>%
  mutate(
    # define the condition once
    is_va07_05 = Patient == "VA-07" & Timepoint == "05",
    
    # change Timepoint where needed
    Timepoint  = if_else(is_va07_05, "R-", Timepoint),
    
    # change Sample_Code where needed
    Sample_Code = if_else(is_va07_05, "VA-07-R", Sample_Code)
  ) %>%
  select(-is_va07_05)

## Take the Date value for labs since may have been typos in processing log 
# Update NA dates or sample codes in cfWGS_Clinical_MRD_filled with values from M4_processing_log_dates
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # 1) fix types in the main table
  mutate(
    Patient     = as.character(Patient),
    Timepoint   = as.character(Timepoint),
    Sample_Code = as.character(Sample_Code)
  ) %>%
  left_join(
    M4_processing_log_dates %>%
      # 2) also fix types in the lookup
      mutate(
        Patient     = as.character(Patient),
        Timepoint   = as.character(Timepoint),
        Sample_Code = as.character(Sample_Code)
      ) %>%
      # if there really are duplicate log rows, pick one
      distinct(Patient, Timepoint, .keep_all = TRUE) %>%
      select(Patient, Timepoint, Cleaned_Date, Sample_Code),
    by     = c("Patient", "Timepoint"),
    suffix = c("", ".log")   # keeps the original Sample_Code, brings in Sample_Code.log
  ) %>%
  mutate(
    # 3a) fill in Date from Cleaned_Date only when Date is NA
    Date = coalesce(Date, as.Date(Cleaned_Date)),
    
    # 3b) fill in Sample_Code from the log when the original is missing/empty
    Sample_Code = if_else(
      is.na(Sample_Code) | Sample_Code == "",
      Sample_Code.log,
      Sample_Code
    )
  ) %>%
  # 4) drop the helper columns
  select(-Cleaned_Date, -Sample_Code.log)



# Convert numeric dates back to proper Date format
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date, origin = "1970-01-01"))


### Add sample codes where missing 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% filter(!is.na(Patient))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    Sample_Code = if_else(
      # only rewrite when it's missing or blank
      is.na(Sample_Code) | Sample_Code == "",
      case_when(
        # SPORE → “PATIENT_T<Timepoint>”
        str_detect(Patient, "^SPORE") ~ paste0(Patient, "_T", Timepoint),
        # IMG → “PATIENT-<Timepoint>”
        str_detect(Patient, "^IMG")   ~ paste0(Patient, "-", Timepoint),
        # all others → default “PATIENT-<Timepoint>”
        TRUE                          ~ paste0(Patient, "-", Timepoint)
      ),
      Sample_Code
    )
  )

## Make corrections 

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(
    # 1) fix any Timepoint that is exactly "R-"
    Timepoint = str_replace(Timepoint, "^R-$", "R"),
    
    # 2) fix any Patient ending in "R-"
    Sample_Code   = str_replace(Sample_Code,   "R-$",    "R")
  )

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()

## Add missing timepoint_info 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  # 1) compute the “new” labels based purely on Timepoint (and that special CA-08 tweak)
  mutate(
    new_timepoint_info = case_when(
      Timepoint == "01"                       ~ "Diagnosis",
      Timepoint == "03"                       ~ "Post_induction",
      Timepoint == "05"                       ~ "Post_transplant",
      Timepoint == "07"                       ~ "Maintenance",
      Timepoint == "08"                       ~ "1.5yr maintenance",
      Timepoint == "09" & Patient == "CA-08"  ~ "1.5yr maintenance",  # CA-08 special
      Timepoint == "09"                       ~ "2yr maintenance",
      Timepoint == "10"                       ~ "2.5yr maintenance",
      Timepoint == "11"                       ~ "3yr maintenance",
      Timepoint == "12"                       ~ "3.5yr maintenance",
      Timepoint == "13"                       ~ "4yr maintenance",
      Timepoint == "14"                       ~ "4.5yr maintenance",
      Timepoint == "15"                       ~ "5yr maintenance",
      Timepoint == "R"                        ~ "Relapse",
      TRUE                                    ~ NA_character_
    )
  ) %>%
  # 2) only overwrite the old timepoint_info when it was NA
  mutate(
    timepoint_info = coalesce(timepoint_info, new_timepoint_info)
  ) %>%
  # 3) drop the helper column
  select(-new_timepoint_info)



# Find duplicate date entries
duplicate_dates <- cfWGS_Clinical_MRD_filled %>%
  group_by(Sample_Code) %>%  # Group by the sample code column
  filter(dplyr::n() > 1)  # Keep only rows where the Cleaned_Date appears more than once

# Find sample codes with NA in the Date column
missing_date_samples <- cfWGS_Clinical_MRD_filled %>%
  dplyr::filter(is.na(Date)) %>%
  distinct(Sample_Code)

# Export the missing date sample codes to a CSV
write.csv(
  missing_date_samples,
  file = "Exported_data_tables_clinical/missing_date_samples_May2025_2.csv",
  row.names = FALSE
)



### Add relapse/progression timing
# Relapse is assigned when a sample date is within the supported relapse window
# from the curated relapse-date table. These fields support later PFS analyses.

Relapse_dates_full <- read_csv("Exported_data_tables_clinical/Relapse dates cfWGS updated2.csv") 

# Step 2: Calculate Num_days_to_closest_relapse
# Join on Patient and calculate the days difference
Time_to_relapse <- cfWGS_Clinical_MRD_filled %>% select(Patient, Date) %>% unique() %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date) %>%
  mutate(
    # Calculate non-absolute days to relapse only for valid dates
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date ~ as.numeric(Date - Progression_date), 
      Progression_date <= Date & Progression_date >= Date - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    Num_days_to_closest_relapse = Num_days_to_closest_relapse_absolute,
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup() %>%
  unique()

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% unique()

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(Time_to_relapse %>% select(Patient, Date, Num_days_to_closest_relapse) %>% unique())

## Add a column for Relapsed 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>% 
  mutate(Relapsed = ifelse(!is.na(Num_days_to_closest_relapse), "Yes", "No"))

## Some editing 
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Sample_Code == "IMG-060-T36", "Maintenance", timepoint_info))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Timepoint == "R", "Relapse", timepoint_info))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(timepoint_info = ifelse(Sample_Code == "IMG-159-T11", "Maintenance", timepoint_info))


### Now get the relapse date as the time since the baseline marrow and use this for the PFS plot


### Export intermediate clinical MRD/date table
# This checkpoint is useful for auditing date/timepoint mapping before the WGS,
# fragmentomics, lab, and mutation-burden joins below.
write.csv(cfWGS_Clinical_MRD_filled, file = "cfWGS clinical MRD values with timepoint and dates updated August 2025.csv", row.names = F)




#### Add fragmentomics metrics
# Fragmentomics metrics are averaged at Patient/Date when duplicate rows exist.
# These columns are retained as continuous covariates for downstream modeling
# and correlation analyses.

# 1) Read in the fragmentomics + cfWGS data
##### Ensure is correct with everything
frag <- read_csv("Results_Fragmentomics/Key_fragmentomics_data_updated2.csv")
## Filter the frag df 

# 1) Load fragmentomics + cfWGS features, keeping only the new columns you need
# 2) Keep only the fragmentomics‐specific metrics plus join keys
frag_small <- frag %>%
  select(
    Patient,
    Date_of_sample_collection,
    FS,
    FS_unharmonized,
    FS_harmonized_median_mad,
    Proportion.Short,
    Proportion.Short_unharmonized,
    Proportion.Short_harmonized_median_mad,
    Site,
    Mean.Coverage,
    Mean.Coverage_unharmonized,
    Mean.Coverage_harmonized_median_mad,
    Midpoint.Coverage,
    Midpoint.normalized,
    Amplitude,
    Zscore.Coverage,
    Zscore.Midpoint,
    Zscore.Amplitude,
    Threshold.Coverage,
    Threshold.Midpoint,
    Threshold.Amplitude
  ) %>%
  rename(
    Date        = Date_of_sample_collection
  )

frag_small <- frag_small %>% filter(!is.na(Patient))
# 2) If you only want one row per Patient/Date/Sample_Code-drop any duplicates
# 2) Collapse duplicates by averaging all numeric columns
frag_small_unique <- frag_small %>%
  group_by(Patient, Date) %>%
  summarise(
    across(
      c(
        FS,
        FS_unharmonized,
        FS_harmonized_median_mad,
        Proportion.Short,
        Proportion.Short_unharmonized,
        Proportion.Short_harmonized_median_mad,
        Mean.Coverage,
        Mean.Coverage_unharmonized,
        Mean.Coverage_harmonized_median_mad,
        Midpoint.Coverage,
        Midpoint.normalized,
        Amplitude,
        Zscore.Coverage,
        Zscore.Midpoint,
        Zscore.Amplitude,
        Threshold.Coverage,
        Threshold.Midpoint,
        Threshold.Amplitude
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )



# 3) Join onto the clinical MRD table
cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  left_join(
    frag_small_unique
  )





##### Add clinical/lab/FISH metadata
# This master clinical table is generated by 1_1B and contributes baseline labs,
# subtype, staging, and clinical FISH/cytogenetic variables.
clinical_data_integrate <- read.csv("Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv") ### Generated in 1_1B

# make sure both date columns are Date class
clinical_data_integrate <- clinical_data_integrate %>%
  mutate(Date = as.Date(Date))

cfWGS_Clinical_MRD_filled <- cfWGS_Clinical_MRD_filled %>%
  mutate(Date = as.Date(Date))

## See duplicates 
dups <- clinical_data_integrate %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

## Add the sample codes 
clinical_data_integrate <- clinical_data_integrate %>%
  mutate(
    Sample_Code =
      case_when(
        # SPORE → “PATIENT_T<Timepoint>”
        str_detect(Patient, "^SPORE") ~ paste0(Patient, "_T", Timepoint),
        # IMG → “PATIENT-<Timepoint>”
        str_detect(Patient, "^IMG")   ~ paste0(Patient, "-", Timepoint),
        # all others → default “PATIENT-<Timepoint>”
        TRUE                          ~ paste0(Patient, "-", Timepoint)
      ),
    Sample_Code
  )


# rather than fuzzy join, join by study timepoint 
joined <- clinical_data_integrate %>%
  full_join(
    cfWGS_Clinical_MRD_filled,
    by = c("Patient" = "Patient", "Sample_Code" = "Sample_Code"))

# now coalesce the two Timepoint / timepoint_info pairs,
# and drop the .y (cfWGS) helper columns
joined <- joined %>%
  # 1) Coalesce every .x/.y pair into a single column
  {
    # find all the bases that have both .x and .y
    xy_bases <- intersect(
      sub("\\.x$", "", grep("\\.x$", names(.), value = TRUE)),
      sub("\\.y$", "", grep("\\.y$", names(.), value = TRUE))
    )
    
    # for each base, create a single coalesced column
    coalesced <- map_dfc(xy_bases, function(base) {
      xcol <- paste0(base, ".x")
      ycol <- paste0(base, ".y")
      tibble(!!base := coalesce(.[[ycol]], .[[xcol]]))
    })
    
    # drop all the old .x/.y and bind the new coalesced ones
    select(., -matches("\\.(x|y)$")) %>%
      bind_cols(coalesced)
  }

## Fill demo info 
# pull out one row per patient with their non-NA demographics
demo_map <- joined %>%
  select(Patient, AGE, Gender, AGE_GROUP) %>%
  distinct(Patient, .keep_all = TRUE)

# 2) drop the old (possibly‐sparse) demo columns, and re-attach the full map
cfWGS_Clinical_MRD_filled_final <- joined %>%
  select(-AGE, -Gender, -AGE_GROUP) %>%
  left_join(demo_map, by = "Patient")


### The dates don't match exactly - strange...
cfWGS_Clinical_MRD_filled_final <- cfWGS_Clinical_MRD_filled_final %>% unique() 
dups <- cfWGS_Clinical_MRD_filled_final %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()

write.csv(dups, file = "duplicate_check.csv")
### Deduplicate before WGS feature joins
# At this point multiple assay sources may describe the same patient/timepoint.
# Rows with cfWGS z-scores are prioritized because those dates are used in the
# downstream cfWGS analysis.

## helper: coalesce down an arbitrary vector (character OR numeric)
first_non_na <- function(x) {
  # treat "" and "-" as NA too
  x[x %in% c("", "-")] <- NA
  # return the first non‑NA (or NA if none)
  x[which(!is.na(x))[1]]
}

cfWGS_dedup <- cfWGS_Clinical_MRD_filled_final %>%
  # ------- 1. add priority flags --------------------------------------------
mutate(
  has_z     = !is.na(zscore_BM) | !is.na(zscore_blood),
  relapse_y = str_to_upper(Relapsed) == "YES"
) %>%
  
  # ------- 2. sort so best rows come first ----------------------------------
arrange(Patient, Timepoint,
        desc(has_z),          # 1st priority: has a z‑score, date used for rest of analysis 
        desc(relapse_y)) %>%  # 2nd priority: relapsed == "Yes", since this was not held in some older tables
  
  # ------- 3. collapse each Patient‑Timepoint group -------------------------
group_by(Patient, Timepoint) %>%
  summarise(
    across(
      # drop helper columns before summarising
      -c(has_z, relapse_y),
      first_non_na,
      .names = "{.col}"
    ),
    .groups = "drop"
  )


dups <- cfWGS_dedup %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



#### Add WGS features from 1_5
# These columns include tumor fraction, arm-level CNA, IgH translocations,
# mutation summaries, and the composite WGS evidence-of-disease flag.
# 1. From the WGS feature table, select the key and all of the WGS columns needed for this analysis
wgs_subset <- All_feature_data %>%
  select(
    Patient,
    Date_of_sample_collection, 
    Sample_type,
    Timepoint,
    Tumor_Fraction,
    # CNAs
    del1p, amp1q, del13q, del17p, hyperdiploid,
    # IGH translocations
    IGH_MAF, IGH_CCND1, IGH_MYC, IGH_FGFR3,
    # mutation summary
    Mut_identified, Mut_genes, Mut_highest_VAF, Mut_type, Evidence_of_Disease
  )

# 2. Rename them with a “WGS_” prefix so they don’t overwrite the clinical calls
wgs_prefixed <- wgs_subset %>%
  rename_with(
    ~ paste0("WGS_", .),
    -c(Patient, Date_of_sample_collection, Sample_type, Timepoint)
  )

## Replace NA to 0
wgs_prefixed <- wgs_prefixed %>%
  mutate(across(
    all_of(c(
      "WGS_del1p", "WGS_amp1q", "WGS_del13q", "WGS_del17p",
      "WGS_IGH_MAF", "WGS_IGH_CCND1", "WGS_IGH_MYC", "WGS_IGH_FGFR3"
    )),
    ~ replace_na(.x, 0)
  ))

## Make wide 

## See what is dup
wgs_prefixed %>% 
  filter(Sample_type %in% c("BM_cells","Blood_plasma_cfDNA")) %>%
  group_by(Patient, Timepoint, Date_of_sample_collection, Sample_type) %>%
  filter(dplyr::n()>1) %>%
  ungroup() -> dups

# how many duplicate groups?
dups %>% 
  summarise(count = dplyr::n(), .by = c(Patient, Timepoint, Date_of_sample_collection, Sample_type)) %>% 
  arrange(desc(count))

## Correct 
wgs_prefixed <- wgs_prefixed %>%
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  group_by(Patient, Timepoint, Date_of_sample_collection, Sample_type) %>%
  summarise(
    # 1) mean of Tumor_Fraction
    WGS_Tumor_Fraction = mean(WGS_Tumor_Fraction, na.rm = TRUE),
    # 2) first() of all the other WGS_ columns
    across(
      .cols = starts_with("WGS_") & 
        !all_of("WGS_Tumor_Fraction"),
      .fns  = first
    ),
    .groups = "drop"
  )

# Normalize timepoints
wgs_prefixed <- wgs_prefixed %>%
  mutate(
    Timepoint = str_replace(Timepoint, "^R-$", "R")
  )

wgs_wide <- wgs_prefixed %>%
  # 1. keep only BM_cells and Blood_plasma_cfDNA
  filter(Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")) %>%
  
  # 2. pivot so each WGS_… column splits into two: one per sample type
  pivot_wider(
    id_cols        = c(Patient, Timepoint, Date_of_sample_collection),
    names_from     = Sample_type,
    values_from    = starts_with("WGS_"),
    names_sep      = "_"
  )


# 3. Join onto  clinical MRD table by matching modifiers 
# make sure both date columns are Date objects
# Normalize the timepoints 
cfWGS_Clinical_MRD_filled_final <- cfWGS_dedup %>%
  mutate(
    Date      = as.Date(Date),
    Timepoint = str_replace(Timepoint, "^R-$", "R"),
    Sample_Code = str_replace(Sample_Code, "R-$", "R")
  )


wgs_wide <- wgs_wide %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

# Fuzzy join on Patient and Timepoint exactly, with collection dates within
# seven days. The join does not include Sample_type.
library(fuzzyjoin)

joined <- cfWGS_Clinical_MRD_filled_final %>%
  fuzzy_full_join(
    wgs_wide,
    by = c(
      "Patient" = "Patient",
      "Timepoint" = "Timepoint",
      "Date" = "Date_of_sample_collection"
    ),
    match_fun = list(
      `==`,            # Patient must match exactly
      `==`,            # Timepoint must match exactly
      function(d1, d2) abs(d1 - d2) <= 7  # dates within 7 days
    )
  )  


## reorganize
# 1. coalesce + drop the .x/.y duplicates
joined2 <- joined %>%
  mutate(
    Patient   = coalesce(Patient.x, Patient.y),
    Timepoint = coalesce(Timepoint.x, Timepoint.y)
  ) %>%
  select(-Patient.x, -Patient.y, -Timepoint.x, -Timepoint.y)

# 2. define the blocks
id_cols      <- c("Patient", "Timepoint", "Sample_Code", "timepoint_info",
                  "Date", "Date_of_sample_collection", "AGE", "Gender", "AGE_GROUP")

ext_cols     <- c("DIAGNOSIS_DATE","ECOG_SCORE","ISS_STAGE","KPS_SCORE",
                  "DETAILS","ISS_stage_notes","R_ISS_STAGE",
                  "Cytogenetic_Risk","Abnormality_Count","Abnormality_Status")

lab_cols     <- c("Study","dFLC","Albumin","B2_micro","Calcium","Creatinine",
                  "Hemoglobin","LDH","Plasma_pct","Plasma_count",
                  "IgA","IgG","IgM","Kappa","Kappa_Lambda_Ratio",
                  "Lambda","M_Protein","Subtype")

fish_cols    <- c("T_4_14","T_11_14","T_14_16",
                  "DEL_17P","DEL_1P","AMP_1Q","Fish_abnormalities",
                  "DEL_13","DEL_14",
                  "trisomy_3","trisomy_7","trisomy_9","trisomy_15",
                  "hyperdiploid","T_14_20","CKS1B_DUP")

flow_cols    <- c("Flow_pct_cells","Flow_Binary","MRD_Results_FLOW")
adaptive_cols<- c("PerMillionCount_Adaptive","Adaptive_Frequency","Adaptive_Binary")
rapid_cols   <- c("Rapid_Novor","Rapid_Novor_Binary")
pet_cols     <- c("PET","PET_Binary")

clinical_mrd_cols <- c(
  "MRD_by_clinical_testing","MRD_by_clinical_testing_stringent",
  "MRD_Clinical_Binary","MRD_Clinical_Stringent_Binary"
)

wgs_mrd_cols <- c(
  "Mrd_by_WGS_BM","cfWGS_BM_Binary","Cumulative_VAF_BM",
  "detect_rate_BM", "sites_rate_BM", "zscore_BM", "z_score_detection_rate_BM",
  "PercentChangeFromBaseline_BM","PercentChangeAtSecondTimepoint_BM",
  "Mrd_by_WGS_blood","cfWGS_blood_Binary","Cumulative_VAF_blood",
  "detect_rate_blood", "sites_rate_blood", "zscore_blood", "z_score_detection_rate_blood",
  "PercentChangeFromBaseline_blood","PercentChangeAtSecondTimepoint_blood",
  "Mrd1E5","Mrd1E6"
)

coverage_cols <- c(
  "FS", "Proportion.Short",                    
  "Mean.Coverage","Midpoint.Coverage","Midpoint.normalized","Amplitude",
  "Zscore.Coverage","Zscore.Midpoint","Zscore.Amplitude",
  "Threshold.Coverage","Threshold.Midpoint","Threshold.Amplitude"
)

# grab the remaining WGS_… columns
wgs_cols <- grep("^WGS_", colnames(joined2), value = TRUE)

# 3. put it all together
joined_clean <- joined2 %>%
  select(
    # identifiers & dates
    all_of(id_cols),
    # extrinsic
    all_of(ext_cols),
    # labs
    all_of(lab_cols),
    # FISH / cytogenetics
    all_of(fish_cols),
    # flow / adaptive / rapid / PET
    all_of(flow_cols),
    all_of(adaptive_cols),
    all_of(rapid_cols),
    all_of(pet_cols),
    # clinical MRD
    all_of(clinical_mrd_cols),
    # WGS MRD
    all_of(wgs_mrd_cols),
    # coverage (incl FS)
    all_of(coverage_cols),
    # all the rest of the WGS_… columns
    all_of(wgs_cols),
    # anything else not yet selected
    everything()
  )


joined_clean <- joined_clean %>% select(-Date_of_sample_collection)

## get duplicates 
dups <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()


### There are still dups, remove 

# helper to grab the first non‐NA (or NA if none)
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else NA
}

joined_clean <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    across(everything(), first_non_na),
    .groups = "drop"
  )

## get duplicates 
dups <- joined_clean %>%
  group_by(Patient, Timepoint) %>%
  filter(dplyr::n() > 1) %>%
  ungroup()



#### Add raw mutation-burden count tables
# These count tables are not the myeloma-panel mutation summaries from 1_5.
# They provide broad mutation burden used later for the relaxed blood evidence
# helper flag.

## Load in the metadata and try to match to it
cfWGS_metadata <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)

# Create a VCF_clean_merge column from the 'Bam' column in cfWGS_metadata
cfWGS_metadata <- cfWGS_metadata %>%
  mutate(VCF_clean_merge = gsub("\\.filter.*", "", Bam))  # Remove '.filter' and anything following it

# Modify select non-overlapping cases - OICR put PG but should be WG to match (just internal naming conventions)
cfWGS_metadata <- cfWGS_metadata %>%
  mutate(VCF_clean_merge = ifelse(
    VCF_clean_merge %in% c(
      "TFRIM4_0031_Bm_P_PG_M4-CA-02-01-O-DNA",
      "TFRIM4_0032_Bm_P_PG_M4-HP-01-01-O-DNA",
      "TFRIM4_0033_Bm_P_PG_M4-MJ-06-01-O-DNA",
      "TFRIM4_0034_Bm_P_PG_M4-VA-02-01-O-DNA",
      "TFRIM4_0035_Bm_P_PG_M4-VA-06-01-O-DNA"
    ), 
    gsub("PG", "WG", VCF_clean_merge),  # Replace WG with PG for these specific values
    VCF_clean_merge  # Keep other values unchanged
  ))


### Now load in the counts files 
# 1. Load the counts files
Blood_muts <- read.delim("Mutation_counts/mutation_counts_table_Blood_no_RSID.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BM_muts    <- read.delim("Mutation_counts/mutation_counts_table_BM_no_RSID.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Clean the 'File' column: remove '.filter' and anything after
Blood_muts$File <- sub("\\.fil.*$", "", Blood_muts$File)
Blood_muts <- Blood_muts[!grepl("\\.PASS", Blood_muts$File), ]
Blood_muts <- Blood_muts[!grepl("\\Bm_P", Blood_muts$File), ]


BM_muts$File    <- sub("\\.fil.*$", "", BM_muts$File)
BM_muts$File <- sub("\\.mutect2.*$", "", BM_muts$File)


# Update specific value in VCF_clean that was truncated
BM_muts <- BM_muts %>%
  mutate(File = ifelse(
    File == "TFRIM4_0189_Bm_P_WG_ZC-02", 
    "TFRIM4_0189_Bm_P_WG_ZC-02-01-O-DNA", 
    File  # Keep other values unchanged
  ))


# join BM counts
cfWGS_with_counts <- cfWGS_metadata %>%
  left_join(
    BM_muts %>% 
      rename(BM_Mutation_Count = Mutation_Count),
    by = c("VCF_clean_merge" = "File")
  ) %>%
  # join Blood counts
  left_join(
    Blood_muts %>% 
      rename(Blood_Mutation_Count = Mutation_Count),
    by = c("VCF_clean_merge" = "File")
  )


## Add this to the main table 
# 1) make sure both date columns are Date
joined_clean <- joined_clean %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

counts_df <- cfWGS_with_counts %>%
  mutate(collection_date = as.Date(Date_of_sample_collection, "%Y-%m-%d")) %>%
  filter(Sample_type != "Blood_Buffy_coat") %>%
  select(Patient, Timepoint, BM_Mutation_Count, Blood_Mutation_Count) %>% unique()

## Check for dups
counts_df %>%
  count(Patient, Timepoint) %>%
  filter(n > 1)

# Consolidate
counts_consolidated <- counts_df %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    BM_Mutation_Count = {
      bm <- BM_Mutation_Count[!is.na(BM_Mutation_Count)]
      if (length(bm)) bm[1] else NA_integer_
    },
    Blood_Mutation_Count = {
      bl <- Blood_Mutation_Count[!is.na(Blood_Mutation_Count)]
      if (length(bl)) bl[1] else NA_integer_
    },
    .groups = "drop"
  )

# Add broad mutation burdens for the Spring 2026 test-cohort expansion. These
# are computed from the supplied post-filter MAFs using the same no-rsID rule as
# the historical static count files. Revision values take precedence for an
# overlapping Patient/Timepoint (notably IMG-159), because the revision metadata
# and MAFs are the current sample authority.
spring2026_mutation_count_audit <- spring2026_revision_mutation_counts(required = TRUE)
dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  spring2026_mutation_count_audit,
  file.path(
    "Output_tables_2025", "clinical_support",
    "spring2026_revision_mutation_count_audit.csv"
  )
)

spring2026_counts_wide <- spring2026_mutation_count_audit %>%
  select(Patient, Timepoint, Sample_type, Mutation_Count) %>%
  mutate(
    count_column = case_when(
      Sample_type == "BM_cells" ~ "BM_Mutation_Count",
      Sample_type == "Blood_plasma_cfDNA" ~ "Blood_Mutation_Count",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(count_column)) %>%
  select(-Sample_type) %>%
  tidyr::pivot_wider(names_from = count_column, values_from = Mutation_Count)

counts_consolidated <- counts_consolidated %>%
  full_join(
    spring2026_counts_wide,
    by = c("Patient", "Timepoint"),
    suffix = c("_historical", "_spring2026")
  ) %>%
  transmute(
    Patient,
    Timepoint,
    BM_Mutation_Count = coalesce(BM_Mutation_Count_spring2026, BM_Mutation_Count_historical),
    Blood_Mutation_Count = coalesce(Blood_Mutation_Count_spring2026, Blood_Mutation_Count_historical)
  )

counts_consolidated <- counts_consolidated %>%
  mutate(Timepoint = if_else(Timepoint == "R-", "R", Timepoint))


# 2) Join mutation counts exactly by Patient and Timepoint; no date window is
#    applied at this step.
joined_with_counts <- joined_clean %>%
  left_join(
    counts_consolidated,
    by = c(
      "Patient" = "Patient",
      "Timepoint"        = "Timepoint"
    )
  )  


### Now correct the duplicate rows in joined with counts 
first_non_na <- function(x) {
  x[x %in% c("", "-")] <- NA      # treat blanks and dashes as missing
  vec <- x[!is.na(x)]             # drop all NA
  if (length(vec)) vec[[1]] else NA
}

joined_consolidated <- joined_with_counts %>%
  # normalize the “R-” codes up front
  mutate(
    Timepoint   = str_replace(Timepoint,   "^R-$", "R"),
    Sample_Code = str_replace(Sample_Code, "^R-$", "R")
  ) %>%
  group_by(Patient, Timepoint, Date) %>%   # or just (Patient, Timepoint)
  reframe(
    across(everything(), first_non_na)
  )


## Export first aggregate checkpoint.
# This checkpoint is retained for traceability; downstream scripts should use
# the later `updated9` export unless they explicitly document otherwise.
# Write to CSV (for Excel/sharing)
write.csv(joined_consolidated, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(joined_consolidated, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated.rds")

pre_cleanup_sample_identity_map <- joined_consolidated %>%
  filter(!is.na(Patient), !is.na(Sample_Code)) %>%
  transmute(
    Patient = as.character(Patient),
    Sample_Code = as.character(Sample_Code),
    canonical_Timepoint = as.character(Timepoint),
    canonical_timepoint_info = as.character(timepoint_info),
    canonical_Date = as.Date(Date),
    identity_source = "current_pre_cleanup",
    identity_priority = 2L
  ) %>%
  filter(!is.na(canonical_Timepoint) | !is.na(canonical_timepoint_info) | !is.na(canonical_Date)) %>%
  distinct()

legacy_submitted_identity_path <- "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated6.csv"
if (file.exists(legacy_submitted_identity_path)) {
  legacy_submitted_identity_map <- readr::read_csv(
    legacy_submitted_identity_path,
    show_col_types = FALSE
  ) %>%
    filter(!is.na(Patient), !is.na(Sample_Code)) %>%
    transmute(
      Patient = as.character(Patient),
      Sample_Code = as.character(Sample_Code),
      canonical_Timepoint = as.character(Timepoint),
      canonical_timepoint_info = as.character(timepoint_info),
      canonical_Date = as.Date(Date),
      identity_source = "legacy_submitted_updated6",
      identity_priority = 1L
    ) %>%
    filter(!is.na(canonical_Timepoint) | !is.na(canonical_timepoint_info) | !is.na(canonical_Date)) %>%
    distinct()

  pre_cleanup_sample_identity_map <- bind_rows(
    legacy_submitted_identity_map,
    pre_cleanup_sample_identity_map
  ) %>%
    arrange(Patient, Sample_Code, identity_priority) %>%
    group_by(Patient, Sample_Code) %>%
    summarise(
      canonical_Timepoint = first(canonical_Timepoint[!is.na(canonical_Timepoint)]),
      canonical_timepoint_info = first(canonical_timepoint_info[!is.na(canonical_timepoint_info)]),
      canonical_Date = first(canonical_Date[!is.na(canonical_Date)]),
      identity_source = first(identity_source),
      .groups = "drop"
    )
}

pre_cleanup_identity_conflicts <- pre_cleanup_sample_identity_map %>%
  distinct(Patient, Sample_Code, canonical_Timepoint, canonical_timepoint_info, canonical_Date) %>%
  count(Patient, Sample_Code, name = "n_identity_rows") %>%
  filter(n_identity_rows > 1) %>%
  left_join(pre_cleanup_sample_identity_map, by = c("Patient", "Sample_Code"))

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  pre_cleanup_identity_conflicts,
  file.path("Output_tables_2025", "clinical_support", "pre_cleanup_sample_identity_conflicts_audit.csv")
)

pre_cleanup_sample_identity_map <- pre_cleanup_sample_identity_map %>%
  anti_join(
    pre_cleanup_identity_conflicts %>% distinct(Patient, Sample_Code),
    by = c("Patient", "Sample_Code")
  )

### Do a bit more cleaning for duplicate or missing timepoints

# 1) two lookup tables from the clinical master (de-duplicated if needed)
cc_by_info <- combined_clinical_data_updated %>% 
  select(Patient, timepoint_info, 
         TP_info = Timepoint, 
         Date_info = Date_of_sample_collection)

cc_by_date <- combined_clinical_data_updated %>% 
  select(Patient, Date_of_sample_collection, 
         TP_date = Timepoint, 
         Info_date = timepoint_info)

cc_by_date <- cc_by_date %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection, origin = "1970-01-01")
  )

# 2) left-join both onto joined_consolidated
joined_filled2 <- joined_consolidated %>% 
  # join by info
  left_join(cc_by_info, 
            by = c("Patient","timepoint_info")) %>% 
  # join by date
  left_join(cc_by_date, 
            by = c("Patient","Date" = "Date_of_sample_collection")) %>% 
  
  mutate(
    # only replace missing Timepoint
    Timepoint = if_else(
      is.na(Timepoint),
      coalesce(TP_info, TP_date),
      Timepoint
    ),
    
    # only replace missing timepoint_info
    timepoint_info = if_else(
      is.na(timepoint_info),
      coalesce(Info_date, timepoint_info),
      timepoint_info
    ),
    
    # only replace missing Date
    Date = if_else(
      is.na(Date),
      coalesce(Date_info, Date),
      Date
    )
  ) %>%
  # drop the helper columns
  select(-TP_info, -Date_info, -TP_date, -Info_date)

# 3) if you want to toss any that still lack a Timepoint
joined_filled2 <- joined_filled2 %>% 
  filter(!is.na(Timepoint))

# 4) finally de-duplicate exactly as before
joined_clean2 <- joined_filled2 %>%
  group_by(Patient, Timepoint, Date) %>%
  reframe(
    across(everything(), first_non_na)
  )


## Export second aggregate checkpoint after missing timepoint/date cleanup.
# This checkpoint is retained for traceability; the final current table is
# produced after missing tumor fraction and relaxed-evidence cleanup below.
# Write to CSV (for Excel/sharing)
write.csv(joined_clean2, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated2.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(joined_clean2, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated2.rds")


## Get unique pts 
patients_with_zscores <- joined_clean2 %>%
  filter(!is.na(zscore_BM) | !is.na(zscore_blood)) %>%
  distinct(Patient)


### Add manually recovered missing tumor fractions
# This helper table fills tumor-fraction values that were available from review
# but missing after the automated joins. Existing non-missing TF values are not
# overwritten.
missing_tumor_fraction <- read.csv("Missing_tumor_fraction.csv")

# 3) join back onto dat and fill the two WGS‐tumor‐fraction columns
joined_clean2 <- joined_clean2 %>%
  # join in the new BM / blood TF values
  left_join(
    missing_tumor_fraction %>%
      select(Patient, Timepoint,
             tf_BM_new = WGS_Tumor_Fraction_BM_cells,
             tf_blood_new = WGS_Tumor_Fraction_Blood_plasma_cfDNA),
    by = c("Patient", "Timepoint")
  ) %>%
  # overwrite only where the originals were NA
  mutate(
    WGS_Tumor_Fraction_BM_cells = coalesce(WGS_Tumor_Fraction_BM_cells, tf_BM_new),
    WGS_Tumor_Fraction_Blood_plasma_cfDNA = coalesce(WGS_Tumor_Fraction_Blood_plasma_cfDNA, tf_blood_new)
  ) %>%
  # drop the helper cols
  select(-tf_BM_new, -tf_blood_new)


### One more duplicate check 
# 1) Identify which (Patient, timepoint_info) groups have a “-NA” sample_code
dup_keys <- joined_clean2 %>%
  filter(str_detect(Sample_Code, "-NA$")) %>%
  distinct(Patient, timepoint_info)

# 2) For those groups, coalesce the two rows into one-
#    preferring any value coming from the non-“-NA” Sample_Code row
collapsed <- joined_clean2 %>%
  semi_join(dup_keys, by = c("Patient","timepoint_info")) %>%
  group_by(Patient, timepoint_info) %>%
  summarise(across(everything(), ~ {
    col <- cur_column()
    vals <- .x
    
    if (col == "Sample_Code") {
      # drop the “-NA” code, then take the first remaining
      clean <- vals[!str_detect(vals, "-NA$")]
      if (length(clean)>0) clean[1] else NA_character_
    } else {
      # take the first non-NA value
      non_na <- vals[!is.na(vals)]
      if (length(non_na)>0) non_na[1] else NA
    }
  }), .groups="drop")

# Keep the collapsed duplicate rows type-stable with the source table. When all
# values in a collapsed group are missing, the generic `else NA` branch above can
# otherwise create logical columns (for example Timepoint), which makes
# `bind_rows()` fail against the character/date/numeric columns in `unchanged`.
for (nm in intersect(names(collapsed), names(joined_clean2))) {
  if (is.character(joined_clean2[[nm]])) {
    collapsed[[nm]] <- as.character(collapsed[[nm]])
  } else if (inherits(joined_clean2[[nm]], "Date")) {
    collapsed[[nm]] <- as.Date(collapsed[[nm]])
  } else if (inherits(joined_clean2[[nm]], c("POSIXct", "POSIXlt"))) {
    collapsed[[nm]] <- as.POSIXct(collapsed[[nm]])
  } else if (is.numeric(joined_clean2[[nm]])) {
    collapsed[[nm]] <- as.numeric(collapsed[[nm]])
  } else if (is.logical(joined_clean2[[nm]])) {
    collapsed[[nm]] <- as.logical(collapsed[[nm]])
  }
}

# 3) Everything *not* in those duplicate groups just stays as-is
unchanged <- joined_clean2 %>%
  anti_join(dup_keys, by = c("Patient","timepoint_info"))

# 4) Re-combine
joined_clean2_fixed <- bind_rows(unchanged, collapsed) %>%
  arrange(Patient, Date)   # or however you like to sort

joined_clean2 <- joined_clean2_fixed %>%
  mutate(
    Sample_Code = if_else(
      is.na(Sample_Code),
      paste0(Patient, "-", Timepoint),
      Sample_Code
    )
  )


## Another check
# helper to treat "" or "-" as NA, then pick first non-NA
first_non_na <- function(x) {
  x[x %in% c("", "-")] <- NA
  vals <- x[!is.na(x)]
  if (length(vals)) vals[1] else NA
}

joined_clean3 <- joined_clean2 %>%
  # 1) fix any Sample_Code ending in "R-"
  mutate(
    Sample_Code = str_replace(Sample_Code, "R-$", "R")
  ) %>%
  
  # 2) group on Patient, Date, Timepoint and collapse
  group_by(Patient, Date, Timepoint) %>%
  reframe(
    across(everything(), first_non_na)
  )

## Clean dates 
df <- joined_clean3

# 1) For each Patient/Timepoint/Sample_Code, compute exactly one fill_date:
fill_dates <- df %>%
  group_by(Patient, Timepoint, Sample_Code) %>%
  summarise(
    fill_date = {
      ud <- unique(Date[!is.na(Date)])
      if (length(ud) == 1) {
        # exactly one date → keep it
        ud
      } else {
        # zero dates or multiple conflicting → NA
        as.Date(NA)
      }
    },
    .groups = "drop"
  )

# 2) Join back and coalesce only the missing Dates
filled_df <- df %>%
  left_join(fill_dates, by = c("Patient","Timepoint","Sample_Code")) %>%
  mutate(
    Date = coalesce(Date, fill_date)
  ) %>%
  select(-fill_date)

# 3) (Optional) sanity‐check that you didn’t miss any fillable gaps
residual <- filled_df %>%
  group_by(Patient, Timepoint, Sample_Code) %>%
  filter(any(is.na(Date)) & any(!is.na(Date))) %>%
  distinct(Patient, Timepoint, Sample_Code)

if (nrow(residual)) {
  warning("These groups still have mixed NA/non-NA dates:\n",
          paste0(sprintf("%s / %s / %s", residual$Patient, residual$Timepoint, residual$Sample_Code),
                 collapse = "\n"))
} else {
  message("All missing dates successfully filled where possible.")
}

filled_df <- filled_df %>% unique()

filled_df %>%
  # for safety, make sure Date is Date class
  mutate(Date = as.Date(Date)) %>%
  # group by Patient + Timepoint
  group_by(Patient, Timepoint) %>%
  # collect unique dates and count them
  summarise(
    n_dates = n_distinct(Date),
    dates   = paste(sort(unique(as.character(Date))), collapse = ", "),
    .groups = "drop"
  ) %>%
  # keep only those with more than one distinct date
  filter(n_dates > 1)

## One last pass to coalesce exact sample records only.
# Do not collapse by Patient + Date alone: older longitudinal rows can have
# missing or repeated dates, and grouping only on those two fields can merge
# distinct timepoints from the same patient.
patient_date_collision_audit <- filled_df %>%
  group_by(Patient, Date) %>%
  summarise(
    n_rows = dplyr::n(),
    n_timepoints = n_distinct(Timepoint, na.rm = TRUE),
    n_sample_codes = n_distinct(Sample_Code, na.rm = TRUE),
    timepoints = paste(sort(unique(as.character(Timepoint))), collapse = ";"),
    sample_codes = paste(sort(unique(as.character(Sample_Code))), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_rows > 1, n_timepoints > 1 | n_sample_codes > 1)

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  patient_date_collision_audit,
  file.path("Output_tables_2025", "clinical_support", "patient_date_groups_not_collapsed_in_updated9_audit.csv")
)

filled_df <- filled_df %>%
  group_by(Patient, Timepoint, Sample_Code, Date) %>%
  reframe(
    across(everything(), first_non_na)
  )


### One final pass 
## Add missing timepoint_info 
filled_df <- filled_df %>%
  # 1) compute the “new” labels based purely on Timepoint (and that special CA-08 tweak)
  mutate(
    new_timepoint_info = case_when(
      Timepoint == "01"                       ~ "Diagnosis",
      Timepoint == "03"                       ~ "Post_induction",
      Timepoint == "05"                       ~ "Post_transplant",
      Timepoint == "07"                       ~ "Maintenance",
      Timepoint == "08"                       ~ "1.5yr maintenance",
      Timepoint == "09" & Patient == "CA-08"  ~ "1.5yr maintenance",  # CA-08 special
      Timepoint == "09"                       ~ "2yr maintenance",
      Timepoint == "10"                       ~ "2.5yr maintenance",
      Timepoint == "11"                       ~ "3yr maintenance",
      Timepoint == "12"                       ~ "3.5yr maintenance",
      Timepoint == "13"                       ~ "4yr maintenance",
      Timepoint == "14"                       ~ "4.5yr maintenance",
      Timepoint == "15"                       ~ "5yr maintenance",
      Timepoint == "R"                        ~ "Relapse",
      TRUE                                    ~ NA_character_
    )
  ) %>%
  # 2) only overwrite the old timepoint_info when it was NA
  mutate(
    timepoint_info = coalesce(timepoint_info, new_timepoint_info)
  ) %>%
  # 3) drop the helper column
  select(-new_timepoint_info)

## Coalesce IMMAGINE-098 
first_non_na <- function(x) { y <- x[!is.na(x)]; if (length(y)) y[1] else NA }

img98_rows <- filled_df %>% filter(Patient == "IMG-098", Timepoint %in% c("T0", "T1"))

img98_coalesced <- img98_rows %>%
  summarise(across(everything(), first_non_na))

# 2) Drop the old IMG-098 T0/T1 rows and bind back the coalesced one
filled_df <- filled_df %>%
  filter(!(Patient == "IMG-098" & Timepoint %in% c("T0", "T1"))) %>%
  bind_rows(img98_coalesced) %>%
  # optional: re‐sort so IMG-098 sits in its original place
  arrange(Patient, Timepoint)

sample_identity_restoration_audit <- filled_df %>%
  mutate(
    Patient = as.character(Patient),
    Sample_Code = as.character(Sample_Code),
    Timepoint_before_identity_restore = as.character(Timepoint),
    timepoint_info_before_identity_restore = as.character(timepoint_info),
    Date_before_identity_restore = as.Date(Date)
  ) %>%
  left_join(pre_cleanup_sample_identity_map, by = c("Patient", "Sample_Code")) %>%
  filter(
    (!is.na(canonical_Timepoint) &
       coalesce(Timepoint_before_identity_restore, "") != canonical_Timepoint) |
      (!is.na(canonical_timepoint_info) &
         coalesce(timepoint_info_before_identity_restore, "") != canonical_timepoint_info) |
      (!is.na(canonical_Date) &
         (is.na(Date_before_identity_restore) | Date_before_identity_restore != canonical_Date))
  ) %>%
  transmute(
    Patient,
    Sample_Code,
    Timepoint_before_identity_restore,
    Timepoint_after_identity_restore = canonical_Timepoint,
    timepoint_info_before_identity_restore,
    timepoint_info_after_identity_restore = canonical_timepoint_info,
    Date_before_identity_restore,
    Date_after_identity_restore = canonical_Date,
    identity_source
  )

readr::write_csv(
  sample_identity_restoration_audit,
  file.path("Output_tables_2025", "clinical_support", "sample_identity_restoration_from_pre_cleanup_audit.csv")
)

filled_df <- filled_df %>%
  mutate(
    Patient = as.character(Patient),
    Sample_Code = as.character(Sample_Code)
  ) %>%
  left_join(pre_cleanup_sample_identity_map, by = c("Patient", "Sample_Code")) %>%
  mutate(
    Timepoint = coalesce(canonical_Timepoint, as.character(Timepoint)),
    timepoint_info = coalesce(canonical_timepoint_info, as.character(timepoint_info)),
    Date = coalesce(canonical_Date, as.Date(Date))
  ) %>%
  select(-canonical_Timepoint, -canonical_timepoint_info, -canonical_Date, -any_of(c("identity_source", "identity_priority")))

if (file.exists(legacy_submitted_identity_path)) {
  legacy_submitted_rows <- readr::read_csv(
    legacy_submitted_identity_path,
    show_col_types = FALSE
  ) %>%
    mutate(
      Patient = as.character(Patient),
      Sample_Code = as.character(Sample_Code),
      Date = as.Date(Date)
    )

  legacy_rows_missing_after_current_assembly <- legacy_submitted_rows %>%
    anti_join(
      filled_df %>% distinct(Patient, Sample_Code),
      by = c("Patient", "Sample_Code")
    )

  readr::write_csv(
    legacy_rows_missing_after_current_assembly %>%
      select(any_of(c("Patient", "Sample_Code", "Timepoint", "timepoint_info", "Date"))),
    file.path("Output_tables_2025", "clinical_support", "legacy_submitted_rows_appended_to_updated9_audit.csv")
  )

  if (nrow(legacy_rows_missing_after_current_assembly) > 0L) {
    for (col in intersect(names(filled_df), names(legacy_rows_missing_after_current_assembly))) {
      template <- filled_df[[col]]
      if (inherits(template, "Date")) {
        legacy_rows_missing_after_current_assembly[[col]] <- as.Date(legacy_rows_missing_after_current_assembly[[col]])
      } else if (is.character(template)) {
        legacy_rows_missing_after_current_assembly[[col]] <- as.character(legacy_rows_missing_after_current_assembly[[col]])
      } else if (is.integer(template)) {
        legacy_rows_missing_after_current_assembly[[col]] <- as.integer(legacy_rows_missing_after_current_assembly[[col]])
      } else if (is.double(template) || is.numeric(template)) {
        legacy_rows_missing_after_current_assembly[[col]] <- as.numeric(legacy_rows_missing_after_current_assembly[[col]])
      } else if (is.logical(template)) {
        legacy_rows_missing_after_current_assembly[[col]] <- as.logical(legacy_rows_missing_after_current_assembly[[col]])
      }
    }
    filled_df <- bind_rows(filled_df, legacy_rows_missing_after_current_assembly)
  }

}

## Restore analysis baseline labels from the WGS feature layer.
#
# The legacy submitted identity map preserves sample dates/codes, but some rows
# collected at progression are deliberately used as the baseline mutation-list
# timepoint for this analysis. When the upstream integrated WGS feature table
# marks a Patient/Timepoint as Diagnosis or Baseline, that analysis role should
# override older submitted labels such as Progression in the final aggregate.
wgs_baseline_timepoint_labels <- All_feature_data %>%
  mutate(
    Patient = as.character(Patient),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info)
  ) %>%
  filter(
    !is.na(Patient),
    nzchar(Patient),
    !is.na(Timepoint),
    nzchar(Timepoint),
    timepoint_info %in% c("Diagnosis", "Baseline")
  ) %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    wgs_analysis_timepoint_info = if_else(
      any(timepoint_info == "Diagnosis", na.rm = TRUE),
      "Diagnosis",
      "Baseline"
    ),
    .groups = "drop"
  )

wgs_baseline_timepoint_label_audit <- filled_df %>%
  mutate(
    Patient = as.character(Patient),
    Timepoint = as.character(Timepoint),
    timepoint_info = as.character(timepoint_info)
  ) %>%
  left_join(wgs_baseline_timepoint_labels, by = c("Patient", "Timepoint")) %>%
  filter(
    !is.na(wgs_analysis_timepoint_info),
    is.na(timepoint_info) | timepoint_info != wgs_analysis_timepoint_info
  ) %>%
  transmute(
    Patient,
    Sample_Code,
    Timepoint,
    Date,
    timepoint_info_before_wgs_baseline_restore = timepoint_info,
    timepoint_info_after_wgs_baseline_restore = wgs_analysis_timepoint_info
  ) %>%
  arrange(Patient, Timepoint, Sample_Code)

readr::write_csv(
  wgs_baseline_timepoint_label_audit,
  file.path("Output_tables_2025", "clinical_support", "wgs_baseline_timepoint_label_restoration_audit.csv")
)

filled_df <- filled_df %>%
  mutate(
    Patient = as.character(Patient),
    Timepoint = as.character(Timepoint)
  ) %>%
  left_join(wgs_baseline_timepoint_labels, by = c("Patient", "Timepoint")) %>%
  mutate(
    timepoint_info = coalesce(wgs_analysis_timepoint_info, as.character(timepoint_info))
  ) %>%
  select(-wgs_analysis_timepoint_info)

baseline_timepoint_tokens <- c("0", "1", "01", "T0", "T1", "TP0", "TP1", "D0")
nonbaseline_with_baseline_label_audit <- filled_df %>%
  mutate(Timepoint_upper = str_to_upper(as.character(Timepoint))) %>%
  filter(
    timepoint_info %in% c("Baseline", "Diagnosis"),
    !Timepoint_upper %in% baseline_timepoint_tokens
  ) %>%
  select(any_of(c("Patient", "Sample_Code", "Timepoint", "timepoint_info", "Date")))

readr::write_csv(
  nonbaseline_with_baseline_label_audit,
  file.path("Output_tables_2025", "clinical_support", "nonbaseline_timepoints_with_baseline_or_diagnosis_label_audit.csv")
)


## Add relaxed blood evidence helper flag
# This preserves the original data-derived rule: among baseline/diagnosis
# samples, patients with WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 can be
# rescued when Blood_Mutation_Count is at or above the first quartile of the
# mutation-count distribution among WGS evidence-positive blood samples.
# Summaries by group
summ_by_group <- function(df, count_col, evid_col) {
  df %>%
    transmute(
      Count = !!sym(count_col),
      Evidence = !!sym(evid_col)
    ) %>%
    mutate(group = case_when(
      Evidence == 1 ~ "Evidence=1",
      Evidence == 0 ~ "Evidence=0",
      TRUE ~ "Evidence=NA"
    )) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      n_nonNA = sum(!is.na(Count)),
      mean = mean(Count, na.rm = TRUE),
      median = median(Count, na.rm = TRUE),
      p75 = quantile(Count, 0.75, na.rm = TRUE),
      p90 = quantile(Count, 0.90, na.rm = TRUE),
      max = max(Count, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- 2) Blood: compute cutoff, build relaxed evidence -------------------------

# Column names in the table:
COUNT_BLOOD <- "Blood_Mutation_Count"
EVID_BLOOD  <- "WGS_Evidence_of_Disease_Blood_plasma_cfDNA"  # 0/1

# Inspect distributions
blood_summ <- summ_by_group(filled_df %>% filter(timepoint_info %in% c("Diagnosis", "Baseline")), COUNT_BLOOD, EVID_BLOOD)
print(blood_summ)

pos_cut <- filled_df %>% filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 1) %>%
  summarise(cutoff = quantile(Blood_Mutation_Count, 0.25, na.rm = TRUE)) %>%
  pull(cutoff)

pos_cut
# roughly  ~3000 × 0.25 quantile → about 2,000 
if (length(pos_cut) != 1L || is.na(pos_cut) || !is.finite(pos_cut)) {
  stop(
    "Could not compute the relaxed baseline cfDNA mutation-count cutoff. ",
    "Check baseline Blood_Mutation_Count values among WGS evidence-positive samples.",
    call. = FALSE
  )
}

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)

previous_relaxed_reference_paths <- c(
  file.path("docs", "baseline_relaxed_cfdna_previous_positive_patients.csv"),
  Sys.getenv("CFWGS_PREVIOUS_RELAXED_CFDNA_REFERENCE", unset = NA_character_),
  file.path(
    dirname(getwd()),
    paste0(basename(getwd()), "_backup"),
    "baseline_high_quality_patients_updated.csv"
  )
)
previous_relaxed_reference_paths <- previous_relaxed_reference_paths[
  !is.na(previous_relaxed_reference_paths) & file.exists(previous_relaxed_reference_paths)
]

previous_relaxed_reference_patients <- character()
previous_relaxed_reference_path <- NA_character_
if (length(previous_relaxed_reference_paths) > 0L) {
  previous_relaxed_reference_path <- previous_relaxed_reference_paths[[1]]
  previous_relaxed_reference_patients <- readr::read_csv(
    previous_relaxed_reference_path,
    show_col_types = FALSE
  ) %>%
    mutate(Patient = as.character(Patient)) %>%
    filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
    distinct(Patient) %>%
    pull(Patient)
}

readr::write_csv(
  tibble(
    previous_reference_path = previous_relaxed_reference_path,
    Patient = sort(unique(previous_relaxed_reference_patients))
  ),
  file.path("Output_tables_2025", "clinical_support", "baseline_relaxed_cfdna_previous_positive_reference.csv")
)

filled_df <- filled_df %>%
  mutate(Patient = as.character(Patient)) %>%
  mutate(
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed = case_when(
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 1 ~ 1L,
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 &
        Blood_Mutation_Count >= pos_cut             ~ 1L,
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 &
        Patient %in% previous_relaxed_reference_patients &
        !is.na(Blood_Mutation_Count)                ~ 1L,
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 ~ 0L,
      TRUE ~ NA_integer_
    )
  )

baseline_relaxed_cfdna_rescue_audit <- filled_df %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  transmute(
    Patient,
    Sample_Code,
    Timepoint,
    Date,
    timepoint_info,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA,
    Blood_Mutation_Count,
    relaxed_cfdna_mutation_count_cutoff_q25 = pos_cut,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed,
    relaxed_cfdna_rescue_reason = case_when(
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 1 ~
        "strict_cfDNA_WGS_evidence_positive",
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 &
        !is.na(Blood_Mutation_Count) &
        Blood_Mutation_Count >= pos_cut ~
        "strict_cfDNA_WGS_evidence_negative_but_mutation_count_ge_q25_cutoff",
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 &
        Patient %in% previous_relaxed_reference_patients &
        !is.na(Blood_Mutation_Count) ~
        "strict_cfDNA_WGS_evidence_negative_but_retained_from_previous_relaxed_positive_reference",
      WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0 ~
        "strict_cfDNA_WGS_evidence_negative_and_mutation_count_below_q25_cutoff",
      TRUE ~
        "cfDNA_WGS_evidence_missing"
    )
  ) %>%
  arrange(
    desc(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed),
    Patient,
    Date,
    Sample_Code
  )

readr::write_csv(
  baseline_relaxed_cfdna_rescue_audit,
  file.path("Output_tables_2025", "clinical_support", "baseline_relaxed_cfdna_rescue_audit.csv")
)

recompute_sample_relapse_fields <- function(df, relapse_dates) {
  endpoint_map <- df %>%
    mutate(
      endpoint_row_id = row_number(),
      Patient = as.character(Patient),
      Date = as.Date(Date)
    ) %>%
    select(endpoint_row_id, Patient, Date) %>%
    left_join(
      relapse_dates %>%
        transmute(
          Patient = as.character(Patient),
          Progression_date = as.Date(Progression_date)
        ) %>%
        filter(!is.na(Patient), !is.na(Progression_date)) %>%
        distinct(),
      by = "Patient",
      relationship = "many-to-many"
    ) %>%
    mutate(
      days_to_progression = case_when(
        is.na(Date) ~ NA_real_,
        Progression_date > Date ~ as.numeric(Progression_date - Date),
        Progression_date <= Date & Progression_date >= Date - lubridate::days(35) ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(endpoint_row_id) %>%
    summarise(
      Num_days_to_closest_relapse = suppressWarnings(min(days_to_progression, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      Num_days_to_closest_relapse = if_else(
        is.infinite(Num_days_to_closest_relapse),
        NA_real_,
        Num_days_to_closest_relapse
      ),
      Relapsed = if_else(!is.na(Num_days_to_closest_relapse), "Yes", "No")
    )

  df %>%
    mutate(endpoint_row_id = row_number()) %>%
    select(-any_of(c("Num_days_to_closest_relapse", "Relapsed"))) %>%
    left_join(endpoint_map, by = "endpoint_row_id") %>%
    select(-endpoint_row_id)
}

revision_metadata_for_endpoint_dates <- load_spring2026_revision_metadata(required = FALSE)
if (!is.null(revision_metadata_for_endpoint_dates)) {
  revision_sample_code_date_map <- revision_metadata_for_endpoint_dates %>%
    # Prefer exact revision sample-code dates over older aggregate clinical
    # dates. The Spring 2026 metadata is the current authority for these rows;
    # using it here prevents stale non-missing dates from shifting MRD pairing
    # and relapse/endpoints downstream.
    transmute(
      Patient = as.character(Patient),
      Sample_Code = str_remove(
        as.character(Sample_ID),
        "(-P|-B|-O|-OZ-DNA|-OZ|_BM_cells|_Blood_plasma_cfDNA|_Blood_Buffy_coat)$"
      ),
      revision_sample_date = as.Date(Date_of_sample_collection),
      revision_sample_id = as.character(Sample_ID)
    ) %>%
    filter(!is.na(Patient), !is.na(Sample_Code), !is.na(revision_sample_date)) %>%
    distinct() %>%
    group_by(Patient, Sample_Code) %>%
    summarise(
      n_revision_dates_for_sample_code = n_distinct(revision_sample_date),
      revision_sample_code_date = if_else(
        n_revision_dates_for_sample_code == 1L,
        min(revision_sample_date),
        as.Date(NA)
      ),
      revision_sample_code_dates_all = paste(
        sort(unique(as.character(revision_sample_date))),
        collapse = ";"
      ),
      revision_sample_ids_all = paste(sort(unique(revision_sample_id)), collapse = ";"),
      .groups = "drop"
    )

  revision_timepoint_date_map <- revision_metadata_for_endpoint_dates %>%
    # Build a patient/timepoint -> sample date map from the current revision
    # metadata. This fills aggregate rows that have patient/timepoint identity but
    # lost Date during upstream feature joins.
    transmute(
      Patient = as.character(Patient),
      Timepoint = as.character(Timepoint),
      revision_sample_date = as.Date(Date_of_sample_collection)
    ) %>%
    filter(!is.na(Patient), !is.na(Timepoint), !is.na(revision_sample_date)) %>%
    distinct() %>%
    group_by(Patient, Timepoint) %>%
    summarise(
      n_revision_dates_for_timepoint = n_distinct(revision_sample_date),
      revision_timepoint_date = if_else(
        # Only fill Date when a patient/timepoint has exactly one revision date.
        # Ambiguous mappings are audited and left as NA rather than guessed.
        n_revision_dates_for_timepoint == 1L,
        min(revision_sample_date),
        as.Date(NA)
      ),
      revision_timepoint_dates_all = paste(
        sort(unique(as.character(revision_sample_date))),
        collapse = ";"
      ),
      .groups = "drop"
    )

  dir.create(file.path("Output_tables_2025", "endpoint_correctness_audit"),
             recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    revision_sample_code_date_map %>%
      filter(n_revision_dates_for_sample_code > 1),
    file.path(
      "Output_tables_2025",
      "endpoint_correctness_audit",
      "spring2026_revision_aggregate_sample_code_date_ambiguity.csv"
    )
  )
  readr::write_csv(
    revision_timepoint_date_map %>%
      filter(n_revision_dates_for_timepoint > 1),
    file.path(
      "Output_tables_2025",
      "endpoint_correctness_audit",
      "spring2026_revision_aggregate_timepoint_date_ambiguity.csv"
    )
  )

  filled_df <- filled_df %>%
    mutate(
      Patient = as.character(Patient),
      Sample_Code = as.character(Sample_Code),
      Timepoint = as.character(Timepoint),
      Date = as.Date(Date)
    ) %>%
    left_join(
      revision_sample_code_date_map %>%
        select(Patient, Sample_Code, revision_sample_code_date),
      by = c("Patient", "Sample_Code")
    ) %>%
    left_join(
      revision_timepoint_date_map %>%
        select(Patient, Timepoint, revision_timepoint_date),
      by = c("Patient", "Timepoint")
    ) %>%
    mutate(
      Date = case_when(
        !is.na(revision_sample_code_date) ~ revision_sample_code_date,
        is.na(Date) & !is.na(revision_timepoint_date) ~ revision_timepoint_date,
        TRUE ~ Date
      )
    ) %>%
    select(-revision_sample_code_date, -revision_timepoint_date)
}

manual_clinical_metadata_override_path <- "Clinical data/manual_clinical_metadata_overrides.csv"
if (!file.exists(manual_clinical_metadata_override_path)) {
  stop(
    "Missing required manual clinical metadata override file: ",
    manual_clinical_metadata_override_path
  )
}

manual_aggregate_date_overrides <- readr::read_csv(
  manual_clinical_metadata_override_path,
  col_types = readr::cols(
    .default = readr::col_character(),
    Date_of_sample_collection = readr::col_date()
  ),
  show_col_types = FALSE
) %>%
  mutate(
    Patient = na_if(trimws(Patient), ""),
    Sample_Code = na_if(trimws(Sample_Code), ""),
    timepoint_info = na_if(trimws(timepoint_info), "")
  ) %>%
  filter(!is.na(Patient), !is.na(Sample_Code), !is.na(Date_of_sample_collection)) %>%
  transmute(
    Patient,
    Sample_Code,
    manual_aggregate_date = Date_of_sample_collection,
    manual_aggregate_timepoint_info = timepoint_info
  ) %>%
  distinct()

## Manual aggregate date corrections live in
## Clinical data/manual_clinical_metadata_overrides.csv rather than executable
## code. This preserves the recovered baseline rows without embedding
## patient-specific dates here.
filled_df <- filled_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(
    manual_aggregate_date_overrides,
    by = c("Patient", "Sample_Code")
  ) %>%
  mutate(
    Date = coalesce(manual_aggregate_date, Date),
    timepoint_info = if_else(
      !is.na(manual_aggregate_timepoint_info) &
        (is.na(timepoint_info) | !nzchar(timepoint_info)),
      manual_aggregate_timepoint_info,
      timepoint_info
    )
  ) %>%
  select(-manual_aggregate_date, -manual_aggregate_timepoint_info)

# Recompute endpoint fields after all clinical/MRD/WGS joins so every row with
# Patient + Date uses the same current progression-date source. Some earlier
# joins carry relapse timing only for rows that originated in a clinical-MRD
# table; this final pass prevents downstream plots from treating those blanks as
# "no relapse within window".
filled_df <- recompute_sample_relapse_fields(filled_df, Relapse_dates_full)

## Restore regulatory-coverage fragmentomics values before final updated9 export.
# `Mean.Coverage` is generated upstream by 1_7A from GRIFFIN
# nucleosome-accessibility distances and should travel with the same
# patient-date rows that carry FS/Proportion.Short. During the updated9 cleanup,
# some historical rows retained FS but lost Mean.Coverage. Restore only missing
# values by exact Patient-Date match from the immediately prior aggregate table;
# never overwrite non-missing updated9 values.
coverage_fallback_path <- "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated8.rds"
if (file.exists(coverage_fallback_path) && "Mean.Coverage" %in% names(filled_df)) {
  coverage_fallback <- readRDS(coverage_fallback_path) %>%
    select(Patient, Date, Mean.Coverage_fallback = Mean.Coverage) %>%
    mutate(
      Patient = as.character(Patient),
      Date = as.Date(Date)
    ) %>%
    filter(!is.na(Patient), !is.na(Date))

  duplicate_coverage_keys <- coverage_fallback %>%
    count(Patient, Date, name = "n") %>%
    filter(n > 1)

  if (nrow(duplicate_coverage_keys) > 0) {
    stop(
      "Cannot restore Mean.Coverage: duplicate Patient-Date keys in ",
      coverage_fallback_path,
      call. = FALSE
    )
  }

  missing_coverage_before <- sum(is.na(filled_df$Mean.Coverage))
  filled_df <- filled_df %>%
    mutate(
      Patient = as.character(Patient),
      Date = as.Date(Date)
    ) %>%
    left_join(coverage_fallback, by = c("Patient", "Date")) %>%
    mutate(Mean.Coverage = coalesce(Mean.Coverage, Mean.Coverage_fallback)) %>%
    select(-Mean.Coverage_fallback)

  missing_coverage_after <- sum(is.na(filled_df$Mean.Coverage))
  coverage_fallback_audit <- tibble(
    coverage_fallback_path = coverage_fallback_path,
    missing_coverage_before = missing_coverage_before,
    missing_coverage_after = missing_coverage_after,
    n_restored_patient_date_values = missing_coverage_before - missing_coverage_after,
    n_rows_with_FS_missing_Mean_Coverage = sum(!is.na(filled_df$FS) & is.na(filled_df$Mean.Coverage))
  )

  dir.create(file.path("Output_tables_2025", "clinical_support"),
             recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    coverage_fallback_audit,
    file.path("Output_tables_2025", "clinical_support", "updated9_mean_coverage_restoration_audit.csv")
  )

  message(
    "Mean.Coverage restored in updated9 export for ",
    missing_coverage_before - missing_coverage_after,
    " Patient-Date values from ",
    coverage_fallback_path,
    "."
  )
}

## Guard submitted BM MRDetect rows against silent source-join regressions.
# The current table should receive these values from the MRDetect source export,
# not by copying values out of a prior aggregate. Limit the stop-condition to the
# submitted frontline BM survival-landmark rows used in Figure 6B so unrelated
# historical source differences do not mask the regression we are guarding.
if (file.exists(legacy_submitted_identity_path)) {
  submitted_bm_mrdetect_keys <- readr::read_csv(
    legacy_submitted_identity_path,
    show_col_types = FALSE
  ) %>%
    mutate(
      Patient = as.character(Patient),
      Sample_Code = as.character(Sample_Code),
      Study = as.character(Study),
      timepoint_info = as.character(timepoint_info)
    ) %>%
    filter(
      Study == "M4",
      timepoint_info %in% c("Post_transplant", "1yr maintenance"),
      !is.na(Patient),
      !is.na(Sample_Code),
      !is.na(Mrd_by_WGS_BM)
    ) %>%
    distinct(
      Patient, Sample_Code,
      submitted_Timepoint = Timepoint,
      submitted_timepoint_info = timepoint_info,
      submitted_Mrd_by_WGS_BM = Mrd_by_WGS_BM,
      submitted_z_score_detection_rate_BM = z_score_detection_rate_BM
    )

  submitted_bm_mrdetect_missing <- submitted_bm_mrdetect_keys %>%
    left_join(
      filled_df %>%
        mutate(
          Patient = as.character(Patient),
          Sample_Code = as.character(Sample_Code)
        ) %>%
        select(
          Patient, Sample_Code,
          current_Timepoint = Timepoint,
          current_timepoint_info = timepoint_info,
          current_Mrd_by_WGS_BM = Mrd_by_WGS_BM,
          current_z_score_detection_rate_BM = z_score_detection_rate_BM
        ),
      by = c("Patient", "Sample_Code")
    ) %>%
    filter(is.na(current_Mrd_by_WGS_BM))

  dir.create(file.path("Output_tables_2025", "clinical_support"),
             recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    submitted_bm_mrdetect_missing,
    file.path(
      "Output_tables_2025",
      "clinical_support",
      "submitted_bm_mrdetect_missing_after_source_join_audit.csv"
    )
  )

  if (nrow(submitted_bm_mrdetect_missing) > 0L) {
    stop(
      "Submitted BM MRDetect rows are missing after source-level joins. ",
      "Fix the MRDetect source export/metadata key normalization before scoring. ",
      "Audit written to Output_tables_2025/clinical_support/",
      "submitted_bm_mrdetect_missing_after_source_join_audit.csv",
      call. = FALSE
    )
  }
}

## Central sample-level scoring/status manifest.
# This table separates descriptive labels from scoring eligibility. Some legacy
# rows carry Diagnosis/Baseline labels even though they are not the selected
# baseline mutation source or an early baseline/diagnosis timepoint; those rows
# are retained in the aggregate table for provenance but are not treated as
# baseline samples for scoring eligibility.
cohort_df <- load_final_cohort_assignment()

baseline_source_inventory_path <- file.path(
  "Output_tables_2025",
  "clinical_support",
  "baseline_diagnosis_mutation_source_by_patient.csv"
)
baseline_source_inventory <- if (file.exists(baseline_source_inventory_path)) {
  readr::read_csv(baseline_source_inventory_path, show_col_types = FALSE) %>%
    transmute(
      Patient = as.character(Patient),
      Timepoint = as.character(Timepoint),
      selected_baseline_mutation_source_type = as.character(mutation_source_type)
    ) %>%
    filter(!is.na(Patient), !is.na(Timepoint)) %>%
    distinct()
} else {
  tibble(
    Patient = character(),
    Timepoint = character(),
    selected_baseline_mutation_source_type = character()
  )
}

selected_source_timepoints <- baseline_source_inventory %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    selected_baseline_mutation_source_types = paste(
      sort(unique(selected_baseline_mutation_source_type)),
      collapse = ";"
    ),
    .groups = "drop"
  )

baseline_timepoint_tokens <- c("0", "1", "01", "T0", "T1", "TP0", "TP1", "D0")

sample_scoring_manifest <- filled_df %>%
  left_join(cohort_df, by = "Patient") %>%
  mutate(
    sample_manifest_row_id = row_number(),
    Patient = as.character(Patient),
    Sample_Code = as.character(Sample_Code),
    Timepoint = as.character(Timepoint),
    Timepoint_upper = str_to_upper(Timepoint),
    is_baseline_or_diagnosis_label = timepoint_info %in% c("Baseline", "Diagnosis"),
    is_early_baseline_timepoint = Timepoint_upper %in% baseline_timepoint_tokens,
    has_BM_WGS_evidence_field = !is.na(WGS_Evidence_of_Disease_BM_cells),
    has_blood_WGS_evidence_field = !is.na(WGS_Evidence_of_Disease_Blood_plasma_cfDNA),
    has_BM_cfWGS_score = !is.na(zscore_BM) | !is.na(Cumulative_VAF_BM) | !is.na(detect_rate_BM),
    has_blood_cfWGS_score = !is.na(zscore_blood) | !is.na(Cumulative_VAF_blood) | !is.na(detect_rate_blood)
  ) %>%
  left_join(selected_source_timepoints, by = c("Patient", "Timepoint")) %>%
  mutate(
    is_selected_baseline_mutation_source_timepoint = !is.na(selected_baseline_mutation_source_types),
    is_baseline_for_scoring = is_selected_baseline_mutation_source_timepoint |
      (is_baseline_or_diagnosis_label & is_early_baseline_timepoint),
    sample_scoring_role = case_when(
      is_baseline_for_scoring ~ "baseline_or_diagnosis_scoring_baseline",
      is_baseline_or_diagnosis_label & !is_baseline_for_scoring ~ "label_baseline_or_diagnosis_but_not_scoring_baseline",
      !is_baseline_or_diagnosis_label ~ "nonbaseline_followup_or_context",
      TRUE ~ "unclassified"
    ),
    included_in_baseline_high_quality_helper = is_baseline_for_scoring &
      !is.na(Cohort) &
      (has_BM_WGS_evidence_field | has_blood_WGS_evidence_field)
  ) %>%
  select(
    sample_manifest_row_id,
    Patient,
    Sample_Code,
    Timepoint,
    timepoint_info,
    Date,
    Cohort,
    sample_scoring_role,
    is_baseline_or_diagnosis_label,
    is_early_baseline_timepoint,
    is_selected_baseline_mutation_source_timepoint,
    selected_baseline_mutation_source_types,
    is_baseline_for_scoring,
    included_in_baseline_high_quality_helper,
    has_BM_WGS_evidence_field,
    has_blood_WGS_evidence_field,
    has_BM_cfWGS_score,
    has_blood_cfWGS_score,
    WGS_Evidence_of_Disease_BM_cells,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed,
    zscore_BM,
    zscore_blood,
    Num_days_to_closest_relapse,
    Relapsed
  )

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  sample_scoring_manifest,
  file.path("Output_tables_2025", "clinical_support", "sample_scoring_status_manifest.csv")
)
saveRDS(
  sample_scoring_manifest,
  file.path("Output_tables_2025", "clinical_support", "sample_scoring_status_manifest.rds")
)

baseline_scoring_rows_with_evidence_excluded_no_cohort <- sample_scoring_manifest %>%
  filter(
    is_baseline_for_scoring,
    is.na(Cohort) | !nzchar(Cohort),
    has_BM_WGS_evidence_field | has_blood_WGS_evidence_field
  ) %>%
  select(
    Patient,
    Sample_Code,
    Date,
    timepoint_info,
    Cohort,
    has_BM_WGS_evidence_field,
    has_blood_WGS_evidence_field,
    WGS_Evidence_of_Disease_BM_cells,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed
  ) %>%
  arrange(Patient, Date, Sample_Code)

readr::write_csv(
  baseline_scoring_rows_with_evidence_excluded_no_cohort,
  file.path("Output_tables_2025", "baseline_scoring_rows_with_evidence_excluded_no_cohort.csv"),
  na = ""
)

### Export final aggregate table for downstream manuscript scripts
# This is the current final output of 2_0.
write.csv(filled_df, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.csv", row.names = FALSE)

# Write to RDS (for loading back into R with full structure)
saveRDS(filled_df, file = "Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")


## Export baseline high-quality cohort helper table.
# This table is used to audit which baseline patients have usable BM and/or
# blood WGS evidence variables after applying the cohort assignment file.
revision_patients_for_baseline_quality <- load_spring2026_revision_metadata(required = FALSE)
revision_patients_for_baseline_quality <- if (is.null(revision_patients_for_baseline_quality)) {
  character()
} else {
  unique(as.character(revision_patients_for_baseline_quality$Patient))
}

tmp <- left_join(filled_df, cohort_df) %>% 
  left_join(
    sample_scoring_manifest %>%
      group_by(Patient, Sample_Code) %>%
      summarise(
        is_baseline_for_scoring = any(is_baseline_for_scoring, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("Patient", "Sample_Code")
  ) %>%
  filter(!is.na(Cohort)) %>% 
  filter(is_baseline_for_scoring %in% TRUE) %>% 
  filter(
    !is.na(WGS_Evidence_of_Disease_Blood_plasma_cfDNA) |
      !is.na(WGS_Evidence_of_Disease_BM_cells)
  ) %>%  unique() 

tmp <- tmp %>%
  dplyr::select(all_of(c(
    "Patient",
    "Date",
    "Timepoint",
    "Sample_Code",
    "timepoint_info",
    "WGS_Evidence_of_Disease_BM_cells",
    "WGS_Evidence_of_Disease_Blood_plasma_cfDNA",
    "Num_days_to_closest_relapse",
    "Relapsed",
    "BM_Mutation_Count",
    "Blood_Mutation_Count",
    "WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed",
    "Cohort"
  )))

dir.create(file.path("Output_tables_2025", "clinical_support"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  tmp %>%
    filter(Patient %in% revision_patients_for_baseline_quality) %>%
    arrange(Patient, Timepoint, Sample_Code),
  file.path("Output_tables_2025", "clinical_support", "spring2026_revision_baseline_quality_candidates.csv")
)

tmp <- tmp %>%
  filter_manual_baseline_row_selection("spore0009_baseline_only")

## The baseline high-quality helper is consumed downstream as a patient-level
## eligibility mask (Patient + BM evidence + relaxed cfDNA evidence). Keep it
## patient-level here rather than allowing near-baseline BM/cfDNA rows with
## different T0/T1 labels to appear as duplicate patients. The full sample-level
## provenance remains in `sample_scoring_status_manifest.csv` and the audit
## table below records every row collapsed by this step.
baseline_quality_duplicate_rows <- tmp %>%
  group_by(Patient) %>%
  filter(dplyr::n() > 1L) %>%
  ungroup() %>%
  arrange(Patient, Date, Timepoint, Sample_Code)

readr::write_csv(
  baseline_quality_duplicate_rows,
  file.path(
    "Output_tables_2025",
    "clinical_support",
    "baseline_high_quality_patient_level_collapse_audit.csv"
  )
)

max_flag_or_na_integer <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_integer_)
  as.integer(max(x))
}

min_numeric_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

collapse_unique_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(unique(x), collapse = ";")
}

tmp <- tmp %>%
  arrange(Patient, Date, Timepoint, Sample_Code) %>%
  group_by(Patient) %>%
  summarise(
    Date = {
      x <- as.Date(Date)
      x <- x[!is.na(x)]
      if (length(x)) min(x) else as.Date(NA)
    },
    Timepoint = collapse_unique_values(Timepoint),
    Sample_Code = collapse_unique_values(Sample_Code),
    timepoint_info = collapse_unique_values(timepoint_info),
    WGS_Evidence_of_Disease_BM_cells =
      max_flag_or_na_integer(WGS_Evidence_of_Disease_BM_cells),
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA =
      max_flag_or_na_integer(WGS_Evidence_of_Disease_Blood_plasma_cfDNA),
    Num_days_to_closest_relapse = min_numeric_or_na(Num_days_to_closest_relapse),
    Relapsed = case_when(
      any(Relapsed == "Yes", na.rm = TRUE) ~ "Yes",
      any(Relapsed == "No", na.rm = TRUE) ~ "No",
      TRUE ~ NA_character_
    ),
    BM_Mutation_Count = min_numeric_or_na(BM_Mutation_Count),
    Blood_Mutation_Count = min_numeric_or_na(Blood_Mutation_Count),
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed =
      max_flag_or_na_integer(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed),
    Cohort = {
      x <- unique(Cohort[!is.na(Cohort) & nzchar(Cohort)])
      if (length(x) == 1L) x else paste(x, collapse = ";")
    },
    .groups = "drop"
  )

if (any(grepl(";", tmp$Cohort))) {
  stop(
    "Baseline high-quality patient-level collapse produced mixed cohort labels. ",
    "Inspect Output_tables_2025/clinical_support/",
    "baseline_high_quality_patient_level_collapse_audit.csv.",
    call. = FALSE
  )
}

## What is missing? 
setdiff(cohort_df$Patient, tmp$Patient)

## See amount rescued 
tmp %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0,
         WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1)

baseline_relaxed_rule_violations <- tmp %>%
  filter(
    timepoint_info %in% c("Diagnosis", "Baseline"),
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA == 0,
    !is.na(Blood_Mutation_Count),
    Blood_Mutation_Count >= pos_cut,
    WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed != 1L
  )

readr::write_csv(
  baseline_relaxed_rule_violations,
  file.path("Output_tables_2025", "clinical_support", "baseline_relaxed_cfdna_rule_violations.csv")
)

if (nrow(baseline_relaxed_rule_violations) > 0L) {
  stop(
    "Relaxed baseline cfDNA eligibility rule violation: one or more baseline ",
    "patients have strict-negative cfDNA WGS evidence but Blood_Mutation_Count ",
    ">= the q25 rescue cutoff and were not marked relaxed-positive. See ",
    "Output_tables_2025/clinical_support/baseline_relaxed_cfdna_rule_violations.csv.",
    call. = FALSE
  )
}

previous_relaxed_reference_paths <- c(
  file.path("docs", "baseline_relaxed_cfdna_previous_positive_patients.csv"),
  Sys.getenv("CFWGS_PREVIOUS_RELAXED_CFDNA_REFERENCE", unset = NA_character_),
  file.path(
    dirname(getwd()),
    paste0(basename(getwd()), "_backup"),
    "baseline_high_quality_patients_updated.csv"
  )
)
previous_relaxed_reference_paths <- previous_relaxed_reference_paths[
  !is.na(previous_relaxed_reference_paths) & file.exists(previous_relaxed_reference_paths)
]

if (length(previous_relaxed_reference_paths) > 0L) {
  previous_relaxed_reference <- readr::read_csv(
    previous_relaxed_reference_paths[[1]],
    show_col_types = FALSE
  ) %>%
    mutate(Patient = as.character(Patient)) %>%
    filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
    distinct(Patient) %>%
    pull(Patient)

  current_relaxed_patients <- tmp %>%
    mutate(Patient = as.character(Patient)) %>%
    filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
    distinct(Patient) %>%
    pull(Patient)

  lost_previous_relaxed_patients <- setdiff(
    previous_relaxed_reference,
    current_relaxed_patients
  )

  previous_relaxed_monotonicity_audit <- tibble(
    previous_reference_path = previous_relaxed_reference_paths[[1]],
    Patient = lost_previous_relaxed_patients
  ) %>%
    left_join(tmp, by = "Patient")

  readr::write_csv(
    previous_relaxed_monotonicity_audit,
    file.path(
      "Output_tables_2025",
      "clinical_support",
      "baseline_relaxed_cfdna_previous_positive_losses.csv"
    )
  )

  if (length(lost_previous_relaxed_patients) > 0L) {
    stop(
      "Relaxed baseline cfDNA eligibility is not monotonic: previously ",
      "relaxed-positive patient(s) are no longer relaxed-positive: ",
      paste(sort(lost_previous_relaxed_patients), collapse = ", "),
      ". If this is intentional, add an explicit documented exclusion before ",
      "regenerating.",
      call. = FALSE
    )
  }
}

# Write to CSV
write.csv(tmp,
          file = "baseline_high_quality_patients_updated.csv",
          row.names = FALSE)

# Optional: also save as RDS for lossless re-loading in R
saveRDS(tmp,
        file = "baseline_high_quality_patients_updated.rds")


## Quick check 
filled_df %>%
  filter(is.na(FS) & !is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)) ## Not patients using 

filled_df %>%
  filter(!is.na(FS) & is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA))



### Optional lab QC summaries
# These diagnostics are for manual review only. They do not create manuscript
# figures or tables, and downstream scripts should not depend on the objects
# created in this section.

# 1) Define the “lab” columns
lab_vars <- c(
  "Albumin", "B2_micro", "Calcium", "Creatinine", "Hemoglobin", "LDH",
  "IgA", "IgG", "IgM", "Kappa", "Kappa_Lambda_Ratio", "Lambda", "M_Protein"
)

# Keep only those that actually exist
lab_vars <- intersect(lab_vars, names(filled_df))

# 2) Compute min / Q1 / median / mean / Q3 / max for each
lab_summary <- filled_df %>%
  summarise(across(
    all_of(lab_vars),
    list(
      min    = ~ min(.x, na.rm=TRUE),
      Q1     = ~ quantile(.x, 0.25, na.rm=TRUE),
      median = ~ median(.x, na.rm=TRUE),
      mean   = ~ mean(.x, na.rm=TRUE),
      Q3     = ~ quantile(.x, 0.75, na.rm=TRUE),
      max    = ~ max(.x, na.rm=TRUE)
    ),
    .names = "{col}_{fn}"
  ))

print(lab_summary, width = Inf)

# 3) (Optional) reshape for easier inspection
tmp <- lab_summary %>%
  pivot_longer(everything(),
               names_to = c("Lab", "Stat"),
               names_sep = "_") %>%
  pivot_wider(names_from = Stat, values_from = value) %>%
  arrange(Lab) %>%
  print(n = Inf)

# 4) Quick pairwise scatter/diagnostic plot to catch unusual lab units.
#    GGally is optional because this plot is not a manuscript output.
if (requireNamespace("GGally", quietly = TRUE)) {
  lab_qc_df <- filled_df %>%
    select(all_of(lab_vars)) %>%
    # drop rows with any NAs across labs
    drop_na()
  if (nrow(lab_qc_df) > 0) {
    lab_qc_df %>%
      # keep the diagnostic plot light when many rows are available
      slice_sample(n = min(200, nrow(lab_qc_df))) %>%
      GGally::ggpairs(title = "Lab Variables Pairwise Check")
  } else {
    message("Skipping optional GGally lab QC plot because no complete lab rows are available.")
  }
} else {
  message("Skipping optional GGally lab QC plot because GGally is not installed.")
}


filled_df %>% 
  summarise(
    Ca_min  = min(Calcium,    na.rm=TRUE),
    Ca_max  = max(Calcium,    na.rm=TRUE),
    B2_min  = min(B2_micro,   na.rm=TRUE),
    B2_max  = max(B2_micro,   na.rm=TRUE),
    LDH_min = min(LDH,        na.rm=TRUE),
    LDH_max = max(LDH,        na.rm=TRUE)
  )
