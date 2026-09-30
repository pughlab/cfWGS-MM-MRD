# ==============================================================================
# 1_1B_Process_clinical_labs.R
#
# Script purpose:
#   Compile and clean laboratory measurements, staging markers, demographic
#   fields, and FISH/cytogenetic annotations across M4, SPORE, and IMMAGINE.
#   The output is an active upstream clinical table used by
#   2_0_Assemble_Table_With_All_Features.R and, downstream of that, Table 1 and
#   clinical-feature summaries.
#
# Inputs:
#   M4 demographics, chemotherapy, marrow-biopsy, staging, and laboratory
#   workbooks; the sample-level clinical table from 1_0; curated IMMAGINE and
#   SPORE baseline clinical workbooks; and the cohort FISH helper tables from
#   1_1A. The complete executable list is validated in `input_files` below.
#
# Active downstream output:
#   • Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv
#     This is not a final manuscript table. It is an intermediate clinical
#     helper table consumed by script 2_0.
# Support/QC output:
#   • Output_tables_2025/clinical_support/lab_value_plausibility_audit.csv
#     Lists parsed laboratory values excluded by the explicit plausibility
#     ranges below; it is a review aid, not a manuscript table.
#
# How to run from the repository root:
#   Rscript 1_1B_Process_clinical_labs.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes clinical laboratory
#     measurements that feed Table 1 and downstream clinical-feature summaries.
#
# Author:    Dory Abelman
# Date:      May 2025
# ==============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table. It creates the cross-cohort clinical/lab
#   helper table consumed by 2_0, which then feeds manuscript figures/tables.
#
# Script outline and analysis unit:
#   1. Build M4 demographic and cytogenetic fields.
#   2. Match M4 biopsy/FISH and laboratory records to clinical timepoints.
#   3. Parse laboratory values, quarantine implausible values, and reshape the
#      selected measurements to patient-timepoint rows.
#   4. Derive a heuristic immunoglobulin/light-chain subtype.
#   5. Harmonize curated SPORE and IMMAGINE baseline clinical/FISH fields.
#   6. Export the cross-cohort clinical/lab helper table used by script 2_0.
#   The working unit is a clinical-event or patient-timepoint staging row.
#   Distinct relapse events and unmatched full-join rows may remain separate;
#   script 2_0 performs the final patient/timepoint consolidation.
#
# Matching and interpretation rules:
#   - Date-based FISH/biopsy and laboratory joins allow records from the same
#     patient within +/-60 days of the clinical sample date.
#   - Laboratory matches retain the smallest absolute date difference. The M4
#     FISH block first consolidates candidate rows by Patient and INTENT using
#     the first available values, so it does not always select the closest date.
#   - Duplicate laboratory values on the same patient/date/test are averaged
#     when all are zero or all are nonzero; when zero and nonzero values coexist,
#     the nonzero maximum is retained.
#   - Values outside the explicit `lab_plausibility_ranges` table are written to
#     the plausibility audit and set to missing in the analysis table.
#   - Immunoglobulin/light-chain subtype labels are deterministic helper
#     annotations derived from available lab values, not adjudicated clinical
#     subtype assignments.
#

# Load only the packages actually used below
library(dplyr)       # data manipulation
library(readxl)      # Excel I/O
library(stringr)     # string processing
library(fuzzyjoin)   # fuzzy / date-slack joins
library(data.table)  # fast grouping & reshaping
library(tidyverse)

ensure_exists <- function(path, producer = NULL) {
  if (file.exists(path)) return(invisible(path))
  msg <- paste0("Missing required input file: ", path)
  if (!is.null(producer)) msg <- paste0(msg, ". Run ", producer, " first.")
  stop(msg, call. = FALSE)
}

combined_clinical_data_path <- "combined_clinical_data_updated_April2025.csv"
input_files <- c(
  combined_clinical_data_path,
  "M4_CMRG_Data/M4_COHORT_DEMO.xlsx",
  "M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx",
  "M4_CMRG_Data/M4_COHORT_BONE_MARROW_BIOPSY.xlsx",
  "M4_CMRG_Data/M4_COHORT_STAGING.xlsx",
  "M4_CMRG_Data/M4_COHORT_LABS.xlsx",
  "Clinical data/IMMAGINE/Clinical data for IMG patients at diagnosis_filled_DA_edited.xlsx",
  "Clinical data/SPORE/Spore_baseline_clinical_demographics_DA_edited.xlsx",
  "Clinical data/SPORE/Spore patient timeline clinical info sent to Toronto.xlsx",
  "Clinical data/SPORE/tidy_fish.csv",
  "Clinical data/Exported clinical data April 2025/IMMAGINE_fish_flags.csv"
)
purrr::walk(input_files, ensure_exists)
message("All required clinical/lab inputs are present.")

write_active_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  message("Active clinical/lab input written: ", path)
  invisible(path)
}

combined_clinical_data_updated <- readr::read_csv(
  combined_clinical_data_path,
  show_col_types = FALSE
) %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

### Upstream clinical demographics fields used later by Table 1
## This script creates the lab/staging/cytogenetic fields. The formatted Table 1
## itself is created downstream in 2_1_Clinical_Demographics_Table.R.

##### First build the M4 clinical/lab table - updated May 2025
#### SPORE and IMMAGINE are harmonized into this structure later in the script.

# Set the path to the data folder
data_path <- "M4_CMRG_Data/"

# Step 1: Load the data files
# Load Demographics Data
demo_file <- paste0(data_path, "M4_COHORT_DEMO.xlsx")
demographics <- read_excel(demo_file)

# Load Chemotherapy Data
chemo_file <- paste0(data_path, "M4_COHORT_CHEMOTHERAPY.xlsx")
chemotherapy <- read_excel(chemo_file)

# Load Bone Marrow Biopsy Data
bmb_file <- paste0(data_path, "M4_COHORT_BONE_MARROW_BIOPSY.xlsx")
bone_marrow <- read_excel(bmb_file)

# Load Staging Data
staging_file <- paste0(data_path, "M4_COHORT_STAGING.xlsx")
staging <- read_excel(staging_file)

# Step 2: Data Refinement and Preprocessing

# 2.1 Convert date columns to Date type in demographics data
date_columns_demo <- c("DOB", "CONSENT_DATE", "DATE_OF_DEATH", "DATE_OF_LAST_FOLLOWUP")
demographics[date_columns_demo] <- lapply(demographics[date_columns_demo], as.Date, format = "%d-%b-%y")

# 2.2 Calculate Age at Consent
demographics <- demographics %>%
  mutate(AGE_AT_CONSENT = floor(as.numeric(difftime(CONSENT_DATE, DOB, units = "days")) / 365.25))

# 2.3 Categorize Age Groups
demographics <- demographics %>%
  mutate(AGE_GROUP = case_when(
    AGE_AT_CONSENT < 50 ~ "<50",
    AGE_AT_CONSENT >= 50 & AGE_AT_CONSENT < 60 ~ "50-59",
    AGE_AT_CONSENT >= 60 & AGE_AT_CONSENT < 70 ~ "60-69",
    AGE_AT_CONSENT >= 70 & AGE_AT_CONSENT < 80 ~ "70-79",
    AGE_AT_CONSENT >= 80 ~ "80+",
    TRUE ~ "Unknown"
  ))

# 2.4 Process Translocation Data from Bone Marrow Biopsy
# Select relevant columns
translocations_old <- bone_marrow %>%
  select(M4_id, INTENT, PROCEDURE_DATE, contains("T_"), contains("DEL_"), contains("AMP_"))  %>%
  select(-study_patient_id) %>%
  #  filter(str_detect(INTENT, regex("diagnosis", ignore_case = TRUE))) %>% # Keep only diagnosis
  mutate(across(
    contains("T_") | contains("DEL_") | contains("AMP_"),
    ~ifelse(. %in% c("Yes", "Positive"), "Positive",
            ifelse(. %in% c("No", "Negative"), "Negative", "Unknown"))
  ))


## Take more
translocations <- bone_marrow %>%
  # 1) pick only the cols needed for this analysis
  select(
    M4_id,
    INTENT,
    PROCEDURE_DATE,
    CHROMOSOME_13_DEL,
    T_4_14, T_11_14, T_14_16,
    DEL_17P, DEL_1P,
    AMP_1Q,
    HAS_NUMERICAL_ABNORMALITY, HAS_TRANSLOCATION, HAS_TRISOMY, HAS_MONOSOMY,
    TRISOMY_CHR3, TRISOMY_CHR7, TRISOMY_CHR9, TRISOMY_CHR15
  ) %>%
  # 2) rename the 13q deletion to match the other DEL_ columns
  rename(
    DEL_13 = CHROMOSOME_13_DEL
  ) %>%
  # 3) recode every “Y/Yes/Positive” → “Positive”, “N/No/Negative” → “Negative”, else “Unknown”
  mutate(
    across(
      c(
        # all of the existing T_, DEL_, AMP_ cols
        starts_with("T_"),
        starts_with("DEL_"),
        starts_with("AMP_"),
        # the HAS_ flags
        starts_with("HAS_"),
        # the four chromosome-trisomy columns
        TRISOMY_CHR3, TRISOMY_CHR7, TRISOMY_CHR9, TRISOMY_CHR15
      ),
      ~ case_when(
        . %in% c("Y", "Yes", "Positive") ~ "Positive",
        . %in% c("N", "No", "Negative")  ~ "Negative",
        TRUE                             ~ "Unknown"
      )
    )
  )

# Combine translocation statuses into a summary column
translocations_annotated <- translocations_old %>%
  mutate(
    TRANSLOCATIONS = paste(
      ifelse(T_4_14 == "Positive", "t(4;14)", ""),
      ifelse(T_11_14 == "Positive", "t(11;14)", ""),
      ifelse(T_14_16 == "Positive", "t(14;16)", ""),
      ifelse(DEL_17P == "Positive", "del(17p)", ""),
      ifelse(DEL_1P == "Positive", "del(1p)", ""),
      ifelse(AMP_1Q == "Positive", "amp(1q)", ""),
      sep = ", "
    )
  ) %>%
  mutate(TRANSLOCATIONS = gsub(", $", "", gsub("^, ", "", TRANSLOCATIONS))) %>%
  select(M4_id, TRANSLOCATIONS)

## Clean up 
# Clean up the TRANSLOCATIONS column
translocations_annotated <- translocations_annotated %>%
  mutate(
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",\\s*", ", "), # Remove extra spaces after commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "^,\\s*", ""),  # Remove leading commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",\\s*$", ""),  # Remove trailing commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ", ,", ","),    # Remove double commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",\\s*,", ","), # Remove empty entries between commas
    TRANSLOCATIONS = ifelse(TRANSLOCATIONS == "", "None", TRANSLOCATIONS) # Replace empty strings with 'None'
  )

translocations_annotated <- translocations_annotated %>%
  mutate(
    TRANSLOCATIONS = str_squish(TRANSLOCATIONS),                        # Remove extra spaces
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "\\s*,\\s*", ", "), # Fix misplaced commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "\\s+", " "),      # Ensure single spaces between terms
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",$", ""),         # Remove trailing commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "^,", ""),         # Remove leading commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "\\s(?=\\w)", ", "), # Add commas between terms missing them
    TRANSLOCATIONS = ifelse(TRANSLOCATIONS == "" | TRANSLOCATIONS == " ", "None", TRANSLOCATIONS) # Replace empty or single space with 'None'
  )


translocations_annotated <- translocations_annotated %>%
  mutate(
    TRANSLOCATIONS = str_squish(TRANSLOCATIONS),                        # Remove extra spaces
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",\\s*,", ", "),   # Replace consecutive commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, "^,", ""),         # Remove leading commas
    TRANSLOCATIONS = str_replace_all(TRANSLOCATIONS, ",$", ""),         # Remove trailing commas
    TRANSLOCATIONS = ifelse(TRANSLOCATIONS == "" | TRANSLOCATIONS == " ", "None", TRANSLOCATIONS) # Replace empty strings with 'None'
  )


### Lastly classify the translocation risk 
# Classify patients based on cytogenetic abnormalities
classified_translocations <- translocations %>%
  mutate(
    Cytogenetic_Risk = case_when(
      # High Risk: At least one high-risk abnormality is "Positive"
      T_4_14 == "Positive" | T_14_16 == "Positive" | DEL_17P == "Positive" ~ "High Risk",
      
      # Standard Risk: All high-risk abnormalities are "Negative"
      T_4_14 == "Negative" & T_14_16 == "Negative" & DEL_17P == "Negative" ~ "Standard Risk",
      
      # Unknown/Missing: Insufficient or unknown data
      TRUE ~ "Unknown/Missing"
    )
  )


## and add counts 
# Step 1: Count the number of "Positive" abnormalities for all relevant columns
classified_translocations <- classified_translocations %>%
  mutate(
    Abnormality_Count = rowSums(select(., T_4_14, T_11_14, T_14_16, DEL_17P, DEL_1P, AMP_1Q) == "Positive", na.rm = TRUE),
    Abnormality_Status = case_when(
      Abnormality_Count == 0 ~ "No cytogenetic abnormality",
      Abnormality_Count == 1 ~ "1 cytogenetic abnormality",
      Abnormality_Count == 2 ~ "2 cytogenetic abnormalities",
      Abnormality_Count == 3 ~ "3 cytogenetic abnormalities",
      Abnormality_Count > 3 ~ ">3 cytogenetic abnormalities",
      TRUE ~ "Unknown/missing"  # Catch any unexpected or insufficient data
    )
  )


classified_translocations <- classified_translocations %>% filter(!is.na(PROCEDURE_DATE))


### Now get the timepoints from the dates 
library(fuzzyjoin)

# 1) Prep the translocation table
ct2 <- classified_translocations %>%
  rename(Patient = M4_id) %>%
  mutate(PROCEDURE_DATE = as.Date(PROCEDURE_DATE))

# 2) Prep the clinical samples table
ccu_small <- combined_clinical_data_updated %>%
  select(
    Patient,
    Date         = Date_of_sample_collection,
    Timepoint,
    timepoint_info
  ) %>%
  mutate(Date = as.Date(Date)) %>% unique()

# 3) Fuzzy join: match Patient exactly and dates within +/-60 days.
joined2 <- fuzzy_left_join(
  ct2, ccu_small,
  by = c(
    "Patient"        = "Patient",
    "PROCEDURE_DATE" = "Date"
  ),
  match_fun = list(
    `==`,                                # exact match on Patient
    function(d1, d2) abs(as.numeric(d1 - d2)) <= 60 # 60-day window permits delayed clinical testing
  )
)

# 4) Preserve the historical consolidation order. The next block first
# collapses candidate matches by Patient + INTENT using the first available
# date/value. Consequently, this section does not guarantee that the closest
# clinical sample is chosen when one INTENT record has multiple +/-60-day
# candidates; `best_match` below only resolves rows remaining after collapse.
joined2 <- joined2 %>%
  rename(Patient = Patient.x) %>%
  select(-Patient.y)

### Clean duplicate rows 
# Deduplicate by Patient + INTENT,
# choosing the non-"Unknown" value for each translocation column when available
trans_cols <- c("T_4_14", "T_11_14", "T_14_16", "DEL_17P", "DEL_1P", "AMP_1Q")

joined2_clean <- joined2 %>%
  group_by(Patient, INTENT) %>%
  summarise(
    # Dates: pick the first non-NA
    PROCEDURE_DATE   = coalesce(first(PROCEDURE_DATE[!is.na(PROCEDURE_DATE)]), first(PROCEDURE_DATE)),
    Date             = coalesce(first(Date[!is.na(Date)]), first(Date)),
    Timepoint        = coalesce(first(Timepoint[!is.na(Timepoint)]), first(Timepoint)),
    timepoint_info   = coalesce(first(timepoint_info[!is.na(timepoint_info)]), first(timepoint_info)),
    
    # For each translocation column, prefer any non-"Unknown" value
    across(all_of(trans_cols),
           ~ {
             vals <- unique(.x)
             known <- vals[vals != "Unknown"]
             if (length(known)) known[1] else vals[1]
           }
    ),
    
    # Cytogenetic_Risk: same logic
    Cytogenetic_Risk = {
      vals <- unique(Cytogenetic_Risk)
      known <- vals[vals != "Unknown/Missing"]
      if (length(known)) known[1] else vals[1]
    },
    
    # Abnormality_Count: take the maximum
    Abnormality_Count = max(Abnormality_Count, na.rm = TRUE),
    
    # Abnormality_Status: prefer anything not "No cytogenetic abnormality"
    Abnormality_Status = {
      vals <- unique(Abnormality_Status)
      abnormal <- vals[vals != "No cytogenetic abnormality"]
      if (length(abnormal)) abnormal[1] else vals[1]
    },
    
    .groups = "drop"
  )

best_match <- joined2_clean %>%
  mutate(days_diff = abs(as.numeric(PROCEDURE_DATE - Date))) %>%
  group_by(Patient, PROCEDURE_DATE) %>%
  slice_min(order_by = days_diff, n = 1, with_ties = FALSE) %>%
  ungroup()


## Clean up the remaining timepoints based on the INTENT info 
best_match <- best_match %>%
  # 1) Extract any numeric visit number (1 or 2 digits) after the word “visit”
  mutate(
    visit_num = str_match(
      INTENT,
      regex("visit[^0-9]*([0-9]{1,2})", ignore_case = TRUE)
    )[,2]
  ) %>%
  # 2) Build a new Timepoint, filling in only where it was NA
  mutate(
    Timepoint = case_when(
      # keep existing if not missing
      !is.na(Timepoint) ~ Timepoint,
      # “Diagnosis” or “first consult” → 01
      str_detect(INTENT, regex("diagnosis|first consult", ignore_case = TRUE)) ~ "01",
      # “relapse” or visit R → R
      str_detect(INTENT, regex("relapse|visit[^0-9]*R", ignore_case = TRUE)) ~ "R-",
      # any numeric visit → zero‐pad to two digits
      !is.na(visit_num) ~ str_pad(visit_num, width = 2, pad = "0"),
      # otherwise leave NA
      TRUE ~ NA_character_
    )
  ) %>%
  select(-visit_num)


## Add timepoint info 
best_match$Study <- "M4"
best_match <- best_match %>%
  mutate(timepoint_info = case_when(
    # M4 study: baseline, induction, transplant, maintenance
    is.na(timepoint_info) & Study == "M4" & Timepoint == "01" ~ "Diagnosis",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "03" ~ "Post_induction",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "05" ~ "Post_transplant",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "07" ~ "Maintenance",
    # Relapse
    is.na(timepoint_info) & Timepoint == "R-"                   ~ "Relapse",
    TRUE                                                      ~ timepoint_info
  )) %>%
  # Special overrides for CA-08
  mutate(timepoint_info = case_when(
    is.na(timepoint_info) & Study == "M4" & Timepoint == "08"                             ~ "1.5yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "09" & Patient == "CA-08"         ~ "1.5yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "09"                             ~ "2yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "10"                             ~ "2.5yr maintenance",
    TRUE                                                                                ~ timepoint_info
  ))

Translocations_FISH_all_tp <- best_match


## Omit counts for now since so many NA values



### Add M-protein, free light chains, and other clinically relevant labs
# These fields contribute to baseline clinical characterization and downstream
# feature/correlation analyses. The selected lab list is intentionally explicit
# so additions from new data exports can be reviewed.
labs_file <- paste0(data_path, "M4_COHORT_LABS.xlsx")
labs <- read_excel(labs_file)

# Define the lab tests of interest for serum/urine protein electrophoresis and free light chains
relevant_labs <- c(
  "Protein Electrophoresis: M- protein",            # Key for M-protein detection
  "Protein Electrophoresis: Total Ur Protein",      # Monitors proteinuria
  "Serum free light chain: Kappa",                  # Measures Kappa light chains
  "Serum free light chain: Lambda",                # Measures Lambda light chains
  "Serum free light chain: Kappa / Lambda ratio",   # Determines clonal excess
  "B2 microglobulin",                               # Part of ISS staging
  "Albumin, Plasma",                                # Part of ISS staging
  "Hemoglobin",                                     # Assesses anemia
  "Calcium Total, Plasma",                          # Indicates hypercalcemia
  "Creatinine, Plasma",                             # Monitors renal function
  "LDH",                                            # Tumor burden marker
  "Plasma cells",                                   # Assesses marrow involvement
  "Absolute number of plasma cells",               # Quantifies plasma cells
  "dFLC (difference between involved FLC and uninvolved FLC)", # Tracks disease burden
  "Immunoglobulin Quantitation: M-protein"
)


# Pre-process the PURPOSE column to remove '#' if present
labs <- labs %>%
  mutate(PURPOSE = gsub("#", "", PURPOSE))  # Remove '#' from PURPOSE

labs <- labs %>%
  mutate(
    Timepoint = case_when(
      grepl("At diagnosis|At Diagnosis", PURPOSE, ignore.case = TRUE) ~ "01",
      grepl("visit 3|Visit 3", PURPOSE, ignore.case = TRUE) ~ "03",
      grepl("visit 2|Visit 2", PURPOSE, ignore.case = TRUE) ~ "02",
      grepl("\\bvisit 1\\b|\\bVisit 1\\b", PURPOSE, ignore.case = TRUE) ~ "01",
      grepl("visit\\s*3", PURPOSE, ignore.case = TRUE) ~ "03",
      grepl("visit\\s*4", PURPOSE, ignore.case = TRUE) ~ "04",
      grepl("visit\\s*5", PURPOSE, ignore.case = TRUE) ~ "05",
      grepl("visit\\s*6", PURPOSE, ignore.case = TRUE) ~ "06",
      grepl("visit\\s*7", PURPOSE, ignore.case = TRUE) ~ "07",
      grepl("visit\\s*8", PURPOSE, ignore.case = TRUE) ~ "08",
      grepl("visit\\s*9", PURPOSE, ignore.case = TRUE) ~ "09",
      grepl("visit\\s*10", PURPOSE, ignore.case = TRUE) ~ "10",
      grepl("visit\\s*11", PURPOSE, ignore.case = TRUE) ~ "11",
      grepl("visit\\s*12", PURPOSE, ignore.case = TRUE) ~ "12",
      grepl("visit\\s*13", PURPOSE, ignore.case = TRUE) ~ "13",
      grepl("visit\\s*14", PURPOSE, ignore.case = TRUE) ~ "14",
      grepl("visit R|Visit R", PURPOSE, ignore.case = TRUE) ~ "R",
      grepl("relapse|Relapse", PURPOSE, ignore.case = TRUE) ~ "R",
      grepl("visit11|Visit11", PURPOSE, ignore.case = TRUE) ~ "11",
      grepl("visit12|Visit12", PURPOSE, ignore.case = TRUE) ~ "12",
      TRUE ~ NA_character_  # Leave as NA for other cases
    )
  )

# Convert LAB_VALUE to numeric, preserving the raw value for auditability.
# The parser is deliberately conservative: it removes explicit units and
# handles whole-field scientific notation, but it must not rewrite internal
# digit sequences. A previous broad "10...<digit>" regex converted values like
# "8.7100000000000009" into 8.7e9 and "2.107 mmol/L" into 2.0e7.
clean_lab_values <- function(lab_values) {
  x <- as.character(lab_values)
  x <- stringr::str_squish(x)
  x <- stringr::str_replace_all(x, "\u00d7", "x")
  x <- stringr::str_replace_all(x, ",", ".")
  x <- stringr::str_replace_all(x, "^\\.", "0.")
  x <- stringr::str_replace_all(x, regex("(?i)^negative$"), "0")
  x <- stringr::str_replace_all(x, regex("(?i)^positive$"), "1")
  x <- stringr::str_replace_all(x, "^[<>]\\s*([+-]?\\d+\\.?\\d*)", "\\1")
  x <- ifelse(
    stringr::str_detect(x, regex("^999(\\b|\\s|$)")),
    NA_character_,
    x
  )

  # Whole-field scientific notation only. This supports forms like
  # "2.7 E+9/L", "2.7 x 10^9/L", and "2.7 x10e9/L" without touching
  # ordinary long decimals such as "8.7100000000000009".
  x <- stringr::str_replace(
    x,
    regex("^([+-]?\\d+(?:\\.\\d+)?)\\s*E\\+?([+-]?\\d+)\\s*/?\\s*[A-Za-z]*\\s*$", ignore_case = TRUE),
    "\\1E\\2"
  )
  x <- stringr::str_replace(
    x,
    regex("^([+-]?\\d+(?:\\.\\d+)?)\\s*x\\s*10\\s*(?:\\^|\\*|e)?\\s*([+-]?\\d+)\\s*/?\\s*[A-Za-z]*\\s*$", ignore_case = TRUE),
    "\\1E\\2"
  )

  # Remove explicit measurement units. Do not remove bare letters such as "m"
  # or "s"; doing so corrupts legitimate numeric strings with long decimals.
  x %>%
    stringr::str_remove_all(regex("\\s*(mg/L|g/L|mmol/L|umol/L|µmol/L|nmol/L|U/L|mIU/L|IU/L|pg/ml|units/mL|mm/h|g/CP|%)\\s*", ignore_case = TRUE)) %>%
    stringr::str_trim()
}

# Use
labs <- labs %>%
  mutate(LAB_VALUE_raw = as.character(LAB_VALUE)) %>%
  mutate(LAB_VALUE_clean = clean_lab_values(LAB_VALUE)) ## doesn't do scientific notation but not needed

## Also now remove 999 as NA

labs <- labs %>%
  mutate(LAB_VALUE = as.numeric(LAB_VALUE_clean))

# Quarantine values that remain biologically/analytically implausible after
# parsing. These rows are exported for manual review rather than silently
# corrected, because values such as "242 mmol/L" or "68 mmol/L" likely reflect
# source data entry issues but the correct decimal placement cannot be inferred
# reproducibly.
lab_plausibility_ranges <- tibble::tribble(
  ~LAB_NAME, ~min_value, ~max_value,
  "Calcium Total, Plasma", 0.5, 5,
  "Albumin, Plasma", 0, 80,
  "Hemoglobin", 0, 250,
  "B2 microglobulin", 0, 10000,
  "Creatinine, Plasma", 0, 5000,
  "LDH", 0, 10000,
  "Plasma cells", 0, 100,
  "Absolute number of plasma cells", 0, 1000000000000,
  "Serum free light chain: Kappa / Lambda ratio", 0, 1000000,
  "Serum free light chain: Kappa", 0, 1000000,
  "Serum free light chain: Lambda", 0, 1000000,
  "Protein Electrophoresis: M- protein", 0, 200,
  "Immunoglobulin Quantitation: M-protein", 0, 200,
  "Protein Electrophoresis: Total Ur Protein", 0, 1000000,
  "dFLC (difference between involved FLC and uninvolved FLC)", 0, 1000000
)

lab_value_plausibility_audit <- labs %>%
  inner_join(lab_plausibility_ranges, by = "LAB_NAME") %>%
  filter(!is.na(LAB_VALUE), LAB_VALUE < min_value | LAB_VALUE > max_value) %>%
  select(
    M4_id,
    study_patient_id,
    LAB_TYPE,
    LAB_DATE,
    PURPOSE,
    LAB_NAME,
    LAB_VALUE_raw,
    LAB_VALUE_clean,
    parsed_value = LAB_VALUE,
    min_value,
    max_value,
    RANGE,
    RESULT,
    DETAILS
  ) %>%
  arrange(LAB_NAME, desc(parsed_value), M4_id, LAB_DATE)

write_active_csv(
  lab_value_plausibility_audit,
  "Output_tables_2025/clinical_support/lab_value_plausibility_audit.csv"
)

labs <- labs %>%
  left_join(lab_plausibility_ranges, by = "LAB_NAME") %>%
  mutate(
    LAB_VALUE = if_else(
      !is.na(min_value) &
        !is.na(LAB_VALUE) &
        (LAB_VALUE < min_value | LAB_VALUE > max_value),
      NA_real_,
      LAB_VALUE
    )
  ) %>%
  select(-min_value, -max_value) %>%
  filter(!is.na(LAB_VALUE))

# Transform selected labs to one row per patient/timepoint where possible.
# Transform the labs table to wide format with special handling for 0 vs non-zero
labs_transformed <- labs %>%
  filter(LAB_NAME %in% relevant_labs) %>%
  # Combine LAB_TYPE and LAB_NAME to create unique column names
  mutate(LAB_COLUMN = paste(LAB_TYPE, LAB_NAME, sep = " - ")) %>%
  # Convert LAB_VALUE to numeric for processing
  mutate(LAB_VALUE = as.numeric(LAB_VALUE)) %>%
  # Group by key columns and LAB_COLUMN
  group_by(M4_id, PURPOSE, Timepoint, LAB_COLUMN, LAB_DATE) %>%
  summarise(
    LAB_VALUE = if (all(LAB_VALUE == 0) || all(LAB_VALUE != 0)) {
      mean(LAB_VALUE, na.rm = TRUE)  # Take average if all are 0 or all are non-zero
    } else {
      max(LAB_VALUE)  # Select the non-zero value when there's both 0 and non-zero
    },
    .groups = "drop"
  ) %>%
  # Pivot to wide format
  pivot_wider(
    names_from = LAB_COLUMN,  # Use LAB_COLUMN for new column names
    values_from = LAB_VALUE   # Fill new columns with LAB_VALUE
  )

# Keep the downstream schema stable even when an optional lab has no usable
# values after parsing/plausibility filtering. Without this, pivot_wider() omits
# all-absent columns and the explicit clinical_all_tp select below can fail.
expected_lab_columns <- paste("Blood labs set", relevant_labs, sep = " - ")
missing_lab_columns <- setdiff(expected_lab_columns, names(labs_transformed))
if (length(missing_lab_columns) > 0) {
  message(
    "Adding absent lab columns as NA_real_: ",
    paste(missing_lab_columns, collapse = ", ")
  )
  labs_transformed[missing_lab_columns] <- NA_real_
}


### Add the timepoint if have from other table 
# 1) Prep the translocation table
ct2 <- labs_transformed %>%
  rename(Patient = M4_id) %>%
  mutate(PROCEDURE_DATE = as.Date(LAB_DATE))

# 3) Fuzzy join: match Patient exactly and dates within +/-60 days.
joined2 <- fuzzy_left_join(
  ct2, ccu_small,
  by = c(
    "Patient"        = "Patient",
    "PROCEDURE_DATE" = "Date"
  ),
  match_fun = list(
    `==`,                                # exact match on Patient
    function(d1, d2) abs(as.numeric(d1 - d2)) <= 60 # 60 day difference in dates permitted in case test done later
  )
)

# 4) Of the (possibly multiple) matches, pick the one with smallest date‐difference
joined2 <- joined2 %>%
  rename(Patient = Patient.x) %>%
  select(-Patient.y)

joined2 <- joined2 %>%
  mutate(days_diff = abs(as.numeric(PROCEDURE_DATE - Date))) %>%
  group_by(Patient, PROCEDURE_DATE) %>%
  slice_min(order_by = days_diff, n = 1, with_ties = FALSE) %>%
  ungroup()

joined2 <- joined2 %>%
  mutate(
    # If Timepoint.x is NA but Timepoint.y is not, use Timepoint.y
    Timepoint.x = coalesce(Timepoint.x, Timepoint.y)
  ) %>%
  # Drop the now-redundant Timepoint.y
  select(-Timepoint.y) %>%
  # Rename Timepoint.x back to Timepoint
  rename(Timepoint = Timepoint.x)

## Change R to R- 
joined2 <- joined2 %>%
  mutate(
    Timepoint = if_else(Timepoint == "R", "R-", Timepoint)
  )

## Remove previous columns 
joined2 <- joined2 %>% select(-days_diff, -Date)

## Add timepoint info 
joined2$Study <- "M4"
joined2 <- joined2 %>%
  mutate(timepoint_info = case_when(
    # M4 study: baseline, induction, transplant, maintenance
    is.na(timepoint_info) & Study == "M4" & Timepoint == "01" ~ "Diagnosis",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "03" ~ "Post_induction",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "05" ~ "Post_transplant",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "07" ~ "Maintenance",
    # Relapse
    is.na(timepoint_info) & Timepoint == "R-"                   ~ "Relapse",
    TRUE                                                      ~ timepoint_info
  )) %>%
  # Special overrides for CA-08
  mutate(timepoint_info = case_when(
    is.na(timepoint_info) & Study == "M4" & Timepoint == "08"                             ~ "1.5yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "09" & Patient == "CA-08"         ~ "1.5yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "09"                             ~ "2yr maintenance",
    is.na(timepoint_info) & Study == "M4" & Timepoint == "10"                             ~ "2.5yr maintenance",
    TRUE                                                                                ~ timepoint_info
  ))

Labs_all_timepoints <- joined2

### Derive immunoglobulin/light-chain subtype from available lab values
# These rules are deterministic heuristics based on the available quantitative
# immunoglobulin and free-light-chain measurements. Ambiguous cases are retained
# as Other or subtype-unknown rather than forced into a category.
# Step 1: Filter Relevant Data for Immunoglobulin and Light Chains
filtered_data <- labs %>%
  filter(grepl("Immunoglobulin Quantitation|free light chain", LAB_NAME)) %>%
  mutate(
    Type = case_when(
      grepl("Ig G", LAB_NAME, ignore.case = TRUE) ~ "IgG",
      grepl("Ig A", LAB_NAME, ignore.case = TRUE) ~ "IgA",
      grepl("Ig M", LAB_NAME, ignore.case = TRUE) ~ "IgM",
      grepl("^Serum free light chain: Kappa$", LAB_NAME, ignore.case = TRUE) ~ "Kappa",
      grepl("free light chain: Lambda", LAB_NAME, ignore.case = TRUE) ~ "Lambda",
      grepl("Kappa / Lambda ratio", LAB_NAME, ignore.case = TRUE) ~ "Kappa/Lambda Ratio",
      grepl("Immunoglobulin Quantitation: M-protein", LAB_NAME, ignore.case = TRUE) ~ "M-protein",
      TRUE ~ NA_character_
    ),
    LAB_VALUE = as.numeric(LAB_VALUE)  # Ensure LAB_VALUE is numeric
  )


# Summarize duplicate rows by averaging values and removing NAs
filtered_data_summary <- filtered_data %>%
  filter(!is.na(Type)) %>%
  group_by(M4_id, Timepoint, Type) %>%  # Group by unique identifier columns
  summarise(
    LAB_VALUE = mean(LAB_VALUE, na.rm = TRUE),  # Take the average, ignoring NAs
    .groups = "drop"  # Drop the grouping after summarization
  )

# Step 2: Reshape the data to wide format
patient_data <- filtered_data_summary %>%
  pivot_wider(
    names_from = Type,
    values_from = LAB_VALUE,
    values_fill = NA  # Fill missing values with NA
  )

# Step 3: Derive exploratory subtypes from the available quantitative values.
# These are deterministic helper annotations, not adjudicated clinical subtype
# labels; missing/discordant values remain subtype-unknown or Other.

## update to also score if Kappa or Lambda
classified_subtypes <- patient_data %>%
  mutate(
    # Ensure Kappa, Lambda, and Kappa/Lambda Ratio are numeric
    Kappa = as.numeric(Kappa),
    Lambda = as.numeric(Lambda),
    `Kappa/Lambda Ratio` = as.numeric(`Kappa/Lambda Ratio`),
    
    # Compute Kappa/Lambda Ratio if missing but Kappa & Lambda are available
    `Kappa/Lambda Ratio` = ifelse(
      is.na(`Kappa/Lambda Ratio`) & !is.na(Kappa) & !is.na(Lambda),
      Kappa / Lambda,
      `Kappa/Lambda Ratio`
    ),
    
    # Define subtypes based on immunoglobulin and light chain restriction
    Subtype = case_when(
      # Light Chain Only (Kappa or Lambda)
      (is.na(IgG) | IgG < 6) &  
        (is.na(IgA) | IgA < 0.8) &  
        (is.na(IgM) | IgM < 0.5) &  
        (!is.na(Kappa) & Kappa > Lambda) &  
        (is.na(`Kappa/Lambda Ratio`) | `Kappa/Lambda Ratio` > 1.65) ~ "Light Chain Only (Kappa)",
      
      (is.na(IgG) | IgG < 6) &  
        (is.na(IgA) | IgA < 0.8) &  
        (is.na(IgM) | IgM < 0.5) &  
        (!is.na(Lambda) & Lambda > Kappa) &  
        (is.na(`Kappa/Lambda Ratio`) | `Kappa/Lambda Ratio` < 0.26) ~ "Light Chain Only (Lambda)",
      
      # IgG Kappa or Lambda
      (!is.na(IgG) & IgG > IgA & IgG > IgM) & 
        (!is.na(Kappa) & Kappa > Lambda) ~ "IgG Kappa",
      
      (!is.na(IgG) & IgG > IgA & IgG > IgM) & 
        (!is.na(Lambda) & Lambda > Kappa) ~ "IgG Lambda",
      
      # IgA Kappa or Lambda
      (!is.na(IgA) & IgA > IgG & IgA > IgM) & 
        (!is.na(Kappa) & Kappa > Lambda) ~ "IgA Kappa",
      
      (!is.na(IgA) & IgA > IgG & IgA > IgM) & 
        (!is.na(Lambda) & Lambda > Kappa) ~ "IgA Lambda",
      
      # IgM Kappa or Lambda (rare in MM)
      (!is.na(IgM) & IgM > IgG & IgM > IgA) & 
        (!is.na(Kappa) & Kappa > Lambda) ~ "IgM Kappa",
      
      (!is.na(IgM) & IgM > IgG & IgM > IgA) & 
        (!is.na(Lambda) & Lambda > Kappa) ~ "IgM Lambda",
      
      # Cases where Kappa or Lambda is NA but IgG or IgA is present
      (!is.na(IgG) & IgG > IgA & IgG > IgM) & (is.na(Kappa) | is.na(Lambda)) ~ "IgG (Subtype Unknown)",
      (!is.na(IgA) & IgA > IgG & IgA > IgM) & (is.na(Kappa) | is.na(Lambda)) ~ "IgA (Subtype Unknown)",
      
      # If cannot classify based on available data
      TRUE ~ "Other"
    )
  )


# Step 4: Summarize Subtypes Across the Cohort
subtype_summary <- classified_subtypes %>%
  group_by(Subtype) %>%
  summarise(
    Count = dplyr::n(),
    Proportion = round(100 * Count / nrow(classified_subtypes), 1)
  ) %>%
  mutate(Summary = paste0(Count, " (", Proportion, "%)"))


## Got the subtypes for all pts

## Join back to labs 
tmp <- classified_subtypes %>%
  rename(Patient = M4_id)

# 2) Join labs + subtypes by Patient & Timepoint
labs_with_subtype <- Labs_all_timepoints %>%
  left_join(
    tmp,
    by = c("Patient", "Timepoint")
  )

# 2.5 Merge Translocations and Staging Data with Demographics
clinical_data <- demographics %>%
  left_join(translocations_annotated, by = c("M4_id")) %>%
  left_join(staging, by = c("M4_id", "study_patient_id")) %>%
  left_join(labs_transformed %>% filter(Timepoint == "01") %>% select(-Timepoint, -PURPOSE)) %>%
  left_join(classified_subtypes %>% filter(Timepoint == "01") %>% select(-Timepoint)) %>% 
  left_join(classified_translocations %>% select(M4_id, Cytogenetic_Risk)) %>% unique()

clinical_data <- clinical_data %>% filter(DISEASE_STATUS == "Initial diagnosis") %>% unique()


## For all timepoints
## -- NO filtering on Timepoint here --------------------------

## Edit staging 
staging2 <- staging %>%
  # Rename M4_id → Patient
  rename(Patient = M4_id) %>%
  # Create a Timepoint column based on DISEASE_STATUS
  mutate(
    Timepoint = case_when(
      DISEASE_STATUS == "Initial diagnosis"     ~ "01",
      DISEASE_STATUS == "Progression/relapse"   ~ "R-",
      TRUE                                       ~ NA_character_
    )
  )

clinical_all_tp <- labs_with_subtype %>%
  full_join(staging2) %>%
  full_join(Translocations_FISH_all_tp %>% select(-PROCEDURE_DATE, -Study))


clinical_all_tp <- clinical_all_tp %>%
  select(
    Patient,
    Timepoint,
    timepoint_info,
    Date = LAB_DATE,
    Study,
    # Key labs
    dFLC              = `Blood immunology labs set - dFLC (difference between involved FLC and uninvolved FLC)`,
    Albumin           = `Blood labs set - Albumin, Plasma`,
    B2_micro          = `Blood labs set - B2 microglobulin`,
    Calcium           = `Blood labs set - Calcium Total, Plasma`,
    Creatinine        = `Blood labs set - Creatinine, Plasma`,
    Hemoglobin        = `Blood labs set - Hemoglobin`,
    LDH               = `Blood labs set - LDH`,
    Plasma_pct        = `Blood labs set - Plasma cells`,
    Plasma_count      = `Blood labs set - Absolute number of plasma cells`,
    # Immuno‐subtype
    IgA,
    IgG,
    IgM,
    Kappa,
    Kappa_Lambda_Ratio = `Kappa/Lambda Ratio`,
    Lambda,
    M_Protein = `M-protein`,
    Subtype,
    # Clinical staging
    DIAGNOSIS_DATE,
    ECOG_SCORE,
    ISS_STAGE,
    KPS_SCORE,
    DETAILS,
    # Cytogenetics
    T_4_14,
    T_11_14,
    T_14_16,
    DEL_17P,
    DEL_1P,
    AMP_1Q,
    Cytogenetic_Risk,
    Abnormality_Count,
    Abnormality_Status
  ) %>%
  arrange(Patient, Timepoint)

clinical_all_tp <- clinical_all_tp %>% filter(!is.na(Timepoint))

### Consolidate 
# 1) Convert to a data.table
dt <- as.data.table(clinical_all_tp)

# 2) Define a helper that returns first non‐NA or a typed NA
first_non_na_or_type <- function(col) {
  vals <- col[!is.na(col)]
  if (length(vals)) {
    # if there is a real value, return it
    return(vals[1])
  }
  # otherwise return NA of the same class
  if (inherits(col, "POSIXct")) {
    return(as.POSIXct(NA))
  }
  if (inherits(col, "Date")) {
    return(as.Date(NA))
  }
  if (is.double(col)) {
    return(NA_real_)
  }
  if (is.integer(col)) {
    return(NA_integer_)
  }
  if (is.logical(col)) {
    return(NA)
  }
  # everything else, treat as character
  return(NA_character_)
}

## Group relapse timepoints within two weeks, else keep seperate. Group other timepoint duplicates
# 3) Compute per-patient relapse windows
dt[, `:=`(
  min_relapse = min(Date), 
  max_relapse = max(Date)
), by = .(Patient, Timepoint, timepoint_info)]

# 4) Now build group_date 
dt[, group_date := {
  # 1) start with an all-NA Date vector
  gd <- rep(as.Date(NA), .N)
  
  # 2) logical mask of relapse rows, force NAs → FALSE
  is_relapse <- timepoint_info == "Relapse"
  is_relapse[is.na(is_relapse)] <- FALSE
  
  # 3) among those relapse rows, are they *all* within 14 days?
  #    use the already-calculated min_relapse/max_relapse
  within_14d <- is_relapse & (max_relapse - min_relapse) <= 14
  
  # 4) find the relapse rows that are *not* within 14 days
  idx_keep_separate <- which(is_relapse & !within_14d)
  
  # 5) for those “far-apart” relapses, preserve their own Date
  gd[idx_keep_separate] <- Date[idx_keep_separate]
  
  # 6) everything else (close relapses, or non-relapses) stays NA
  gd
}]


# 5) Collapse by Patient+Timepoint+group_date
# collapse into a helper then drop columns in a separate step:
tmp <- dt[
  , lapply(.SD, first_non_na_or_type),
  by = .(Patient, Timepoint, group_date)
]

clinical_consolidated <- tmp[
  , .SD, .SDcols = setdiff(names(tmp), c("min_relapse","max_relapse","group_date"))
]


# 3) Verify for duplicates
dup_check <- clinical_consolidated[
  , .N, by = .(Patient, Timepoint)
][N > 1]

if (nrow(dup_check)) {
  warning("Still have duplicates for:\n", paste(dup_check$Patient, dup_check$Timepoint, sep=":", collapse="\n"))
} else {
  message("All Patient+Timepoint combinations are now unique.")
}

## Just for relapse cases - this is ok


## Correct ECOG score 
clinical_consolidated <- clinical_consolidated %>%
  mutate(
    ECOG_SCORE = case_when(
      str_starts(ECOG_SCORE, "0") ~ 0,
      str_starts(ECOG_SCORE, "1") ~ 1,
      str_starts(ECOG_SCORE, "2") ~ 2,
      str_starts(ECOG_SCORE, "3") ~ 3,
      str_starts(ECOG_SCORE, "4") ~ 4,
      ECOG_SCORE == "Unknown"     ~ NA_real_,
      TRUE                        ~ NA_real_  # catch-all for unexpected values
    )
  )




### Harmonize curated IMMAGINE and SPORE baseline clinical sheets
# These curated sheets fill cross-cohort baseline clinical fields that are not
# available in the M4 CMRG raw-workbook structure.
clinical_data_IMMAGINE <- read_excel("Clinical data/IMMAGINE/Clinical data for IMG patients at diagnosis_filled_DA_edited.xlsx")
clinical_data_SPORE <- read_excel("Clinical data/SPORE/Spore_baseline_clinical_demographics_DA_edited.xlsx")

clinical_data_IMMAGINE <- clinical_data_IMMAGINE %>%
  mutate(GENDER = recode(GENDER, "M" = "Male", "F" = "Female"))

clinical_data_add <- bind_rows(clinical_data_IMMAGINE, clinical_data_SPORE)

# 2.3 Categorize Age Groups
clinical_data_add <- clinical_data_add %>%
  mutate(AGE_GROUP = case_when(
    AGE_AT_CONSENT < 50 ~ "<50",
    AGE_AT_CONSENT >= 50 & AGE_AT_CONSENT < 60 ~ "50-59",
    AGE_AT_CONSENT >= 60 & AGE_AT_CONSENT < 70 ~ "60-69",
    AGE_AT_CONSENT >= 70 & AGE_AT_CONSENT < 80 ~ "70-79",
    AGE_AT_CONSENT >= 80 ~ "80+",
    TRUE ~ "Unknown"
  ))

# Get numerals to match M4 
clinical_data_add <- clinical_data_add %>%
  mutate(
    ISS_STAGE = case_when(
      str_detect(ISS_STAGE, regex("stage 1", ignore_case = TRUE)) ~ "Stage I",
      str_detect(ISS_STAGE, regex("stage 2", ignore_case = TRUE)) ~ "Stage II",
      str_detect(ISS_STAGE, regex("stage 3", ignore_case = TRUE)) ~ "Stage III",
      TRUE ~ ISS_STAGE  # Keep original if unmatched
    )
  )

### Add in new table from Esteban and other info if have from David 
SPORE_new_info <- read_excel("Clinical data/SPORE/Spore patient timeline clinical info sent to Toronto.xlsx")

# helper to clean one vector of character values into numeric
clean_flc <- function(x) {
  x2 <- x %>%
    str_replace_all(regex("negative", ignore_case = TRUE), "0") %>%
    str_replace_all("[<>]", "") %>%
    str_trim() %>%
    na_if("") %>%
    na_if("NA")
  num <- as.numeric(x2)
  # in case any parsed < 0, bump up to zero
  num[num < 0] <- 0
  num
}

SPORE_new_info <- SPORE_new_info %>%
  mutate(
    `M-Spike`    = clean_flc(`M-Spike`),
    `Kappa FLC`  = clean_flc(`Kappa FLC`),
    `Lambda FLC` = clean_flc(`Lambda FLC`),
    `Ratio FLC`  = clean_flc(`Ratio FLC`)
  )

## Add it in 
# 1) prep the SPORE fragmentomics info
spore_labs <- SPORE_new_info %>%
  rename(
    M_Protein    = `M-Spike`,
    Kappa  = `Kappa FLC`,
    Lambda = `Lambda FLC`,
    Kappa_Lambda_Ratio  = `Ratio FLC`
  ) %>%
  # pull out the number after "_T" in Sample_ID
  mutate(
    Timepoint = str_extract(Sample_ID, "(?<=_T)\\d+")
  ) %>%
  select(Patient, Timepoint, M_Protein, Kappa, Lambda, Kappa_Lambda_Ratio, Date_of_sample_collection)


### Now add the FISH to this
### Generated from the processing clinical data script
spore_fish <- read_csv("Clinical data/SPORE/tidy_fish.csv") %>%
  rename(
    Patient            = patient,
    fish_date          = fish_date,
    Fish_abnormalities = fish_abnormalities,
    DEL_1P             = del_1p,
    DEL_13             = del_13,
    DEL_14             = del_14,
    DEL_17P            = del_17p,
    AMP_1Q             = gain_1q,
    trisomy_3          = trisomy_3,
    trisomy_7          = trisomy_7,
    trisomy_9          = trisomy_9,
    trisomy_15         = trisomy_15,
    hyperdiploid       = hyperdiploid,
    T_11_14            = t_11_14,
    T_4_14             = t_4_14,
    T_14_16            = t_14_16,
    T_14_20            = t_14_20,
    CKS1B_DUP          = cks1b_dup
  ) %>%
  mutate(fish_date = as.Date(fish_date))

# make sure both date‐columns are Date (not POSIXct)
spore_fish2 <- spore_fish %>%
  mutate(fish_date = as.Date(fish_date))

spore_labs2 <- spore_labs %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

# 2) fuzzy full join on Patient + within ±60 days
joined2 <- fuzzy_full_join(
  spore_fish2, spore_labs2,
  by = c(
    "Patient"                   = "Patient",
    "fish_date"                 = "Date_of_sample_collection"
  ),
  match_fun = list(
    `==`,                                     # exact match on Patient
    function(d1, d2) abs(as.numeric(d1 - d2)) <= 60 # Up to 60 days apart 
  )
)

joined2 <- joined2 %>%
  mutate(
    # choose fish vs cfWGS date
    Date = case_when(
      !is.na(DEL_1P)                      ~ fish_date,
      TRUE                                ~ Date_of_sample_collection
    ),
    # merge the Patient columns as before
    Patient = coalesce(Patient.x, Patient.y)
  ) %>%
  select(
    Patient, Date, everything(),
    -Patient.x, -Patient.y,
    -Date_of_sample_collection, -fish_date
  )

SPORE_labs_joined <- joined2 %>% unique()

## Get it to match same format
SPORE_labs_joined <- SPORE_labs_joined %>%
  mutate(
    across(
      where(is.logical),
      ~ case_when(
        is.na(.)      ~ NA_character_,
        .             ~ "Positive",
        TRUE          ~ "Negative"
      ),
      .names = "{.col}"
    )
  )

### Add the SPORE info to it 
clinical_data_SPORE <- clinical_data_add %>%
  filter(str_detect(Patient_ID, "SPORE"))

## Get dates - taking first since was at dx
# 1. Extract first date per patient from SPORE_labs_joined
spore_first_dates <- SPORE_labs_joined %>%
  group_by(Patient) %>%
  summarise(Date = min(Date, na.rm = TRUE), .groups = "drop")

# 2. Join the diagnosis date to clinical_data_SPORE
clinical_data_SPORE <- clinical_data_SPORE %>%
  left_join(spore_first_dates, by = c("Patient_ID" = "Patient"))


clinical_data_SPORE <- clinical_data_SPORE %>%
  rename(
    Patient            = Patient_ID,
    ECOG_SCORE         = ECOG_score,
    Date     = Date,
    Fish_abnormalities = FISH_abnormalities,
    AGE                = AGE_AT_CONSENT,
    Gender             = GENDER,
    AMP_1Q             = TRANSLOCATIONS,      # Placeholder: only rename if this column encodes 1q gain
    Cytogenetic_Risk   = Cytogenetic_Risk,    # Keep name, already aligned
    Subtype            = Subtype,
    ISS_STAGE          = ISS_STAGE,           # Already matches
    R_ISS_STAGE        = R_ISS_stage,
    AGE_GROUP          = AGE_GROUP,
  )

clinical_data_SPORE <- clinical_data_SPORE %>% select(-Fish_abnormalities, -AMP_1Q, -Translocations)

SPORE_labs_joined <- full_join(SPORE_labs_joined, clinical_data_SPORE)

SPORE_labs_joined$ECOG_SCORE <- as.numeric(SPORE_labs_joined$ECOG_SCORE)

# 2) join onto the consolidated clinical table
clinical_consolidated <- clinical_consolidated %>%
  full_join(SPORE_labs_joined)



### Add IMMAGINE FISH helper fields
# This file is generated in 1_1A. Because 1_1A carries an audit warning for the
# IMMAGINE FISH regex export, downstream analyses should treat these flags as
# reviewable clinical annotations rather than immutable raw data.
IMMAGINE_fish <- read.csv("Clinical data/Exported clinical data April 2025/IMMAGINE_fish_flags.csv")

## Get it to match same format
IMMAGINE_fish <- IMMAGINE_fish %>%
  mutate(
    across(
      where(is.logical),
      ~ case_when(
        is.na(.)      ~ NA_character_,
        .             ~ "Positive",
        TRUE          ~ "Negative"
      ),
      .names = "{.col}"
    )
  )

IMMAGINE_fish <- IMMAGINE_fish %>%
  rename(
    Patient           = id,
    Fish_abnormalities = fish_text,
    DEL_13            = del_13q,
    DEL_17P           = del_17p,
    DEL_1P            = loss_1p,
    AMP_1Q            = gain_1q,
    T_11_14           = t_11_14,
    T_4_14            = t_4_14,
    T_14_16           = t_14_16,
    T_14_20           = t_14_20,
    CKS1B_DUP         = cks1b_dup
    # trisomy_* columns already match
  )

IMMAGINE_fish$timepoint_info <- "Diagnosis"


### Add M-protein or ISS stage 
### Add the SPORE info to it 
clinical_data_IMG <- clinical_data_add %>%
  filter(str_detect(Patient_ID, "IMG"))

clinical_data_IMG <- clinical_data_IMG %>%
  rename(
    Patient            = Patient_ID,
    ECOG_SCORE         = ECOG_score,
    Fish_abnormalities = FISH_abnormalities,
    AGE                = AGE_AT_CONSENT,
    Gender             = GENDER,
    AMP_1Q             = TRANSLOCATIONS,      # Placeholder: only rename if this column encodes 1q gain
    Cytogenetic_Risk   = Cytogenetic_Risk,    # Keep name, already aligned
    Subtype            = Subtype,
    ISS_STAGE          = ISS_STAGE,           # Already matches
    R_ISS_STAGE        = R_ISS_stage,
    AGE_GROUP          = AGE_GROUP,
  )

clinical_data_IMG <- clinical_data_IMG %>% select(-Fish_abnormalities, -AMP_1Q, -Translocations)

clinical_data_IMG$timepoint_info <- "Diagnosis"
IMMAGINE_labs_joined <- full_join(IMMAGINE_fish, clinical_data_IMG)

IMMAGINE_labs_joined$ECOG_SCORE <- as.numeric(IMMAGINE_labs_joined$ECOG_SCORE)


## Add back 
clinical_consolidated <- clinical_consolidated %>%
  full_join(IMMAGINE_labs_joined)


### Add some age and other data for M4 
# 1) Prepare demographics with the same “Patient” key
demog_clean <- demographics %>%
  rename(
    Patient    = M4_id,
    AGE.demo   = AGE_AT_CONSENT,
    Gender.demo= GENDER,
    AGE_GROUP.demo = AGE_GROUP
  ) %>%
  select(Patient, AGE.demo, Gender.demo, AGE_GROUP.demo)

# 2) Left-join them onto the clinical table, bringing in the “.demo” columns
clinical_filled <- clinical_consolidated %>%
  left_join(demog_clean, by = "Patient")

# 3) Coalesce so that if the main columns are NA, they take the .demo value,
#    otherwise they keep their original
clinical_filled <- clinical_filled %>%
  mutate(
    AGE       = coalesce(AGE, AGE.demo),
    Gender    = coalesce(Gender, Gender.demo),
    AGE_GROUP = coalesce(AGE_GROUP, AGE_GROUP.demo)
  ) %>%
  # 4) Drop the helper .demo columns
  select(-ends_with(".demo"))


### Final export: cross-cohort master clinical/lab table
# Downstream scripts use this table for clinical demographics, baseline
# covariates, cytogenetic/FISH comparisons, and clinical-feature summaries. The
# intended key is Patient + Timepoint. The duplicate warning above evaluates the
# consolidated M4 table before the later SPORE/IMMAGINE full joins, so a refreshed
# final export also requires a downstream key check in script 2_0.
write_active_csv(
  clinical_filled,
  "Clinical data/Master_clinical_data_table_all_projects_May2025_updated2.csv"
)
