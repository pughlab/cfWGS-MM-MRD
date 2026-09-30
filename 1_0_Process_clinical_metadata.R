# ──────────────────────────────────────────────────────────────────────────────
#  1_0_Process_clinical_metadata.R
#
#  Purpose:
#    • Ingest and harmonize multiple clinical/extract metadata sources (SPORE, M4, IMMAGINE, MyC).
#    • Clean and normalize sample identifiers, sample types, and collection dates.
#    • Annotate each sample with a unified timepoint_info (Diagnosis, Baseline, Post-induction,
#      Post-transplant, Maintenance, Relapse/Progression, etc.).
#    • Integrate FISH flags, lab timepoints, relapse/progression dates, and OS follow-up data.
#    • Compute time-to-relapse (absolute and non-absolute) and maximum follow-up durations for PFS analyses.
#    • Generate summary tables (counts by study, sample type, timepoint) and export QC plots.
#    • Export intermediary and final “combined_clinical_data_updated” master table along with
#      derived tables: baseline dates, PFS_days, latest_dates_per_patient, patient_counts, BAM lists, etc.
#
# How to run from the repository root:
#   Rscript 1_0_Process_clinical_metadata.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script harmonizes clinical metadata across
#     cohorts for downstream Table 1, swim-plot, model, and survival scripts.
#
# Pipeline role:
#   This script creates the clinical metadata inputs used by later
#   scripts. Historical root-level outputs that are read downstream are kept in
#   place for compatibility. Support-only QC summaries, temporary checkpoints,
#   and BAM review lists are written under
#   `Output_tables_2025/clinical_metadata_support/` so they are not confused
#   with manuscript figure/table outputs.
#
# Analysis unit and order:
#   The primary working unit is one sequenced or inventoried sample nested
#   within a patient. Patient-level progression, censoring, baseline, and
#   follow-up tables are derived later in the script. The major phases are:
#   source ingestion; sample/timepoint harmonization; curated corrections;
#   progression and follow-up integration; PFS endpoint construction; final
#   sample-level relapse annotation; and availability/QC exports.
#
# Main outputs used downstream:
#   - combined_clinical_data_updated_April2025.csv: sample-level clinical table.
#   - Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds:
#     current patient-level baseline/event/censor table.
#   - Exported_data_tables_clinical/Relapse_dates_full_updated.rds: retained
#     inventory of distinct progression dates, including later events.
#   - Exported_data_tables_clinical/patient_followup_dates_updated.rds:
#     patient-level last-known follow-up for prospective evaluability windows.
#   Earlier similarly named files are intermediate or compatibility exports;
#   comments beside each write identify their role.
#
#  Author:        Dory Abelman
#  Last updated:  2025-06-06
#
#  Notes:
#    • All file paths are assumed relative to the project root.
#    • Before running, ensure these inputs exist at the project-relative paths
#      used below:
#         – cfWGS_MS_Integrated_Metadata.xlsx
#         – SPORE: Spore_bams_with_collected_dates.txt,
#                  Extracted_clinical_data_cfDNA_project.xlsx,
#                  Spore_timepoint_info_updated.xlsx,
#                  SPORE_OS_info.xlsx
#         – M4: M4_COHORT_STAGING.xlsx,
#               M4_COHORT_STEM_CELL_TRANSPLANT.xlsx,
#               M4_COHORT_LABS.xlsx,
#               Relapse_dates.xlsx,
#               M4 Status and best response update Dec 2 24 just dates.xlsx,
#               M4_COHORT_DEMO.xlsx,
#               M4_COHORT_DIAGNOSTIC_BIOPSY.xlsx,
#               M4_COHORT_CHEMOTHERAPY.xlsx
#         – Root-level processing log: TFRIM4_Processing Log_Nov2024.xlsx
#         – IMMAGINE: IMG_request_20241009 (2).xlsx,
#                     Cleaned transplant dates just dates.xlsx,
#                     Extracted_clinical_MRD_data.xlsx,
#                     Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv,
#                     Additional_relapse_dates_IMG_from_Esther.xlsx,
#                     NAs in timepoint info - IMMAGINE additional timepoint info.csv
#         – Curated project inputs: cohort_assignment_table_updated.rds and
#                   Clinical data/manual_clinical_metadata_overrides.csv
#         – Optional endpoint refreshes, used only when present:
#                   M4_CMRG_Data/March 2026/M4_COHORT_CHEMOTHERAPY.csv,
#                   M4_CMRG_Data/March 2026/M4_COHORT_DEMO.csv, and
#                   New OICR Submissions/derived_metadata/
#                     oicr_revision_repo_style_metadata.csv
#
#    • Outputs:
#         – “combined_clinical_data_updated_April2025.csv” (master metadata table)
#         – Various exported tables under “Exported_data_tables_clinical” and “Oct2024_exported_data”
#         – Support/QC summaries and BAM review lists under
#           “Output_tables_2025/clinical_metadata_support”
#
#   This is a controlled upstream rebuild. Review source versions, manual
#   override tables, and endpoint audit outputs before refreshing downstream
#   manuscript analyses.
# ──────────────────────────────────────────────────────────────────────────────
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#



# ─── 1.  Load libraries ─────────────────────────────────────────────────────────
library(readxl)       # reading Excel files
library(data.table)   # fast fread/fwrite if needed
library(tidyverse)    # dplyr, tidyr, ggplot2, etc.
library(lubridate)    # date parsing

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

clinical_support_dir <- file.path("Output_tables_2025", "clinical_metadata_support")
clinical_support_plot_dir <- file.path(clinical_support_dir, "plots")
dir.create(clinical_support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_support_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Type-stable date normalization used before joins and patient-specific date
# corrections. Some source spreadsheets provide R Date objects, some provide
# POSIX datetimes, and some become numeric after mixed-column reshaping. Numeric
# values above 30000 are treated as Excel serial dates; smaller numeric values
# are treated as R Date offsets.
normalize_clinical_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))

  if (is.numeric(x)) {
    out <- rep(as.Date(NA), length(x))
    excel_serial <- !is.na(x) & x > 30000
    r_serial <- !is.na(x) & !excel_serial
    out[excel_serial] <- as.Date(x[excel_serial], origin = "1899-12-30")
    out[r_serial] <- as.Date(x[r_serial], origin = "1970-01-01")
    return(out)
  }

  x_chr <- as.character(x)
  numeric_text <- suppressWarnings(as.numeric(x_chr))
  parsed <- rep(as.Date(NA), length(x_chr))
  excel_serial_text <- !is.na(numeric_text) & numeric_text > 30000
  r_serial_text <- !is.na(numeric_text) & !excel_serial_text
  parsed[excel_serial_text] <- as.Date(numeric_text[excel_serial_text], origin = "1899-12-30")
  parsed[r_serial_text] <- as.Date(numeric_text[r_serial_text], origin = "1970-01-01")

  needs_text_parse <- is.na(parsed)
  parsed[needs_text_parse] <- suppressWarnings(ymd(x_chr[needs_text_parse]))
  still_missing <- is.na(parsed) & !is.na(x) & nzchar(as.character(x))
  parsed[still_missing] <- suppressWarnings(mdy(x_chr[still_missing]))
  still_missing <- is.na(parsed) & !is.na(x) & nzchar(as.character(x))
  parsed[still_missing] <- suppressWarnings(dmy(x_chr[still_missing]))
  as.Date(parsed)
}

max_clinical_date_or_na <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(as.Date(NA))
  max(x)
}

normalize_m4_endpoint_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
  if (is.numeric(x)) return(normalize_clinical_date(x))

  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "N/A", "Unknown", "NULL")] <- NA_character_

  parsed <- rep(as.Date(NA), length(x_chr))
  compact_day_month_year <- !is.na(x_chr) &
    stringr::str_detect(x_chr, "^[0-9]{1,2}[A-Za-z]{3}[0-9]{4}$")
  parsed[compact_day_month_year] <- suppressWarnings(
    lubridate::dmy(x_chr[compact_day_month_year])
  )

  still_missing <- is.na(parsed) & !is.na(x_chr) & nzchar(x_chr)
  parsed[still_missing] <- normalize_clinical_date(x_chr[still_missing])
  still_missing <- is.na(parsed) & !is.na(x_chr) & nzchar(x_chr)
  parsed[still_missing] <- suppressWarnings(lubridate::my(x_chr[still_missing]))
  as.Date(parsed)
}

read_march2026_csv_if_present <- function(filename) {
  path <- file.path("M4_CMRG_Data", "March 2026", filename)
  if (!file.exists(path)) {
    return(tibble())
  }
  readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
}

# ─── 2.  Configuration & helper functions (optional) ───────────────────────────
# [ Single place to change paths ]
#
# CFG <- list(
#   root        = here::here(),
#   out_tables  = "out/tables",
#   out_plots   = "out/plots",
#   match_window = 60   # use when matching SPORE sample dates to clinical dates
# )
# 
# dir.create(CFG$out_tables, showWarnings = FALSE, recursive = TRUE)
# dir.create(CFG$out_plots,  showWarnings = FALSE, recursive = TRUE)
# 
# pth <- function(...) fs::path(CFG$root, ...)
# write_out <- function(df, name) {
#   write_csv(df, pth(CFG$out_tables, paste0(name, ".csv")))
# }


# ─── 3.  Load raw clinical files ────────────────────────────────────────────────
# This block intentionally keeps each source file visible at the top of the
# script. If a new clinical export is received, update the path here and rerun
# from the command line before regenerating downstream figures/tables.
##  3.1  cf/WGS integrated metadata (M4, MyC etc.)
clinical_data <- read_excel("cfWGS_MS_Integrated_Metadata.xlsx", sheet = 1)

##  3.2  SPORE “bams with collected dates” spreadsheet
spore_data <- read_tsv("Clinical data/SPORE/Spore_bams_with_collected_dates.txt")

##  3.3  M4 staging (diagnosis dates)
M4_diagnosis_date      <- read_excel("Clinical data/M4/M4_COHORT_STAGING.xlsx")

##  3.4  SPORE relapse/extracted info (long list of diagnosis/relapse dates)
spore_data_extracted   <- read_excel("Clinical data/SPORE/Extracted_clinical_data_cfDNA_project.xlsx")

##  3.5  SPORE timepoint info (dates of interest)
spore_timepoint_info   <- read_excel("Clinical data/SPORE/Spore_timepoint_info_updated.xlsx")

##  3.6  SPORE overall survival info (for later OS annotations)
spore_OS_info          <- read_excel("Clinical data/SPORE/SPORE_OS_info.xlsx")

cohort_df <- readRDS("cohort_assignment_table_updated.rds")


# ─── 4.  Tidy SPORE “bams” table ────────────────────────────────────────────────
##  4.1  Create Sample_type from “Bams have” patterns
spore_data <- spore_data %>%
  mutate(
    Sample_type = case_when(
      grepl("Bm_T", `Bams have`) ~ "BM_cells",
      grepl("Pl_T", `Bams have`) ~ "Blood_plasma_cfDNA",
      grepl("Pb_R", `Bams have`) ~ "Blood_Buffy_coat",
      TRUE ~ NA_character_
    )
  )

##  4.2  Sort by Patient & date to compute “Timepoint” (≥14 days → new timepoint)
# RATIONALE: SPORE samples do not carry pre-assigned timepoint codes like M4.
# New timepoints are inferred: if a subsequent sample was collected >14 days
# after the prior one it is treated as a new clinical visit. 14 days was chosen
# as the minimum expected inter-visit gap; samples within the same 14-day window
# (e.g. same-day BM + blood draws) receive the same timepoint number. A shorter
# gap would over-split same-visit aliquots; a longer gap would merge distinct
# clinical assessments.
spore_data <- spore_data %>%
  arrange(Patient, Date_of_sample_collection) %>%
  group_by(Patient) %>%
  mutate(
    Timepoint = {
      date_diff      <- c(0, as.numeric(diff(Date_of_sample_collection)))
      timepoint_grp  <- cumsum(date_diff > 14) + 1
      as.character(timepoint_grp)
    }
  ) %>%
  ungroup()

##  4.3  Add “Study = SPORE” + unified “Sample_ID” and rename “Bams have” → “Bam”
spore_data <- spore_data %>%
  mutate(
    Study     = "SPORE",
    Sample_ID = paste0(Patient, "_T", Timepoint, "_", Sample_type)
  ) %>%
  rename(Bam = `Bams have`)


# ─── 5.  Tidy cf/WGS “clinical_data” (M4 + MyC) ─────────────────────────────────
clinical_data <- clinical_data %>%
  mutate(
    Sample_type = case_when(
      Sample_type == "P" ~ "Blood_plasma_cfDNA",
      Sample_type == "O" ~ "BM_cells",
      Sample_type == "B" ~ "Blood_Buffy_coat",
      TRUE ~ Sample_type
    )
  )

# (Note: at this point clinical_data already has its own Study column 
#  or you’ll add one manually, e.g., mutate(Study="M4") or mutate(Study="MyC").)


# ─── 6.  Combine SPORE + clinical_data into one “combined_clinical_data” ───────
combined_clinical_data <- bind_rows(spore_data, clinical_data) %>%
  filter(!is.na(Timepoint))   # drop any rows that lack a Timepoint


# ─── 7.  Annotate M4 timepoint_info (diagnosis, post_induction, etc.) ──────────
# M4 study protocol uses numeric codes to label each planned visit:
#   "01" = At diagnosis / pre-treatment (Timepoint_01_Prior_TX in the labs table)
#   "03" = After induction chemotherapy (post-induction BM / blood)
#   "05" = After autologous stem-cell transplant (post-ASCT)
#   "07" = ~1-year maintenance (the primary MRD assessment timepoint)
#   "08" = ~1.5-year maintenance  (see section 10 below for extended timepoints)
#   "09" = ~2-year maintenance
#   "R-" = Relapse / progression (visit triggered by clinical events, not protocol)
# These codes originate from the M4 study and are hard-coded in
# the TFRIM4 processing log; they should not change unless the protocol is amended.
combined_clinical_data <- combined_clinical_data %>%
  mutate(
    timepoint_info = case_when(
      Study == "M4" & Timepoint == "01" ~ "Diagnosis",
      Study == "M4" & Timepoint == "03" ~ "Post_induction",
      Study == "M4" & Timepoint == "05" ~ "Post_transplant",
      Study == "M4" & Timepoint == "07" ~ "Maintenance",
      Timepoint == "R-"                 ~ "Relapse",
      TRUE                               ~ NA_character_
    )
  )


# ─── 8.  Build SPORE timepoints long → integrate with combined_clinical_data ─────
##  8.1  Ensure SPORE timepoint info dates are Date class
spore_timepoint_info <- spore_timepoint_info %>%
  mutate(`Timepoint of interest` = as.Date(`Timepoint of interest`))

##  8.2  Transform “spore_data_extracted” to long: Actual_Diagnosis → Diagnosis,
##       Relapse_date* → Progression + add Status_at_timepoint = NA
spore_data_long <- spore_data_extracted %>%
  select(Patient, Actual_Diagnosis_Date, Relapse_date, Relapse_date_2, Relapse_date_3,
         Relapse_date_4, Relapse_date_5, Relapse_date_6) %>%
  pivot_longer(
    cols      = c(Actual_Diagnosis_Date, Relapse_date, Relapse_date_2, Relapse_date_3,
                  Relapse_date_4, Relapse_date_5, Relapse_date_6),
    names_to  = "timepoint_info",
    values_to = "Timepoint of interest"
  ) %>%
  filter(!is.na(`Timepoint of interest`)) %>%
  mutate(
    timepoint_info = case_when(
      timepoint_info == "Actual_Diagnosis_Date" ~ "Diagnosis",
      grepl("Relapse_date", timepoint_info)    ~ "Progression",
      TRUE                                      ~ timepoint_info
    ),
    Status_at_timepoint = NA_character_
  )

##  8.3  Bind the two SPORE sources & remove duplicates
combined_spore_timepoint_info <- bind_rows(
  spore_timepoint_info,
  spore_data_long
) %>%
  select(Patient, `Timepoint of interest`, timepoint_info) %>%
  mutate(timepoint_info = gsub("Relapse", "Progression", timepoint_info)) %>%
  distinct()


# ─── 9.  Join combined_clinical_data (only SPORE rows) → timepoint_info by Patient ─
##      and pick the row with smallest date_diff ≤ 60 days; else “Awaiting_integration”
combined_clinical_data <- combined_clinical_data %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

temp <- combined_clinical_data %>%
  filter(Study == "SPORE") %>%
  select(-timepoint_info) %>%
  left_join(combined_spore_timepoint_info, by = "Patient") %>%
  mutate(
    date_diff = abs(as.numeric(
      difftime(Date_of_sample_collection, `Timepoint of interest`, units = "days")
    ))
  ) %>%
  mutate(
    timepoint_info = ifelse(!is.na(date_diff) & date_diff <= 60, timepoint_info, "Awaiting_integration")
  ) %>%
  filter(timepoint_info != "Awaiting_integration") %>%
  group_by(Bam) %>%
  slice_min(date_diff, with_ties = FALSE) %>%
  select(-`Timepoint of interest`, -date_diff) %>%
  distinct() %>%
  ungroup()

##  9.1  All SPORE rows without a match get “Awaiting_integration”
combined_spore_data <- combined_clinical_data %>%
  filter(Study == "SPORE") %>%
  select(-timepoint_info) %>%
  left_join(temp, by = c("Bam","Patient","Sample_type","Timepoint")) %>%
  mutate(timepoint_info = coalesce(timepoint_info, "Awaiting_integration"))

##  9.2  Re‐bind SPORE + non‐SPORE (to update combined_clinical_data_updated)
combined_clinical_data_updated <- combined_clinical_data %>%
  filter(Study != "SPORE") %>%
  bind_rows(combined_spore_data) %>%
  mutate(timepoint_info = if_else(is.na(timepoint_info), "Awaiting_integration", timepoint_info))


# ─── 10.  Fix special M4 cases (CA-08, IMG-060, IMG-098, etc.) ───────────────────
##   (e.g. “CA-08” mislabeled “09” → “08”, IMG-060 T0 → Diagnosis, etc.)
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    timepoint_info = case_when(
      Study == "M4" & Timepoint == "08" ~ "1.5yr maintenance",
      Study == "M4" & Timepoint == "09" & Patient == "CA-08" ~ "1.5yr maintenance",
      Study == "M4" & Timepoint == "09"              ~ "2yr maintenance",
      Study == "M4" & Timepoint == "10"              ~ "2.5yr maintenance",
      Patient == "IMG-060" & Timepoint == "T0"       ~ "Diagnosis",
      Patient == "IMG-098" & (Timepoint == "T0" | Timepoint == "T1") ~ "Diagnosis",
      grepl("IMG", Patient) & Timepoint == "T0"      ~ "Diagnosis",
      TRUE                                           ~ timepoint_info
    )
  )

manual_clinical_metadata_override_path <- "Clinical data/manual_clinical_metadata_overrides.csv"
if (!file.exists(manual_clinical_metadata_override_path)) {
  stop(
    "Missing required manual clinical metadata override file: ",
    manual_clinical_metadata_override_path
  )
}

manual_clinical_metadata_overrides <- readr::read_csv(
  manual_clinical_metadata_override_path,
  col_types = readr::cols(
    .default = readr::col_character(),
    Date_of_sample_collection = readr::col_date(),
    include_in_recovered_test_cohort = readr::col_logical()
  ),
  show_col_types = FALSE
) %>%
  mutate(
    across(
      c(Patient, Sample_ID, Sample_Code, Timepoint, timepoint_info, Cohort),
      ~ na_if(trimws(.x), "")
    ),
    include_in_recovered_test_cohort = coalesce(
      include_in_recovered_test_cohort,
      FALSE
    )
  )

required_manual_override_cols <- c(
  "Patient",
  "Sample_ID",
  "Sample_Code",
  "Timepoint",
  "Date_of_sample_collection",
  "timepoint_info",
  "Cohort",
  "include_in_recovered_test_cohort"
)
missing_manual_override_cols <- setdiff(
  required_manual_override_cols,
  names(manual_clinical_metadata_overrides)
)
if (length(missing_manual_override_cols) > 0) {
  stop(
    "Manual clinical metadata override file is missing required columns: ",
    paste(missing_manual_override_cols, collapse = ", ")
  )
}
if (any(is.na(manual_clinical_metadata_overrides$Patient))) {
  stop("Manual clinical metadata override file contains rows without Patient.")
}

manual_sample_id_date_overrides <- manual_clinical_metadata_overrides %>%
  filter(!is.na(Date_of_sample_collection), !is.na(Sample_ID)) %>%
  transmute(
    Patient,
    Sample_ID,
    manual_sample_id_date = Date_of_sample_collection,
    manual_sample_id_timepoint_info = timepoint_info
  ) %>%
  distinct()

manual_patient_timepoint_date_overrides <- manual_clinical_metadata_overrides %>%
  filter(
    !is.na(Date_of_sample_collection),
    is.na(Sample_ID),
    !is.na(Timepoint)
  ) %>%
  transmute(
    Patient,
    Timepoint,
    manual_timepoint_date = Date_of_sample_collection,
    manual_timepoint_info = timepoint_info
  ) %>%
  distinct()

## Patient-specific sample-date corrections live in
## Clinical data/manual_clinical_metadata_overrides.csv rather than executable
## code, so reruns remain reproducible without embedding PHI-like dates here.
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(
    manual_sample_id_date_overrides,
    by = c("Patient", "Sample_ID")
  ) %>%
  left_join(
    manual_patient_timepoint_date_overrides,
    by = c("Patient", "Timepoint")
  ) %>%
  mutate(
    Date_of_sample_collection = coalesce(
      manual_sample_id_date,
      manual_timepoint_date,
      normalize_clinical_date(Date_of_sample_collection)
    ),
    timepoint_info = coalesce(
      manual_sample_id_timepoint_info,
      manual_timepoint_info,
      timepoint_info
    )
  ) %>%
  select(
    -manual_sample_id_date,
    -manual_sample_id_timepoint_info,
    -manual_timepoint_date,
    -manual_timepoint_info
  ) %>%
  mutate(
    timepoint_info = case_when(
      Patient == "IMG-181" & Timepoint %in% c("T6", "T12", "T15") ~ "Maintenance",
      TRUE ~ timepoint_info
    ),
    Timepoint = ifelse(Patient == "IMG-181" & Timepoint == "T6", "07", Timepoint)
  )

recoverable_img_wgs_sample_dates <- manual_clinical_metadata_overrides %>%
  filter(
    include_in_recovered_test_cohort,
    !is.na(Date_of_sample_collection),
    !is.na(Sample_ID)
  ) %>%
  transmute(Patient, Sample_ID, Date_of_sample_collection) %>%
  distinct()

# ─── 11.  Write out the first combined clinical-data intermediate ──────────────
# This temporary checkpoint is useful for auditing timepoint harmonization before
# later relapse, PFS, and patient-specific corrections are applied. It is not the
# final clinical metadata table used by downstream manuscript scripts.
write.csv(
  combined_clinical_data_updated,
  file.path(clinical_support_dir, "combined_clinical_data_updated_tmp.csv"),
  row.names = FALSE
)


# ─── 12.  QC: get counts & write summary tables ────────────────────────────────
# These counts are command-line QC outputs. They help detect unexpected shifts
# when clinical exports change or new samples are added; they are not formatted
# as manuscript tables.
##  12.1  Number of unique patients
num_patients <- combined_clinical_data_updated %>%
  summarise(num_patients = n_distinct(Patient))

##  12.2  Samples by Sample_type
samples_by_type      <- combined_clinical_data_updated %>%
  group_by(Sample_type) %>%
  summarise(num_samples = dplyr::n())

##  12.3  Samples by timepoint_info
samples_by_timepoint <- combined_clinical_data_updated %>%
  group_by(timepoint_info) %>%
  summarise(num_samples = dplyr::n())

##  12.4  Samples by Study
samples_by_study     <- combined_clinical_data_updated %>%
  group_by(Study) %>%
  summarise(num_samples = dplyr::n())

##  12.5  Print to screen & write CSVs
print(num_patients)
print(samples_by_type)
print(samples_by_timepoint)
print(samples_by_study)

# Export support-only count summaries used to audit clinical metadata updates.
write.csv(num_patients, file.path(clinical_support_dir, "num_patients.csv"), row.names = FALSE)
write.csv(samples_by_type, file.path(clinical_support_dir, "samples_by_type.csv"), row.names = FALSE)
write.csv(samples_by_timepoint, file.path(clinical_support_dir, "samples_by_timepoint.csv"), row.names = FALSE)
write.csv(samples_by_study, file.path(clinical_support_dir, "samples_by_study.csv"), row.names = FALSE)


# ─── 13.  Make simple QC barplots (Sample_type, timepoint_info, Study) ─────────
# These plots are quick data-integrity checks for the clinical metadata build.
# They are not mapped manuscript figure panels.
output_dir <- clinical_support_plot_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##  13.1  Samples by Sample_type
plot_samples_by_type <- ggplot(samples_by_type,
                               aes(x = Sample_type, y = num_samples, fill = Sample_type)) +
  geom_col(stat = "identity") +
  labs(title = "Number of Samples by Sample Type",
       x = "Sample Type", y = "Number of Samples") +
  theme_minimal()
ggsave(file.path(output_dir, "samples_by_type.png"), width = 6, height = 4)

##  13.2  Samples by Timepoint (factor levels in order)
samples_by_timepoint$timepoint_info <- factor(
  samples_by_timepoint$timepoint_info,
  levels = c("Diagnosis","Baseline","Post_induction","Post_transplant",
             "Maintenance","1.5yr maintenance","2yr maintenance",
             "Progression","Relapse","Awaiting_integration","Unknown")
)
plot_samples_by_timepoint <- ggplot(samples_by_timepoint,
                                    aes(x = timepoint_info, y = num_samples, fill = timepoint_info)) +
  geom_col(stat = "identity") +
  geom_text(aes(label = paste0("n=", num_samples)), vjust = -0.5, size = 3.5) +
  labs(title = "Number of Samples by Timepoint",
       x = "Timepoint", y = "Number of Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "samples_by_timepoint.png"), width = 6, height = 4)

##  13.3  Stacked bar: Samples by Study & Sample_type
samples_by_study_and_type <- combined_clinical_data %>%
  group_by(Study, Sample_type) %>%
  summarise(num_samples = dplyr::n(), .groups = "drop")

plot_samples_by_study_and_type <- ggplot(samples_by_study_and_type,
                                         aes(x = Study, y = num_samples, fill = Sample_type)) +
  geom_col(stat = "identity", position = "stack") +
  labs(title = "Samples by Study and Sample Type",
       x = "Study", y = "Number of Samples") +
  theme_minimal()
ggsave(file.path(output_dir, "samples_by_study_and_type.png"), width = 6, height = 4)


##  13.4  Number of Patients by Study
patients_by_study <- combined_clinical_data %>%
  group_by(Study) %>%
  summarise(num_patients = n_distinct(Patient), .groups = "drop")

plot_patients_by_study <- ggplot(patients_by_study,
                                 aes(x = Study, y = num_patients, fill = Study)) +
  geom_col(stat = "identity") +
  geom_text(aes(label = paste0("n=", num_patients)), vjust = -0.5, size = 5) +
  labs(title = "Number of Patients by Study",
       x = "Study", y = "Number of Patients") +
  theme_minimal()
ggsave(file.path(output_dir, "patients_by_study.png"), width = 6, height = 4)


# ─── 14.  Tidy M4 labs & extract timepoints → combine back into main table ──────
##  14.1  Read M4 lab workbook
M4_Labs <- read_excel("M4_CMRG_Data/M4_COHORT_LABS.xlsx")

##  14.2  Select columns
M4_Labs <- M4_Labs %>%
  select(M4_id, study_patient_id, LAB_TYPE, LAB_DATE, PURPOSE, LAB_NAME, LAB_VALUE)

##  14.3  Define extract_timepoint() (purpose → “01”, “03”, … “R”)
# Define the function to extract the timepoint
extract_timepoint <- function(purpose) {
  if (is.na(purpose) || purpose == "") return(NA)
  
  # Convert to lower case for case-insensitive matching
  purpose_lower <- tolower(purpose)
  
  # Remove leading and trailing whitespace
  purpose_lower <- trimws(purpose_lower)
  
  # Initialize timepoint as NA
  timepoint <- NA
  
  # Check for 'Diagnosis' related terms
  if (grepl("diagnosis|at diagnosis|baseline|at first consult|baseline line|dara baseline", purpose_lower)) {
    timepoint <- "01"
  } 
  # Check for 'Relapse'
  else if (grepl("relapse", purpose_lower)) {
    timepoint <- "R"
  } 
  else {
    # Try to extract visit number patterns
    match <- regexpr("v(isit)?[ .-]*#?[ ]*(\\d+|r)", purpose_lower, perl=TRUE)
    if (match[1] != -1) {
      visit_str <- regmatches(purpose_lower, match)
      visit_num <- sub("v(isit)?[ .-]*#?[ ]*", "", visit_str, perl=TRUE)
      visit_num <- toupper(visit_num)
      if (visit_num == "R") {
        timepoint <- "R"
      } else if (grepl("^\\d+$", visit_num)) {
        timepoint <- sprintf("%02d", as.integer(visit_num))
      }
    }
  }
  return(timepoint)
}

##  14.4  Apply extract_timepoint() → TIMEPOINT + map to descriptions
description_mapping <- c(
  "01" = "Diagnosis", "03" = "Post-induction", "05" = "Post-transplant",
  "07" = "1yr maintenance", "08" = "1.5yr maintenance", "09" = "2yr maintenance",
  "10" = "2.5yr maintenance", "11" = "3yr maintenance", "12" = "3.5yr maintenance",
  "13" = "4yr maintenance", "14" = "4.5yr maintenance", "15" = "5yr maintenance",
  "R"  = "Relapse"
)

M4_Labs <- M4_Labs %>%
  mutate(
    Patient         = M4_id,
    TIMEPOINT       = sapply(PURPOSE, extract_timepoint),
    TIMEPOINT_DESC  = description_mapping[TIMEPOINT]
  )

write.csv(M4_Labs, file = "M4_labs_cleaned.csv", row.names = FALSE)

##  14.5  Make M4_dates_cmrg with unique M4_id, LAB_DATE, TIMEPOINT, TIMEPOINT_DESC
M4_dates_cmrg <- M4_Labs %>%
  select(M4_id, LAB_DATE, PURPOSE, TIMEPOINT, TIMEPOINT_DESC) %>%
  distinct()


# ─── 15.  (Legacy) “M4 Status & best response” Excel → long ─────────────────────
data <- read_excel("Clinical data/M4/M4 Status and best response update Dec 2 24 just dates.xlsx")
names(data) <- make.names(names(data))
data <- data %>% select(-Best.Response.CR.or.VGPR, -Date)

colnames(data) <- c(
  "SITE", "PATIENT", "Current.Status", "CONSENT",
  "Timepoint_01_Prior_TX", "Timepoint_03_Post_Chemo",
  "Timepoint_05_Post_ASCT_100_or_1-3_DYS", "Timepoint_07_Post_Maintenance_Therapy_12_MONTHS",
  "Timepoint_08_6_Months_Post_07", "Timepoint_09_Every_6_Months",
  "Timepoint_10_Every_6_Months", "Timepoint_11_Every_6_Months",
  "Timepoint_12_Every_6_Months", "Timepoint_13_Every_6_Months",
  "Timepoint_R_Relapse"
)

# Transform data to long format and clean it to standardize timepoint names
data_long <- data %>%
  pivot_longer(
    cols      = -c(SITE, PATIENT, Current.Status),
    names_to  = "Timepoint",
    values_to = "Date"
  ) %>%
  mutate(
    Timepoint = str_trim(Timepoint),
    Timepoint = gsub("\\s+", " ", Timepoint),
    Date      = sapply(Date, function(x) {
      parsed <- parse_date_time(as.character(x), orders = c("dmy","mdy","ymd"), quiet = TRUE)
      as_date(parsed)
    }),
    Timepoint_Code = case_when(
      Timepoint == "Timepoint_01_Prior_TX"                  ~ "01",
      Timepoint == "Timepoint_03_Post_Chemo"                ~ "03",
      Timepoint == "Timepoint_05_Post_ASCT_100_or_1-3_DYS"   ~ "05",
      Timepoint == "Timepoint_07_Post_Maintenance_Therapy_12_MONTHS" ~ "07",
      Timepoint == "Timepoint_08_6_Months_Post_07"           ~ "08",
      Timepoint == "Timepoint_09_Every_6_Months"             ~ "09",
      Timepoint == "Timepoint_10_Every_6_Months"             ~ "10",
      Timepoint == "Timepoint_11_Every_6_Months"             ~ "11",
      Timepoint == "Timepoint_12_Every_6_Months"             ~ "12",
      Timepoint == "Timepoint_13_Every_6_Months"             ~ "13",
      Timepoint == "Timepoint_R_Relapse"                     ~ "R",
      Timepoint == "CONSENT"                                 ~ "CONSENT",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(PATIENT = str_pad(as.character(PATIENT), width = 2, pad = "0")) %>%
  filter(!is.na(Date), Timepoint_Code != "CONSENT") %>%
  rename(Patient_Code = PATIENT) %>%
  mutate(Patient = paste(SITE, Patient_Code, sep = "-")) %>%
  distinct()

M4_dates <- data_long
rm(data_long)

##  15.1  Convert “R” → “R-” in M4_dates
M4_dates <- M4_dates %>%
  mutate(Timepoint_Code = ifelse(Timepoint_Code == "R", "R-", Timepoint_Code))


# ─── 16.  Merge M4_dates into combined_clinical_data_updated ────────────────────
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(
    M4_dates %>% select(Patient, Timepoint_Code, Date),
    by = c("Patient", "Timepoint" = "Timepoint_Code")
  ) %>%
  mutate(
    Date_of_sample_collection = coalesce(
      normalize_clinical_date(Date),
      normalize_clinical_date(Date_of_sample_collection)
    )
  ) %>%
  select(-Date)

##  16.1  Write intermediate file + M4_dates to disk)
write.csv(M4_dates,
          "Exported_data_tables_clinical/M4_dates_Nov2024.csv",
          row.names = FALSE)


# ─── 17.  Use M4_Labs_earliest to fill any remaining missing M4 dates ───────────
M4_Labs_filtered <- M4_Labs %>%
  select(M4_id, PURPOSE, LAB_DATE, TIMEPOINT, TIMEPOINT_DESC) %>%
  filter(!grepl("Baseline|baseline", PURPOSE)) %>%
  distinct()

M4_Labs_filtered <- M4_Labs_filtered %>%
  mutate(Timepoint_Code = ifelse(TIMEPOINT == "R", "R-", TIMEPOINT)) %>%
  filter(!is.na(TIMEPOINT))

M4_Labs_earliest <- M4_Labs_filtered %>%
  group_by(M4_id, TIMEPOINT) %>%
  filter(LAB_DATE == min(LAB_DATE, na.rm = TRUE)) %>%
  ungroup()

##  17.1  Fix CA-08 mislabeled “08”
M4_Labs_earliest <- M4_Labs_earliest %>%
  mutate(Timepoint_Code = if_else(M4_id == "CA-08" & Timepoint_Code == "08", "09", Timepoint_Code))

##  17.2  Read additional relapse dates from Excel → Relapse_dates_M4_clean
Relapse_dates_M4 <- read_excel("M4_CMRG_Data/Relapse_dates.xlsx")

Relapse_dates_M4_clean <- Relapse_dates_M4 %>%
  select(SITE, PATIENT, `Current Status`, `CMRG Relapse Date`, `M4 Relapse Visit`) %>%
  rename(
    Site            = SITE,
    Patient         = PATIENT,
    Status          = `Current Status`,
    CMRG_Relapse_Date = `CMRG Relapse Date`,
    M4_Relapse_Visit  = `M4 Relapse Visit`
  ) %>%
  mutate(
    CMRG_Relapse_Date = as.numeric(CMRG_Relapse_Date),
    M4_Relapse_Visit  = as.numeric(M4_Relapse_Visit),
    CMRG_Relapse_Date = as.Date(CMRG_Relapse_Date, origin = "1899-12-30"),
    M4_Relapse_Visit  = as.Date(M4_Relapse_Visit,  origin = "1899-12-30")
  )

##  17.3  Correct error in sample IDs (diagnosis_bam_values, unknown_timepoint_bam)
diagnosis_bam_values   <- c(
  "TFRIM4_0166_Bm_P_WG_CA-06.filter.deduped.recalibrated.bam",
  "TFRIM4_0168_Bm_P_WG_CA-15.filter.deduped.recalibrated.bam",
  "TFRIM4_0169_Bm_P_WG_CA-16.filter.deduped.recalibrated.bam",
  "TFRIM4_0175_Bm_P_WG_EK-09.filter.deduped.recalibrated.bam"
)
unknown_timepoint_bam  <- "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam"

# Update Timepoint, timepoint_info, and Date_of_sample_collection accordingly
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    Timepoint = case_when(
      Bam %in% diagnosis_bam_values ~ "01",
      Bam == unknown_timepoint_bam ~ NA_character_,
      TRUE ~ Timepoint
    ),
    timepoint_info = ifelse(Bam %in% diagnosis_bam_values, "Diagnosis", timepoint_info),
    Date_of_sample_collection = if_else(
      Bam %in% c(diagnosis_bam_values, unknown_timepoint_bam),
      as.Date(NA),
      as.Date(Date_of_sample_collection)
    )
  )

## Edit the IDs to match
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(Sample_ID = ifelse(Bam %in% diagnosis_bam_values, gsub("-R-O", "-01-O", Sample_ID), Sample_ID))


## Now add to previous dataframe 
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(M4_dates %>% select(Patient, Timepoint_Code, Date), 
            by = c("Patient", "Timepoint" = "Timepoint_Code")) %>%
  mutate(
    Date_of_sample_collection = coalesce(
      normalize_clinical_date(Date),
      normalize_clinical_date(Date_of_sample_collection)
    )
  ) %>%
  select(-Date)  # Remove the Date column after updating

## If still NA, go off of labs
## Merge 
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(M4_Labs_earliest %>% dplyr::select(M4_id, Timepoint_Code, LAB_DATE),
            by = c("Patient" = "M4_id", "Timepoint" = "Timepoint_Code")) %>%
  # Update Date_of_sample_collection only if it is NA, using LAB_DATE where available.
  # Coerce both branches to Date to avoid RStudio/session-dependent Date/datetime
  # coercion when dplyr checks strict column types.
  mutate(Date_of_sample_collection = coalesce(
    normalize_clinical_date(Date_of_sample_collection),
    normalize_clinical_date(LAB_DATE)
  )) %>%
  # Drop the temporary LAB_DATE column added from the join
  select(-LAB_DATE) %>% 
  unique()

##  17.4  Merge M4_dates back (again) + then M4_Labs_earliest to fill any remaining NAs
## Get all unique dates for PFS
# Step 1: Separate SITE and Patient_Code in M4_Labs_earliest
tmp <- M4_Labs_earliest %>%
  separate(M4_id, into = c("SITE", "Patient_Code"), sep = "-")

# Step 2: Identify new combinations in M4_Labs_earliest not present in M4_dates
new_combinations <- tmp %>%
  rename(Date = LAB_DATE) %>%
  mutate(Date = normalize_clinical_date(Date)) %>%
  select(SITE, Patient_Code, Timepoint_Code, Date) %>%
  anti_join(M4_dates %>% select(SITE, Patient_Code, Timepoint_Code), by = c("SITE", "Patient_Code", "Timepoint_Code"))

# Step 3: Combine new combinations with M4_dates, prioritizing dates in M4_dates
combined_dates <- M4_dates %>%
  select(SITE, Patient_Code, Timepoint_Code, Date) %>%
  mutate(Date = normalize_clinical_date(Date)) %>%
  bind_rows(new_combinations)

# Step 4: Calculate the date range for each patient
patient_date_range <- combined_dates %>%
  group_by(SITE, Patient_Code) %>%
  summarise(
    Earliest_Date = min(Date, na.rm = TRUE),
    Latest_Date = max(Date, na.rm = TRUE),
    Days_Between = as.numeric(Latest_Date - Earliest_Date)
  ) %>%
  ungroup()

# Create the Patient column by concatenating SITE and Patient_Code
patient_date_range <- patient_date_range %>%
  mutate(Patient = paste(SITE, Patient_Code, sep = "-")) %>%
  select(Patient, everything())  # Optional: Move Patient column to the front

# Export the dataframe to a CSV file
write.csv(patient_date_range, "Oct2024_exported_data/patient_date_range_M4_PFS.csv", row.names = FALSE)


##  17.5 Update the VA-05 with updated info
va07_recode <- manual_clinical_correction_rows("va07_relapse_sample_recode")
target_bam <- va07_recode$Bam[[1]]

# Update Date_of_sample_collection, Timepoint, timepoint_info, and modify Sample_ID
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(Date_of_sample_collection = normalize_clinical_date(Date_of_sample_collection))

combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    Date_of_sample_collection = if_else(
      Bam == target_bam,
      va07_recode$Date_of_sample_collection[[1]],
      as.Date(Date_of_sample_collection)
    ),
    Timepoint = ifelse(Bam == target_bam, va07_recode$Timepoint[[1]], Timepoint),
    timepoint_info = ifelse(Bam == target_bam, va07_recode$timepoint_info[[1]], timepoint_info),
    Sample_ID = ifelse(
      Bam == target_bam,
      gsub(va07_recode$Sample_ID_find[[1]], va07_recode$Sample_ID_replace[[1]], Sample_ID),
      Sample_ID
    )
  )

combined_clinical_data_updated$Date_of_sample_collection <- as.Date(combined_clinical_data_updated$Date_of_sample_collection)
combined_clinical_data_updated$Date_of_sample_collection <- normalize_clinical_date(combined_clinical_data_updated$Date_of_sample_collection)


## 17.6 Correct IMG cases
# Set Plot_timepoint_numeric to 6 for the specific BAM entry
combined_clinical_data_updated$Timepoint[combined_clinical_data_updated$Bam == "MYC_0181_Pl_T_PG_MyP-181-T2-P.filter.deduped.recalibrated.bam"] <- "T6"

## 17.7 Fix FZ-07 BM sample in combined_clinical_data_updated
source_row <- combined_clinical_data_updated %>%
  filter(Bam == "TFRIM4_0178_Cf_P_PG_FZ-07-01-P-DNA.filter.deduped.recalibrated.bam")

new_row <- source_row
new_row$Bam         <- "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam"
new_row$Sample_ID   <- "FZ-07-01-O"
new_row$Sample_type <- "BM_cells"

combined_clinical_data_updated <- combined_clinical_data_updated %>%
  filter(Bam != "TFRIM4_0178_Bm_P_WG_FZ-07.filter.deduped.recalibrated.bam") %>%
  bind_rows(new_row)




# ─── 18.  Get progression info  ─────────────────

# 18.1 Calculate the progression info 
# Add Sample_ID column with Site and zero-padded Patient
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  mutate(Sample_ID = paste0(Site, "-", sprintf("%02d", Patient)))

### Add in missing dates of progressive disease from CMRG
CMRG_progressive_disease <- read_excel("M4_CMRG_Data/M4_COHORT_CHEMOTHERAPY.xlsx") %>%
  transmute(
    M4_id = as.character(M4_id),
    PROGRESSION_DATE = normalize_m4_endpoint_date(PROGRESSION_DATE),
    endpoint_source = "M4_COHORT_CHEMOTHERAPY.xlsx"
  )

march2026_chemo_progression <- read_march2026_csv_if_present("M4_COHORT_CHEMOTHERAPY.csv") %>%
  transmute(
    M4_id = as.character(M4_id),
    PROGRESSION_DATE = normalize_m4_endpoint_date(PROGRESSION_DATE),
    endpoint_source = "March 2026/M4_COHORT_CHEMOTHERAPY.csv"
  )

CMRG_progressive_disease <- bind_rows(
  CMRG_progressive_disease,
  march2026_chemo_progression
) %>%
  filter(
    !is.na(M4_id),
    stringr::str_detect(M4_id, "^[A-Z]{2}-[0-9]{2}$"),
    !is.na(PROGRESSION_DATE)
  ) %>%
  distinct(M4_id, PROGRESSION_DATE, .keep_all = TRUE)

readr::write_csv(
  CMRG_progressive_disease %>%
    arrange(M4_id, PROGRESSION_DATE),
  file.path(clinical_support_dir, "m4_progression_dates_with_march2026_augmentation.csv"),
  na = ""
)

# Perform a left join to bring PROGRESSION_DATE into Relapse_dates_M4_clean based on Sample_ID and M4_id
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  left_join(
    CMRG_progressive_disease %>%
      dplyr::select(M4_id, PROGRESSION_DATE),
    by = c("Sample_ID" = "M4_id")
  )

### Add another column called integrated progression date which puts the value of M4_Relapse Visit if Progression is NA
# Create the Integrated_Progression_Date column
Relapse_dates_M4_clean <- Relapse_dates_M4_clean %>%
  mutate(
    Integrated_Progression_Date = as.Date(ifelse(
      is.na(PROGRESSION_DATE),
      as.Date(M4_Relapse_Visit),
      as.Date(PROGRESSION_DATE))
    )
  )

## Add everything together 
# Step 1: Filter combined_spore_timepoint_info for rows with timepoint_info == "Progression"
spore_progression <- combined_spore_timepoint_info %>%
  filter(timepoint_info == "Progression") %>%
  rename(Progression_date = `Timepoint of interest`)

## Add additional dates 
# Step 1: Modify spore_OS_info
tmp <- spore_OS_info %>%
  mutate(
# Legacy fallback retained from the historical workflow: a deceased patient's
# status-last-follow-up date is appended as a candidate progression date. This
# is not an independently adjudicated progression event. Readers interpreting
# PFS should review this source alongside the later curated endpoint layers and
# patient-level endpoint audit outputs.
    timepoint_info = ifelse(`Current status` == "Deceased", "Progression", NA),
    Progression_date = `Status last follow up`
  ) %>%
  dplyr::select(Patient, Progression_date, timepoint_info) %>%
  filter(!is.na(timepoint_info))

# Step 2: Bind rows with spore_progression
spore_progression <- bind_rows(spore_progression, tmp)


# Step 2: Select Patient and PROGRESSION_DATE columns from Relapse_dates_M4_clean
m4_progression <- Relapse_dates_M4_clean %>%
  select(Sample_ID, Integrated_Progression_Date) %>%
  rename(Patient = Sample_ID, Progression_date = Integrated_Progression_Date)

# Step 3: Combine the two dataframes into Relapse_dates_full
Relapse_dates_full <- bind_rows(spore_progression, m4_progression) %>%
  filter(!is.na(Progression_date)) %>% select(-timepoint_info)

### Add the IMMAGINE dates here 
IMMAGINE_progression <- read_excel("Clinical data/IMMAGINE/Extracted_clinical_MRD_data.xlsx", sheet = 1)

# Fix dates
IMMAGINE_progression <- IMMAGINE_progression %>%
  # Split the 'Relapse Dates' column into separate rows
  separate_rows(`Relapse Dates`, sep = ", ") %>%
  # Rename columns to the desired names
  rename(Diagnosis_date = `Diagnosis Date`, Relapse_date = `Relapse Dates`) %>%
  # Convert 'Diagnosis_date' and 'Relapse_date' columns to Date format
  mutate(
    Diagnosis_date = as.Date(Diagnosis_date),
    Relapse_date = as.Date(Relapse_date, format = "%Y-%m-%d")
  )

## Join to relapse_dates 
temp <- IMMAGINE_progression %>% 
  select(ID, Relapse_date) %>%
  rename(Patient = ID, Progression_date = Relapse_date) %>% 
  filter(!is.na(Progression_date)) %>% 
  unique()

Relapse_dates_full <- bind_rows(Relapse_dates_full, temp)
Relapse_dates_full <- bind_rows(Relapse_dates_full, tmp %>% select(Patient, Progression_date))
## Now calculate the time to closest progression for each sample in the combined_clinical_info

# Ensure Date_of_sample_collection and Progression_date are in Date format
combined_clinical_data_updated$Date_of_sample_collection <- as.Date(combined_clinical_data_updated$Date_of_sample_collection)
Relapse_dates_full$Progression_date <- as.Date(Relapse_dates_full$Progression_date)

# Step 1: Add 'Relapsed' column
## Need to edit for those that relapsed before sample collected
test <- combined_clinical_data_updated %>%
  mutate(
    Relapsed = ifelse(Patient %in% Relapse_dates_full$Patient, "Y", "N")
  )

## Add the other dates got from David
new_dates <- manual_clinical_correction_rows(
  c("img060_additional_progression", "img098_additional_progression")
) %>%
  transmute(Patient, Progression_date)

# Add the new dates to the existing data frame
Relapse_dates_full <- bind_rows(Relapse_dates_full, new_dates) %>% unique()

## Add other dates in 
# load the additional dates from Esther
new_dates_esther <- read_excel(
  "Clinical data/IMMAGINE/Additional_relapse_dates_IMG_from_Esther.xlsx",
  sheet = 1,
  col_types = c("text", "text")
) %>%
  transmute(
    Patient          = `IMMAGINE ID`,
    Progression_date = PROGRESSION)

new_dates_esther <- new_dates_esther %>%
  mutate(
    Progression_date = as.Date(
      as.numeric(Progression_date),
      origin = "1899-12-30"
    )
  )

# 3) drop any existing entries for these patients
updated_relapse <- Relapse_dates_full %>%
  filter(!Patient %in% new_dates_esther$Patient)

# 4) append the new dates
updated_relapse <- bind_rows(updated_relapse, new_dates_esther) %>%
  arrange(Patient)

## Fix back to original name for consistency 
Relapse_dates_full <- updated_relapse

## Add additional dates from Sarah
new_rows <- manual_clinical_correction_rows(
  c("va24_additional_progression", "hp05_additional_progression")
) %>%
  transmute(Patient, Progression_date)

Relapse_dates_full <- bind_rows(Relapse_dates_full, new_rows)

oicr_revision_endpoint_path <- file.path(
  # The revision metadata contains endpoint/follow-up fields for Spring 2026
  # patients that may not exist in the older clinical workbooks. This file is
  # used here only to augment endpoint availability; sample-level integration is
  # handled centrally by read_combined_clinical_metadata_with_revision().
  "New OICR Submissions",
  "derived_metadata",
  "oicr_revision_repo_style_metadata.csv"
)

oicr_revision_endpoints <- if (file.exists(oicr_revision_endpoint_path)) {
  readr::read_csv(oicr_revision_endpoint_path, show_col_types = FALSE) %>%
    transmute(
      # Normalize all candidate endpoint dates immediately so downstream relapse
      # and follow-up logic works with Date objects rather than mixed Excel/text
      # formats.
      Patient = as.character(Patient),
      sample_collection_date = normalize_clinical_date(Date_of_sample_collection),
      date_diagnosis = normalize_clinical_date(date_diagnosis),
      first_progression_date = normalize_clinical_date(first_progression_date),
      latest_progression_date = normalize_clinical_date(latest_progression_date),
      relapse_or_censor_date = normalize_clinical_date(relapse_or_censor_date),
      relapse_or_censor_status = as.character(relapse_or_censor_status)
    ) %>%
    distinct()
} else {
  tibble(
    Patient = character(),
    sample_collection_date = as.Date(character()),
    date_diagnosis = as.Date(character()),
    first_progression_date = as.Date(character()),
    latest_progression_date = as.Date(character()),
    relapse_or_censor_date = as.Date(character()),
    relapse_or_censor_status = character()
  )
}

oicr_revision_progression_dates <- oicr_revision_endpoints %>%
  mutate(
    relapse_or_censor_progression_date = if_else(
      # relapse_or_censor_date is only a progression date when the status says
      # relapse/progression. Censored dates are follow-up endpoints, not events.
      str_detect(str_to_lower(coalesce(relapse_or_censor_status, "")), "relapse|progress"),
      relapse_or_censor_date,
      as.Date(NA)
    )
  ) %>%
  select(
    Patient,
    first_progression_date,
    latest_progression_date,
    relapse_or_censor_progression_date
  ) %>%
  pivot_longer(
    # Long format preserves the source of each appended progression date in the
    # audit output, which is important when first/latest/censor-derived dates
    # disagree.
    cols = -Patient,
    names_to = "revision_progression_date_source",
    values_to = "Progression_date"
  ) %>%
  filter(!is.na(Patient), !is.na(Progression_date)) %>%
  distinct(Patient, Progression_date, .keep_all = TRUE)

dir.create(clinical_support_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  oicr_revision_endpoints,
  file.path(clinical_support_dir, "oicr_revision_endpoint_dates_used.csv"),
  row.names = FALSE
)
write.csv(
  oicr_revision_progression_dates,
  file.path(clinical_support_dir, "oicr_revision_progression_dates_appended.csv"),
  row.names = FALSE
)

Relapse_dates_full <- bind_rows(
  Relapse_dates_full,
  oicr_revision_progression_dates %>% select(Patient, Progression_date)
) %>%
  distinct(Patient, Progression_date, .keep_all = TRUE)

# Recompute patient-level relapse availability after all manually curated
# progression dates have been appended. The first `test` object above is built
# before later David/Esther/Sarah date additions, so using it here would leave
# newly added relapse patients with stale sample-level `Relapsed = "N"` flags.
test <- combined_clinical_data_updated %>%
  mutate(
    Relapsed = ifelse(Patient %in% Relapse_dates_full$Patient, "Y", "N")
  )

## Export relapse dates 
write.csv(Relapse_dates_full, file = "Relapse dates cfWGS updated.csv", row.names = F)
write.csv(Relapse_dates_M4_clean, "Relapse_dates_M4_clean.csv", row.names = FALSE)


# 18.2  Find the maximum number of days available per patient 
# Join on Patient and calculate the days difference
# Only relapse dates AFTER the sample was collected, unless within 35 days
## Then count that as a relapse timepoint
Time_to_relapse2 <- test %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date_of_sample_collection) %>%
  mutate(
    # Historical sign convention: future progression is stored as
    # sample_date - progression_date and is therefore negative. The primary
    # Num_days_to_closest_relapse value exported below uses its absolute value.
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date_of_sample_collection ~ as.numeric(Date_of_sample_collection - Progression_date), 
      Progression_date <= Date_of_sample_collection & Progression_date >= Date_of_sample_collection - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    Relapsed = first(Relapsed),
    
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup()

## Now add this back 
tmp <- combined_clinical_data_updated
combined_clinical_data_updated <- combined_clinical_data_updated %>% 
  dplyr::select(-any_of(c(
    "Num_days_to_closest_relapse",
    "Relapsed",
    "Num_days_to_closest_relapse_absolute",
    "Num_days_to_closest_relapse_non_absolute"
  ))) %>%
  left_join(Time_to_relapse2)

## Add other timepoint_info
tmp <- read.csv("Clinical data/IMMAGINE/NAs in timepoint info - IMMAGINE additional timepoint info.csv")
tmp <- tmp %>% 
  rename(timepoint_info_updated = timepoint_info) %>% 
  select(-Date_of_sample_collection)

# Merge combined_clinical_data_updated with specimen_collection, replacing `timepoint_info` with `timepoint_info_updated` wherever there is a match
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(tmp, by = c("Patient", "Sample_ID")) %>%
  dplyr::mutate(
    # Replace `timepoint_info` with `timepoint_info_updated` wherever there is a match
    timepoint_info = if_else(!is.na(timepoint_info_updated), timepoint_info_updated, timepoint_info)
  ) %>%
  dplyr::select(-timepoint_info_updated)  # Drop `timepoint_info_updated` after merging

# Replace "MRD" with "Maintenance" in the `timepoint_info` column
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  mutate(
    timepoint_info = if_else(timepoint_info == "MRD", "Maintenance", timepoint_info)
  )


### Take the date of the relapse from the baseline marrow 
# Filter baseline bone marrow samples with "Diagnosis" or "baseline" in timepoint_info
## Edit the relapse dates of additional SPORE patients 
baseline_BM_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "BM_cells", timepoint_info %in% c("Diagnosis", "Baseline"))

baseline_with_relapse <- baseline_BM_samples

## Export baseline BM relapse-timing helper table for downstream QC.
dir.create("Oct2024_exported_data", showWarnings = FALSE, recursive = TRUE)
write.csv(baseline_with_relapse, file = "Oct2024_exported_data/Days to relapse baseline BM cell samples updated.txt", row.names = F)

combined_clinical_data_updated$Num_days_to_closest_relapse <- combined_clinical_data_updated$Num_days_to_closest_relapse_absolute

# Update Relapsed to "No" only if it is NA
combined_clinical_data_updated$Relapsed[is.na(combined_clinical_data_updated$Relapsed)] <- "N"

# View the result to confirm changes
table(combined_clinical_data_updated$Relapsed)


### Calculate longest follwowup time have on patients 
### Need additional data on immagine for this to see how long we have
# Step 1: Identify baseline marrow sample date for each patient
baseline_dates <- combined_clinical_data_updated %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  group_by(Patient) %>%
  summarise(Baseline_Date = min(Date_of_sample_collection, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Join baseline dates with main data
combined_with_baseline <- combined_clinical_data_updated %>% filter(Sample_type == "Blood_plasma_cfDNA") %>%
  left_join(baseline_dates, by = "Patient") %>%
  # Step 3: Calculate follow-up days only for dates after baseline
  mutate(
    Followup_Days = ifelse(Date_of_sample_collection > Baseline_Date,
                           as.numeric(Date_of_sample_collection - Baseline_Date),
                           NA)
  )

# Step 4: Summarize maximum follow-up days for each patient
patient_followup <- combined_with_baseline %>%
  group_by(Patient) %>%
  summarise(Max_Followup_Days = max(Followup_Days, na.rm = TRUE)) %>%
  ungroup()

patient_followup <- patient_followup %>% filter(Max_Followup_Days > 0)

## Now join with relapse info to get the PFS
Time_to_relapse2 <- combined_clinical_data_updated %>%
  filter(timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  left_join(
    Relapse_dates_full,
    by = c("Patient")
  ) %>%
  group_by(Patient, Date_of_sample_collection) %>%
  mutate(
    # Calculate non-absolute days to relapse only for valid dates
    days_to_relapse_non_absolute = case_when(
      Progression_date > Date_of_sample_collection ~ as.numeric(Date_of_sample_collection - Progression_date), 
      Progression_date <= Date_of_sample_collection & Progression_date >= Date_of_sample_collection - 35 ~ 0,
      TRUE ~ NA_real_  # Set to NA for dates that do not meet conditions
    )
  ) %>%
  filter(!is.na(days_to_relapse_non_absolute)) %>%
  summarise(
    Relapsed = first(Relapsed),
    
    # Calculate the absolute days to relapse for valid dates
    Num_days_to_closest_relapse_absolute = min(abs(days_to_relapse_non_absolute), na.rm = TRUE),
    
    # Get the closest relapse date that came strictly after the sample collection
    Num_days_to_closest_relapse_non_absolute = days_to_relapse_non_absolute[which.min(abs(days_to_relapse_non_absolute))]
  ) %>%
  ungroup()

# Step to keep the highest Num_days_to_closest_relapse_absolute for each patient
Time_to_relapse2 <- Time_to_relapse2 %>% 
  group_by(Patient) %>%
  slice_max(Num_days_to_closest_relapse_absolute, with_ties = FALSE) %>%
  ungroup()

## Now add this back 
patient_followup <- left_join(patient_followup, Time_to_relapse2)

### Have to think about multiple relapse cases, what to do for that for the kaplan myer curves 
# Keep only the max num_days for each patient
max_num_days_per_patient <- patient_followup %>%
  group_by(Patient) %>%
  summarise(
    max_num_days = max(Max_Followup_Days, na.rm = TRUE),
    max_num_weeks = max_num_days / 7
  ) %>%
  ungroup()


## Export
write.table(max_num_days_per_patient, file = "Oct2024_exported_data/max_num_days_per_patient_Jan28_updated.txt", sep = "\t", quote = F, row.names = F)



### 18.3 Clean processing log to get additional dates and export 
file_path <- "TFRIM4_Processing Log_Nov2024.xlsx"  # Replace with the input file path
processing_log <- read_excel(file_path)

# Step 3: Select only the first three columns
processing_log_cleaned <- processing_log %>%
  select(
    `Plasma Pugh Lab Identifier`,
    `Red is on 35-Patient list`,
    `Blood Drawn Date [YYYY/MM/DD]`
  )

# Ensure naming is consistent
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(`Plasma Pugh Lab Identifier` = str_replace_all(`Plasma Pugh Lab Identifier`, "_", "-"))


# Step 4: Rename and process Sample_Code
processing_log_cleaned <- processing_log_cleaned %>%
  rename(Sample_Code = `Plasma Pugh Lab Identifier`) %>%  # Rename the column
  mutate(
    Sample_Code = str_remove(Sample_Code, "M4-"),  # Remove "M4-" prefix
    Sample_Code = str_remove_all(Sample_Code, "(_[0-9]|-[0-9])$"),  # Remove endings like -1, -2, _1, _2
    # Extract the true M4 patient code from the sample identifier. A historical
    # processing-log typo (`MA-VA-25-01-2`) previously became the invalid patient
    # `MA-VA` when the first five characters were used. Regex extraction keeps
    # valid IDs such as VA-06 and recovers VA-25 from that malformed sample code.
    Patient = str_extract(Sample_Code, "[A-Z]{2}-[0-9]{2}"),
    Timepoint = str_extract(Sample_Code, "[^-]+$")  # Extract everything after the last "-"
  )

## Make date consistent
# Define a function to process mixed date formats
process_dates <- function(date_column) {
  date_column <- gsub("\\?", NA, date_column) # Replace question marks with NA
  date_column <- trimws(date_column)  # Remove leading/trailing whitespace
  
  processed_dates <- sapply(date_column, function(date) {
    if (is.na(date) || date == "") {
      return(NA)
    } else if (grepl("^\\d{5}$", date)) {
      # Convert Excel serial numbers (Julian dates)
      as.character(as.Date(as.numeric(date), origin = "1899-12-30"))
    } else if (grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) {
      # ISO 8601 format: YYYY-MM-DD
      as.character(as.Date(date))
    } else if (grepl("^\\d{2}/\\d{1,2}/\\d{4}$", date)) {
      # Format: DD/MM/YYYY or D/M/YYYY
      as.character(dmy(date))
    } else if (grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", date)) {
      # Format: D-M-YYYY or DD-MM-YYYY
      as.character(dmy(date))
    } else if (grepl("^\\d{4}/\\d{1,2}/\\d{1,2}$", date)) {
      # Format: YYYY/MM/DD or YYYY/M/D
      as.character(ymd(date))
    } else {
      # Invalid or unknown format
      return(NA)
    }
  })
  
  return(processed_dates)
}

# Apply the function to clean the dates
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(
    Cleaned_Date = process_dates(`Blood Drawn Date [YYYY/MM/DD]`),
    Cleaned_Date = as.Date(Cleaned_Date)  # Convert to Date type
  )

## Fix +1 Notation 
processing_log_cleaned <- processing_log_cleaned %>%
  mutate(Timepoint = str_replace_all(Timepoint, "\\+.*", ""))

## Extract relevant columns 
# Process the data
processing_log_cleaned <- processing_log_cleaned %>%
  select(Patient, Timepoint, Cleaned_Date) %>%  # Keep only the specified columns
  distinct() %>%  # Remove duplicate rows
  mutate(Sample_type = "Blood_plasma_cfDNA") %>% # Add the new column with the specified value 
  filter(!is.na(Cleaned_Date))

### Export
# Ensure the folder exists
dir.create("Exported_data_tables_clinical", showWarnings = FALSE)

# Export the processed data as a CSV file
write.csv(processing_log_cleaned, file = "Exported_data_tables_clinical/M4_processing_log_dates.csv", row.names = FALSE)


### 18.4 Get fill dates table for PFS

## Make M4_dates_full to use for the PFS table 
### Join all of the ones available 
# 1) Start with the values in M4_dates

# Replace 'R-' with 'R' in the Timepoint column
M4_dates <- M4_dates %>%
  mutate(Timepoint_Code = gsub("^R-", "R", Timepoint_Code))

df1 <- M4_dates %>%
  select(
    Patient = Patient, 
    Timepoint = Timepoint_Code, 
    Date = Date
  ) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Ensure Date is of Date class
df1$Date <- as.Date(df1$Date)

# 2) Add additional timepoints or patients from M4_dates_cmrg not included in M4_dates
df2 <- M4_dates_cmrg %>%
  select(
    Patient = M4_id, 
    Timepoint = TIMEPOINT, 
    Date = LAB_DATE
  ) %>%
  mutate(
    Date = as.Date(Date), 
    Timepoint = as.character(Timepoint),
    Timepoint = unname(Timepoint)  # Remove names attribute from Timepoint
  ) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Find entries in df2 not already in df1
df2_new <- anti_join(df2, df1, by = c("Patient", "Timepoint"))

# Combine df1 and df2_new
df_combined <- bind_rows(df1, df2_new)

# 3) Add remaining entries from processing_log_cleaned not included in previous data
df3 <- processing_log_cleaned %>%
  select(
    Patient = Patient, 
    Timepoint = Timepoint, 
    Date = Cleaned_Date
  ) %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Patient, Timepoint, Date) %>%
  group_by(Patient, Timepoint) %>%
  slice(1) %>%  # Take the first date if multiple exist
  ungroup()

# Find entries in df3 not already in df_combined
df3_new <- anti_join(df3, df_combined, by = c("Patient", "Timepoint"))

# Combine all data into the super table
df_combined <- bind_rows(df_combined, df3_new)

df_combined <- df_combined %>% filter(!is.na(Timepoint))

## Now make a table with the longest dates have per pt and then the time to relapse for PGS

## Use this for the PFS curves
# Calculate the number of days from the baseline (Timepoint == "01") for each sample
df_combined <- df_combined %>%
  group_by(Patient) %>%
  mutate(Days_from_baseline = as.numeric(Date - min(Date[Timepoint == "01"], na.rm = TRUE))) %>%
  ungroup()

# Determine the maximum number of days to a sample for each patient
max_days <- df_combined %>%
  group_by(Patient) %>%
  summarize(Max_Days = max(Days_from_baseline, na.rm = TRUE)) %>%
  ungroup()

# Identify if the patient relapsed and the number of days to the relapse sample (Timepoint == "R")
relapse_info <- df_combined %>%
  filter(Timepoint == "R", !is.na(Days_from_baseline)) %>%
  group_by(Patient) %>%
  summarize(Relapse_Days = first(Days_from_baseline), Relapsed = 1) %>%
  ungroup()


# Merge the maximum days and relapse information
PFS_days <- max_days %>%
  left_join(relapse_info, by = "Patient") %>%
  mutate(
    Relapse_Days = ifelse(is.na(Relapse_Days), 0, Relapse_Days),
    Relapsed = ifelse(is.na(Relapsed), 0, Relapsed),
    Final_Days = ifelse(Relapsed == 1, Relapse_Days, Max_Days)  # Add a new column for the final days
  )

# Export the PFS_dates table to a CSV file
# No longer used
write.csv(PFS_days, file = "Exported_data_tables_clinical/PFS_days.csv", row.names = FALSE)

### Get the latest date per patient 
latest_dates <- df_combined %>%
  group_by(Patient) %>%
  summarise(latest_date = max_clinical_date_or_na(Date), .groups = "drop")

#### Add SPORE to this and IMG 
# Extract the Patient and 'Status last follow up' from spore_OS_info,
# renaming it to latest_date
tmp <- spore_OS_info %>%
  select(Patient, latest_date = `Status last follow up`)

## Add additional IMG last-follow-up information from Esther.
##
## Important distinction:
##   The PFS table uses `censor_date` as the first progression/censoring date
##   for survival analysis. That is not always the same as the last date a
##   relapsed patient was known to be followed clinically. For prospective
##   time-window analyses in the non-frontline/test cohort, downstream scripts
##   need the latter: a patient-level "known through" date. Esther's workbook
##   provides that date explicitly for updated IMMAGINE cases such as IMG-098.
tmp_IMG <- read_excel(
  "Clinical data/IMMAGINE/Additional_relapse_dates_IMG_from_Esther.xlsx",
  sheet = 2,
  col_types = c("text","text")
) %>%
  transmute(
    Patient     = `IMMAGINE ID`,
    latest_date = `Date of last Follow-up`
  )

tmp_IMG <- tmp_IMG %>%
  mutate(
    latest_date = normalize_clinical_date(latest_date),
    followup_source = "Clinical data/IMMAGINE/Additional_relapse_dates_IMG_from_Esther.xlsx:Date_of_last_followup"
  )

## Additional cleaned IMMAGINE follow-up table used for recovered WGS-baseline
## cases not represented in Esther's relapse workbook. This is deliberately
## follow-up/censor information only; relapse events are handled separately.
tmp_IMG_cleaned_followup <- readr::read_csv(
  "Clinical data/IMMAGINE/Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv",
  show_col_types = FALSE
) %>%
  transmute(
    Patient = Patient_ID,
    latest_date = normalize_clinical_date(Last_Followup_Date),
    followup_source = "Clinical data/IMMAGINE/Cleaned_Patient_Follow-Up_Table_IMMAGINE.csv:Last_Followup_Date"
  ) %>%
  filter(!is.na(Patient), !is.na(latest_date))

tmp_oicr_revision_followup <- oicr_revision_endpoints %>%
  rowwise() %>%
  mutate(
    # For prospective time-window analyses, use the latest known clinical
    # endpoint among censor date and progression dates. This is a "known through"
    # date, not necessarily a PFS event date.
    latest_date = max_clinical_date_or_na(c(
      relapse_or_censor_date,
      latest_progression_date,
      first_progression_date,
      sample_collection_date
    ))
  ) %>%
  ungroup() %>%
  transmute(
    Patient,
    latest_date,
    followup_source = case_when(
      # Record which revision field supported the follow-up date so downstream
      # audits can distinguish censored follow-up from progression-derived
      # follow-up.
      !is.na(sample_collection_date) & sample_collection_date == latest_date ~
        "Spring 2026 OICR revision metadata:Date_of_sample_collection",
      str_detect(str_to_lower(coalesce(relapse_or_censor_status, "")), "censor") ~
        "Spring 2026 OICR revision metadata:relapse_or_censor_date",
      !is.na(latest_progression_date) ~
        "Spring 2026 OICR revision metadata:latest_progression_date",
      !is.na(first_progression_date) ~
        "Spring 2026 OICR revision metadata:first_progression_date",
      TRUE ~ "Spring 2026 OICR revision metadata:no endpoint date"
    )
  ) %>%
  filter(!is.na(Patient))


# Combine the two tables
latest_dates <- bind_rows(
  latest_dates %>% mutate(followup_source = "latest clinical/sample date available in combined metadata"),
  tmp_IMG,
  tmp_IMG_cleaned_followup,
  tmp_oicr_revision_followup
)

latest_dates <- latest_dates %>%
  group_by(Patient) %>%
  summarise(
    latest_date = max_clinical_date_or_na(latest_date),
    followup_source = paste(sort(unique(followup_source[!is.na(latest_date)])), collapse = "; "),
    .groups = "drop"
  )

## Export this 
write.csv(latest_dates, file = "Exported_data_tables_clinical/latest_dates_per_patient.csv", row.names = FALSE)

## Export an explicitly named follow-up table for downstream scripts that need
## evaluability windows rather than PFS event/censor dates. Keeping this separate
## prevents relapsed patients' PFS event dates from being overwritten by later
## last-contact dates.
patient_followup_dates_updated <- latest_dates %>%
  transmute(
    Patient,
    followup_end_date = latest_date,
    followup_source
  )
write.csv(
  patient_followup_dates_updated,
  file = "Exported_data_tables_clinical/patient_followup_dates_updated.csv",
  row.names = FALSE
)
saveRDS(
  patient_followup_dates_updated,
  file = "Exported_data_tables_clinical/patient_followup_dates_updated.rds"
)




#### Re-do with IMG and update method
# Helper to extract one source
extract_latest <- function(df, pt_col, tp_col, date_col){
  df %>% 
    transmute(Patient = .data[[pt_col]],
              Timepoint = as.character(.data[[tp_col]]),
              Date = as.Date(.data[[date_col]])) %>% 
    filter(!is.na(Timepoint)) %>% 
    group_by(Patient, Timepoint) %>% 
    summarise(Date = max_clinical_date_or_na(Date), .groups="drop")
}

# 1) pull from each source
d1 <- extract_latest(M4_dates,       "Patient",       "Timepoint_Code", "Date")
d2 <- extract_latest(M4_dates_cmrg,  "M4_id",         "TIMEPOINT",      "LAB_DATE")
d3 <- extract_latest(processing_log_cleaned, "Patient","Timepoint",      "Cleaned_Date")
d4 <- extract_latest(combined_clinical_data_updated, "Patient","Timepoint",      "Date_of_sample_collection")

### Add the other CMRG table to try get latest date.
M4_DEMO_current <- read_excel("M4_CMRG_Data/M4_COHORT_DEMO.xlsx") %>%
  transmute(
    Patient = as.character(M4_id),
    Date = normalize_m4_endpoint_date(DATE_OF_LAST_FOLLOWUP),
    endpoint_source = "M4_COHORT_DEMO.xlsx"
  )

M4_DEMO_march2026 <- read_march2026_csv_if_present("M4_COHORT_DEMO.csv") %>%
  transmute(
    Patient = as.character(M4_id),
    Date = normalize_m4_endpoint_date(DATE_OF_LAST_FOLLOWUP),
    endpoint_source = "March 2026/M4_COHORT_DEMO.csv"
  )

M4_DEMO <- bind_rows(M4_DEMO_current, M4_DEMO_march2026) %>%
  filter(
    !is.na(Patient),
    stringr::str_detect(Patient, "^[A-Z]{2}-[0-9]{2}$"),
    !is.na(Date)
  )

M4_DEMO_followup_march2026_comparison <- full_join(
  M4_DEMO_current %>%
    filter(!is.na(Patient), stringr::str_detect(Patient, "^[A-Z]{2}-[0-9]{2}$")) %>%
    select(Patient, current_demo_followup = Date),
  M4_DEMO_march2026 %>%
    filter(!is.na(Patient), stringr::str_detect(Patient, "^[A-Z]{2}-[0-9]{2}$")) %>%
    select(Patient, march2026_demo_followup = Date),
  by = "Patient"
) %>%
  mutate(
    march2026_followup_change = case_when(
      is.na(current_demo_followup) & !is.na(march2026_demo_followup) ~ "new_march_followup",
      !is.na(current_demo_followup) & is.na(march2026_demo_followup) ~ "missing_in_march_demo",
      !is.na(current_demo_followup) & !is.na(march2026_demo_followup) &
        march2026_demo_followup > current_demo_followup ~ "march_followup_later",
      !is.na(current_demo_followup) & !is.na(march2026_demo_followup) &
        march2026_demo_followup < current_demo_followup ~ "march_followup_earlier",
      !is.na(current_demo_followup) & !is.na(march2026_demo_followup) &
        march2026_demo_followup == current_demo_followup ~ "same_followup",
      TRUE ~ "no_followup_in_either_demo"
    ),
    followup_days_added = as.numeric(march2026_demo_followup - current_demo_followup)
  ) %>%
  arrange(desc(march2026_followup_change == "march_followup_later"), Patient)

readr::write_csv(
  M4_DEMO_followup_march2026_comparison,
  file.path(clinical_support_dir, "m4_demo_march2026_followup_comparison.csv"),
  na = ""
)

readr::write_csv(
  M4_DEMO %>%
    arrange(Patient, Date, endpoint_source),
  file.path(clinical_support_dir, "m4_last_followup_dates_with_march2026_augmentation.csv"),
  na = ""
)

M4_DEMO <- M4_DEMO %>%
  select(Patient, Date) %>%
  group_by(Patient) %>%
  summarise(Date = max_clinical_date_or_na(Date), .groups = "drop") %>%
  ungroup() %>%
  distinct()

## Check if later in DEMO file
other_dates <- bind_rows(d1, d2, d3, d4) %>%
  group_by(Patient) %>%
  summarise(latest_other_date = max_clinical_date_or_na(Date), .groups = "drop")

# 2. Compare with M4_DEMO
more_recent_in_demo <- M4_DEMO %>%
  left_join(other_dates, by = "Patient") %>%
  filter(!is.na(latest_other_date), Date > latest_other_date)

# Print patients whose demographic-table follow-up is later than all other M4
# date sources. This is console QC; no endpoint is changed in this statement.
more_recent_in_demo


# 2) combine and get per-patient max
latest_dates <- bind_rows(d1, d2, d3, d4, M4_DEMO) %>% 
  group_by(Patient) %>% 
  summarise(latest_date = max_clinical_date_or_na(Date), .groups="drop")

# 3) incorporate SPORE follow-up
spore_dates <- spore_OS_info %>% 
  transmute(Patient, status_date = as.Date(`Status last follow up`))

latest_dates <- latest_dates %>% 
  left_join(spore_dates, by="Patient") %>% 
  bind_rows(tmp_IMG, tmp_IMG_cleaned_followup, tmp_oicr_revision_followup) %>% 
  rowwise() %>%
  mutate(latest_date = max_clinical_date_or_na(c(latest_date, status_date))) %>%
  ungroup() %>%
  select(Patient, latest_date) %>%
  group_by(Patient) %>%
  summarise(latest_date = max_clinical_date_or_na(latest_date), .groups = "drop") %>%
  distinct()


# 4) write out
write.csv(latest_dates,
          "Exported_data_tables_clinical/latest_dates_per_patient_updated.csv",
          row.names = FALSE)


# ─── 19.  Filter/compare patient sets for cfDNA vs. BM longitudinal → export counts ─

df <- combined_clinical_data_updated

# 19.1  Identify who has ≥2 distinct cfDNA timepoints (“longitudinal cfDNA”)
longitudinal_cfDNA <- df %>%
  filter(Sample_type=="Blood_plasma_cfDNA") %>%
  group_by(Patient) %>%
  summarise(n_timepoints=n_distinct(Timepoint)) %>%
  filter(n_timepoints>=2) %>%
  pull(Patient)

# 19.2  Baseline/Diagnosis BM & cfDNA sets, Progression/Relapse BM & cfDNA sets
baseline_BM       <- df %>% filter(Sample_type=="BM_cells", timepoint_info=="Diagnosis")         %>% distinct(Patient) %>% pull(Patient)
baseline_cfDNA    <- df %>% filter(Sample_type=="Blood_plasma_cfDNA", timepoint_info=="Diagnosis") %>% distinct(Patient) %>% pull(Patient)
progression_BM    <- df %>% filter(Sample_type=="BM_cells", timepoint_info %in% c("Progression","Relapse")) %>% distinct(Patient) %>% pull(Patient)
progression_cfDNA <- df %>% filter(Sample_type=="Blood_plasma_cfDNA", timepoint_info %in% c("Progression","Relapse")) %>% distinct(Patient) %>% pull(Patient)

# 19.3  Helper to count intersections
count_intersect <- function(A,B) length(intersect(A,B))

n_a <- count_intersect(baseline_BM, longitudinal_cfDNA)
n_b <- count_intersect(baseline_cfDNA, longitudinal_cfDNA)
n_c <- count_intersect(progression_BM, longitudinal_cfDNA)
n_d <- count_intersect(progression_cfDNA, longitudinal_cfDNA)

no_baseline_BM    <- setdiff(progression_BM, baseline_BM)
n_e <- count_intersect(no_baseline_BM, longitudinal_cfDNA)

no_baseline_cfDNA <- setdiff(progression_cfDNA, baseline_cfDNA)
n_f <- count_intersect(no_baseline_cfDNA, longitudinal_cfDNA)

no_baseline_BM_cfDNA <- setdiff(baseline_cfDNA, baseline_BM)
n_g <- count_intersect(no_baseline_BM_cfDNA, longitudinal_cfDNA)

# 19.4  Build results + detailed lists → write CSVs
results <- tibble(
  description = c(
    "Baseline/diagnosis BM + longitudinal cfDNA",
    "Baseline/diagnosis cfDNA + longitudinal cfDNA",
    "Progression/relapse BM + longitudinal cfDNA",
    "Progression/relapse cfDNA + longitudinal cfDNA",
    "Prog/relapse BM BUT NOT baseline BM + longitudinal cfDNA",
    "Prog/relapse cfDNA BUT NOT baseline cfDNA + longitudinal cfDNA",
    "Baseline cfDNA BUT NOT baseline BM + longitudinal cfDNA"
  ),
  count = c(n_a, n_b, n_c, n_d, n_e, n_f, n_g)
)

find_intersect <- function(A,B) intersect(A,B)

patients_a <- find_intersect(baseline_BM, longitudinal_cfDNA)
patients_b <- find_intersect(baseline_cfDNA, longitudinal_cfDNA)
patients_c <- find_intersect(progression_BM, longitudinal_cfDNA)
patients_d <- find_intersect(progression_cfDNA, longitudinal_cfDNA)
patients_e <- find_intersect(setdiff(progression_BM, baseline_BM), longitudinal_cfDNA)
patients_f <- find_intersect(setdiff(progression_cfDNA, baseline_cfDNA), longitudinal_cfDNA)
patients_g <- find_intersect(no_baseline_BM_cfDNA, longitudinal_cfDNA)

detailed_results <- tibble(
  description = rep(results$description, times = c(length(patients_a), length(patients_b),
                                                   length(patients_c), length(patients_d),
                                                   length(patients_e), length(patients_f),
                                                   length(patients_g))),
  patient     = c(patients_a, patients_b, patients_c, patients_d, patients_e, patients_f, patients_g)
)

write.csv(results,          file.path(clinical_support_dir, "patient_counts_summary.csv"),  row.names = FALSE)
write.csv(detailed_results, file.path(clinical_support_dir, "patient_counts_detailed.csv"),  row.names = FALSE)


# ─── 20.  Filter BAM lists (IMG only, SPORE only) → write out CSVs ─────────────
bams_filtered_IMG <- combined_clinical_data_updated %>%
  filter(grepl("IMG", Patient, ignore.case = TRUE)) %>%
  select(Bam, Patient, Sample_ID)
write.csv(bams_filtered_IMG, file.path(clinical_support_dir, "bams_with_IMG_patients.csv"), row.names = FALSE)

bams_filtered_SPORE <- combined_clinical_data_updated %>%
  filter(grepl("SPORE", Patient, ignore.case = TRUE)) %>%
  select(Bam, Patient, Sample_ID)
write.csv(bams_filtered_SPORE, file.path(clinical_support_dir, "SPORE_bams.csv"), row.names = FALSE)

cat("Filtered BAM lists exported under ", clinical_support_dir, "\n", sep = "")


# ─── 21.  Final tweaks & write combined_clinical_data_updated_April2025.csv ────
# This is the primary clinical-metadata output consumed by downstream scripts.
# When adding new samples, regenerate this file before rerunning feature
# integration, clinical-demographics tables, swim plots, and survival analyses.
combined_clinical_data_updated <- combined_clinical_data_updated %>%
  left_join(
    manual_clinical_correction_rows("spore0007_timepoint_info") %>%
      transmute(
        Patient,
        Date_of_sample_collection = Date_of_sample_collection,
        manual_timepoint_info = timepoint_info
      ),
    by = c("Patient", "Date_of_sample_collection")
  ) %>%
  mutate(timepoint_info = coalesce(manual_timepoint_info, timepoint_info)) %>%
  select(-manual_timepoint_info)

write.csv(combined_clinical_data_updated,
          "combined_clinical_data_updated_April2025.csv",
          row.names = FALSE)


# ─── 22.  Build “baseline_dates_all” (SPORE, M4, IMMAGINE) ───────────────────────
###  22.1  SPORE: pick the first BM “Baseline” for SPORE_0009, else first “Diagnosis”
baseline_0009 <- combined_spore_timepoint_info %>%
  filter(Patient=="SPORE_0009", timepoint_info=="Baseline") %>%
  slice_min(`Timepoint of interest`, with_ties=FALSE)

diagnosis_others <- combined_spore_timepoint_info %>%
  filter(Patient!="SPORE_0009", timepoint_info=="Diagnosis") %>%
  group_by(Patient) %>%
  slice_min(`Timepoint of interest`, with_ties=FALSE) %>%
  ungroup()

spore_baseline_dates <- bind_rows(baseline_0009, diagnosis_others)

###  22.2  M4: from M4_diagnosis_date (rename M4_id→Patient)
diagnosis_dates_m4 <- M4_diagnosis_date %>%
  select(Patient = M4_id, Diagnosis_date = DIAGNOSIS_DATE) %>%
  filter(!is.na(Diagnosis_date)) %>%
  distinct()

###  22.2b  M4 fallback: use Timepoint_01/Prior_TX when staging lacks a patient
m4_timepoint01_baseline_dates <- M4_dates %>%
  filter(Timepoint_Code == "01", !is.na(Date)) %>%
  transmute(Patient, Diagnosis_date = as.Date(Date)) %>%
  arrange(Patient, Diagnosis_date) %>%
  group_by(Patient) %>%
  slice(1) %>%
  ungroup() %>%
  anti_join(diagnosis_dates_m4 %>% distinct(Patient), by = "Patient")

###  22.3  IMMAGINE: from IMMAGINE_progression → Diagnosis_date
diagnosis_dates_immagine <- IMMAGINE_progression %>%
  select(Patient = ID, Diagnosis_date) %>%
  filter(!is.na(Diagnosis_date)) %>%
  distinct()

###  22.4  Bind all three sources, drop NA
spore_baseline_clean <- spore_baseline_dates %>%
  select(Patient, Diagnosis_date = `Timepoint of interest`) %>%
  distinct()

baseline_dates_all <- bind_rows(
  spore_baseline_clean,
  diagnosis_dates_m4,
  m4_timepoint01_baseline_dates,
  diagnosis_dates_immagine,
  oicr_revision_endpoints %>%
    select(Patient, Diagnosis_date = date_diagnosis) %>%
    filter(!is.na(Diagnosis_date))
) %>%
  distinct() %>%
  filter(!is.na(Diagnosis_date))

n_distinct(baseline_dates_all$Patient)

baseline_dates_all$Baseline_Date <- baseline_dates_all$Diagnosis_date

## WGS-baseline overrides for recovered IMG cases. These patients have historical
## disease courses before the WGS sample; using original diagnosis dates would
## make survival baselines scientifically wrong for the WGS MRD comparison.
img_wgs_baseline_overrides <- recoverable_img_wgs_sample_dates %>%
  transmute(Patient, Diagnosis_date = Date_of_sample_collection)

baseline_dates_all <- baseline_dates_all %>%
  filter(!Patient %in% img_wgs_baseline_overrides$Patient) %>%
  bind_rows(img_wgs_baseline_overrides) %>%
  arrange(Patient, Diagnosis_date) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  mutate(Baseline_Date = Diagnosis_date)

dir.create("exported_clinical_data_April2025", showWarnings = FALSE)
write.csv(baseline_dates_all,
          file.path("exported_clinical_data_April2025", "all_baseline_dates.csv"),
          row.names = FALSE)


### Now re-do the PFS table 
# 1) Relapsed patients

## Remove EK-09 since progression was from smouldering MM
Relapse_dates_full <- Relapse_dates_full %>% filter(Patient != "EK-09")

prog_info <- Relapse_dates_full %>%
  # 1) coerce to Date
  mutate(Progression_date = as.Date(Progression_date)) %>%
  # 2) bring in each patient’s baseline
  left_join(baseline_dates_all %>% select(Patient, Baseline_Date),
            by = "Patient") %>%
  # 3) drop any prog that happened before that patient’s baseline
  filter(is.na(Baseline_Date) | Progression_date >= Baseline_Date) %>%
  # 4) now group & collect
  group_by(Patient) %>%
  summarise(
    all_prog     = list(sort(Progression_date)),
    .groups = "drop"
  ) %>%
  # 5) re-join baseline so you can pick first post-baseline event
  left_join(baseline_dates_all, by = "Patient") %>%
  mutate(
    prog_after_baseline = all_prog,          # already filtered above
    Censor_date             = map_dbl(prog_after_baseline, ~ .x[1]),
    Other_Progression_Dates = map(prog_after_baseline, 
                                  ~ if(length(.x)>1) .x[-1] else as.Date(NA)),
    Relapsed      = 1,
    High_Priority = as.integer(Patient %in% cohort_df$Patient)
  ) %>%
  select(Patient, Censor_date, Other_Progression_Dates, Relapsed, High_Priority)

# 2) Non-relapsed patients
nonprog_info_pre_baseline_censor_audit <- latest_dates %>%
  # drop anyone who has at least one progression event
  filter(! Patient %in% prog_info$Patient) %>%
  left_join(
    baseline_dates_all %>% select(Patient, Baseline_Date),
    by = "Patient"
  ) %>%
  rename(Censor_date = latest_date) %>%
  mutate(
    Censor_date = as.Date(Censor_date),
    pre_baseline_nonprogress_censor = !is.na(Baseline_Date) &
      !is.na(Censor_date) &
      Censor_date < Baseline_Date
  )

readr::write_csv(
  nonprog_info_pre_baseline_censor_audit %>%
    filter(pre_baseline_nonprogress_censor) %>%
    arrange(Patient),
  file.path(clinical_support_dir, "nonrelapsed_pre_baseline_censor_dates_removed.csv"),
  na = ""
)

nonprog_info <- nonprog_info_pre_baseline_censor_audit %>%
  mutate(
    # A non-relapsed/censored endpoint before the analysis baseline creates
    # negative PFS and is not scientifically interpretable for baseline WGS
    # survival analyses. Keep the patient visible, but remove that invalid
    # censor date unless another later source is available upstream.
    Censor_date = if_else(
      pre_baseline_nonprogress_censor,
      as.Date(NA),
      Censor_date
    ),
    Other_Progression_Dates = as.Date(NA),    # no “other” dates
    Relapsed                = 0,
    High_Priority           = as.integer(Patient %in% cohort_df$Patient)
  ) %>%
  select(Patient, Censor_date, Other_Progression_Dates, Relapsed, High_Priority)


# 3) Combine
# If prog_info$Censor_date came out numeric, convert it back:
prog_info <- prog_info %>%
  mutate(
    Censor_date = as.Date(Censor_date, origin = "1970-01-01")
  )

# Ensure nonprog_info$Censor_date is Date as well
nonprog_info <- nonprog_info %>%
  mutate(
    Censor_date = as.Date(Censor_date)
  )

# Convert nonprog_info's single Date into a list-column of length-1 Date vectors:
nonprog_info <- nonprog_info %>%
  mutate(
    # wrap each scalar Censor_date in its own list (of length 1)
    Other_Progression_Dates = map(Censor_date, ~ as.Date(NA))
  )

# Now bind:
final_tbl <- bind_rows(prog_info, nonprog_info)

## Add baseline
## Edit to be correct for EK-09 and EK-06 since current dates are for smouldering 
manual_baseline_date_corrections <- manual_clinical_correction_rows(
  c("ek09_baseline_date", "ek06_baseline_date")
) %>%
  transmute(Patient, manual_diagnosis_date = Diagnosis_date, manual_baseline_date = Baseline_Date)

baseline_dates_all <- baseline_dates_all %>%
  left_join(manual_baseline_date_corrections, by = "Patient") %>%
  mutate(
    Diagnosis_date = coalesce(manual_diagnosis_date, Diagnosis_date),
    Baseline_Date = coalesce(manual_baseline_date, Baseline_Date)
  ) %>%
  select(-manual_diagnosis_date, -manual_baseline_date)

final_tbl <-final_tbl %>%
  # 1) join baseline
  left_join(baseline_dates_all, by = "Patient") 

# Reorder
final_tbl <- final_tbl %>%
  mutate(
    # `MA-VA` is a parser artifact, not a real patient ID; map it back to the
    # corresponding VA-series participant before writing PFS support tables.
    Patient = case_when(Patient == "MA-VA" ~ "VA-25", TRUE ~ Patient)
  ) %>%
  select(Patient, Baseline_Date, Censor_date, Relapsed, Other_Progression_Dates, High_Priority)

# Collapse
tmp <- final_tbl %>%
  mutate(
    Other_Progression_Dates = map_chr(Other_Progression_Dates, function(dts) {
      # if it’s all NA or empty, leave blank
      if (all(is.na(dts)) || length(dts)==0) {
        ""
      } else {
        # convert Dates to "YYYY-MM-DD" and join with semicolons
        paste(as.character(dts), collapse = ";")
      }
    })
  )

write.csv(
  tmp,
  "Exported_data_tables_clinical/Censore_dates_per_patient_for_PFS_updated.csv",
  row.names = FALSE
)

tmp <- tmp %>% filter(!grepl("^(SPORE|IMG)", Patient))
tmp <- tmp %>% filter(!is.na(Baseline_Date))
write.csv(
  tmp,
  "Exported_data_tables_clinical/Censore_dates_per_patient_for_PFS_just_M4.csv",
  row.names = FALSE
)

saveRDS(final_tbl, "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS.rds")

normalise_patient_id <- function(patient) {
  dplyr::case_when(
    # `MA-VA` is a parser artifact, not a real patient ID; map it back to the
    # corresponding VA-series participant before writing PFS support tables.
    patient == "MA-VA" ~ "VA-25",
    TRUE ~ patient
  )
}

final_tbl_for_update <- final_tbl %>%
  rename(
    baseline_date = Baseline_Date,
    censor_date = Censor_date,
    relapsed = Relapsed
  ) %>%
  mutate(Patient = normalise_patient_id(Patient)) %>%
  arrange(Patient, baseline_date, censor_date) %>%
  distinct(Patient, .keep_all = TRUE)

max_date_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(as.Date(NA))
  max(x)
}

first_date_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(as.Date(NA))
  x[1]
}

max_integer_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_integer_)
  as.integer(max(x))
}


### Add updated dates from Sarah 
censor_dates_updated <- read.csv("Clinical data/M4/Censore_dates_per_patient_for_PFS_just_M4 SB.csv")

sarah_endpoint_status <- censor_dates_updated %>%
  select(Patient, Has.the.patient.Relapsed) %>%
  mutate(
    Patient = as.character(Patient),
    sarah_relapse_status = case_when(
      tolower(str_trim(Has.the.patient.Relapsed)) == "yes" ~ "yes",
      tolower(str_trim(Has.the.patient.Relapsed)) == "no" ~ "no",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Patient, sarah_relapse_status)

censor_tbl <- censor_dates_updated %>% select(-CMRG.ID) %>%
  mutate(across(where(is.character), str_trim)) %>%      # remove stray blanks
  mutate(
    # Preserve the historical Sarah B. update behavior: anything other than an
    # explicit "yes" is treated as not relapsed within this update table. Newer
    # March 2026 endpoint exports are applied below as a separate source layer.
    relapsed    = if_else(tolower(Has.the.patient.Relapsed) == "yes", 1L, 0L),
    relapse_dt  = dmy(If.yes..First.Relapse.Date , quiet = TRUE),
    followup_dt = dmy(If.no..Date.of.last.follow.up, quiet = TRUE),
    censor_date = coalesce(relapse_dt, followup_dt)       # relapse > last‑FU
  ) %>% 
  select(Patient, censor_date, relapsed)

## join to existing baseline dates
final_like_tbl <- censor_tbl %>% 
  left_join(final_tbl_for_update |> select(Patient, baseline_date), by = "Patient") %>% 
  relocate(baseline_date, .after = Patient)              # keep column order

# See difference 
# 3. For Patients in both, show which fields differ
diff_by_patient <- full_join(final_tbl_for_update,
                             final_like_tbl,
                             by     = "Patient",
                             suffix = c(".final", ".new")) %>%
  filter(baseline_date.final != baseline_date.new |
           censor_date.final   != censor_date.new   |
           relapsed.final      != relapsed.new) %>%
  select(Patient,
         baseline_date.final, baseline_date.new,
         censor_date.final,   censor_date.new,
         relapsed.final,      relapsed.new)

# Command-line replacement for the original interactive RStudio inspection.
# This table documents exactly which PFS/censor fields changed after applying
# the manually curated Sarah B. update file.
write.csv(
  diff_by_patient,
  "Exported_data_tables_clinical/PFS_censor_date_manual_update_differences.csv",
  row.names = FALSE
)

## Coalesce
# 1. Drop EK-09 from the diff list
diff_by_patient_flt <- diff_by_patient %>%
  filter(Patient != "EK-09")

# 2. Build the patch: keep baseline_date.final. For patients with any event
#    evidence, keep the earliest event date; for non-relapsed/censored
#    patients, keep the latest follow-up date. This prevents a later
#    progression/follow-up row from replacing the first PFS event.
diff_update <- diff_by_patient_flt %>%
  mutate(
    any_relapsed = relapsed.final == 1L | relapsed.new == 1L
  ) %>%
  rowwise() %>%
  mutate(
    censor_date = if (isTRUE(any_relapsed)) {
      first_date_or_na(sort(c(censor_date.final, censor_date.new), na.last = NA))
    } else {
      max_date_or_na(c(censor_date.final, censor_date.new))
    },
    relapsed = as.integer(any_relapsed)
  ) %>%
  ungroup() %>%
  transmute(
    Patient,
    baseline_date = baseline_date.final,
    censor_date,
    relapsed
  ) %>%
  group_by(Patient) %>%
  summarise(
    baseline_date = first_date_or_na(baseline_date),
    censor_date = max_date_or_na(censor_date),
    relapsed = max_integer_or_na(relapsed),
    .groups = "drop"
  )

final_tbl_updated <- final_tbl_for_update %>%
  rows_update(diff_update, by = "Patient")

## Apply March 2026 endpoint updates on top of the existing workflow.
#
# This is intentionally narrow:
#   * March chemotherapy progression dates are added as newer progression
#     evidence. March stem-cell-transplant progression dates are deliberately
#     not used for this endpoint update.
#   * March last-follow-up dates extend censoring only for patients with no
#     progression date, and only when the March date is later than the current
#     censor/follow-up date.
march2026_progression_by_patient <- march2026_chemo_progression %>%
  transmute(
    Patient = as.character(M4_id),
    march_progression_date = as.Date(PROGRESSION_DATE),
    march_progression_source = endpoint_source
  ) %>%
  filter(
    !is.na(Patient),
    stringr::str_detect(Patient, "^[A-Z]{2}-[0-9]{2}$"),
    !is.na(march_progression_date)
  ) %>%
  distinct(Patient, march_progression_date, .keep_all = TRUE) %>%
  left_join(final_tbl_updated %>% select(Patient, baseline_date), by = "Patient") %>%
  filter(is.na(baseline_date) | march_progression_date >= baseline_date) %>%
  group_by(Patient) %>%
  summarise(
    march_first_progression = min(march_progression_date),
    march_progression_dates = list(sort(unique(march_progression_date))),
    march_progression_sources = paste(sort(unique(march_progression_source)), collapse = "; "),
    .groups = "drop"
  )

march2026_followup_by_patient <- M4_DEMO_march2026 %>%
  transmute(
    Patient = as.character(Patient),
    march_last_followup = as.Date(Date)
  ) %>%
  filter(
    !is.na(Patient),
    stringr::str_detect(Patient, "^[A-Z]{2}-[0-9]{2}$"),
    !is.na(march_last_followup)
  ) %>%
  group_by(Patient) %>%
  summarise(march_last_followup = max_clinical_date_or_na(march_last_followup), .groups = "drop")

date_max_or_first <- function(...) {
  vals <- as.Date(c(...), origin = "1970-01-01")
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return(as.Date(NA))
  max(vals)
}

march2026_endpoint_update_audit <- final_tbl_updated %>%
  select(Patient, baseline_date, censor_date, relapsed, Other_Progression_Dates) %>%
  left_join(sarah_endpoint_status, by = "Patient") %>%
  left_join(march2026_progression_by_patient, by = "Patient") %>%
  left_join(march2026_followup_by_patient, by = "Patient") %>%
  mutate(
    apply_march_progression = !is.na(march_first_progression),
    march_progression_changes_endpoint = apply_march_progression &
      (as.integer(relapsed) != 1L |
         is.na(censor_date) |
         march_first_progression < as.Date(censor_date)),
    apply_march_followup = relapsed == 0L &
      !apply_march_progression &
      !is.na(march_last_followup) &
      (is.na(censor_date) | march_last_followup > censor_date),
    updated_relapsed = if_else(
      march_progression_changes_endpoint,
      1L,
      as.integer(relapsed)
    ),
    updated_censor_date = case_when(
      march_progression_changes_endpoint ~ march_first_progression,
      apply_march_followup ~ march_last_followup,
      TRUE ~ as.Date(censor_date)
    ),
    march_update_reason = case_when(
      march_progression_changes_endpoint ~ "march_progression_updates_event",
      apply_march_progression ~ "march_progression_confirmed_or_later_than_existing_event",
      apply_march_followup ~ "march_followup_extended_nonrelapsed_censor",
      TRUE ~ "no_march_endpoint_change"
    )
  )

readr::write_csv(
  march2026_endpoint_update_audit %>%
    filter(march_update_reason != "no_march_endpoint_change") %>%
    mutate(
      march_progression_dates = purrr::map_chr(
        march_progression_dates,
        ~ paste(as.character(.x), collapse = ";")
      )
    ) %>%
    select(
      Patient, baseline_date, censor_date, relapsed,
      sarah_relapse_status,
      march_first_progression,
      march_progression_dates,
      march_last_followup, updated_censor_date, updated_relapsed,
      march_update_reason
    ),
  file.path(clinical_support_dir, "march2026_endpoint_update_audit.csv"),
  na = ""
)

march2026_endpoint_patch <- march2026_endpoint_update_audit %>%
  filter(march_progression_changes_endpoint | apply_march_followup) %>%
  transmute(
    Patient,
    censor_date = updated_censor_date,
    relapsed = updated_relapsed
  )

if (nrow(march2026_endpoint_patch) > 0) {
  final_tbl_updated <- final_tbl_updated %>%
    rows_update(march2026_endpoint_patch, by = "Patient", unmatched = "ignore")
}

pfs_progression_events_for_other_dates <- Relapse_dates_full %>%
  transmute(
    Patient = as.character(Patient),
    Progression_date = as.Date(Progression_date)
  ) %>%
  bind_rows(
    final_tbl_updated %>%
      filter(relapsed == 1L) %>%
      transmute(
        Patient = as.character(Patient),
        Progression_date = as.Date(censor_date)
      )
  ) %>%
  filter(!is.na(Patient), !is.na(Progression_date)) %>%
  distinct()

final_tbl_updated <- final_tbl_updated %>%
  left_join(
    pfs_progression_events_for_other_dates %>%
      group_by(Patient) %>%
      summarise(
        all_pfs_progression_dates = list(sort(unique(Progression_date))),
        .groups = "drop"
      ),
    by = "Patient"
  ) %>%
  mutate(
    Other_Progression_Dates = purrr::pmap(
      list(all_pfs_progression_dates, censor_date, relapsed),
      function(events, censor, relapsed_flag) {
        if (!isTRUE(as.integer(relapsed_flag) == 1L) ||
            is.null(events) ||
            length(events) == 0 ||
            all(is.na(events)) ||
            is.na(censor)) {
          return(as.Date(NA))
        }
        later_events <- events[!is.na(events) & events > as.Date(censor)]
        if (length(later_events) == 0) {
          as.Date(NA)
        } else {
          later_events
        }
      }
    )
  ) %>%
  select(-all_pfs_progression_dates)

## Now export 
final_tbl_updated_csv <- final_tbl_updated %>%
  rename(
    Baseline_Date = baseline_date,
    Censor_date = censor_date,
    Relapsed = relapsed
  ) %>%
  mutate(
    Other_Progression_Dates = map_chr(Other_Progression_Dates, function(dts) {
      if (all(is.na(dts)) || length(dts) == 0) {
        ""
      } else {
        paste(as.character(dts), collapse = ";")
      }
    })
  )

write.csv(
  final_tbl_updated_csv,
  "Exported_data_tables_clinical/Censore_dates_per_patient_for_PFS_updated.csv",
  row.names = FALSE
)
saveRDS(final_tbl_updated, "Exported_data_tables_clinical/Censor_dates_per_patient_for_PFS_updated.rds")

## Add this to the previous relapse dates 
# Keep original for comparison
Relapse_dates_full_old <- Relapse_dates_full

# Extract relapse rows from final_tbl_updated
new_relapses <- final_tbl_updated %>%
  filter(relapsed == 1) %>%
  select(Patient, Progression_date = censor_date)

# Add finalized patient-level PFS events to the relapse-date table without
# deleting historical progression evidence. This table is an event inventory;
# the current PFS event/censor per patient is controlled by `final_tbl_updated`
# above, not by pruning this historical table.
Relapse_dates_full <- Relapse_dates_full %>%
  mutate(
    Patient = as.character(Patient),
    Progression_date = as.Date(Progression_date)
  ) %>%
  select(Patient, Progression_date) %>%
  bind_rows(new_relapses) %>%
  filter(!is.na(Patient), !is.na(Progression_date)) %>%
  distinct(Patient, Progression_date, .keep_all = TRUE)

# Check if any changes occurred
changes <- anti_join(Relapse_dates_full, Relapse_dates_full_old) %>%
  bind_rows(anti_join(Relapse_dates_full_old, Relapse_dates_full))

if (nrow(changes) > 0) {
  message("Changes detected:")
  print(changes)
} else {
  message("No changes detected.")
}

saveRDS(Relapse_dates_full, "Exported_data_tables_clinical/Relapse_dates_full_updated.rds")
write.csv(Relapse_dates_full, "Exported_data_tables_clinical/Relapse dates cfWGS updated2.csv",   row.names = FALSE)

# Write the final combined clinical CSV after all endpoint sources have been
# integrated. Sample-level relapse fields are derived from the finalized
# progression-date inventory so downstream tables start from the same endpoint
# finalized progression dates as the patient-level PFS table.
recompute_sample_relapse_fields <- function(clinical_df, relapse_dates_tbl) {
  relapse_dates_tbl <- relapse_dates_tbl %>%
    transmute(
      Patient = as.character(Patient),
      Progression_date = as.Date(Progression_date)
    ) %>%
    filter(!is.na(Patient), !is.na(Progression_date)) %>%
    distinct()

  clinical_base <- clinical_df %>%
    mutate(
      .clinical_row_id = row_number(),
      Patient = as.character(Patient),
      Date_of_sample_collection = as.Date(Date_of_sample_collection)
    ) %>%
    select(-any_of(c(
      "Num_days_to_closest_relapse",
      "Relapsed",
      "Num_days_to_closest_relapse_absolute",
      "Num_days_to_closest_relapse_non_absolute"
    )))

  sample_relapse <- clinical_base %>%
    select(.clinical_row_id, Patient, Date_of_sample_collection) %>%
    left_join(relapse_dates_tbl, by = "Patient", relationship = "many-to-many") %>%
    mutate(
      days_to_relapse_non_absolute = case_when(
        is.na(Date_of_sample_collection) | is.na(Progression_date) ~ NA_real_,
        Progression_date > Date_of_sample_collection ~
          as.numeric(Date_of_sample_collection - Progression_date),
        Progression_date <= Date_of_sample_collection &
          Progression_date >= Date_of_sample_collection - 35 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(.clinical_row_id) %>%
    summarise(
      Num_days_to_closest_relapse_absolute = {
        valid_days <- days_to_relapse_non_absolute[!is.na(days_to_relapse_non_absolute)]
        if (length(valid_days) == 0) NA_real_ else min(abs(valid_days))
      },
      Num_days_to_closest_relapse_non_absolute = {
        valid_days <- days_to_relapse_non_absolute[!is.na(days_to_relapse_non_absolute)]
        if (length(valid_days) == 0) NA_real_ else valid_days[which.min(abs(valid_days))]
      },
      .groups = "drop"
    ) %>%
    mutate(
      Num_days_to_closest_relapse = Num_days_to_closest_relapse_absolute,
      Relapsed = if_else(!is.na(Num_days_to_closest_relapse_absolute), "Y", "N")
    )

  clinical_base %>%
    left_join(sample_relapse, by = ".clinical_row_id") %>%
    select(-.clinical_row_id)
}

combined_clinical_data_updated <- recompute_sample_relapse_fields(
  clinical_df = combined_clinical_data_updated,
  relapse_dates_tbl = Relapse_dates_full
)

write.csv(
  combined_clinical_data_updated,
  "combined_clinical_data_updated_April2025.csv",
  row.names = FALSE
)

# ─── 23.  Get the patient info, number of patients with each feature ───────────────────────

# Define baseline and monitoring groups
baseline_timepoints <- c("Diagnosis", "Baseline")
monitoring_timepoints <- c("Post_transplant", "Post_induction", "Maintenance", "1.5yr maintenance", "2yr maintenance")

# Filter for baseline cfDNA samples
baseline_cfDNA <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% baseline_timepoints) %>%
  dplyr::select(Patient) %>%
  distinct()

# Count unique patients with baseline cfDNA samples
num_baseline_patients <- n_distinct(baseline_cfDNA$Patient)

# Filter for monitoring cfDNA samples
monitoring_cfDNA <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% monitoring_timepoints) %>%
  dplyr::select(Patient) %>%
  distinct()

# Count unique patients with monitoring cfDNA samples
num_monitoring_patients <- n_distinct(monitoring_cfDNA$Patient)

# Find patients with both baseline and monitoring samples
patients_with_both <- intersect(baseline_cfDNA$Patient, monitoring_cfDNA$Patient)
num_patients_with_both <- length(patients_with_both)

# Print results
cat("Number of patients with baseline cfDNA samples:", num_baseline_patients, "\n")
cat("Number of patients with monitoring cfDNA samples:", num_monitoring_patients, "\n")
cat("Number of patients with both baseline and monitoring cfDNA samples:", num_patients_with_both, "\n")

# Availability summary only: retain patients with both baseline and monitoring
# cfDNA regardless of baseline signal. Evidence-of-disease filtering belongs in
# downstream scientific analyses, not in this clinical inventory step.

# Number of cfDNA samples for patients with both baseline and monitoring cfDNA types
num_samples_with_both <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", Patient %in% patients_with_both) %>%
  nrow()

cat("Number of cfDNA samples for patients with both baseline and monitoring cfDNA types:", num_samples_with_both, "\n")

write.csv(patients_with_both, file = "patients_with_cfDNA_at_baseline_and_monitoring.csv")

samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", Patient %in% patients_with_both)

# Total unique patients with cfDNA samples
total_unique_patients <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  dplyr::select(Patient) %>%
  distinct() %>%
  nrow()

# Total number of cfDNA samples
total_cfDNA_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA") %>%
  nrow()

# Print results
cat("Total number of unique patients with cfDNA samples:", total_unique_patients, "\n")
cat("Total number of cfDNA samples:", total_cfDNA_samples, "\n")

# Number of baseline samples and patients
baseline_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% baseline_timepoints)

num_baseline_samples <- nrow(baseline_samples)
num_baseline_patients <- n_distinct(baseline_samples$Patient)

# Number of monitoring samples and patients
monitoring_samples <- combined_clinical_data_updated %>%
  filter(Sample_type == "Blood_plasma_cfDNA", timepoint_info %in% monitoring_timepoints)

num_monitoring_samples <- nrow(monitoring_samples)
num_monitoring_patients <- n_distinct(monitoring_samples$Patient)

# Print results
cat("Number of baseline cfDNA samples:", num_baseline_samples, "\n")
cat("Number of unique patients with baseline cfDNA samples:", num_baseline_patients, "\n")
cat("Number of monitoring cfDNA samples:", num_monitoring_samples, "\n")
cat("Number of unique patients with monitoring cfDNA samples:", num_monitoring_patients, "\n")

# Combine patients from baseline and monitoring groups
patients_baseline <- baseline_samples$Patient
patients_monitoring <- monitoring_samples$Patient

# Union of baseline and monitoring patients
patients_with_either <- union(patients_baseline, patients_monitoring)

# Number of unique patients with either baseline or monitoring samples
num_patients_with_either <- length(patients_with_either)

# Print the result
cat("Number of unique patients with either baseline or monitoring cfDNA samples:", num_patients_with_either, "\n")

# Extract site identifiers from patients_with_either
site_identifiers <- sapply(patients_with_either, function(patient) {
  strsplit(patient, "[-_]")[[1]][1]  # Split by hyphen or underscore and take the first part
})

# Exclude IMG from the count
filtered_sites <- unique(site_identifiers[!grepl("^IMG", site_identifiers)])

# Number of unique sites
num_sites <- length(filtered_sites)

# Print the result
cat("Number of unique sites (excluding IMG):", num_sites, "\n")

write.csv(patients_with_either, file = "patients_with_either_cfDNA_at_baseline_and_monitoring.csv", row.names = FALSE)



#### Check against CMRG baseline dates

Baseline_dates_CMRG <- read_excel("M4_CMRG_Data/M4_COHORT_DIAGNOSTIC_BIOPSY.xlsx") %>%
  mutate(
    # if the biopsy-date column contains Excel serial dates:
    Baseline_Date = as.Date(as.numeric(BIOPSY_DATE), origin = "1899-12-30"), 
    Patient = M4_id
  )

## Get diagnosis only 
# Keep only rows where PURPOSE indicates a diagnostic baseline
diagnostic_baselines <- Baseline_dates_CMRG %>%
  filter(
    str_detect(PURPOSE, regex("^At diagnosis",      ignore_case = TRUE)) |
      str_detect(PURPOSE, regex("^Diagnosis",         ignore_case = TRUE)) |
      str_detect(PURPOSE, regex("^Diagnostic\\b",     ignore_case = TRUE))
  )

# Inspect
table(diagnostic_baselines$PURPOSE)
# 2) Get the official baseline_date from the curated PFS table written above.
baseline_official <- final_tbl_updated %>%
  select(Patient, baseline_date)

# 3) Compute the difference
baseline_comparison <- baseline_official %>%
  left_join(diagnostic_baselines %>% select(Patient, Baseline_Date, PURPOSE), by = "Patient") %>%
  mutate(
    diff_days   = as.numeric(baseline_date - Baseline_Date),
    diff_months = diff_days / 30.44
  )

baseline_comparison <- baseline_comparison %>%
  left_join(cohort_df, by = "Patient")

baseline_comparison <- baseline_comparison %>%
  rename(
    Baseline_Date_CMRG_BM_Biopsy     = Baseline_Date,
    Baseline_Date_current_info       = baseline_date
  )

### Now see which blood are not at baseline? 
# 1) Filter to timepoint “01” (baseline draws)  
processing_baseline <- processing_log_cleaned %>%
  filter(Timepoint == "01") %>%
  rename(sample_date = Cleaned_Date) %>%
  select(Patient, sample_date)

## See mismatch
baseline_mismatch <- processing_baseline %>%
  left_join(
    baseline_comparison %>%
      select(Patient,
             Baseline_Date_CMRG_BM_Biopsy,
             Baseline_Date_current_info,
             Cohort),
    by = "Patient"
  ) %>%
  mutate(
    diff_vs_CMRG     = as.numeric(sample_date - Baseline_Date_CMRG_BM_Biopsy),
    diff_vs_current  = as.numeric(sample_date - Baseline_Date_current_info),
    mismatch_CMRG    = diff_vs_CMRG    != 0,
    mismatch_current = diff_vs_current != 0
  ) %>%
  filter(mismatch_CMRG | mismatch_current) %>% 
  distinct()

# Inspect
print(baseline_mismatch)

# Export baseline-date mismatch QC tables for command-line review.
outdir <- "Output_tables_2025/detection_progression"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  baseline_mismatch,
  file.path(outdir, "baseline_date_mismatches_both.csv")
)
# 4) Write out to CSV for reporting  
write_csv(
  baseline_mismatch,
  file.path(outdir, "baseline_date_mismatches.csv")
)


# This block is audit-only: baseline mismatches are exported for manual review.
# It does not automatically replace an official baseline date with a second draw.
