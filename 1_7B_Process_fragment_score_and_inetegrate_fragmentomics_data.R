# ──────────────────────────────────────────────────────────────────────────────
# 1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R
#
# Purpose:
#   • Load and process “insert_size_summary.tsv” and “fragment_scores.tsv” for cfWGS samples 
#     (MM baseline/monitoring) and PON (healthy controls).
#   • Harmonize sample names, merge with “combined_clinical_data_updated_April2025.csv” 
#   • Define comparison groups:
#        – Diagnosis/Baseline vs Monitoring vs Healthy (insert size)
#   • For insert‐size (“Proportion.Short”):
#        – Test normality (Shapiro–Wilk), run ANOVA or Kruskal–Wallis
#        – Run all pairwise t‐tests (or Wilcoxon) and save results to CSV.
#   • For fragment‐score (“FS”):
#        – Merge Fragment Score to clinical data
#        – Deduplicate by sample and sequencing platform
#   • Load MM‐DARs chromatin activation metrics (MM_DARs_chromatin_activation_data.csv),
#     join into the fragment‐score data frame, and export the unified
#     `Key_fragmentomics_data_updated2` table.
#
# Inputs:
#   • Historical and revision insert_size_summary.tsv / fragment_scores.tsv
#     tables resolved under Fragmentomics_data and the revision data root
#   • Platform-specific healthy-control summaries
#   • combined_clinical_data_updated_April2025.csv
#   • Results_Fragmentomics/MM_DARs_chromatin_activation_data.rds or .csv
#   • Required legacy recovery file (read unconditionally by the current code):
#     Results_Fragmentomics/MM_Griffin_all_relevant_sites_data_updated_SPORE.csv
#   • patients_with_either_cfDNA_at_baseline_and_monitoring.csv and
#     patients_with_cfDNA_at_baseline_and_monitoring.csv are read for historical
#     compatibility but are not referenced by the active calculations.
#   • Scripts_2025/Final_Scripts/helpers.R
#
# Outputs:
#   • Pairwise_t_Test_Results_fragment_scores.csv
#   • Results_Fragmentomics/Insert_size/FS_cutoffs_table_healthy_controls.csv
#   • Results_Fragmentomics/Insert_size/HC_Ranges_Selected_MM_DARs_Metrics.csv
#   • Results_Fragmentomics/fragmentomics_platform_harmonization_reference_parameters.csv
#   • Results_Fragmentomics/fragmentomics_platform_harmonization_sample_audit.csv
#   • Results_Fragmentomics/fragmentomics_duplicate_sample_key_audit.csv
#   • Results_Fragmentomics/fragmentomics_MyP060_identity_audit.csv
#   • Results_Fragmentomics/Key_fragmentomics_data_updated2.csv
#   • Results_Fragmentomics/Key_fragmentomics_data_updated2.rds
#
# Analysis units and limitations:
#   The intended integrated-output key is fragmentomics sample + platform +
#   sample role. The final metadata reattachment can retain more than one row for
#   that key when the FS/clinical and MM-DAR source dates disagree; downstream
#   consumers must audit key uniqueness before treating rows as samples.
#   Group-comparison tests operate on sample rows, not independent patients;
#   repeated longitudinal samples can therefore contribute more than once.
#   The adaptive Shapiro/test-selection results and unadjusted pairwise p-values
#   are exploratory.
#
# Notes:
#   – Plotting will happen later in a separate script.
# 
# How to run:
#   Rscript Scripts_2025/Final_Scripts/1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes fragment-size scores and
#     integrates fragmentomics features used by model training, longitudinal
#     summaries, and dilution-series comparisons.
# ──────────────────────────────────────────────────────────────────────────────
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#
# Dependencies:
#   BoutrosLab.plotting.general, ggplot2, conflicted, survival, survminer,
#   dplyr, readr, stringr, rstatix, and tidyr.
#

### PREPARE SESSION ################################################################################
library(BoutrosLab.plotting.general)
library(ggplot2)        # used only for tibble/data‐frame handling (no plots)
library(conflicted)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)
library(rstatix)        # for pairwise_wilcox_test, if needed
library(tidyr)

conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::mutate,
  dplyr::summarise,
  dplyr::arrange,
  dplyr::rename,
  dplyr::lag,
  .quiet = TRUE
)

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  x[[1]]
}

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)


### SET PATHS #######################################################################################
input.dir <- "Fragmentomics_data"
pon.dir   <- file.path("Fragmentomics_data", "Normals")
out.dir   <- file.path("Results_Fragmentomics", "Insert_size")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)


### STEP 1: PROCESS INSERT‐SIZE (“Proportion.Short”) ###############################################
# "insert_size_summary.tsv" contains the cfDNA fragment-length distribution
# per sample. "Proportion.Short" is the short-fragment score: the fraction of
# all plasma cfDNA fragment lengths that fall below a short-fragment threshold
# (< 150 bp). In cancer patients, nucleosome-free tumour-derived
# cfDNA leads to an elevated proportion of short fragments relative to healthy
# donors. Higher Proportion.Short therefore correlates with higher tumour
# fraction and is a tumour-naive MRD feature.
historical_insert_files <- list.files(
  path = input.dir,
  pattern = "insert_size_summary.tsv",
  full.names = TRUE
)
revision_insert_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-25_cfWGS_MM_fragmentomics_Revisions_Spring2026_insert_size_summary[.]tsv$"
)
pon.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "insert_size_summary.tsv",
             full.names = TRUE)
))
xplus_pon_insert_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_insert_size_summary[.]tsv$"
)

read_fragmentomics_summary_batch <- function(paths, platform, role) {
  if (!length(paths)) return(tibble())
  bind_rows(lapply(paths, function(path) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(
        fragmentomics_sequencing_platform = platform,
        fragmentomics_sample_role = role,
        fragmentomics_source_file = basename(path)
      )
  }))
}

filter_revision_exclusions <- function(data) {
  data %>%
    filter(
      !(.data$fragmentomics_sample_role == "patient_revision" &
          is_spring2026_revision_primary_analysis_excluded(.data$Sample))
    )
}

cfWGS_insert <- bind_rows(
  read_fragmentomics_summary_batch(historical_insert_files, "NovaSeq 6000", "patient"),
  read_fragmentomics_summary_batch(revision_insert_files, "NovaSeq XPlus", "patient_revision")
) %>%
  filter_revision_exclusions()

if (!length(pon.files) || !length(xplus_pon_insert_files)) {
  stop("Both historical and XPlus insert-size healthy-control inputs are required.", call. = FALSE)
}
pon_insert <- bind_rows(
  read_fragmentomics_summary_batch(pon.files[[1]], "NovaSeq 6000", "healthy_control") %>%
    filter_fragmentomics_original_historical_controls(),
  read_fragmentomics_summary_batch(xplus_pon_insert_files[[1]], "NovaSeq XPlus", "healthy_control") %>%
    filter_fragmentomics_shared_healthy_controls()
)

# Load clinical tables
combined_clinical <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)

# Construct a “Merge_ID” for joining
tmp.clin <- combined_clinical %>%
  mutate(
    Merge_ID = if_else(
      str_starts(Bam, "SPORE"),
      str_extract(Bam, "^[^\\.]+"),
      str_replace(Sample_ID, "_Blood_plasma_cfDNA", "-P")
    )
  )
spring2026_pwgval_metadata <- load_spring2026_pwgval_dilution_metadata(required = FALSE)
if (!is.null(spring2026_pwgval_metadata)) {
  tmp.clin <- bind_rows(
    tmp.clin,
    spring2026_pwgval_metadata %>%
      transmute(
        Merge_ID = paste0(.data$Group_ID, "-P"),
        Patient = .data$Patient,
        Sample_ID = .data$Sample_ID,
        Bam = .data$Bam,
        Date_of_sample_collection = as.Date(NA),
        Timepoint = as.character(.data$Timepoint),
        timepoint_info = .data$timepoint_info,
        Study = "PWGVAL_M4CHIP"
      )
  )
}
tmp.clin <- bind_rows(
  tmp.clin,
  tmp.clin %>% mutate(Merge_ID = str_replace(Merge_ID, "^IMG", "MyP"))
) %>%
  distinct(Merge_ID, Patient, Sample_ID, .keep_all = TRUE)

tmp.clin <- tmp.clin %>%
  arrange(.data$Merge_ID, desc(!is.na(.data$Date_of_sample_collection))) %>%
  distinct(.data$Merge_ID, .keep_all = TRUE)

# Standardize cfWGS_insert sample names, then join to clinical keys
cfWGS_insert <- cfWGS_insert %>%
  mutate(
    Sample = str_replace(Sample, "_Blood_plasma_cfDNA", "-P"),
    Sample = if_else(
      !str_starts(Sample, "SPORE"),
      str_replace_all(Sample, "_", "-"),
      Sample
    ),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(tmp.clin, by = c("Sample" = "Merge_ID"), relationship = "many-to-one")

# Mark PON samples as “Healthy”
classify_fragmentomics_material <- function(sample) {
  case_when(
    str_detect(sample, "(^|[-_])Cf([-_]|$)|_Cf_|-Cf-") ~ "cfDNA",
    str_detect(sample, "(^|[-_])Pb([-_]|$)|_Pb_|-Pb-") ~ "PB_cellular_or_buffy",
    str_detect(sample, "(^|[-_])Ct([-_]|$)|_Ct_|-Ct-") ~ "cellular_or_tissue",
    TRUE ~ "unknown"
  )
}

pon_insert <- pon_insert %>%
  mutate(
    Group = "Healthy",
    fragmentomics_material = classify_fragmentomics_material(.data$Sample)
  )

# Define “Monitoring” vs “Diagnosis/Baseline” groups
# Clinical timepoint groupings used for group comparisons:
# group_1 = on-treatment / post-transplant monitoring timepoints
#   (short-fragment score expected to be reduced or normalised vs. Baseline).
# group_2 = pre-treatment diagnostic/baseline timepoints
#   (typically highest tumour fraction, therefore highest Proportion.Short).
group_1 <- c("Post_transplant", "Post_induction", "Maintenance",
             "1.5yr maintenance", "2yr maintenance")
group_2 <- c("Diagnosis", "Baseline")

data_insert <- cfWGS_insert %>%
  mutate(Group = ifelse(timepoint_info %in% group_1,
                        "Monitoring", "Diagnosis/Baseline"))

# Read patient lists
patients_with_either <- read_csv("patients_with_either_cfDNA_at_baseline_and_monitoring.csv")
patients_with_both   <- read_csv("patients_with_cfDNA_at_baseline_and_monitoring.csv")

# Combine cfWGS + PON (dropping one duplicate PON sample if present)
pon_insert <- pon_insert %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

# Prepare pon_insert with the same columns (as NA)
pon_insert2 <- pon_insert %>%
  mutate(
    Patient        = NA_character_,   # same type as data_insert$Patient
    timepoint_info = NA_character_    # same type as data_insert$timepoint_info
  )

# Bind together
all_insert <- bind_rows(
  data_insert %>%
    select(Sample, Patient, timepoint_info, Group, Proportion.Short,
           fragmentomics_sequencing_platform, fragmentomics_sample_role),
  
  pon_insert2 %>%
    select(Sample, Patient, timepoint_info, Group, Proportion.Short,
           fragmentomics_sequencing_platform, fragmentomics_sample_role)
)

# Preserve the original NovaSeq 6000-only group-comparison analysis. XPlus
# patient and control rows remain in `all_insert` for feature integration and
# harmonization, but must not alter the historical descriptive statistics.
all_insert_stats <- all_insert %>%
  filter(.data$fragmentomics_sequencing_platform == "NovaSeq 6000")


# TEST NORMALITY (Shapiro–Wilk per group)
# If all groups are normally distributed (p > 0.05), use parametric ANOVA
# and pairwise t-tests. If any group fails normality, fall back to
# non-parametric Kruskal-Wallis global test and Wilcoxon pairwise tests.
# This adaptive branching preserves valid inference regardless of distribution.
shapiro_results <- all_insert_stats %>%
  group_by(Group) %>%
  summarise(p_value = shapiro.test(Proportion.Short)$p.value)

print(shapiro_results)

# DECIDE whether to use ANOVA or Kruskal–Wallis based on shapiro_results
if (all(shapiro_results$p_value > 0.05)) {
  # parametric ANOVA
  test_insert_global <- aov(Proportion.Short ~ Group, data = all_insert_stats)
  print(summary(test_insert_global))
} else {
  # nonparametric Kruskal–Wallis
  test_insert_global <- kruskal.test(Proportion.Short ~ Group, data = all_insert_stats)
  print(test_insert_global)
}

# PAIRWISE comparisons – use t.test if normal, else Wilcoxon
# Pairwise comparisons across all three groups (Healthy / Diagnosis+Baseline /
# Monitoring). Results capture mean Proportion.Short per group, the difference,
# and the p-value from either a t-test (normal) or Wilcoxon (non-normal).
# Saved to CSV for downstream reporting in manuscript Table / Figure.
pairwise_t_results <- list()
groups <- unique(all_insert_stats$Group)

for (i in seq_len(length(groups) - 1)) {
  for (j in seq((i + 1), length(groups))) {
    g1 <- groups[i]
    g2 <- groups[j]
    vec1 <- all_insert_stats %>%
      filter(Group == g1) %>%
      pull(Proportion.Short)
    vec2 <- all_insert_stats %>%
      filter(Group == g2) %>%
      pull(Proportion.Short)
    
    if (all(shapiro_results$p_value > 0.05)) {
      ttest <- t.test(vec1, vec2, na.rm = TRUE)
      pval <- ttest$p.value
      method_used <- "t.test"
    } else {
      wt <- wilcox.test(vec1, vec2, na.rm = TRUE)
      pval <- wt$p.value
      method_used <- "wilcox.test"
    }
    
    pairwise_t_results[[paste(g1, "vs", g2)]] <- data.frame(
      Group1 = g1,
      Group2 = g2,
      Mean_Group1 = mean(vec1, na.rm = TRUE),
      Mean_Group2 = mean(vec2, na.rm = TRUE),
      Mean_Difference = mean(vec1, na.rm = TRUE) - mean(vec2, na.rm = TRUE),
      P_Value = pval,
      CI_Lower = if (method_used == "t.test") ttest$conf.int[1] else NA,
      CI_Upper = if (method_used == "t.test") ttest$conf.int[2] else NA,
      Method = method_used,
      stringsAsFactors = FALSE
    )
  }
}

pairwise_t_results_df <- do.call(rbind, pairwise_t_results)

write_csv(
  pairwise_t_results_df,
  file = "Pairwise_t_Test_Results_fragment_scores.csv"
)


### STEP 2: PROCESS FRAGMENT SCORE (“FS”) ############################################################
# "fragment_scores.tsv" contains the CHARM Fragment Score (FS) per sample.
# FS is a genome-wide summary of short-fragment enrichment derived by
# comparing cfDNA coverage across CHARM-defined bins to a healthy-donor
# reference panel. It integrates fragment-length information over many loci
# simultaneously, making it more robust than a simple short-fragment ratio.
# Higher FS → more tumour-like fragment landscape.
historical_fs_files <- list.files(
  path = input.dir,
  pattern = "fragment_scores.tsv",
  full.names = TRUE
)
revision_fs_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-25_cfWGS_MM_fragmentomics_Revisions_Spring2026_fragment_scores[.]tsv$"
)
pon.fs.files <- rev(sort(
  list.files(path = pon.dir,
             pattern = "fragment_scores.tsv",
             full.names = TRUE)
))
xplus_pon_fs_files <- spring2026_revision_files(
  "Fragmentomics_Pipeeline_Suite_all_outputs",
  "^2026-06-30_cfWGS_MM_fragmentomics_Revisions_Spring2026_CHARM_Xplus_HC_fragment_scores[.]tsv$"
)

cfWGS_fs <- bind_rows(
  read_fragmentomics_summary_batch(historical_fs_files, "NovaSeq 6000", "patient"),
  read_fragmentomics_summary_batch(revision_fs_files, "NovaSeq XPlus", "patient_revision")
) %>%
  filter_revision_exclusions()

# Drop any “TGL49” rows if they snuck in
if (any(grepl("TGL49", cfWGS_fs$Sample))) {
  cfWGS_fs <- cfWGS_fs[!grepl("TGL49", cfWGS_fs$Sample), ]
}

if (!length(pon.fs.files) || !length(xplus_pon_fs_files)) {
  stop("Both historical and XPlus fragment-score healthy-control inputs are required.", call. = FALSE)
}
pon_fs <- bind_rows(
  read_fragmentomics_summary_batch(pon.fs.files[[1]], "NovaSeq 6000", "healthy_control") %>%
    filter_fragmentomics_original_historical_controls(),
  read_fragmentomics_summary_batch(xplus_pon_fs_files[[1]], "NovaSeq XPlus", "healthy_control") %>%
    filter_fragmentomics_shared_healthy_controls()
) %>%
  mutate(fragmentomics_material = classify_fragmentomics_material(.data$Sample))

# Standardize sample names & join clinical keys just as above
cfWGS_fs <- cfWGS_fs %>%
  mutate(
    Sample = str_replace(Sample, "_Blood_plasma_cfDNA", "-P"),
    Sample = if_else(
      !str_starts(Sample, "SPORE"),
      str_replace_all(Sample, "_", "-"),
      Sample
    ),
    Sample = str_replace(Sample, "IMG", "MyP")
  ) %>%
  left_join(
    tmp.clin %>% dplyr::select(Patient, Timepoint, timepoint_info, Merge_ID),
    by = c("Sample" = "Merge_ID"),
    relationship = "many-to-one"
  ) %>%
  mutate(Cohort = "cfWGS")

cfWGS_fs <- cfWGS_fs %>% unique()

pon_fs <- pon_fs %>%
  mutate(Cohort = "HBC") %>%
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample))

# Combine them
combined_fs <- bind_rows(cfWGS_fs, pon_fs)

# Join MRD data and drop unnecessary columns
combined_fs <- combined_fs %>%
 # dplyr::select(-Relapsed, -Num_days_to_closest_relapse, -Num_days_to_closest_relapse_absolute, -Num_days_to_closest_relapse_non_absolute) %>%
  left_join(combined_clinical %>% select(-Bam),
            by = c("Patient", "Timepoint", "timepoint_info"))

combined_fs <- combined_fs %>% dplyr::select(-dplyr::any_of("...1"))
  
# Deduplicate numeric fields (take mean) and non‐numeric (first non‐NA)
combined_fs_dedup <- combined_fs %>%
  group_by(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ !is.numeric(.x)), ~ first(.x)),
    .groups = "drop"
  )




# STEP 2b: DEFINE HEALTHY CONTROL CUTOFF FOR FRAGMENT SCORE ("FS")
# The upper cutoff anchors "normal" FS variation from a CHARM healthy-donor
# panel (TGL49 samples). Samples above this cutoff are flagged as having
# an elevated tumour-like fragment signature. If the healthy FS distribution
# is Gaussian, use mean ± 1.96·SD (a central 95% reference interval, whose
# upper boundary is approximately the 97.5th percentile); otherwise use the
# empirical 2.5th and 97.5th percentiles.

safe_platform_cutoffs <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 3L) {
    return(tibble(Method = "insufficient controls", Lower_Cutoff = NA_real_, Upper_Cutoff = NA_real_))
  }
  shp_p <- shapiro.test(values)$p.value
  if (shp_p > 0.05) {
    tibble(
      Method = "mean +/- 1.96 SD",
      Lower_Cutoff = mean(values) - 1.96 * sd(values),
      Upper_Cutoff = mean(values) + 1.96 * sd(values)
    )
  } else {
    qs <- quantile(values, probs = c(0.025, 0.975), na.rm = TRUE)
    tibble(
      Method = "2.5th / 97.5th percentiles",
      Lower_Cutoff = qs[[1]],
      Upper_Cutoff = qs[[2]]
    )
  }
}

fs_cutoffs_tbl <- pon_fs %>%
  filter(.data$Cohort == "HBC") %>%
  group_by(.data$fragmentomics_sequencing_platform) %>%
  group_modify(~ safe_platform_cutoffs(.x$FS)) %>%
  ungroup()

# 2) (Optional) print to console
print(fs_cutoffs_tbl)

# 3) Write out to CSV in the output directory
readr::write_csv(
  fs_cutoffs_tbl,
  file.path(out.dir, "FS_cutoffs_table_healthy_controls.csv")
)





### STEP 3: MERGE MM-DARs CHROMATIN ACTIVATION METRICS ###############################################
# Loads the per-sample GRIFFIN metrics for the MM_DARs_chromatin_activation
# site, pre-computed by script 1_7A. Contains Mean.Coverage, Midpoint.Coverage,
# Midpoint.normalized, Amplitude, z-scores vs. healthy, and binary
# Threshold flags. These are joined with FS and Proportion.Short below
# to form the unified fragmentomics feature table.
mm_dars_rds_path <- "Results_Fragmentomics/MM_DARs_chromatin_activation_data.rds"
mm_dars_csv_path <- "Results_Fragmentomics/MM_DARs_chromatin_activation_data.csv"
if (
  file.exists(mm_dars_csv_path) &&
    (!file.exists(mm_dars_rds_path) ||
       file.info(mm_dars_csv_path)$mtime > file.info(mm_dars_rds_path)$mtime)
) {
  warning(
    "MM_DARs RDS is missing or older than the CSV; reading the CSV to avoid stale fragmentomics features."
  )
  mm_dars <- read_csv(mm_dars_csv_path, show_col_types = FALSE)
} else {
  mm_dars <- readRDS(mm_dars_rds_path)
}

## In case run in two batches (optional)
mm_dars2 <- read.csv("Results_Fragmentomics/MM_Griffin_all_relevant_sites_data_updated_SPORE.csv")

mm_dars2_small <- mm_dars2 %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection),
    Date = as.Date(Date),
    fragmentomics_sequencing_platform = "NovaSeq 6000",
    fragmentomics_sample_role = if_else(str_detect(.data$Sample, "^TGL"), "healthy_control", "patient"),
    n_healthy_reference = NA_integer_
  ) %>%
  select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized,
    Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    fragmentomics_sequencing_platform, fragmentomics_sample_role,
    n_healthy_reference,
    Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct() %>%
  filter(
    Site == "MM_DARs_chromatin_activation",
    !str_detect(.data$Sample, "^TGL49[-_]")
  )

# Keep only the “MM_DARs_chromatin_activation” site rows
mm_dars_small <- mm_dars %>%
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection, format = "%Y-%m-%d"),
    Date = if ("Date" %in% names(.)) as.Date(.data$Date, format = "%Y-%m-%d") else Date_of_sample_collection
  ) %>%
  select(
    Sample, Site,
    Mean.Coverage, Midpoint.Coverage, Midpoint.normalized,
    Amplitude, Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
    Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
    fragmentomics_sequencing_platform, fragmentomics_sample_role,
    n_healthy_reference,
    Bam, Patient, Date_of_sample_collection
  ) %>%
  distinct() %>%
  filter(Site == "MM_DARs_chromatin_activation")

# Combine only genuinely missing legacy patient rows. Healthy controls must
# come exclusively from the newly regenerated, platform-filtered PON table;
# the optional legacy batch contains unmatched controls as well as differently
# punctuated duplicates of the current matched-control rows.
mm_dars2_missing <- mm_dars2_small %>%
  anti_join(mm_dars_small %>% distinct(.data$Sample), by = "Sample")
mm_dars_combined <- bind_rows(mm_dars_small, mm_dars2_missing) %>%
  # normalize the date columns
  mutate(
    Date_of_sample_collection = as.Date(Date_of_sample_collection)
  ) %>%
  # pick first occurrence per key
  group_by(Sample, Site, Bam, Date_of_sample_collection) %>%
  slice_head(n = 1) %>%
  ungroup()

# For consistency
mm_dars_small <- mm_dars_combined %>% unique()

# Merge into “combined_fs_dedup” by keys (Sample, Patient, Date_of_sample_collection)
# (First, restore combined_fs_dedup into a data frame we can join onto)
fragdata_for_merge <- combined_fs_dedup

fragdata_for_merge <- fragdata_for_merge %>%
  dplyr::select(
    Sample, Patient, Date_of_sample_collection, FS,
    fragmentomics_sequencing_platform, fragmentomics_sample_role
  )


merged_data <- fragdata_for_merge %>%
  full_join(
    mm_dars_small %>%
      dplyr::select(
        Sample, Patient, Date_of_sample_collection, Site,
        Mean.Coverage, Midpoint.Coverage, Midpoint.normalized, Amplitude,
        Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude,
        Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
        fragmentomics_sequencing_platform, fragmentomics_sample_role,
        n_healthy_reference
      ),
    by = c(
      "Sample", "Patient", "Date_of_sample_collection",
      "fragmentomics_sequencing_platform", "fragmentomics_sample_role"
    )
  )


### If minor difference in dates due to different versions of log being used 
# Date-windowing: cfDNA samples collected within 7 days of each other for the
# same patient are treated as a single clinical timepoint. This handles cases
# where the FS pipeline and the GRIFFIN pipeline slightly disagree on the
# sample date (e.g. due to different pipeline log versions). Metrics are
# averaged within each window; the earliest date in the window is retained.
merged_data_grouped <- merged_data %>%
  arrange(Patient, Date_of_sample_collection) %>%
  group_by(Patient, fragmentomics_sequencing_platform, fragmentomics_sample_role) %>%
  mutate(
    date_diff = as.numeric(Date_of_sample_collection - lag(Date_of_sample_collection)),
    new_group = is.na(date_diff) | date_diff > 7,
    group_id = cumsum(new_group)
  ) %>%
  group_by(Patient, fragmentomics_sequencing_platform, fragmentomics_sample_role, group_id) %>%
  summarise(
    Sample = first(Sample),  # Pick one sample name (can be adapted to paste0 if needed)
    Date_of_sample_collection = min(Date_of_sample_collection),
    Site = first_nonmissing(Site),
    Threshold.Coverage = first_nonmissing(Threshold.Coverage),
    Threshold.Midpoint = first_nonmissing(Threshold.Midpoint),
    Threshold.Amplitude = first_nonmissing(Threshold.Amplitude),
    n_healthy_reference = first_nonmissing(n_healthy_reference),
    Mean.Coverage = mean(Mean.Coverage, na.rm = TRUE),
    Midpoint.Coverage = mean(Midpoint.Coverage, na.rm = TRUE),
    Midpoint.normalized = mean(Midpoint.normalized, na.rm = TRUE),
    Amplitude = mean(Amplitude, na.rm = TRUE),
    Zscore.Coverage = mean(Zscore.Coverage, na.rm = TRUE),
    Zscore.Midpoint = mean(Zscore.Midpoint, na.rm = TRUE),
    Zscore.Amplitude = mean(Zscore.Amplitude, na.rm = TRUE),
    FS = mean(FS, na.rm = T),
    .groups = "drop"
  )

merged_data_grouped <- merged_data_grouped %>% select(-group_id) %>% unique()

# Identify any missing (should be only healthy controls)
missing_samples <- setdiff(
  fragdata_for_merge$Sample,
  mm_dars_small$Sample
)
cat("Samples missing from MM_DARs data (likely healthy controls):\n")
print(missing_samples)


## Add missing data for VA-07 - changed when 05 timepoint became R in script 1_0 due to relapse at this timepoint update
va07_recode <- manual_clinical_correction_rows("va07_relapse_sample_recode")
merged_data_grouped <- merged_data_grouped %>%
  mutate(
    Patient = case_when(
      Sample == "VA-07-05-P" ~ "VA-07",
      TRUE ~ Patient
    ),
    Date_of_sample_collection = case_when(
      Sample == "VA-07-05-P" ~ va07_recode$Date_of_sample_collection[[1]],
      TRUE ~ Date_of_sample_collection
    )
  )

## Add insert size data 
# 1) Deduplicate all_insert so there's only one Proportion.Short per Sample
all_insert_unique <- all_insert %>%
  distinct(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role, Proportion.Short)

# 2) Left‐join onto merged_data_grouped by Sample
merged_with_short <- merged_data_grouped %>%
  left_join(
    all_insert_unique,
    by = c("Sample", "fragmentomics_sequencing_platform", "fragmentomics_sample_role")
  )

# 3) Quick check: how many got filled?
merged_with_short %>%
  summarize(
    total       = dplyr::n(),
    has_short   = sum(!is.na(Proportion.Short)),
    missing     = sum(is.na(Proportion.Short))
  )

# Correct naming
merged_with_short <- merged_with_short %>%
  mutate(
    Sample = if_else(
      str_detect(Sample, "^TGL49"),
      str_replace_all(Sample, "-", "_"),
      Sample
    )
  )
## Consolidate 
# 1) The list of metric columns want to merge
metrics <- c(
  "Mean.Coverage",       "Midpoint.Coverage",    "Midpoint.normalized",
  "Amplitude",           "Zscore.Coverage",      "Zscore.Midpoint",
  "Zscore.Amplitude",    "FS",                   "Proportion.Short"
)

# 2) A helper that turns NaN → NA then picks the first non-NA
# Required because the date-windowing summarise() step can produce NaN
# (mean of all-NA vectors) for samples where only one of the two pipeline
# outputs (GRIFFIN or CHARM) was available. first_non_na_numeric returns
# the first real value, preferring GRIFFIN-derived metrics where present.
first_non_na_numeric <- function(x) {
  # replace NaN with NA
  x <- replace(x, is.nan(x), NA_real_)
  # drop NA
  y <- x[!is.na(x)]
  if (length(y)) y[[1]] else NA_real_
}

# 3) Summarise per Sample
merged_final_metrics <- merged_with_short %>%
  group_by(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role) %>%
  summarise(
    across(
      all_of(metrics),
      first_non_na_numeric
    ),
    .groups = "drop"
  )

# 4) Add columns back 
merged_final_metrics <- merged_final_metrics %>%
  left_join(
    merged_with_short %>%
      select(
        Patient, Sample, Date_of_sample_collection, Site,
        Threshold.Coverage, Threshold.Midpoint, Threshold.Amplitude,
        fragmentomics_sequencing_platform, fragmentomics_sample_role,
        n_healthy_reference
      ) %>%
      distinct(),
    by = c("Sample", "fragmentomics_sequencing_platform", "fragmentomics_sample_role")
  )

# The metadata reattachment above can create two rows for one sample/platform
# when the FS/clinical date and MM-DAR date disagree. Retain both historical
# rows because choosing either date could change which Patient/Date row receives
# the continuous fragmentomics values in 2_0. Export the conflict for review;
# this audit does not alter the active fragmentomics table.
fragmentomics_duplicate_key_audit <- merged_final_metrics %>%
  group_by(Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role) %>%
  filter(dplyr::n() > 1L) %>%
  mutate(
    duplicate_key_n = dplyr::n(),
    site_bearing_candidate_n = sum(!is.na(Site))
  ) %>%
  ungroup()

readr::write_csv(
  fragmentomics_duplicate_key_audit,
  "Results_Fragmentomics/fragmentomics_duplicate_sample_key_audit.csv"
)

# Remove healthy duplicates
# columns to blank for TGL* samples
threshold_cols <- c(
  "Threshold.Coverage",
  "Threshold.Midpoint",
  "Threshold.Amplitude",
  "Site"
)

merged_final_metrics_clean <- merged_final_metrics %>%
  mutate(
    across(
      all_of(threshold_cols),
      ~ replace(.x, str_detect(Sample, "^TGL"), NA)    # NA here is logical NA
    )
  ) %>%
  distinct()   # drop exact duplicates

# Harmonize XPlus model inputs onto the frozen NovaSeq 6000 training scale.
# Platform-matched z-scores remain in the Zscore.* columns; this location/scale
# mapping is only for the raw predictors consumed by the preserved classifier.
build_reference_parameters <- function(data, feature) {
  data %>%
    filter(.data$fragmentomics_sample_role == "healthy_control") %>%
    distinct(.data$Sample, .data$fragmentomics_sequencing_platform, .data[[feature]]) %>%
    group_by(.data$fragmentomics_sequencing_platform) %>%
    summarise(
      feature = feature,
      n = sum(is.finite(.data[[feature]])),
      mean = mean(.data[[feature]], na.rm = TRUE),
      sd = sd(.data[[feature]], na.rm = TRUE),
      median = median(.data[[feature]], na.rm = TRUE),
      mad = mad(.data[[feature]], na.rm = TRUE),
      .groups = "drop"
    )
}

mean_coverage_reference <- mm_dars_small %>%
  filter(.data$fragmentomics_sample_role == "healthy_control") %>%
  transmute(
    Sample, fragmentomics_sequencing_platform, fragmentomics_sample_role,
    Mean.Coverage
  )
matched_pon_fs <- filter_fragmentomics_shared_healthy_controls(pon_fs)
matched_pon_insert <- filter_fragmentomics_shared_healthy_controls(pon_insert)
matched_mean_coverage_reference <- filter_fragmentomics_shared_healthy_controls(
  mean_coverage_reference
)
harmonization_reference <- bind_rows(
  build_reference_parameters(matched_pon_fs, "FS"),
  build_reference_parameters(matched_pon_insert, "Proportion.Short"),
  build_reference_parameters(matched_mean_coverage_reference, "Mean.Coverage")
)

if (any(harmonization_reference$n < 3L) ||
    any(!is.finite(harmonization_reference$sd)) ||
    any(harmonization_reference$sd <= 0) ||
    any(!is.finite(harmonization_reference$mad)) ||
    any(harmonization_reference$mad <= 0)) {
  stop("Invalid fragmentomics platform harmonization reference parameters.", call. = FALSE)
}

harmonize_to_historical_scale <- function(data, feature, parameters) {
  hist <- parameters %>%
    filter(.data$feature == .env$feature, .data$fragmentomics_sequencing_platform == "NovaSeq 6000")
  xplus <- parameters %>%
    filter(.data$feature == .env$feature, .data$fragmentomics_sequencing_platform == "NovaSeq XPlus")
  if (nrow(hist) != 1L || nrow(xplus) != 1L) {
    stop("Missing one-to-one harmonization parameters for feature: ", feature, call. = FALSE)
  }
  raw_name <- paste0(feature, "_unharmonized")
  data[[raw_name]] <- data[[feature]]
  idx <- data$fragmentomics_sequencing_platform == "NovaSeq XPlus" & is.finite(data[[feature]])
  data[[feature]][idx] <-
    ((data[[feature]][idx] - xplus$mean[[1]]) / xplus$sd[[1]]) * hist$sd[[1]] + hist$mean[[1]]
  data
}

for (feature in c("FS", "Proportion.Short", "Mean.Coverage")) {
  merged_final_metrics_clean <- harmonize_to_historical_scale(
    merged_final_metrics_clean,
    feature,
    harmonization_reference
  )
}

add_robust_harmonization <- function(data, feature, parameters) {
  historical <- parameters %>%
    filter(.data$feature == .env$feature, .data$fragmentomics_sequencing_platform == "NovaSeq 6000")
  xplus <- parameters %>%
    filter(.data$feature == .env$feature, .data$fragmentomics_sequencing_platform == "NovaSeq XPlus")
  if (nrow(historical) != 1L || nrow(xplus) != 1L ||
      !is.finite(historical$mad[[1]]) || !is.finite(xplus$mad[[1]]) || xplus$mad[[1]] <= 0) {
    stop("Invalid robust harmonization parameters for feature: ", feature, call. = FALSE)
  }
  raw_name <- paste0(feature, "_unharmonized")
  robust_name <- paste0(feature, "_harmonized_median_mad")
  data[[robust_name]] <- data[[raw_name]]
  idx <- data$fragmentomics_sequencing_platform == "NovaSeq XPlus" & is.finite(data[[raw_name]])
  data[[robust_name]][idx] <-
    ((data[[raw_name]][idx] - xplus$median[[1]]) / xplus$mad[[1]]) *
    historical$mad[[1]] + historical$median[[1]]
  data
}

for (feature in c("FS", "Proportion.Short", "Mean.Coverage")) {
  merged_final_metrics_clean <- add_robust_harmonization(
    merged_final_metrics_clean,
    feature,
    harmonization_reference
  )
}
merged_final_metrics_clean <- merged_final_metrics_clean %>%
  mutate(
    fragmentomics_harmonization_method = if_else(
      .data$fragmentomics_sequencing_platform == "NovaSeq XPlus",
      "XPlus healthy-control mean/SD mapped to NovaSeq 6000 healthy-control scale",
      "No transform; native NovaSeq 6000 scale"
    )
  )

write_csv(
  harmonization_reference,
  "Results_Fragmentomics/fragmentomics_platform_harmonization_reference_parameters.csv"
)
write_csv(
  merged_final_metrics_clean %>%
    select(
      Sample, Patient, Date_of_sample_collection,
      fragmentomics_sequencing_platform, fragmentomics_sample_role,
      fragmentomics_harmonization_method,
      FS_unharmonized, FS,
      FS_harmonized_median_mad,
      Proportion.Short_unharmonized, Proportion.Short,
      Proportion.Short_harmonized_median_mad,
      Mean.Coverage_unharmonized, Mean.Coverage,
      Mean.Coverage_harmonized_median_mad,
      Zscore.Coverage, Zscore.Midpoint, Zscore.Amplitude
    ),
  "Results_Fragmentomics/fragmentomics_platform_harmonization_sample_audit.csv"
)

# Patient-level identity audit for the confirmed MyP-060/IMG-060 alias.
# Do not overwrite the active Patient field here: T9 and T11 share the same
# recorded date but retain distinct historical fragmentomics values. Assigning
# both to IMG-060 before the Patient/Date collapse in 2_0 would average those
# values and change manuscript inputs. The submitted IMG-060-T11 row is restored
# later by 2_0 from its frozen legacy aggregate, so it must not be reconstructed
# here from the current raw-feature rows. `Resolved_Patient` records the
# confirmed identity for review without changing the submitted data flow.
myp060_identity_audit <- merged_final_metrics_clean %>%
  filter(str_starts(Sample, "MyP-060-")) %>%
  mutate(
    Current_Patient = Patient,
    Resolved_Patient = "IMG-060",
    active_patient_field_changed = FALSE
  ) %>%
  select(
    Sample, Current_Patient, Resolved_Patient,
    active_patient_field_changed, Date_of_sample_collection,
    fragmentomics_sequencing_platform, fragmentomics_sample_role,
    FS, Proportion.Short, Mean.Coverage
  ) %>%
  arrange(Date_of_sample_collection, Sample)

if (nrow(myp060_identity_audit) == 0L ||
    any(myp060_identity_audit$Resolved_Patient != "IMG-060") ||
    any(myp060_identity_audit$active_patient_field_changed)) {
  stop(
    "MyP-060 identity audit must preserve active values and resolve to IMG-060.",
    call. = FALSE
  )
}

write_csv(
  myp060_identity_audit,
  "Results_Fragmentomics/fragmentomics_MyP060_identity_audit.csv"
)


# Save “Key_fragmentomics_data.csv” 
write_csv(
  merged_final_metrics_clean,
  file = "Results_Fragmentomics/Key_fragmentomics_data_updated2.csv"
)

saveRDS(
  merged_final_metrics_clean,
  file = "Results_Fragmentomics/Key_fragmentomics_data_updated2.rds"
)



### Now get range of healthy controls 
# Define the specific metrics of interest
selected_metrics <- c(
  "Mean.Coverage",
  "Midpoint.Coverage",
  "Midpoint.normalized",
  "Amplitude"
)

# Compute healthy control ranges for those metrics
healthy_dars <- mm_dars_small %>% 
  filter(!grepl("TGL49_0267_Pb_U_PE_428", Sample)) %>% 
  filter(
    grepl("TGL49", Sample),
    .data$fragmentomics_sequencing_platform == "NovaSeq 6000"
  )

# Test normality 
# 1) Gather into long form
metric_vals <- healthy_dars %>%
  select(all_of(selected_metrics)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value")

# 2) Run Shapiro–Wilk by metric
normality <- metric_vals %>%
  group_by(metric) %>%
  summarise(
    shapiro_p = shapiro.test(value)$p.value,
    .groups   = "drop"
  )

## Assess cutoff
# Compute reference ranges for each MM-DARs coverage metric across
# healthy controls (TGL49 panel). Both Gaussian (mean ± 1.96·SD) and
# empirical (2.5 – 97.5th percentile) bounds are output so downstream
# scripts can choose the appropriate reference interval per metric
# based on the normality test above.
hc_ranges_selected <- healthy_dars %>%
  select(all_of(selected_metrics)) %>%
  pivot_longer(
    cols     = everything(),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  summarise(
    mean_value      = mean(value, na.rm = TRUE),
    sd_value        = sd(value,   na.rm = TRUE),
    lower_gaussian  = mean_value - 1.96 * sd_value,
    upper_gaussian  = mean_value + 1.96 * sd_value,
    lower_empirical = quantile(value, 0.025, na.rm = TRUE),
    upper_empirical = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# View the table
print(hc_ranges_selected)

# Save to CSV in the output directory
readr::write_csv(
  hc_ranges_selected,
  file.path(out.dir, "HC_Ranges_Selected_MM_DARs_Metrics.csv")
)
