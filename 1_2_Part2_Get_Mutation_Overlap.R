# Script: 1_2_Part2_Get_Mutation_Overlap.R
# ──────────────────────────────────────────────────────────────────────────────
# Title       : Mutation Overlap Analysis (BM vs. PB cfDNA)  
# Author      : Dory Abelman  
# Date        : 2025-06-30  
#  
# Description :  
#   1. Reads the annotated bone-marrow and blood cfDNA mutation RDS tables
#      produced by 1_2_Process_Mutation_Data.R
#   2. Merges variant calls, creates per-mutation unique IDs  
#   3. Builds wide table for presence/absence by sample type  
#   4. For each eligible patient:
#        • draws a Venn diagram of BM vs. cfDNA calls  
#        • computes the Jaccard overlap percentage and plots a summary
#  
# Usage from the repository root:
#   Rscript 1_2_Part2_Get_Mutation_Overlap.R
#   - or -  
#   source("1_2_Part2_Get_Mutation_Overlap.R")  
#  
# Dependencies:  
#   library(maftools)    # read.maf()  
#   library(dplyr)       # data wrangling  
#   library(tidyr)       # pivot_wider()  
#   library(ggplot2)     # barplot  
#   library(scales)      # optional formatting  
#   library(grid)        # grid graphics  
#   library(VennDiagram) # draw.pairwise.venn()  
#  
# Input files :
#   - combined_maf_bm_dx.rds
#   - combined_maf_blood_all_muts_updated.rds
#   - combined_clinical_data_updated_April2025.csv
#   - cohort_assignment_table_updated.rds
#   - id_map.rds
#   - Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds
#  
# Output files:
#   Used in the manuscript (all evaluable baseline pairs):
#     Final Tables and Figures/Fig3C_mutation_overlap_lollipop_all_evaluable_baseline.png
#     Final Tables and Figures/Extended_Data_Figure_2G_all_evaluable_baseline_mutation_overlap_*.csv
#   Cohort-specific support version (Frontline only):
#   - Final Tables and Figures/Fig3C_mutation_overlap_lollipop2_updated.png
#   - Final Tables and Figures/Extended_Data_Figure_2G_mutation_overlap_source_data.csv
#   - Final Tables and Figures/Extended_Data_Figure_2G_mutation_overlap_summary.csv
#
#   Support-only QA/review:
#   - Final Tables and Figures/mutation_overlap_support/patient_venn_diagrams/
#     VennDiagram_<Patient>_Sep2025_updatedcolors.png
#   - Final Tables and Figures/mutation_overlap_support/percent_overlap_barplot_Sep2025.png
#   - Final Tables and Figures/mutation_overlap_support/*summary*.csv
#   - Final Tables and Figures/mutation_overlap_support/
#     baseline_bm_cfdna_pairing_union_exact_or_within_30d_audit.csv
#  
# Notes:
#   - The script is command-line runnable from the project root. Clinical
#     metadata are read from disk; no preloaded RStudio workspace object is
#     required.
#   - Analysis unit: patient. Within each eligible patient, unique variant
#     coordinates are compared between BM and cfDNA.
#   - `Percent_Overlap` is 100 * |BM intersection cfDNA| / |BM union cfDNA|.
#     It is a symmetric Jaccard overlap, not the percentage of BM variants
#     recovered in cfDNA (which would instead use the BM count as denominator).
#   - No `t_depth` threshold is applied here. The depth column is carried into
#     the wide table with `max`, but overlap uses mutation presence only.
#   - The Extended Data Figure 2G panel in the final figure PDF is the
#     all-evaluable version and includes both assigned cohorts. The separately
#     retained 28-patient panel contains the Frontline cohort only.
#   - `All_feature_data_Sep2025_updated2.rds`, produced later by 1_5, is loaded
#     only for a printed high-quality-subset summary. Because the file is still
#     an unconditional input, a clean run must have it available already or
#     rerun this script after 1_5.
#  
# How to run from the repository root:
#   Rscript 1_2_Part2_Get_Mutation_Overlap.R
#
# Manuscript outputs created/updated:
#   - Extended Data Figure 2G: all-evaluable patient-level lollipop plot showing
#     the Jaccard percentage overlap between baseline/diagnosis BM and cfDNA
#     mutation sets. A Frontline-only support version is also exported.
#
# Pipeline role:
#   This analysis asks whether baseline tumour mutations seen in diagnostic
#   bone marrow overlap calls in plasma cfDNA. As a conservative filtering rule,
#   calls carrying a dbSNP rsID are removed before overlap calculations. An rsID
#   is not itself proof that a variant is common or germline, so this exclusion
#   should be reported as an operational filter rather than germline adjudication.
# ──────────────────────────────────────────────────────────────────────────────
# Pipeline status:
#   Active in the command-line pipeline. This script creates or stages the
#   manuscript output(s) listed above into final_manuscript_objects/ when the
#   required upstream inputs are available.
#


### Load libraries 
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(stringr)
library(grid)
library(VennDiagram)
library(viridis)

# Shared manuscript-output helpers.
# The manuscript panel produced by this script is copied into
# final_manuscript_objects/ immediately after it is saved below.
.manuscript_helper <- file.path("Scripts_2025", "Final_Scripts", "manuscript_output_helpers.R")
if (!file.exists(.manuscript_helper)) {
  .manuscript_helper <- "manuscript_output_helpers.R"
}
source(.manuscript_helper)
rm(.manuscript_helper)

outdir <- "Final Tables and Figures/"
mutation_overlap_support_dir <- file.path(outdir, "mutation_overlap_support")
venn_output_dir <- file.path(mutation_overlap_support_dir, "patient_venn_diagrams")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(mutation_overlap_support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(venn_output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the MAF file using read.maf
df_bm <- readRDS("combined_maf_bm_dx.rds")
df_blood <- readRDS("combined_maf_blood_all_muts_updated.rds")

## Remove rsids to match new filtering 
# Conservatively exclude variants with an rsID ("rs...") in the dbSNP column to
# match the VCF filtering used by MRDetect. dbSNP also contains rare and
# disease-associated variants, so this rule reduces likely polymorphism carryover
# but does not, by itself, establish somatic or germline status.
# filter out all rows with an RSID
df_blood <- df_blood %>%
  filter(is.na(dbSNP_RS) | dbSNP_RS == "" | dbSNP_RS == "." |
           !grepl("^rs", dbSNP_RS, ignore.case = TRUE))

df_bm <- df_bm %>%
  filter(is.na(dbSNP_RS) | dbSNP_RS == "" | dbSNP_RS == "." |
           !grepl("^rs", dbSNP_RS, ignore.case = TRUE))

# now combine if you like
combined_maf <- bind_rows(df_bm, df_blood)

rm(df_bm)
rm(df_blood)

# Support-only QA: per-patient Venn diagrams between diagnostic BM and cfDNA
# mutations. These help inspect individual overlap patterns but are not mapped
# to a final manuscript panel; the lollipop summary used in the manuscript is
# summary saved later as Extended Data Figure 2G.
# Create a column to indicate presence (1) of the mutation
combined_maf <- combined_maf %>%
  mutate(Presence = 1)

## Re-join with info in case issues 
# Load in the patient info 
.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

metada_df_mutation_comparison <- read_combined_clinical_metadata_with_revision(
  "combined_clinical_data_updated_April2025.csv"
)

# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Tumor_Sample_Barcode = Bam %>%
           # Remove _PG or _WG
           str_remove_all("_PG|_WG") %>%
           # Remove anything after ".filter", ".ded", or ".recalibrate"
           str_replace_all("\\.filter.*|\\.ded.*|\\.recalibrate.*", ""))

# Add the Bam column to combined_maf
combined_maf <- combined_maf %>%
  mutate(Bam = paste0(Tumor_Sample_Barcode, ".filter.deduped.recalibrated.bam"))

# Modify the specific Tumor_Sample_Barcode in combined_maf with error
combined_maf <- combined_maf %>%
  mutate(Tumor_Sample_Barcode = ifelse(Tumor_Sample_Barcode == "TFRIM4_0189_Bm_P_ZC-02", 
                                       paste0(Tumor_Sample_Barcode, "-01-O-DNA"), 
                                       Tumor_Sample_Barcode))

combined_maf <- combined_maf %>%
  select(
    -any_of(c(
      "Patient",
      "Date_of_sample_collection",
      "Sample_type",
      "Timepoint",
      "Study",
      "Sample_ID",
      "timepoint_info",
      "Relapsed",
      "Num_days_to_closest_relapse_absolute",
      "Num_days_to_closest_relapse_non_absolute",
      "Num_days_to_closest_relapse"
    ))
  )

# Join with the metadata dataframe
combined_maf <- left_join(combined_maf %>% select(-Bam), metada_df_mutation_comparison, by = "Tumor_Sample_Barcode")

## First filter to just diagnosis timepoints
combined_maf <- combined_maf %>% filter(timepoint_info %in% c("Diagnosis", "Baseline"))

# Create mutation identifiers only after the metadata join has supplied the
# patient ID used in the analysis. Some legacy MAF rows (including baseline cfDNA for
# IMG-159/IMG-05) have a missing raw Patient field. Building a patient-aware key
# before this join incorrectly makes matched BM/cfDNA variants non-identical.
if (any(is.na(combined_maf$Patient) | !nzchar(combined_maf$Patient))) {
  missing_key_barcodes <- combined_maf %>%
    filter(is.na(Patient) | !nzchar(Patient)) %>%
    distinct(Tumor_Sample_Barcode) %>%
    pull(Tumor_Sample_Barcode)
  stop(
    "Cannot build patient-specific mutation keys after metadata harmonization; ",
    "missing Patient for: ",
    paste(missing_key_barcodes, collapse = ", "),
    call. = FALSE
  )
}

# Two ID fields are built:
#   Unique_ID_patient = chr:start:end:ref:alt:Patient  -> compare BM vs cfDNA
#     within the same harmonized patient.
#   Mutation_ID       = chr:start:end:ref:alt          -> coordinate-level key
#     independent of patient.
combined_maf <- combined_maf %>%
  mutate(
    Unique_ID_patient = paste(
      Chromosome, Start_Position, End_Position, Reference_Allele,
      Tumor_Seq_Allele2, Patient, sep = ":"
    ),
    Mutation_ID = paste(
      Chromosome, Start_Position, End_Position, Reference_Allele,
      Tumor_Seq_Allele2, sep = ":"
    )
  )

## Add cohort df 
cohort_df <- load_final_cohort_assignment()

# Find cohort-assigned patients with no baseline/diagnosis mutation rows after
# the metadata join and timepoint filter. This is diagnostic only; it does not
# change the plotted overlap data.
missing_patients <- setdiff(cohort_df$Patient, combined_maf$Patient)

if (length(missing_patients) == 0) {
  message("✔ All patients in cohort_df are represented in combined_maf.")
} else {
  message("⚠️ These patients are missing from combined_maf: ", 
          paste(missing_patients, collapse = ", "))
}

# Pivot the dataframe to wide format
df_wide <- pivot_wider(
  combined_maf,
  id_cols = c("Hugo_Symbol", 'Patient', "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Unique_ID_patient"),
  names_from = Sample_type,
  values_from = c(VAF, t_depth, Presence),
  values_fn = {max}
)

## Select baseline BM-cfDNA pairs without discarding the historical exact-
## timepoint cohort. The original implementation required identical Timepoint
## labels. That retained valid historical pairs but excluded IMG-142 and
## IMG-235 even though their baseline BM and cfDNA specimens were collected
## only 5 and 7 days apart, respectively. We therefore use a union rule:
##
##   1. retain every patient with an existing exact Patient + Timepoint pair;
##   2. for patients without an exact pair, add the nearest dated
##      baseline/diagnosis BM-cfDNA pair when it is within 30 days.
##
## This is intentionally additive: no historically evaluable patient is
## removed by the date rule. Patients with no callable mutation set remain in
## the downstream audit with Percent_Overlap = NA.
##
## Important scope detail: this block selects eligible patients, not a single
## barcode pair for the mutation calculation. After eligibility is established,
## `df_wide` below contains the union of all Diagnosis/Baseline BM calls and all
## Diagnosis/Baseline cfDNA calls for that patient. The selected nearest dates
## are retained only in the pairing audit; the selected barcode columns are not
## carried into the current audit export.
baseline_pair_metadata <- metada_df_mutation_comparison %>%
  filter(
    timepoint_info %in% c("Baseline", "Diagnosis"),
    Sample_type %in% c("BM_cells", "Blood_plasma_cfDNA")
  ) %>%
  mutate(Date_of_sample_collection = as.Date(Date_of_sample_collection))

exact_timepoint_pairs <- baseline_pair_metadata %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    has_BM = any(Sample_type == "BM_cells"),
    has_cfDNA = any(Sample_type == "Blood_plasma_cfDNA"),
    .groups = "drop"
  ) %>%
  filter(has_BM, has_cfDNA) %>%
  transmute(
    Patient,
    pairing_rule = "exact_timepoint_label",
    BM_Timepoint = as.character(Timepoint),
    cfDNA_Timepoint = as.character(Timepoint),
    BM_date = as.Date(NA),
    cfDNA_date = as.Date(NA),
    abs_days = NA_integer_
  )

nearest_dated_pairs <- baseline_pair_metadata %>%
  filter(!is.na(Date_of_sample_collection)) %>%
  select(
    Patient,
    Sample_type,
    Timepoint,
    Date_of_sample_collection,
    Tumor_Sample_Barcode
  ) %>%
  {
    inner_join(
      filter(., Sample_type == "BM_cells") %>%
        transmute(
          Patient,
          BM_Timepoint = as.character(Timepoint),
          BM_date = Date_of_sample_collection,
          BM_Tumor_Sample_Barcode = Tumor_Sample_Barcode
        ),
      filter(., Sample_type == "Blood_plasma_cfDNA") %>%
        transmute(
          Patient,
          cfDNA_Timepoint = as.character(Timepoint),
          cfDNA_date = Date_of_sample_collection,
          cfDNA_Tumor_Sample_Barcode = Tumor_Sample_Barcode
        ),
      by = "Patient"
    )
  } %>%
  mutate(abs_days = abs(as.integer(BM_date - cfDNA_date))) %>%
  filter(abs_days <= 30L) %>%
  anti_join(exact_timepoint_pairs %>% distinct(Patient), by = "Patient") %>%
  arrange(Patient, abs_days, BM_date, cfDNA_date, BM_Timepoint, cfDNA_Timepoint) %>%
  group_by(Patient) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  mutate(pairing_rule = "nearest_baseline_pair_within_30_days")

pairing_union_audit <- bind_rows(
  exact_timepoint_pairs,
  nearest_dated_pairs %>%
    select(
      Patient,
      pairing_rule,
      BM_Timepoint,
      cfDNA_Timepoint,
      BM_date,
      cfDNA_date,
      abs_days
    )
) %>%
  arrange(Patient, pairing_rule)

expected_date_added_patients <- c("IMG-142", "IMG-235")
missing_expected_date_additions <- setdiff(
  expected_date_added_patients,
  nearest_dated_pairs$Patient
)
if (length(missing_expected_date_additions)) {
  stop(
    "Expected <=30-day baseline pair(s) were not recovered: ",
    paste(missing_expected_date_additions, collapse = ", "),
    call. = FALSE
  )
}

readr::write_csv(
  pairing_union_audit,
  file.path(
    mutation_overlap_support_dir,
    "baseline_bm_cfdna_pairing_union_exact_or_within_30d_audit.csv"
  )
)

patients_with_required_samples <- pairing_union_audit %>%
  distinct(Patient)

# Join the filtered patients with the original df_wide to get the final filtered dataframe
df_wide <- df_wide %>%
  inner_join(patients_with_required_samples, by = "Patient")


## Make plot
patients <- levels(as.factor(patients_with_required_samples$Patient))

## Now create a venn diagram for each patient 
# Loop over all patients
library(scales)
library(grid)
library(VennDiagram)

# Loop through each patient
for (i in 1:length(patients)) {
  
  # Subset the data for the current patient
  df_patient <- df_wide %>% filter(Patient == patients[i])
  
  # Get the mutations for each set
  bm_mutations <- df_patient %>% filter(Presence_BM_cells == 1) %>% pull(Unique_ID_patient) %>% unique() %>% na.omit()
  cf_mutations <- df_patient %>% filter(Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient) %>% unique() %>% na.omit()
  both_mutations <- df_patient %>% filter(Presence_BM_cells == 1 & Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient) %>% unique() %>% na.omit()

  if (length(bm_mutations) == 0 || length(cf_mutations) == 0) {
    message("Skipping Venn diagram for ", patients[i], ": one or both mutation sets are empty.")
    next
  }
  
  # Create a new plot with a white background
  png(
    filename = file.path(venn_output_dir, paste0("VennDiagram_", patients[i], "_Sep2025_updatedcolors.png")),
    width = 2400,
    height = 2400,
    res = 300
  )
  
  # Clear the plot and set background to white
  grid.newpage()
  pushViewport(viewport(gp = gpar(fill = "white")))
  
  # Generate the Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(bm_mutations),
    area2 = length(cf_mutations),
    cross.area = length(intersect(bm_mutations, cf_mutations)),
    category = c("BM cells", "Blood cfDNA"),
    fill = c("#008080", "#FF6347"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("#008080", "#FF6347"),
    margin = 0.1
  )
  
  # Add title using grid graphics
  grid.text(paste0("Patient: ", patients[i]), x = 0.5, y = 0.95, gp = gpar(fontsize = 20))
  grid.draw(venn.plot)
  dev.off()
}


## Make barplot showing percent overlap 

# Create a dataframe to store the overlap percentages
overlap_data <- data.frame(
  Patient = character(),
  BM_Mutation_Count = integer(),
  cfDNA_Mutation_Count = integer(),
  Shared_Mutation_Count = integer(),
  Union_Mutation_Count = integer(),
  Percent_Overlap = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each patient
for (i in 1:length(patients)) {
  
  # Subset the data for the current patient
  df_patient <- df_wide %>% filter(Patient == patients[i])
  
  # Get the mutations for each set
  bm_mutations <- df_patient %>% filter(Presence_BM_cells == 1) %>% pull(Unique_ID_patient)
  cf_mutations <- df_patient %>% filter(Presence_Blood_plasma_cfDNA == 1) %>% pull(Unique_ID_patient)
  
  # Calculate the number of overlapping mutations
  overlap_count <- length(intersect(bm_mutations, cf_mutations))
  
  # Calculate the total number of unique mutations
  total_unique_mutations <- length(unique(c(bm_mutations, cf_mutations)))
  
# Calculate the Jaccard overlap percentage:
#   100 * number of variants called in both compartments /
#         number of variants called in either compartment.
# This is symmetric and is not a BM-recovery/sensitivity denominator.
# If both mutation sets are empty, retain the
  # patient in the audit table but mark the overlap as missing rather than
  # creating NaN/Inf summary statistics.
  percent_overlap <- if (total_unique_mutations > 0) {
    (overlap_count / total_unique_mutations) * 100
  } else {
    NA_real_
  }
  
  # Store the results in the dataframe
  overlap_data <- rbind(
    overlap_data,
    data.frame(
      Patient = patients[i],
      BM_Mutation_Count = length(unique(bm_mutations)),
      cfDNA_Mutation_Count = length(unique(cf_mutations)),
      Shared_Mutation_Count = overlap_count,
      Union_Mutation_Count = total_unique_mutations,
      Percent_Overlap = percent_overlap
    )
  )
}

# Reorder the Patient factor levels based on Percent_Overlap in descending order
overlap_data <- overlap_data %>%
  arrange(desc(Percent_Overlap)) %>%
  mutate(Patient = factor(Patient, levels = unique(Patient)))

# The support overview plot should show only patients with a defined overlap
# percentage. Patients with no callable mutation set remain in the audit CSV.
overlap_plot_data <- overlap_data %>%
  filter(!is.na(Percent_Overlap))

# Create the bar plot
overlap_plot <- ggplot(overlap_plot_data, aes(x = Patient, y = Percent_Overlap)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Percent Overlap Between BM and PB cfDNA mutation calls",
    x = "Patient",
    y = "Percent Overlap"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

# Display the plot
print(overlap_plot)

# Save support-only overview barplot. The final manuscript panel below uses a
# protected-ID lollipop layout instead.
ggsave(
  file.path(mutation_overlap_support_dir, "percent_overlap_barplot_Sep2025.png"),
  plot = overlap_plot,
  width = 6,
  height = 5,
  dpi = 500
)


### Now get the overlap table and summarise stats 

# 1) Left-join the cohort assignments onto the overlap table
overlap_with_cohort <- overlap_data %>%
  left_join(cohort_df, by = "Patient") %>%
  mutate(
    Percent_Overlap = if_else(
      is.nan(Percent_Overlap) | is.infinite(Percent_Overlap),
      NA_real_,
      Percent_Overlap
    )
  )

# Rows with a finite percent overlap are valid for numerical summaries. Patients
# with no callable baseline mutation set are retained in overlap_with_cohort for
# audit but excluded from plotted/summary percentages.
overlap_for_stats <- overlap_with_cohort %>%
  filter(!is.na(Percent_Overlap))

overlap_for_cohort_stats <- overlap_for_stats %>%
  filter(!is.na(Cohort))

# 2) How many patients have both baseline BM and cfDNA samples?
n_patients_overlap <- n_distinct(overlap_with_cohort$Patient)
message("Number of patients with both baseline BM and cfDNA: ", n_patients_overlap)

# 3) Overall summary of Percent_Overlap
overall_stats <- overlap_for_stats %>%
  summarise(
    mean_overlap   = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap    = IQR(Percent_Overlap, na.rm = TRUE),
    min_overlap    = min(Percent_Overlap, na.rm = TRUE),
    max_overlap    = max(Percent_Overlap, na.rm = TRUE)
  )
print(overall_stats)

# 4) Summary by cohort (replace 'Cohort' with the actual cohort column name)
stats_by_cohort <- overlap_for_cohort_stats %>%
  group_by(Cohort) %>%
  summarise(
    n               = dplyr::n(),
    mean_overlap    = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap  = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap     = IQR(Percent_Overlap, na.rm = TRUE),
    min_overlap     = min(Percent_Overlap, na.rm = TRUE),
    max_overlap     = max(Percent_Overlap, na.rm = TRUE)
  )
print(stats_by_cohort)


library(glue)
# Create descriptive sentences for each cohort
cohort_sentences <- stats_by_cohort %>%
  mutate(
    sentence = glue(
      "In the {Cohort} cohort, there were {n} matched baseline samples. ",
      "The percent overlap had a mean of {round(mean_overlap, 1)}% ",
      "(range {round(min_overlap, 1)}%–{round(max_overlap, 1)}%, ",
      "IQR {round(IQR_overlap, 1)}%)."
    )
  ) %>%
  pull(sentence)

# Print them
cat(cohort_sentences, sep = "\n")


### Now make figure 3C 
# 1. Sort overlap_data and add position for plotting
plot_df <- overlap_with_cohort %>%
  filter(Cohort == "Frontline") %>%
  filter(!is.na(Percent_Overlap)) %>%
  arrange(Percent_Overlap) %>%
  mutate(
    Patient = factor(Patient, levels = Patient),      # keep order
    pos     = row_number()                            # y-position
  )

## edit IDs 
# load ID map (Patient -> New_ID)
id_map <- readRDS("id_map.rds") %>% distinct(Patient, New_ID)

# 1) Apply ID map, keep order for plotting
plot_df <- overlap_with_cohort %>%
  filter(Cohort == "Frontline") %>%
  filter(!is.na(Percent_Overlap)) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(Patient = coalesce(New_ID, Patient)) %>%   # swap to New_ID when available
  select(-New_ID) %>%
  arrange(Percent_Overlap) %>%
  mutate(
    Patient = factor(Patient, levels = Patient),    # lock order for y-axis
    pos     = row_number()
  )

## Get stats 
plot_df %>%
  summarise(
    mean_overlap   = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap    = IQR(Percent_Overlap, na.rm = TRUE),
    min    = min(Percent_Overlap, na.rm = TRUE),
    max    = max(Percent_Overlap, na.rm = TRUE)
  )

### see the high quality ones 
All_feature_data <- readRDS("Jan2025_exported_data/All_feature_data_Sep2025_updated2.rds") ## generated in later script

BM_good_pts <- All_feature_data %>%
  filter(Sample_type == "BM_cells",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>%
  unique()

cfDNA_good_pts <- All_feature_data %>%
  filter(Sample_type == "Blood_plasma_cfDNA",
         Evidence_of_Disease == 1,
         timepoint_info %in% c("Diagnosis","Baseline")) %>%
  pull(Patient) %>%
  unique()

# convert BM/cfDNA patient IDs to New_ID first
BM_good_pts_mapped <- BM_good_pts %>%
  tibble(Patient = .) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(Patient = coalesce(New_ID, Patient)) %>%
  pull(Patient) %>%
  unique()

cfDNA_good_pts_mapped <- cfDNA_good_pts %>%
  tibble(Patient = .) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(Patient = coalesce(New_ID, Patient)) %>%
  pull(Patient) %>%
  unique()

# intersect after mapping
common_pts <- intersect(BM_good_pts_mapped, cfDNA_good_pts_mapped)

# 1) Apply ID map, keep order for plotting, restrict to common patients
plot_df %>% filter(Patient %in% common_pts) %>%
  summarise(
    mean_overlap   = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap    = IQR(Percent_Overlap, na.rm = TRUE),
    min    = min(Percent_Overlap, na.rm = TRUE),
    max    = max(Percent_Overlap, na.rm = TRUE)
  )

# 2) Summary stats on the restricted set
plot_df %>%
  summarise(
    n              = n(),
    mean_overlap   = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap    = IQR(Percent_Overlap, na.rm = TRUE),
    IQR_lower      = quantile(Percent_Overlap, 0.25, na.rm = TRUE),
    IQR_upper      = quantile(Percent_Overlap, 0.75, na.rm = TRUE),
    min_overlap    = min(Percent_Overlap, na.rm = TRUE),
    max_overlap    = max(Percent_Overlap, na.rm = TRUE)
  )

# 2. Calculate overall statistics
med  <- median(plot_df$Percent_Overlap, na.rm = TRUE)
iqrL <- quantile(plot_df$Percent_Overlap, 0.25, na.rm = TRUE)
iqrU <- quantile(plot_df$Percent_Overlap, 0.75, na.rm = TRUE)

# 3. Build lollipop plot
p_overlap <- ggplot(plot_df, aes(x = Percent_Overlap, y = Patient)) +
 # annotate(                         # IQR ribbon
 #   "rect", xmin = iqrL, xmax = iqrU,
 #   ymin = 0.5, ymax = nrow(plot_df) + 0.5,
 #   fill = "grey90"
 # ) +
  geom_vline(xintercept = med, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0, xend = Percent_Overlap, yend = Patient),
               colour = "grey65", linewidth = 0.35) +
  geom_point(aes(colour = Percent_Overlap), size = 2) +
  scale_colour_viridis_c(
    option = "D", end = 0.9, name = "% overlap",
    guide  = guide_colourbar(barwidth = 0.4, barheight = 3)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Mutation overlap: BM ∩ cfDNA (%)",
    y = NULL,
    title    = "Patient-level overlap of baseline mutations (BM vs. cfDNA)",
    subtitle = glue::glue(
      "Median = {round(med,1)}%   |   IQR = {round(iqrL,1)}–{round(iqrU,1)}%"
    )
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title      = element_text(face = "bold", size = 9),
    plot.subtitle   = element_text(size = 7, margin = margin(b = 6)),
    axis.text.y     = element_text(size = 6),
    panel.grid.major.y = element_blank(),
    legend.position     = c(0.97, 0.05),     # inside bottom-right
    legend.justification = c(1, 0),
    legend.background    = element_rect(
      fill = scales::alpha("white", 0.7), colour = NA
    ),
    legend.key.width  = unit(0.3, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.title      = element_text(size = 6),
    legend.text       = element_text(size = 6)
  )

# 4. Save the Frontline-only support panel.
# The final Extended Data Figure 2G PDF uses the all-evaluable panel below.
# Historical filename: Fig3C_mutation_overlap_lollipop2_updated.png.
edfig2g_path <- file.path(outdir, "Fig3C_mutation_overlap_lollipop2_updated.png")
ggsave(
  filename = edfig2g_path,
  plot     = p_overlap,
  width    = 4,   # a bit narrower
  height   = 5,
  dpi      = 600
)

# SUPPORT OUTPUT: Frontline-cohort version of Extended Data Figure 2G.
# The artifact ID is retained for compatibility with the existing output map.
ms_copy_artifact(
  source_path = edfig2g_path,
  artifact_id = "EDFIG2G",
  role = "figure_panel_png",
  description = "Extended Data Figure 2G: baseline patient-level BM/cfDNA mutation-overlap lollipop plot.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)


# Export source/helper tables for the Frontline-cohort support version.
edfig2g_source_data_path <- file.path(outdir, "Extended_Data_Figure_2G_mutation_overlap_source_data.csv")
edfig2g_summary_path <- file.path(outdir, "Extended_Data_Figure_2G_mutation_overlap_summary.csv")

edfig2g_summary <- plot_df %>%
  summarise(
    n = n(),
    mean_overlap = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap = IQR(Percent_Overlap, na.rm = TRUE),
    IQR_lower = quantile(Percent_Overlap, 0.25, na.rm = TRUE),
    IQR_upper = quantile(Percent_Overlap, 0.75, na.rm = TRUE),
    min_overlap = min(Percent_Overlap, na.rm = TRUE),
    max_overlap = max(Percent_Overlap, na.rm = TRUE)
  )

write.csv(plot_df, edfig2g_source_data_path, row.names = FALSE)
write.csv(edfig2g_summary, edfig2g_summary_path, row.names = FALSE)

ms_copy_artifact(
  source_path = edfig2g_source_data_path,
  artifact_id = "EDFIG2G",
  role = "source_data_csv",
  description = "Extended Data Figure 2G source data: patient-level baseline BM/cfDNA mutation overlap used for the lollipop plot.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)
ms_copy_artifact(
  source_path = edfig2g_summary_path,
  artifact_id = "EDFIG2G",
  role = "summary_csv",
  description = "Extended Data Figure 2G summary statistics for baseline BM/cfDNA mutation overlap.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)

# Extended Data Figure 2G: all cohort-assigned evaluable diagnosis/baseline
# BM-cfDNA pairs. The final figure PDF uses this version. It keeps every matched
# baseline/diagnosis patient pair with defined mutation overlap, including the
# test cohort, after the same upstream metadata joins and pair filters.
plot_df_all_evaluable <- overlap_with_cohort %>%
  filter(!is.na(Cohort), !is.na(Percent_Overlap)) %>%
  left_join(id_map, by = "Patient") %>%
  mutate(
    Patient = coalesce(New_ID, Patient),
    Cohort = case_when(
      Cohort == "Frontline" ~ "Training Cohort",
      Cohort == "Non-frontline" ~ "Test Cohort",
      TRUE ~ Cohort
    )
  ) %>%
  select(-New_ID) %>%
  arrange(Percent_Overlap) %>%
  mutate(
    Patient = factor(Patient, levels = Patient),
    pos = row_number()
  )

all_eval_med <- median(plot_df_all_evaluable$Percent_Overlap, na.rm = TRUE)
all_eval_iqr_l <- quantile(plot_df_all_evaluable$Percent_Overlap, 0.25, na.rm = TRUE)
all_eval_iqr_u <- quantile(plot_df_all_evaluable$Percent_Overlap, 0.75, na.rm = TRUE)

p_overlap_all_evaluable <- ggplot(
  plot_df_all_evaluable,
  aes(x = Percent_Overlap, y = Patient)
) +
  geom_vline(xintercept = all_eval_med, linetype = "dashed", colour = "grey40") +
  geom_segment(
    aes(x = 0, xend = Percent_Overlap, yend = Patient),
    colour = "grey65",
    linewidth = 0.35
  ) +
  geom_point(aes(colour = Percent_Overlap), size = 2) +
  scale_colour_viridis_c(
    option = "D",
    end = 0.9,
    name = "% overlap",
    guide = guide_colourbar(barwidth = 0.4, barheight = 3)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Mutation overlap: BM ∩ cfDNA (%)",
    y = NULL,
    title = "Patient-level overlap of baseline mutations (BM vs. cfDNA)",
    subtitle = glue::glue(
      "All evaluable baseline pairs (n = {nrow(plot_df_all_evaluable)}); median {round(all_eval_med, 1)}% (IQR {round(all_eval_iqr_l, 1)}-{round(all_eval_iqr_u, 1)}%)"
    )
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7, margin = margin(b = 6)),
    axis.text.y = element_text(size = 6),
    panel.grid.major.y = element_blank(),
    legend.position = c(0.97, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.7),
      colour = NA
    ),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )

edfig2g_all_evaluable_path <- file.path(
  outdir,
  "Fig3C_mutation_overlap_lollipop_all_evaluable_baseline.png"
)
ggsave(
  filename = edfig2g_all_evaluable_path,
  plot = p_overlap_all_evaluable,
  width = 4.6,
  height = 5,
  dpi = 600
)

edfig2g_all_evaluable_source_data_path <- file.path(
  outdir,
  "Extended_Data_Figure_2G_all_evaluable_baseline_mutation_overlap_source_data.csv"
)
edfig2g_all_evaluable_summary_path <- file.path(
  outdir,
  "Extended_Data_Figure_2G_all_evaluable_baseline_mutation_overlap_summary.csv"
)

edfig2g_all_evaluable_summary <- plot_df_all_evaluable %>%
  summarise(
    n = n(),
    n_training_cohort = sum(Cohort == "Training Cohort", na.rm = TRUE),
    n_test_cohort = sum(Cohort == "Test Cohort", na.rm = TRUE),
    mean_overlap = mean(Percent_Overlap, na.rm = TRUE),
    median_overlap = median(Percent_Overlap, na.rm = TRUE),
    IQR_overlap = IQR(Percent_Overlap, na.rm = TRUE),
    IQR_lower = quantile(Percent_Overlap, 0.25, na.rm = TRUE),
    IQR_upper = quantile(Percent_Overlap, 0.75, na.rm = TRUE),
    min_overlap = min(Percent_Overlap, na.rm = TRUE),
    max_overlap = max(Percent_Overlap, na.rm = TRUE)
  )

write.csv(
  plot_df_all_evaluable,
  edfig2g_all_evaluable_source_data_path,
  row.names = FALSE
)
write.csv(
  edfig2g_all_evaluable_summary,
  edfig2g_all_evaluable_summary_path,
  row.names = FALSE
)

# The role and description strings below are retained historical manifest
# labels. Despite the word "alternate" in those strings, the all-evaluable
# panel above is the version used as Extended Data Figure 2G; the Frontline-only
# panel is a cohort-specific support version.
ms_copy_artifact(
  source_path = edfig2g_all_evaluable_path,
  artifact_id = "EDFIG2G",
  role = "all_evaluable_figure_panel_png",
  description = "Alternate all-evaluable baseline/diagnosis BM/cfDNA mutation-overlap lollipop plot for Extended Data Figure 2G.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)
ms_copy_artifact(
  source_path = edfig2g_all_evaluable_source_data_path,
  artifact_id = "EDFIG2G",
  role = "all_evaluable_source_data_csv",
  description = "Source data for the alternate all-evaluable baseline/diagnosis BM/cfDNA mutation-overlap lollipop plot.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)
ms_copy_artifact(
  source_path = edfig2g_all_evaluable_summary_path,
  artifact_id = "EDFIG2G",
  role = "all_evaluable_summary_csv",
  description = "Summary statistics for the alternate all-evaluable baseline/diagnosis BM/cfDNA mutation-overlap lollipop plot.",
  script_name = "1_2_Part2_Get_Mutation_Overlap.R"
)

# Export broader support tables for audit/review. These are not direct
# manuscript panel files and are intentionally kept out of the project root.
write.csv(
  overlap_with_cohort,
  file.path(mutation_overlap_support_dir, "mutation_overlap_with_cohort.csv"),
  row.names = FALSE
)
write.csv(
  overall_stats,
  file.path(mutation_overlap_support_dir, "overall_percent_overlap_stats.csv"),
  row.names = FALSE
)
write.csv(
  stats_by_cohort,
  file.path(mutation_overlap_support_dir, "percent_overlap_stats_by_cohort.csv"),
  row.names = FALSE
)
