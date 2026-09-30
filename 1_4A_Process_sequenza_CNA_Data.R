# =============================================================================
# Script: 1_4A_Process_sequenza_CNA_Data.R
#
# Description:
#   End-to-end arm-level CNA caller for Sequenza segments with **ploidy-aware**
#   labels and the project-specific **hyperdiploidy** definition.
#
# Goal:
#   Convert Sequenza segment calls into ploidy-aware CNA feature tables that
#   complement the ichorCNA outputs from 1_4. This script creates arm-level
#   calls for del1p, amp1q, del13q, del17p, hyperdiploidy, FISH-probe-level
#   Sequenza calls, cytoband-expanded FISH calls, and sample-level Sequenza
#   purity/ploidy estimates.
#
#   What it does:
#     1) Reads all Sequenza “*_segments.txt(.gz)” in `seg_dir`
#        – Requires columns: chromosome, start.pos, end.pos, CNt, A, B
#     2) Reads Sequenza confints files for purity/ploidy and estimates a
#        length-weighted autosomal CNt-mode fallback when confints ploidy is
#        unavailable; saves ploidy/purity audit tables
#     3) Converts each segment to a **baseline-aware categorical call**
#        (LOSS / GAIN / AMP / HLAMP; plus HOMD/HETD and CNLOH at baseline with B=0)
#     4) Merges all samples on (chr, start, end) into `combined_seg_data`
#     5) Maps segments to chromosomal **arms** (1p/1q, …) via cytobands using
#        GRanges and assigns the arm with **maximum bp overlap**
#     6) Precomputes **fixed denominators**:
#        – `arm_lengths` for each p/q arm (bp)
#        – `chr_lengths` for whole chromosomes (bp)
#     7) Computes **length-weighted altered fraction per arm** using fixed arm
#        lengths as denominators, then binarizes:
#        – del1p, amp1q, del13q, del17p (default threshold: > 1/3 of arm length)
#     8) Calls **hyperdiploidy** using relaxed myeloma-appropriate criteria:
#        – A chromosome is considered gained if > `gain_thresh_chr` (default 0.65)
#          of its length is GAIN/HLAMP
#        – A sample is hyperdiploid if ≥ `k_trisomies` (default 5) of
#          {3,5,7,9,11,15,19,21} are gained
#     9) Calls clinical FISH-probe regions both from direct probe coordinates
#        and cytoband-expanded windows for concordance/audit analyses
#    10) Exports:
#        – Combined per-segment calls: CSV/RDS
#        – Per-sample ploidy estimates: CSV
#        – Final Sample × feature matrix (`cna_data`) with
#          {del1p, amp1q, del13q, del17p, hyperdiploid}: RDS/TXT
#        – FISH-style Sequenza call matrices: RDS/TXT
#    11) Prints cohort-level frequencies for quick QC
#
# How to run from the repository root:
#   Rscript 1_4A_Process_sequenza_CNA_Data.R
#
# Manuscript outputs created/updated:
#   - None directly. This upstream CNA-processing script creates Sequenza CNA
#     feature tables for WGS feature integration, CNA/FISH concordance,
#     Supplementary Table 2, and downstream subclonal-evolution support.
#
# Inputs:
#   • Scripts_2025/Final_Scripts/helpers.R
#   • Oct 2024 data/Sequenza/All_Segments_400/*_segments.txt(.gz)
#   • cytoband.txt = hg38 cytobands with columns chr, start, end, band, stain;
#                     hg38 is required because the FISH probe coordinates are hg38
#   • Oct 2024 data/Sequenza/All_confints_400/*_confints_CP.txt
#   • combined_clinical_data_updated_April2025.csv
#   • Clinical data/FISH probe locations.xlsx
#   • Optional Spring 2026 Sequenza_segments.tsv and
#       Sequenza_ploidy_purity.tsv inputs under Data_Spring_2026_Revisions/
#
# Analysis units:
#   The segment layer contains one genomic segment per row and one sample per
#   call column. Final CNA, probe, and ploidy tables contain one row per sample.
#   Arm/chromosome fractions are base-pair weighted against fixed hg38 lengths;
#   they are not proportions of abnormal cells.
#
# Active downstream outputs (gamma = 400 run):
#   • Jan2025_exported_data/cna_data_from_sequenza_400_updated.rds/.txt
#       - consumed by 1_5_Integrate_WGS_Feature_Data.R
#   • Jan2025_exported_data/FISH_data_from_sequenza_400_updated.rds/.txt
#       - consumed by 1_5_Integrate_WGS_Feature_Data.R
#   • Jan2025_exported_data/Sample_ploidy_from_sequenza_400.rds
#       - chosen ploidy baseline/source joined to patient metadata; consumed by
#         2_3_Feature_Concordance_And_Mutation_Counts.R
#   • Jan2025_exported_data/Sample_ploidy_from_sequenza.txt
#       - readable companion for the active ploidy RDS
#
# Support/QC outputs:
#   • Output_tables_2025/sequenza_cna_processing_support/Sep_2025_combined_sequenza_calls_400_updated.csv/.rds
#   • Output_tables_2025/sequenza_cna_processing_support/Sep_2025_FISH_probe_calls_400_updated.csv/.rds
#   • Output_tables_2025/sequenza_cna_processing_support/Sep_2025_sequenza_ploidy_estimates_updated.csv
#   • Output_tables_2025/sequenza_cna_processing_support/sequenza_purity_ploidy_estimates.txt
#   • Output_tables_2025/sequenza_cna_processing_support/FISH_data_from_sequenza_400_by_cytoband_updated.rds/.txt
#   • Output_tables_2025/sequenza_cna_processing_support/cna_data_summary_400_updated.rds/.txt
#   • Console:  one-row tibble with cohort proportions for each feature
#
# Key parameters (tunable):
#   • Arm binary threshold:        > 1/3 of fixed hg38 arm length
#   • Hyperdiploid per-chr gain:   gain_thresh_chr = 0.65
#   • Hyperdiploid sample rule:    k_trisomies = 5 of 8 chromosomes
#
# Assumptions & notes:
#   • Calls are **baseline-aware** per sample (neutral ≈ round(ploidy)).
#     Arm/probe binning uses `gain_labels = {GAIN, AMP, HLAMP}` and
#     `loss_labels = {LOSS, HETD, HOMD, CNLOH}`.
#   • Denominators are **fixed** (full arm/chrom lengths) to avoid NA inflation.
#     The current arm and chromosome joins replace a wholly missing proportion
#     with 0, so an unevaluable arm/chromosome is represented as not altered.
#   • At FISH probes, an unavailable categorical call remains NA but `%in%`
#     converts the corresponding binary alteration indicator to 0.
#   • Genome build of `cb_frame` must match BAM/Sequenza build.
#   • This run uses Sequenza gamma = 400, encoded in its input/output paths.
#
# Dependencies:
#   library(tidyverse)  # dplyr, tidyr, readr, purrr
#   library(GenomicRanges); library(IRanges); library(S4Vectors); library(readxl)
#
# Usage from an interactive R session opened at the repository root:
#   source("1_4A_Process_sequenza_CNA_Data.R")
#   # Produces: combined_seg_data, sample_ploidy, results, cna_data, and exports
#
# Author: Dory Abelman
# Date:   2025-10-06
# =============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#
# Reproducibility note:
#   The preferred end-to-end path starts from raw Sequenza segment files in
#   Oct 2024 data/Sequenza/All_Segments_400/ plus matching confints files in
#   Oct 2024 data/Sequenza/All_confints_400/. Non-consumed segment, probe, and
#   ploidy audit files are written to a support folder so they cannot be
#   mistaken for final manuscript tables.
#



# Load Libraries:
library(tidyverse)       # dplyr, purrr, tidyr, readr, etc.
library(purrr)           # for reduce()
library(tidyr)           # for pivot_longer()
library(GenomicRanges)
library(readxl)
library(dplyr)           # reattach after Bioconductor packages so bare verbs are tidyverse verbs

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

## Load clinical info
# Load in the patient info 
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

metada_df_mutation_comparison <- metada_df_mutation_comparison %>%
  mutate(Bam_clean_tmp = gsub(".bam$", "", Bam))  # Remove the '.bam' suffix


seg_dir <- "Oct 2024 data/Sequenza/All_Segments_400/"

export_dir <- "Jan2025_exported_data"
support_table_dir <- file.path("Output_tables_2025", "sequenza_cna_processing_support")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(support_table_dir, recursive = TRUE, showWarnings = FALSE)

support_file <- function(filename) {
  file.path(support_table_dir, filename)
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

# ---- helper: derive ichor-like categorical calls from Sequenza fields
call_from_CNt_AB <- function(CNt, A, B) {
  case_when(
    is.na(CNt)               ~ NA_character_,
    CNt == 0                 ~ "HOMD",
    CNt == 1                 ~ "HETD",
    CNt == 2 & !is.na(B) & B == 0 ~ "CNLOH",        # retained for CNA audit tables
    CNt == 3                 ~ "GAIN",
    CNt >= 4 & CNt < 6       ~ "AMP",
    CNt >= 6                 ~ "HLAMP",
    TRUE                     ~ "NEUT"
  )
}

# Reference helper for a simplified baseline-aware caller. The executed segment
# loops below use their own inline case_when(), which additionally distinguishes
# AMP from GAIN; this function is currently not called.
call_from_CNt_baseline <- function(CNt, A, B, baseline) {
  case_when(
    is.na(CNt)               ~ NA_character_,
    CNt == 0                 ~ "HOMD",
    CNt == 1                 ~ "HETD",
    CNt >= baseline + 3      ~ "HLAMP",
    CNt >  baseline          ~ "GAIN",
    CNt <  baseline          ~ "LOSS",
    # optional: CNLOH when total equals baseline but B==0
    CNt == baseline & !is.na(B) & B == 0 ~ "CNLOH",
    TRUE                      ~ "NEUT"
  )
}

# length-weighted mode for discrete values
wmode <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  t <- tapply(w[ok], x[ok], sum, na.rm = TRUE)
  as.numeric(names(t)[which.max(t)])
}


# ---- read Sequenza *_segments.txt (optionally .gz), rename cols, make calls
seg_files <- list.files(seg_dir, pattern = "_segments\\.txt(\\.gz)?$", full.names = TRUE)

seg_data_list <- vector("list", length(seg_files))
names(seg_data_list) <- basename(seg_files)

ploidy_list <- vector("list", length(seg_files))   


### Get the ploidy from confints
### Lastly, get the purity and ploidy estimates
confints_dir <- "Oct 2024 data/Sequenza/All_confints_400/"

# Discover files that end with _confints_CP.txt
confints <- list.files(
  confints_dir,
  pattern = "_confints_CP\\.txt$",
  full.names = TRUE
)

if (length(confints) == 0L) {
  stop("No *_confints_CP.txt files found in: ", confints_dir)
}

# Sequenza confints tables contain lower, central, and upper candidate solutions;
# row 2 is the central estimate used here. Files with fewer than two data rows
# cannot supply that estimate and remain explicitly missing.
read_pp_safe <- function(path) {
  df <- suppressMessages(read_tsv(path, show_col_types = FALSE, progress = FALSE))
  if (nrow(df) < 2L) {
    warning("File has <2 rows, skipping values: ", basename(path))
    return(tibble(
      Sample = NA_character_,
      Purity = NA_real_,
      Ploidy = NA_real_,
      File   = path
    ))
  }
  
  # Some sequenza outputs use 'cellularity' and 'ploidy.estimate'
  # Normalize possible name variants just in case.
  nm <- names(df)
  nm <- nm |>
    str_replace("^cellularity$", "cellularity") |>
    str_replace("^ploidy\\.estimate$", "ploidy.estimate")
  names(df) <- nm

  row2 <- df |> dplyr::slice(2)
  
  # Extract clean sample ID from filename
  sample_id <- basename(path) |> str_remove("_confints_CP\\.txt$")
  
  tibble(
    Sample = sample_id,
    Purity = suppressWarnings(as.numeric(row2$cellularity)),
    Ploidy = suppressWarnings(as.numeric(row2$`ploidy.estimate`)),
    File   = path
  )
}

pp_table <- map_dfr(confints, read_pp_safe) |>
  # Keep only labeled columns for downstream use
  dplyr::select(Sample, Purity, Ploidy)

# Basic sanity checks + labeling polish
pp_table <- pp_table |>
  mutate(
    Purity = round(Purity, 4),
    Ploidy = round(Ploidy, 3)
  ) |>
  arrange(Sample)

spring2026_ploidy_files <- spring2026_revision_files(
  # Optional Sequenza purity/ploidy summary from the Spring 2026 DNA pipeline.
  # These rows extend pp_table for new revision samples so segment calls can be
  # interpreted relative to sample-specific ploidy instead of a generic diploid
  # baseline.
  "DNA_pipeline_suite_Sequenza_outputs",
  "Sequenza_ploidy_purity[.]tsv$"
)
if (length(spring2026_ploidy_files) > 0L) {
  spring2026_pp <- purrr::map_dfr(
    spring2026_ploidy_files,
    function(path) {
      readr::read_tsv(path, show_col_types = FALSE) %>%
        dplyr::mutate(.sequenza_source_file = basename(path))
    }
  )

  required_spring_pp_cols <- c("Sample", "cellularity", "ploidy")
  missing_spring_pp_cols <- setdiff(required_spring_pp_cols, names(spring2026_pp))
  if (length(missing_spring_pp_cols)) {
    stop(
      "Spring 2026 Sequenza purity/ploidy tables are missing columns: ",
      paste(missing_spring_pp_cols, collapse = ", "),
      call. = FALSE
    )
  }

  duplicate_spring_pp_samples <- spring2026_pp %>%
    dplyr::filter(!is.na(.data$Sample), nzchar(as.character(.data$Sample))) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      n_source_files = dplyr::n_distinct(.data$.sequenza_source_file),
      source_files = paste(sort(unique(.data$.sequenza_source_file)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_source_files > 1L)
  if (nrow(duplicate_spring_pp_samples) > 0L) {
    stop(
      "Spring 2026 Sequenza purity/ploidy samples occur in multiple input files; ",
      "resolve the ambiguous source before integration: ",
      paste(
        paste0(duplicate_spring_pp_samples$Sample, " [", duplicate_spring_pp_samples$source_files, "]"),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  spring2026_pp <- spring2026_pp %>%
    transmute(
      Sample = as.character(Sample),
      Purity = round(suppressWarnings(as.numeric(cellularity)), 4),
      Ploidy = round(suppressWarnings(as.numeric(ploidy)), 3)
    )
  pp_table <- bind_rows(pp_table, spring2026_pp) %>%
    # Historical rows win if the same Sample already exists. The Spring file is
    # an extension path, not an override of previously reviewed Sequenza outputs.
    distinct(Sample, .keep_all = TRUE) %>%
    arrange(Sample)
}

write.table(
  pp_table,
  support_file("sequenza_purity_ploidy_estimates.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



ploidy_map <- setNames(pp_table$Ploidy, pp_table$Sample)

ploidy_list   <- vector("list", length(seg_files))
seg_data_list <- vector("list", length(seg_files))

for (i in seq_along(seg_files)) {
  f  <- seg_files[i]
  df <- readr::read_tsv(f, show_col_types = FALSE)
  
  # required columns in Sequenza segments
  req <- c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(f), ": ", paste(miss, collapse = ", "))
  
  # sample id must match the id used in pp_table$Sample
  sample_name <- sub("_segments\\.txt(\\.gz)?$", "", basename(f))
  
  # ---- core numeric frame
  seg_core <- df %>%
    dplyr::transmute(
      chr   = stringr::str_remove(chromosome, "^chr"),
      start = suppressWarnings(as.numeric(start.pos)),
      end   = suppressWarnings(as.numeric(end.pos)),
      CNt   = suppressWarnings(as.numeric(CNt)),
      A     = suppressWarnings(as.numeric(A)),
      B     = suppressWarnings(as.numeric(B))
    ) %>%
    dplyr::mutate(seg_len = pmax(end - start + 1, 1))
  
  # ---- try to use confints ploidy
  ploidy_conf <- unname(ploidy_map[sample_name])
  
  # If confints missing, fall back to autosome-weighted mode of rounded CNt
  if (!is.finite(ploidy_conf)) {
    df_auto <- seg_core %>%
      dplyr::filter(chr %in% as.character(1:22)) %>%
      dplyr::mutate(CNt_round = round(CNt))
    est_ploidy <- wmode(df_auto$CNt_round, df_auto$seg_len)  # assumes wmode() defined upstream
    baseline_int <- ifelse(is.na(est_ploidy), NA_real_, pmax(1, round(est_ploidy)))
    source_tag   <- "segments_fallback"
  } else {
    # Use Sequenza confints ploidy; snap to nearest integer for labeling
    baseline_int <- pmax(1, round(ploidy_conf))
    source_tag   <- "confints"
  }
  
  # record ploidy sources for QC
  ploidy_list[[i]] <- tibble::tibble(
    Sample        = sample_name,
    ploidy_conf   = ifelse(is.finite(ploidy_conf), ploidy_conf, NA_real_),
    baseline_int  = baseline_int,
    baseline_src  = source_tag
  )
  
  # ---- build per-segment categorical calls (baseline-aware; fallback to diploid if NA)
  if (is.na(baseline_int)) {
    seg_calls <- seg_core %>%
      dplyr::mutate(Call = call_from_CNt_AB(CNt, A, B))  # the fallback copy-number rule
  } else {
    seg_calls <- seg_core %>%
      dplyr::mutate(
        Call = dplyr::case_when(
          !is.finite(CNt)                    ~ NA_character_,
          CNt == 0                           ~ "HOMD",
          CNt == 1                           ~ "HETD",
          CNt >= baseline_int + 3            ~ "HLAMP",
          CNt >  baseline_int + 1            ~ "AMP",
          CNt >  baseline_int                ~ "GAIN",
          CNt <  baseline_int                ~ "LOSS",
          CNt == baseline_int & !is.na(B) & B == 0 ~ "CNLOH",
          CNt == baseline_int                ~ "NEUT",
          TRUE                               ~ NA_character_
        )
      )
  }
  
  seg_df <- seg_calls %>% dplyr::select(chr, start, end, Call)
  colnames(seg_df)[4] <- sample_name
  seg_data_list[[i]] <- seg_df
}

spring2026_segment_files <- spring2026_revision_files(
  # Optional combined Spring 2026 Sequenza segment table. It is parsed in addition
  # to per-sample historical segment files and only contributes samples not
  # already present in seg_data_list.
  "DNA_pipeline_suite_Sequenza_outputs",
  "Sequenza_segments[.]tsv$"
)
if (length(spring2026_segment_files) > 0L) {
  spring2026_segments <- purrr::map_dfr(
    spring2026_segment_files,
    function(path) {
      readr::read_tsv(path, show_col_types = FALSE) %>%
        dplyr::mutate(.sequenza_source_file = basename(path))
    }
  )
  required_spring_cols <- c("Sample", "chromosome", "start.pos", "end.pos", "CNt", "A", "B")
  missing_spring_cols <- setdiff(required_spring_cols, names(spring2026_segments))
  if (length(missing_spring_cols)) {
    stop(
      "Spring 2026 Sequenza segment table is missing columns: ",
      paste(missing_spring_cols, collapse = ", "),
      call. = FALSE
    )
  }

  duplicate_spring_segment_samples <- spring2026_segments %>%
    dplyr::filter(!is.na(.data$Sample), nzchar(as.character(.data$Sample))) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      n_source_files = dplyr::n_distinct(.data$.sequenza_source_file),
      source_files = paste(sort(unique(.data$.sequenza_source_file)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_source_files > 1L)
  if (nrow(duplicate_spring_segment_samples) > 0L) {
    stop(
      "Spring 2026 Sequenza segment samples occur in multiple input files; ",
      "resolve the ambiguous source before integration: ",
      paste(
        paste0(
          duplicate_spring_segment_samples$Sample,
          " [", duplicate_spring_segment_samples$source_files, "]"
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  spring_samples <- sort(unique(spring2026_segments$Sample))
  existing_samples <- vapply(
    seg_data_list,
    function(x) {
      sample_columns <- setdiff(names(x), c("chr", "start", "end"))
      if (length(sample_columns) != 1L) {
        stop("Each historical Sequenza segment table must contain exactly one sample column.")
      }
      sample_columns
    },
    character(1)
  )
  for (sample_name in setdiff(spring_samples, existing_samples)) {
    # Process each revision sample with the same CNt/A/B -> categorical-call
    # logic used for historical Sequenza segments, so downstream CNA summaries do
    # not depend on which input file format supplied the sample.
    df <- spring2026_segments %>% filter(Sample == sample_name)
    seg_core <- df %>%
      dplyr::transmute(
        chr   = stringr::str_remove(chromosome, "^chr"),
        start = suppressWarnings(as.numeric(start.pos)),
        end   = suppressWarnings(as.numeric(end.pos)),
        CNt   = suppressWarnings(as.numeric(CNt)),
        A     = suppressWarnings(as.numeric(A)),
        B     = suppressWarnings(as.numeric(B))
      ) %>%
      dplyr::mutate(seg_len = pmax(end - start + 1, 1))

    ploidy_conf <- unname(ploidy_map[sample_name])
    if (!is.finite(ploidy_conf)) {
      # If Spring purity/ploidy is missing for a sample, fall back to an
      # autosome-weighted modal CNt estimate and record that source in the QC
      # ploidy table.
      df_auto <- seg_core %>%
        dplyr::filter(chr %in% as.character(1:22)) %>%
        dplyr::mutate(CNt_round = round(CNt))
      est_ploidy <- wmode(df_auto$CNt_round, df_auto$seg_len)
      baseline_int <- ifelse(is.na(est_ploidy), NA_real_, pmax(1, round(est_ploidy)))
      source_tag <- "spring2026_segments_fallback"
    } else {
      baseline_int <- pmax(1, round(ploidy_conf))
      source_tag <- "spring2026_ploidy_purity"
    }

    ploidy_list[[length(ploidy_list) + 1L]] <- tibble::tibble(
      Sample = sample_name,
      ploidy_conf = ifelse(is.finite(ploidy_conf), ploidy_conf, NA_real_),
      baseline_int = baseline_int,
      baseline_src = source_tag
    )

    if (is.na(baseline_int)) {
      seg_calls <- seg_core %>%
        dplyr::mutate(Call = call_from_CNt_AB(CNt, A, B))
    } else {
      seg_calls <- seg_core %>%
        dplyr::mutate(
          Call = dplyr::case_when(
            !is.finite(CNt)                    ~ NA_character_,
            CNt == 0                           ~ "HOMD",
            CNt == 1                           ~ "HETD",
            CNt >= baseline_int + 3            ~ "HLAMP",
            CNt >  baseline_int + 1            ~ "AMP",
            CNt >  baseline_int                ~ "GAIN",
            CNt <  baseline_int                ~ "LOSS",
            CNt == baseline_int & !is.na(B) & B == 0 ~ "CNLOH",
            CNt == baseline_int                ~ "NEUT",
            TRUE                               ~ NA_character_
          )
        )
    }

    seg_df <- seg_calls %>% dplyr::select(chr, start, end, Call)
    colnames(seg_df)[4] <- sample_name
    seg_data_list[[sample_name]] <- seg_df
  }
}

combined_seg_data <- purrr::reduce(seg_data_list, dplyr::full_join, by = c("chr","start","end")) %>%
  dplyr::arrange(factor(chr, levels = c(as.character(1:22), "X", "Y")), start)

sample_ploidy <- dplyr::bind_rows(ploidy_list) %>%
  dplyr::arrange(Sample)


print(sample_ploidy, n = 50)
write.csv(sample_ploidy, support_file("Sep_2025_sequenza_ploidy_estimates_updated.csv"), row.names = FALSE)


# Harmonize chr labels in segments
combined_seg_data <- combined_seg_data %>%
  dplyr::mutate(chr = toupper(as.character(chr))) %>%
  dplyr::mutate(chr = gsub("^CHR", "", chr)) %>%
  dplyr::mutate(chr = dplyr::recode(chr, `23`="X", `24`="Y"))

valid_chr <- c(as.character(1:22), "X", "Y")
combined_seg_data <- combined_seg_data %>% dplyr::filter(chr %in% valid_chr)

cytoband_txt <- "cytoband.txt"
if (!file.exists(cytoband_txt)) {
  stop("Missing required hg38 cytoband file: ", cytoband_txt,
       ". This file is needed to map Sequenza CNA segments to chromosome arms.")
}
cb_frame <- readr::read_tsv(
  cytoband_txt,
  col_names = c("chr", "start", "end", "band", "stain"),
  show_col_types = FALSE
)

# Ensure cytoband column is named 'band' (some tables use 'name')
if (!"band" %in% names(cb_frame) && "name" %in% names(cb_frame)) {
  cb_frame <- cb_frame %>% dplyr::rename(band = name)
}
stopifnot(all(c("chr","start","end","band") %in% names(cb_frame)))

# Collapse cytobands to one p- and one q-range per chromosome
cb_arms <- cb_frame %>%
  dplyr::mutate(chr = gsub("^chr", "", chr)) %>%
  dplyr::filter(chr %in% valid_chr) %>%
  dplyr::mutate(arm_letter = substr(band, 1, 1)) %>%
  dplyr::filter(arm_letter %in% c("p","q")) %>%
  dplyr::group_by(chr, arm_letter) %>%
  dplyr::summarise(start = min(start), end = max(end), .groups = "drop") %>%
  dplyr::mutate(arm = paste0(chr, arm_letter))

# Arm lengths (bp)
arm_lengths <- cb_arms %>%
  dplyr::mutate(arm_len = end - start + 1) %>%
  dplyr::rename(arm_start = start, arm_end = end) %>%
  dplyr::select(chr, arm, arm_start, arm_end, arm_len)

# Chromosome lengths (p+q)
chr_lengths <- cb_arms %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(chr_len = max(end) - min(start) + 1, .groups = "drop")

# Build GRanges
arms_gr <- GRanges(
  seqnames = cb_arms$chr,
  ranges   = IRanges(cb_arms$start, cb_arms$end),
  arm      = cb_arms$arm
)
seg_gr <- GRanges(
  seqnames = combined_seg_data$chr,
  ranges   = IRanges(combined_seg_data$start, combined_seg_data$end)
)

# Overlap and assign arm with maximum basepair overlap
hits <- findOverlaps(seg_gr, arms_gr, ignore.strand = TRUE)
if (length(hits) == 0L) stop("No overlaps between segments and cytoband arms. Check genome build / chr labels.")

ov_w <- width(pintersect(ranges(seg_gr)[queryHits(hits)],
                         ranges(arms_gr)[subjectHits(hits)]))

best_by_row <- tapply(
  X     = seq_along(ov_w),
  INDEX = queryHits(hits),
  FUN   = function(idx) {
    arm_ids <- subjectHits(hits)[idx]
    mcols(arms_gr)$arm[ arm_ids[ which.max(ov_w[idx]) ] ]
  }
)

arm_vec <- rep(NA_character_, nrow(combined_seg_data))
arm_vec[as.integer(names(best_by_row))] <- unname(unlist(best_by_row))

combined_seg_data <- combined_seg_data %>%
  mutate(arm = arm_vec) %>%
  relocate(arm, .after = end)

combined_seg_data <- combined_seg_data %>%
  left_join(arm_lengths, by = c("chr","arm")) %>%
  mutate(
    seg_len_in_arm = pmax(pmin(end, arm_end) - pmax(start, arm_start) + 1, 0),
    seg_len_in_arm = ifelse(is.na(seg_len_in_arm), 0, seg_len_in_arm)
  )



### Get FISH probe overlap
# Tunables
pad_bp  <- 150000L    # ±150 kb padding around each vendor span
min_bp  <- 1000L     # require at least 1 kb overlap 

# 1) Load probe table (if not already loaded)
if (!exists("fish_probe_locations")) {
  fish_probe_locations <- read_excel("Clinical data/FISH probe locations.xlsx")
}
stopifnot(all(c("Chromosome","Location in hg38") %in% names(fish_probe_locations)))

# 2) Map Chromosome -> feature and parse coordinates
feature_map <- c("1p" = "del1p", "1q" = "amp1q", "17p" = "del17p", "13q" = "del13q")

probe_windows <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    loc     = gsub(",", "", `Location in hg38`)
  ) %>%
  tidyr::extract(loc, into = c("chr","start","end"),
                 regex = "^chr([^:]+):(\\d+)-(\\d+)$", remove = TRUE) %>%
  mutate(
    chr   = gsub("^chr","", chr),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(!is.na(feature), !is.na(chr), !is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

# 3) Build GRanges: segments and padded probe windows
sample_cols_seg <- setdiff(
  names(combined_seg_data),
  c("chr","start","end","arm","arm_start","arm_end","arm_len","seg_len_in_arm")
)

seg_gr <- GRanges(
  seqnames = as.character(combined_seg_data$chr),
  ranges   = IRanges(as.integer(combined_seg_data$start),
                     as.integer(combined_seg_data$end))
)

# Ensure calls are uppercase characters
calls_df <- combined_seg_data[, sample_cols_seg] |>
  dplyr::mutate(across(everything(), ~ toupper(as.character(.))))

# If any duplicate sample names, make them unique (prevents subsetting errors)
if (anyDuplicated(names(calls_df))) {
  names(calls_df) <- make.unique(names(calls_df), sep = "_dup")
}

# Correct: one column per sample, no name mangling
mcols(seg_gr) <- S4Vectors::DataFrame(as.list(calls_df), check.names = FALSE)

# Use the actual names from mcols downstream
samples <- colnames(mcols(seg_gr))

probe_gr <- GRanges(
  seqnames = probe_windows$chr,
  ranges   = IRanges(probe_windows$start_pad, probe_windows$end_pad),
  feature  = probe_windows$feature
)

# 4) Overlap probes with segments; compute bp overlap
hits <- findOverlaps(probe_gr, seg_gr, ignore.strand = TRUE)
if (length(hits) == 0L) {
  warning("No overlaps between padded probe windows and segments. Check genome build / chr labels.")
}

ovl_bp <- width(pintersect(ranges(probe_gr)[queryHits(hits)],
                           ranges(seg_gr)[subjectHits(hits)]))

hits_df <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits)),
  seg_idx   = as.integer(subjectHits(hits)),
  ovl       = as.integer(ovl_bp)
) %>%
  dplyr::filter(ovl >= min_bp) %>%                # apply minimum overlap filter
  dplyr::arrange(probe_idx, dplyr::desc(ovl))


# 5) For each probe and each sample, take the first non-NA call among overlapping
#    segments ordered by decreasing overlap.
probe_ids <- seq_along(probe_gr)

final_calls_mat <- matrix(NA_character_, nrow = length(probe_ids), ncol = length(samples),
                          dimnames = list(NULL, samples))

for (p in probe_ids) {
  segs <- hits_df %>% dplyr::filter(probe_idx == p) %>% dplyr::pull(seg_idx)
  if (length(segs) == 0L) next
  # Subset the mcols matrix for these segments
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  # For each sample, pick first non-NA down the rows
  for (j in seq_along(samples)) {
    col_vals <- seg_calls[, j]
    nn <- which(!is.na(col_vals) & nzchar(col_vals))
    if (length(nn)) final_calls_mat[p, j] <- col_vals[nn[1]]
  }
}

# 6) Tidy to long and bin by expected direction
probe_meta <- tibble::tibble(
  probe_idx = probe_ids,
  feature   = as.character(mcols(probe_gr)$feature)
)

probe_calls_long <- as_tibble(final_calls_mat, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids) %>%
  tidyr::pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  left_join(probe_meta, by = "probe_idx") %>%
  select(Sample, feature, probe_call)

gain_labels <- c("GAIN","AMP","HLAMP")
loss_labels <- c("HOMD","HETD","LOSS","CNLOH")
dir_map     <- c(amp1q = "gain", del1p = "loss", del13q = "loss", del17p = "loss")

probe_calls_bin <- probe_calls_long %>%
  mutate(direction = unname(dir_map[feature]),
         is_altered_at_probe = dplyr::if_else(
           direction == "gain",
           as.integer(probe_call %in% gain_labels),
           as.integer(probe_call %in% loss_labels)
         )) %>%
  select(-direction) %>%
  tidyr::pivot_wider(names_from = feature,
                     values_from = c(probe_call, is_altered_at_probe),
                     names_sep = "_")
# =======================================================================

### Now redo, but using the cytoband instead 
# ---- Cytoband-based concordance ----

# Expected cytoband file format (UCSC): chrom  start  end  band  gieStain
# Example row:                          chr1   0      2300000  p36.33  gneg
cytoband_file <- "cytoband.txt"

cyto <- readr::read_tsv(
  cytoband_file,
  col_names = c("chrom", "start", "end", "band", "gieStain"),
  col_types = readr::cols(
    chrom = readr::col_character(),
    start = readr::col_double(),
    end = readr::col_double(),
    band = readr::col_character(),
    gieStain = readr::col_character()
  ),
  progress = FALSE
)

# Clean up chromosome column (remove "chr" prefix if present)
cyto <- cyto %>%
  mutate(
    chrom = gsub("^chr", "", chrom),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  arrange(chrom, start)


# Be forgiving about column names
nm <- names(cyto)
stopifnot(
  any(grepl("^chrom", nm, ignore.case = TRUE)),
  any(grepl("^start$", nm, ignore.case = TRUE)),
  any(grepl("^end$",   nm, ignore.case = TRUE)),
  any(grepl("^band",   nm, ignore.case = TRUE))
)

# Normalize column names we need
names(cyto)[grepl("^chrom", names(cyto), ignore.case = TRUE)] <- "chrom"
names(cyto)[grepl("^start$", names(cyto), ignore.case = TRUE)] <- "start"
names(cyto)[grepl("^end$",   names(cyto), ignore.case = TRUE)] <- "end"
names(cyto)[grepl("^band",   names(cyto), ignore.case = TRUE)]  <- "band"

cyto <- cyto %>%
  mutate(
    chrom = gsub("^chr", "", chrom),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  arrange(chrom, start)

# Helper: expand a band tag like "q21" or "q21.3" to the span covering all matching rows
.band_span <- function(cyto_chr, arm, band_tag) {
  tag <- paste0(arm, band_tag)            # e.g., "q21"
  rows <- which(startsWith(cyto_chr$band, tag))
  if (!length(rows)) return(NULL)
  tibble::tibble(
    start = min(cyto_chr$start[rows]),
    end   = max(cyto_chr$end[rows])
  )
}

# Parse Target like:
#   "17p13.1"            -> chr=17, arm1=p, band1=13.1 (single)
#   "1q21-q22"           -> chr=1,  arm1=q, band1=21 ; arm2=q, band2=22
#   "1q21.1-q22.2"       -> range with decimals
.parse_target <- function(s) {
  s <- gsub("\\s", "", s)
  m <- regexec("^([0-9XYM]+)([pq])([0-9]+(?:\\.[0-9]+)?)(?:-([pq])?([0-9]+(?:\\.[0-9]+)?))?$", s)
  g <- regmatches(s, m)[[1]]
  if (!length(g)) return(NULL)
  chr   <- g[2]
  arm1  <- g[3]; band1 <- g[4]
  arm2  <- g[5]; band2 <- g[6]
  if (is.na(arm2) || identical(arm2, "")) arm2 <- arm1
  if (is.na(band2) || identical(band2, "")) {
    list(chr = chr, arm1 = arm1, band1 = band1, arm2 = arm1, band2 = band1, is_range = FALSE)
  } else {
    list(chr = chr, arm1 = arm1, band1 = band1, arm2 = arm2, band2 = band2, is_range = TRUE)
  }
}

# Map each probe row to a cytoband-derived window
probe_windows_cytoband <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    Target  = Target
  ) %>%
  rowwise() %>%
  mutate(
    .pt = list(.parse_target(Target))
  ) %>%
  ungroup() %>%
  filter(!purrr::map_lgl(.pt, is.null)) %>%
  mutate(
    chr  = purrr::map_chr(.pt, ~ .x$chr),
    arm1 = purrr::map_chr(.pt, ~ .x$arm1),
    b1   = purrr::map_chr(.pt, ~ .x$band1),
    arm2 = purrr::map_chr(.pt, ~ .x$arm2),
    b2   = purrr::map_chr(.pt, ~ .x$band2)
  ) %>%
  select(-.pt) %>%
  group_by(feature, chr, arm1, b1, arm2, b2) %>%
  # Resolve to hg38 span using cytobands
  reframe({
    cy_chr <- dplyr::filter(cyto, chrom == chr) %>% arrange(start)
    sspan  <- .band_span(cy_chr, arm1, b1)
    espan  <- .band_span(cy_chr, arm2, b2)
    if (is.null(sspan) || is.null(espan)) {
      tibble::tibble(start = NA_integer_, end = NA_integer_)
    } else {
      tibble::tibble(
        start = min(sspan$start, espan$start),
        end   = max(sspan$end,   espan$end)
      )
    }
  }) %>%
  ungroup() %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

# Quick sanity check of the tricky one (1q21-q22): should span all q21.* through q22.* on chr1
# print(dplyr::filter(probe_windows_cytoband, grepl("^amp1q$", feature)))

# Build GRanges for cytoband windows
probe_gr_cyto <- GRanges(
  seqnames = probe_windows_cytoband$chr,
  ranges   = IRanges(probe_windows_cytoband$start_pad, probe_windows_cytoband$end_pad),
  feature  = probe_windows_cytoband$feature
)

# Overlaps and bp widths
hits_cyto <- findOverlaps(probe_gr_cyto, seg_gr, ignore.strand = TRUE)
if (length(hits_cyto) == 0L) {
  warning("No overlaps between cytoband windows and segments. Check build/labels.")
}

ovl_bp_cyto <- width(pintersect(
  ranges(probe_gr_cyto)[queryHits(hits_cyto)],
  ranges(seg_gr)[subjectHits(hits_cyto)]
))

hits_df_cyto <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits_cyto)),
  seg_idx   = as.integer(subjectHits(hits_cyto)),
  ovl       = as.integer(ovl_bp_cyto)
) %>%
  dplyr::filter(ovl >= min_bp) %>%
  dplyr::arrange(probe_idx, dplyr::desc(ovl))


# Senestivity/overlap chooser
probe_ids_cyto <- seq_along(probe_gr_cyto)
final_calls_mat_cyto <- matrix(
  NA_character_, nrow = length(probe_ids_cyto), ncol = length(samples),
  dimnames = list(NULL, samples)
)

# Severity ladders (edit if the labels differ)
gain_priority <- c("HLAMP","AMP","GAIN")
loss_priority <- c("HOMD","HETD","LOSS","CNLOH")

# Feature per cytoband probe
feature_vec_cyto <- as.character(mcols(probe_gr_cyto)$feature)

for (p in probe_ids_cyto) {
  seg_rows <- hits_df_cyto %>% dplyr::filter(probe_idx == p)
  segs     <- seg_rows$seg_idx
  if (!length(segs)) next
  
  # Segment calls for just these overlaps (rows = segments, cols = samples)
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  
  # Overlap widths aligned to 'segs'
  ovl_here <- seg_rows$ovl
  
  # Direction-specific priority for this probe
  feat <- feature_vec_cyto[p]
  dir  <- unname(dir_map[feat])   # "gain" or "loss" from the existing dir_map
  pri  <- if (identical(dir, "gain")) gain_priority else loss_priority
  
  # Precompute severity rank per cell (lower is more severe)
  sev_rank <- apply(seg_calls, 2, function(col) {
    match(toupper(col), pri, nomatch = length(pri) + 1L)
  })
  
  # Pick, per sample: most severe; tie-break by larger bp overlap
  for (j in seq_along(samples)) {
    calls_j <- seg_calls[, j]
    valid   <- which(!is.na(calls_j) & nzchar(calls_j))
    if (!length(valid)) next
    
    rnk  <- sev_rank[valid, j]
    best <- which(rnk == min(rnk))
    if (length(best) > 1L) {
      best <- best[which.max(ovl_here[valid][best])]
    }
    final_calls_mat_cyto[p, j] <- calls_j[valid][best]
  }
}

# Tidy long and bin calls by expected direction, mirroring probe_calls_bin
probe_meta_cyto <- tibble::tibble(
  probe_idx = probe_ids_cyto,
  feature   = as.character(mcols(probe_gr_cyto)$feature)
)

probe_calls_long_cyto <- as_tibble(final_calls_mat_cyto, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids_cyto) %>%
  tidyr::pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  dplyr::left_join(probe_meta_cyto, by = "probe_idx") %>%
  dplyr::select(Sample, feature, probe_call)

probe_calls_bin_cytoband <- probe_calls_long_cyto %>%
  mutate(
    direction = unname(dir_map[feature]),
    is_altered_at_probe = dplyr::if_else(
      direction == "gain",
      as.integer(probe_call %in% gain_labels),
      as.integer(probe_call %in% loss_labels)
    )
  ) %>%
  select(-direction) %>%
  tidyr::pivot_wider(
    names_from  = feature,
    values_from = c(probe_call, is_altered_at_probe),
    names_sep   = "_"
  )





# (Optional) quick sanity checks
# combined_seg_data %>% summarize(prop_assigned = mean(!is.na(arm)))
# combined_seg_data %>% filter(is.na(arm)) %>% count(chr, sort = TRUE)

# ---- Save combined Sequenza segment and probe-call intermediates ----
# These files support checks and reuse of the Sequenza layer; the final
# feature integration consumes the cleaned exports written near the end.
write.csv(combined_seg_data, support_file("Sep_2025_combined_sequenza_calls_400_updated.csv"), row.names = FALSE)
saveRDS(combined_seg_data, support_file("Sep_2025_combined_sequenza_calls_400_updated.rds"))

write.csv(probe_calls_bin, support_file("Sep_2025_FISH_probe_calls_400_updated.csv"), row.names = FALSE)
saveRDS(probe_calls_bin, support_file("Sep_2025_FISH_probe_calls_400_updated.rds"))
message("Support Sequenza CNA/probe caches written under: ", support_table_dir)

# ---- long format with segment length (for length-weighting)
sample_cols <- setdiff(
  names(combined_seg_data),
  c("chr","start","end","arm","arm_start","arm_end","arm_len","seg_len_in_arm")
)

long_data <- combined_seg_data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Value") %>%
  mutate(
    Value   = toupper(as.character(Value)),
    seg_len = seg_len_in_arm
  )

# ---- arms of interest
myeloma_arms_df <- tibble::tribble(
  ~Arm,     ~chr,  ~arm,
  "del1p",  "1",   "1p",
  "amp1q",  "1",   "1q",
  "del13q", "13",  "13q",
  "del17p", "17",  "17p"
)

gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD", "LOSS", "CNLOH")

samples <- unique(long_data$Sample)
results <- tibble(Sample = samples)

# ---- length-weighted % altered per arm
for (i in seq_len(nrow(myeloma_arms_df))) {
  arm_name  <- myeloma_arms_df$Arm[i]
  chr_value <- myeloma_arms_df$chr[i]
  arm_value <- myeloma_arms_df$arm[i]
  
  denom <- arm_lengths %>%
    filter(chr == chr_value, arm == arm_value) %>%
    pull(arm_len) %>% unique()
  if (length(denom) != 1) stop("Arm length lookup failed for ", arm_value)
  
  data_arm <- long_data %>%
    filter(chr == chr_value, arm == arm_value, !is.na(Value)) %>%  # key!
    mutate(IsAltered = case_when(
      startsWith(arm_name, "amp") ~ as.integer(Value %in% gain_labels),
      startsWith(arm_name, "del") ~ as.integer(Value %in% loss_labels),
      TRUE ~ 0L
    ))
  
  sample_props <- data_arm %>%
    group_by(Sample) %>%
    summarise(PropAltered = sum(IsAltered * seg_len, na.rm = TRUE) / denom,
              .groups = "drop")
  
  sample_props <- full_join(tibble(Sample = samples), sample_props, by = "Sample") %>%
    mutate(PropAltered = replace_na(PropAltered, 0))
  
  results <- results %>%
    left_join(sample_props, by = "Sample") %>%
    mutate(!!arm_name := as.integer(PropAltered > (1/3))) %>%
    select(-PropAltered)
}

# ---- Hyperdiploidy: length-weighted full-gain on 3,5,7,9,11,15,19,21
hyperdiploid_chrs <- c("3","5","7","9","11","15","19","21")

# Implemented hyperdiploidy parameters.
gain_thresh_chr <- 0.65   # chromosome is gained when >65% of its fixed length is gained
k_trisomies     <- 5      # hyperdiploid when at least 5 of the 8 chromosomes are gained

# compute PropGained per chr
gain_prop_tbl <- purrr::map_dfr(hyperdiploid_chrs, function(chr_value){
  # fixed denominator (whole chromosome length) to avoid NA-driven inflation
  denom_chr <- chr_lengths %>%
    filter(chr == chr_value) %>% pull(chr_len)
  if (length(denom_chr) != 1) stop("Chromosome length lookup failed for chr", chr_value)
  
  long_data %>%
    dplyr::filter(chr == chr_value, !is.na(Value)) %>%
    mutate(IsGain = as.integer(Value %in% gain_labels)) %>%
    group_by(Sample) %>%
    summarise(PropGained = sum(IsGain * seg_len, na.rm = TRUE) / denom_chr,
              .groups = "drop") %>%
    mutate(chr = chr_value)
}) %>%
  tidyr::pivot_wider(names_from = chr, values_from = PropGained,
                     names_prefix = "chr") %>%
  # ensure all samples are present
  right_join(tibble(Sample = unique(long_data$Sample)), by = "Sample") %>%
  mutate(across(starts_with("chr"), ~ tidyr::replace_na(.x, 0)))

# binary per-chr gains with relaxed threshold
hd_bin <- gain_prop_tbl %>%
  mutate(across(starts_with("chr"),
                ~ as.integer(.x > gain_thresh_chr),
                .names = "hyperdiploid_chr{substr(.col,4,6)}"))

# count trisomies and call hyperdiploid
hd_count <- hd_bin %>%
  mutate(hd_count = rowSums(select(., starts_with("hyperdiploid_chr")), na.rm = TRUE),
         hyperdiploid = as.integer(hd_count >= k_trisomies)) %>%
  select(Sample, starts_with("hyperdiploid_chr"), hd_count, hyperdiploid)

# merge into results (drop older strict cols if present)
results <- results %>%
  select(-dplyr::any_of(c("hyperdiploid")),
         -dplyr::starts_with("hyperdiploid_chr")) %>%
  left_join(hd_count, by = "Sample")

## Old way - undercalled with Sequenza
# for (chr_value in hyperdiploid_chrs) {
#   colname <- paste0("hyperdiploid_chr", chr_value)
#   
#   denom_chr <- chr_lengths %>%
#     filter(chr == chr_value) %>% pull(chr_len)
#   if (length(denom_chr) != 1) stop("Chromosome length lookup failed for chr", chr_value)
#   
#   data_chr <- long_data %>%
#     filter(chr == chr_value, !is.na(Value)) %>%          # key!
#     mutate(IsGain = as.integer(Value %in% gain_labels))
#   
#   sample_gain <- data_chr %>%
#     group_by(Sample) %>%
#     summarise(PropGained = sum(IsGain * seg_len, na.rm = TRUE) / denom_chr,
#               .groups = "drop")
#   
#   sample_gain <- full_join(tibble(Sample = samples), sample_gain, by = "Sample") %>%
#     mutate(PropGained = replace_na(PropGained, 0))
#   
#   results <- results %>%
#     left_join(sample_gain, by = "Sample") %>%
#     mutate(!!colname := as.integer(PropGained > 0.80)) %>%
#     select(-PropGained)
# # }
# 
# results <- results %>%
#   mutate(hyperdiploid = as.integer(rowSums(select(., starts_with("hyperdiploid_chr"))) == 8))

myeloma_CNA_matrix_with_HRD <- results

# ---- Build final Sequenza arm-level CNA table ----
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate(across(-Sample, as.character))

## Check if what have is expected
# Convert all alteration columns to numeric
cna_data_summary <- cna_data %>%
  mutate(across(-Sample, as.numeric)) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    prop_del1p       = mean(del1p == 1, na.rm = TRUE),
    prop_amp1q       = mean(amp1q == 1, na.rm = TRUE),
    prop_del13q      = mean(del13q == 1, na.rm = TRUE),
    prop_del17p      = mean(del17p == 1, na.rm = TRUE),
    prop_hyperdiploid= mean(hyperdiploid == 1, na.rm = TRUE)
  )

print(cna_data_summary)


### Edit the same so in same format as other tools expect 
# Add a Tumor_Sample_Barcode column to metada_df_mutation_comparison
# Clean up Sample names in CNA data to match metadata format
cna_data_cleaned <- cna_data %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

cna_data_cleaned <- cna_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join to metadata by the best available Sequenza sample key.
# Historical Sequenza inputs use BAM-derived sample names that match
# Tumor_Sample_Barcode after removing _PG/_WG. The Spring 2026 pipeline-suite
# export instead uses repo-style Sample_ID values such as "IMG-081-T0-OZ".
# Resolve both formats to Bam_clean_tmp, because downstream integration uses
# Bam_clean_tmp as the canonical CNA/translocation sample key.
sequenza_metadata_lookup <- metada_df_mutation_comparison %>%
  mutate(
    Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
    Sample_ID = as.character(Sample_ID),
    Bam_clean_tmp = as.character(Bam_clean_tmp)
  )

sequenza_tsb_lookup <- sequenza_metadata_lookup %>%
  filter(!is.na(Tumor_Sample_Barcode), Tumor_Sample_Barcode != "") %>%
  distinct(Tumor_Sample_Barcode, Bam_clean_tmp)

sequenza_sample_id_lookup <- sequenza_metadata_lookup %>%
  filter(!is.na(Sample_ID), Sample_ID != "", !is.na(Bam_clean_tmp), Bam_clean_tmp != "") %>%
  semi_join(
    cna_data_cleaned %>% distinct(Sample),
    by = c("Sample_ID" = "Sample")
  ) %>%
  distinct(Sample_ID, Bam_clean_tmp)

duplicated_sample_id_lookup <- sequenza_sample_id_lookup %>%
  group_by(Sample_ID) %>%
  summarise(n_bam_keys = n_distinct(Bam_clean_tmp), .groups = "drop") %>%
  filter(n_bam_keys > 1)

if (nrow(duplicated_sample_id_lookup) > 0L) {
  stop(
    "Cannot map Sequenza Sample_ID to Bam_clean_tmp; non-unique Sample_ID values: ",
    paste(duplicated_sample_id_lookup$Sample_ID, collapse = ", "),
    call. = FALSE
  )
}

cna_data_merged <- cna_data_cleaned %>%
  left_join(
    sequenza_tsb_lookup %>% dplyr::rename(Bam_clean_tmp_from_tsb = Bam_clean_tmp),
    by = c("Sample" = "Tumor_Sample_Barcode")
  ) %>%
  left_join(
    sequenza_sample_id_lookup %>% dplyr::rename(Bam_clean_tmp_from_sample_id = Bam_clean_tmp),
    by = c("Sample" = "Sample_ID")
  ) %>%
  mutate(
    Bam_clean_tmp = coalesce(Bam_clean_tmp_from_tsb, Bam_clean_tmp_from_sample_id)
  ) %>%
  select(-Bam_clean_tmp_from_tsb, -Bam_clean_tmp_from_sample_id)

unmapped_spring2026_sequenza <- cna_data_merged %>%
  filter(str_detect(Sample, "^IMG-"), is.na(Bam_clean_tmp))

if (nrow(unmapped_spring2026_sequenza) > 0L) {
  stop(
    "Spring 2026 Sequenza rows could not be mapped to Bam_clean_tmp: ",
    paste(unmapped_spring2026_sequenza$Sample, collapse = ", "),
    call. = FALSE
  )
}

# Check results
message("Merged Sequenza CNA data with metadata: added Bam_clean_tmp column.")

# Reuse the validated two-key Sequenza mapping (historical BAM-derived barcode
# first, revision Sample_ID second) for every downstream Sequenza export. This
# prevents revision samples from mapping in the CNA table but remaining blank in
# the FISH-probe and purity/ploidy helpers.
sequenza_sample_to_bam <- cna_data_merged %>%
  dplyr::select(Sample, Bam_clean_tmp) %>%
  dplyr::distinct()

sequenza_bam_metadata <- metada_df_mutation_comparison %>%
  dplyr::select(Bam_clean_tmp, Patient, Timepoint) %>%
  dplyr::filter(!is.na(Bam_clean_tmp), Bam_clean_tmp != "") %>%
  dplyr::distinct()

ambiguous_sequenza_bam_metadata <- sequenza_bam_metadata %>%
  dplyr::count(Bam_clean_tmp, name = "n_metadata_rows") %>%
  dplyr::filter(n_metadata_rows > 1L)
if (nrow(ambiguous_sequenza_bam_metadata) > 0L) {
  stop(
    "Sequenza BAM keys map to multiple patient/timepoint rows: ",
    paste(ambiguous_sequenza_bam_metadata$Bam_clean_tmp, collapse = ", "),
    call. = FALSE
  )
}


## ---- Export Sequenza CNA results ----
## The number in filenames corresponds to the gamma parameter used by Sequenza (gamma = 400).
## These are upstream helper files consumed by 1_5 and downstream concordance scripts.

# Primary CNA matrix
saveRDS(cna_data_merged, file = file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"))
write.table(cna_data_merged,
            file = file.path(export_dir, "cna_data_from_sequenza_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Active Sequenza CNA feature table written: ",
        file.path(export_dir, "cna_data_from_sequenza_400_updated.rds"))

# Summary statistics (proportion of samples per alteration). This is support/QC,
# not a downstream manuscript input.
saveRDS(cna_data_summary, file = support_file("cna_data_summary_400_updated.rds"))
write.table(cna_data_summary,
            file = support_file("cna_data_summary_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Console confirmation
message("Sequenza CNA export complete.")


### Export FISH-probe-level Sequenza calls based on direct hg38 coordinates
FISH_data_cleaned <- probe_calls_bin %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

FISH_data_cleaned <- FISH_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join through the validated Sequenza sample-to-BAM map.
FISH_data_cleaned <- FISH_data_cleaned %>%
  left_join(sequenza_sample_to_bam, by = "Sample") %>%
  left_join(
    sequenza_bam_metadata %>% dplyr::select(Bam_clean_tmp, Patient),
    by = "Bam_clean_tmp"
  )


# Export
saveRDS(FISH_data_cleaned, file = file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"))
write.table(FISH_data_cleaned,
            file = file.path(export_dir, "FISH_data_from_sequenza_400_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Active Sequenza FISH helper written: ",
        file.path(export_dir, "FISH_data_from_sequenza_400_updated.rds"))


## Export FISH-probe-level Sequenza calls based on cytoband-expanded windows
FISH_data_cleaned <- probe_calls_bin_cytoband %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

FISH_data_cleaned <- FISH_data_cleaned %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join through the validated Sequenza sample-to-BAM map.
FISH_data_cleaned <- FISH_data_cleaned %>%
  left_join(sequenza_sample_to_bam, by = "Sample") %>%
  left_join(
    sequenza_bam_metadata %>% dplyr::select(Bam_clean_tmp, Patient),
    by = "Bam_clean_tmp"
  )


# Export cytoband-expanded probe calls as support/QC. The active downstream
# integration script uses the direct-coordinate FISH helper above.
saveRDS(FISH_data_cleaned, file = support_file("FISH_data_from_sequenza_400_by_cytoband_updated.rds"))
write.table(FISH_data_cleaned,
            file = support_file("FISH_data_from_sequenza_400_by_cytoband_updated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


## Export Sequenza purity/ploidy estimates joined to metadata
sample_ploidy <- sample_ploidy %>%
  mutate(Sample = str_remove_all(Sample, "_PG|_WG"))

sample_ploidy <- sample_ploidy %>%
  mutate(Sample = ifelse(
    Sample == "TFRIM4_0189_Bm_P_ZC-02", 
    "TFRIM4_0189_Bm_P_ZC-02-01-O-DNA", 
    Sample  # Keep other values unchanged
  ))


# Join through the validated Sequenza sample-to-BAM map.
sample_ploidy <- sample_ploidy %>%
  left_join(sequenza_sample_to_bam, by = "Sample") %>%
  left_join(sequenza_bam_metadata, by = "Bam_clean_tmp")


## Export 
saveRDS(sample_ploidy, file = file.path(export_dir, "Sample_ploidy_from_sequenza_400.rds"))
write.table(sample_ploidy,
            file = file.path(export_dir, "Sample_ploidy_from_sequenza.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Active Sequenza ploidy helper written: ",
        file.path(export_dir, "Sample_ploidy_from_sequenza_400.rds"))




## Optional cross-check against downstream clinical/FISH summaries
# `filled_df`, `cohort_df`, and `fish_summary` are created in later scripts or
# interactive review sessions. Keep this block non-blocking so the command-line
# CNA processing step can complete from a fresh R session.
if (exists("filled_df") && exists("cohort_df") && exists("fish_summary")) {
  fish_summary2 <- filled_df %>%
    dplyr::select(Patient, DEL_17P, DEL_1P, AMP_1Q) %>%
    dplyr::left_join(cohort_df, by = "Patient") %>%
    dplyr::filter(!is.na(Cohort)) %>%
    dplyr::filter(
      DEL_17P == "Positive" |
        DEL_1P  == "Positive" |
        AMP_1Q  == "Positive"
    )
  
  fish_summary2 <- fish_summary %>%
    left_join(FISH_data_cleaned)
} else {
  message("Skipping optional FISH cross-check: filled_df, cohort_df, or fish_summary is not available in this session.")
}
