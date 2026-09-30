# =============================================================================
# Script: 1_4_Process_CNA_Data.R
#
# Description:
#   This script imports, harmonizes, and summarizes copy-number segment (SEG)
#   calls from ichorCNA outputs, mapping them to chromosomal arms of interest
#   and generating a binary matrix of arm-level alterations for myeloma analyses.
#   Steps include:
#     1. Recursively list and read all “*.seg” files in seg_dir
#        – Extracts chr, start, end, and the sample’s “Corrected_Call” column
#        – Renames each call column to the sample name (basename without .seg)
#     2. Merges every SEG table by (chr, start, end) via full joins into
#        `combined_seg_data`
#     3. Loads cytoband definitions (cb_frame) and derives an “arm” label
#        (e.g. “1p”, “1q”, “13q”, “17p”) by overlapping each segment
#     4. Saves the wide segment cache under
#        Output_tables_2025/ichor_cna_processing_support/
#     5. Pivots to long format and uppercases calls for consistency
#     6. Defines a list of myeloma-relevant arms (del1p, amp1q,
#        del13q, del17p) and for each:
#        – Flags segments as altered if “Corrected_Call” matches gain/loss labels
#        – Computes the proportion of altered bins per sample
#        – Sets a binary indicator (1 if >33% of that arm is altered)
#     7. Adds hyperdiploidy status by checking chromosome-level gains (>65%)
#        across chromosomes 3,5,7,9,11,15,19,21 and requiring >=5/8 gains
#     8. Outputs `myeloma_CNA_matrix_with_HRD`, a Sample × feature binary table
#
# Goal:
#   Convert ichorCNA segment-level copy-number calls into compact sample-level
#   CNA feature tables for downstream WGS feature integration, baseline genomic
#   landscape summaries, BM-vs-cfDNA concordance analyses, FISH-vs-WGS
#   concordance, and Supplementary Table 2 support.
#
# How to run from the repository root:
#   Rscript 1_4_Process_CNA_Data.R
#
# Inputs:
#   • Scripts_2025/Final_Scripts/helpers.R
#   • seg_dir          = "Oct 2024 data/Ichor_CNA"          # ichorCNA *.seg files
#   • Output_tables_2025/ichor_cna_processing_support/
#       Oct_2024_combined_corrected_calls.rds                 # fallback cache
#   • combined_clinical_data_updated_April2025.csv
#     This table is currently loaded to construct metadata sample keys, but
#     those keys are not subsequently joined into the exported CNA tables.
#   • Clinical data/FISH probe locations.xlsx
#   • cytoband.txt      = UCSC-style hg38 cytoband reference
#   • Optional Spring 2026 ichorCNA SEG files and metadata under
#       Data_Spring_2026_Revisions/
#
# Analysis units:
#   The raw/cache layer contains one genomic segment per row and one sample per
#   call column. Final feature tables contain one row per sample. Arm and whole-
#   chromosome proportions are fractions of evaluated ichorCNA segments; they
#   are not fractions of genomic base pairs or abnormal cells.
#
# Active downstream outputs:
#   • Jan2025_exported_data/FISH_probe_calls_bin_cytoband_ichorCNA.rds
#   • Jan2025_exported_data/FISH_probe_calls_bin_cytoband_ichorCNA.txt
#   • Jan2025_exported_data/cna_data_ichorCNA.rds
#   • Jan2025_exported_data/cna_data_ichorCNA.txt
#
# Support/QC outputs:
#   • Output_tables_2025/ichor_cna_processing_support/Oct_2024_combined_corrected_calls.rds
#   • Output_tables_2025/ichor_cna_processing_support/Oct_2024_combined_corrected_calls.csv
#   • Output_tables_2025/ichor_cna_processing_support/cna_data_compact_backup.txt
#   • Output_tables_2025/ichor_cna_processing_support/cna_data_hyperdiploid_chromosome_qc.txt
#   • Output_tables_2025/ichor_cna_processing_support/cna_arm_call_qc.tsv
#   • R object: myeloma_CNA_matrix_with_HRD (binary arm-level matrix)
#
# Key assumptions and audit notes:
#   • Arm-level del/gain calls use >1/3 of evaluated segments altered for del1p, amp1q,
#     del13q, and del17p.
#   • Hyperdiploidy uses >65% gained segments per canonical chromosome and
#     calls a sample hyperdiploid when >=5/8 canonical chromosomes are gained.
#   • FISH-probe calls use cytoband-derived probe windows with +/-150 kb padding
#     and a nearest-segment fallback up to 10 Mb when no segment overlaps.
#     Only loci present in the workbook are exported; the current workbook has
#     1p, 1q, and 17p probes, while del13q is available only in the arm-level
#     CNA table.
#   • Arm-level calls preserve samples with no evaluated segments as NA. In
#     contrast, the current probe-level `%in%` binning and hyperdiploidy logic
#     treat missing segment calls as not altered/not gained (binary 0).
#   • These thresholds are manuscript logic. Do not change them without
#     explicit scientific review.
#
# Dependencies:
#   library(tidyverse); library(purrr)
#   library(GenomicRanges); library(IRanges); library(S4Vectors); library(readxl)
#
# Usage from an interactive R session opened at the repository root:
#   source("1_4_Process_CNA_Data.R")
#   # creates combined_seg_data and myeloma_CNA_matrix_with_HRD in memory
#
# Manuscript outputs created/updated:
#   - None directly. This upstream script processes ichorCNA CNA calls for WGS
#     feature integration, baseline heatmaps, FISH/WGS concordance, and
#     subclonal-evolution support outputs.
#
# Author: Dory Abelman
# Date:   2025-05-26
# =============================================================================
# Pipeline status:
#   Active upstream dependency. This script does not directly create a named
#   final manuscript figure/table, but downstream scripts depend on its cleaned
#   outputs for figure, table, or model generation.
#
# Reproducibility note:
#   The preferred end-to-end path starts from raw ichorCNA *.seg files in
#   Oct 2024 data/Ichor_CNA/. If those raw segment files are not distributed in a
#   lightweight review bundle, the script can use the preserved combined segment
#   cache in Output_tables_2025/ichor_cna_processing_support/ to regenerate the
#   same downstream CNA feature tables without requiring RStudio state.
#


# Load Libraries:
library(tidyverse)       # dplyr, purrr, tidyr, readr, etc.
library(purrr)           # for reduce()
library(tidyr)           # for pivot_longer()
library(GenomicRanges)
library(IRanges)
library(S4Vectors)

.helpers_path <- file.path("Scripts_2025", "Final_Scripts", "helpers.R")
if (!file.exists(.helpers_path)) {
  .helpers_path <- "helpers.R"
}
source(.helpers_path)
rm(.helpers_path)

export_dir <- "Jan2025_exported_data"
support_table_dir <- file.path("Output_tables_2025", "ichor_cna_processing_support")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(support_table_dir, recursive = TRUE, showWarnings = FALSE)

support_file <- function(filename) {
  file.path(support_table_dir, filename)
}

# Normalize chromosome identifiers before any cross-sample join. Legacy SEG
# files use identifiers such as "1", whereas Spring 2026 revision SEG files use
# "chr1". Keeping both forms creates separate genomic keys and can silently
# turn valid revision calls into missing values downstream.
normalize_chromosome_id <- function(x) {
  normalized <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_remove(stringr::regex("^chr", ignore_case = TRUE)) %>%
    toupper()
  normalized[normalized == "M"] <- "MT"
  normalized[normalized == ""] <- NA_character_
  normalized
}

first_existing_file <- function(paths, description) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop(
      "Missing required ", description, ". Checked:\n  ",
      paste(paths, collapse = "\n  "),
      call. = FALSE
    )
  }
  hit
}


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



#### Now get the ichorCNA data loaded 
# ichorCNA writes one *.cna.seg file per sample containing genome-wide
# copy-number segments. Each row is a contiguous segment with a
# Corrected_Call label (GAIN/AMP/HETD/HOMD/NEUT) derived from the
# ichorCNA HMM after correcting for GC bias, mappability, and the
# estimated tumor fraction in that sample.
# Set the directory containing the SEG files
seg_dir <- "Oct 2024 data/Ichor_CNA"

# Get a list of all SEG files in the directory. The explicit "\\.seg$" pattern
# avoids accidentally reading non-segment files that happen to contain "seg" in
# their names.
seg_files <- if (dir.exists(seg_dir)) {
  list.files(seg_dir, pattern = "\\.seg$", full.names = TRUE)
} else {
  character()
}
spring2026_seg_files <- spring2026_revision_files(
  # Optional Spring 2026 ichorCNA segment files. These extend the main patient
  # CNA table when revision CNA outputs are present.
  "iChorCNA_All_CNA_Seg",
  "[.]seg$"
)
# M4CHIP files are dilution-series controls/samples and should not enter the
# main patient CNA feature matrix.
spring2026_seg_files <- spring2026_seg_files[!grepl("^M4CHIP_", basename(spring2026_seg_files))]
seg_files <- unique(c(seg_files, spring2026_seg_files))

# A small number of Spring 2026 files share a sample basename with a legacy
# copy. A many-file join cannot safely retain two columns with the same sample
# name, so choose one deterministically. Because revision files are appended
# after legacy files above, retaining the last occurrence gives precedence to
# the current Spring 2026 source.
seg_sample_names <- gsub("\\.cna\\.seg$", "", basename(seg_files))
duplicate_seg_sample_names <- unique(
  seg_sample_names[duplicated(seg_sample_names) | duplicated(seg_sample_names, fromLast = TRUE)]
)
if (length(duplicate_seg_sample_names) > 0L) {
  message(
    "Duplicate ichorCNA sample basenames detected; retaining the current ",
    "Spring 2026 occurrence for: ",
    paste(duplicate_seg_sample_names, collapse = ", ")
  )
  keep_seg_file <- !duplicated(seg_sample_names, fromLast = TRUE)
  seg_files <- seg_files[keep_seg_file]
  seg_sample_names <- seg_sample_names[keep_seg_file]
}

combined_seg_cache_rds <- support_file("Oct_2024_combined_corrected_calls.rds")
combined_seg_cache_csv <- support_file("Oct_2024_combined_corrected_calls.csv")

if (length(seg_files) > 0L) {
  message("Reading ", length(seg_files), " raw ichorCNA SEG files from ", seg_dir)
  
  # Loop over each SEG file and keep only the genomic interval and corrected
  # copy-number state. The result is one sample-specific call column per file.
  seg_data_list <- lapply(seg_files, function(file) {
    seg_data <- read.delim(file, header = TRUE)
    required_cols <- c("chr", "start", "end")
    if (!all(required_cols %in% names(seg_data))) {
      stop("SEG file is missing required columns chr/start/end: ", file, call. = FALSE)
    }
    corrected_call_cols <- grep("Corrected_Call$", names(seg_data), value = TRUE)
    if (length(corrected_call_cols) != 1L) {
      stop(
        "Expected exactly one Corrected_Call column in ", file,
        "; found ", length(corrected_call_cols),
        call. = FALSE
      )
    }
    
    sample_name <- gsub("\\.cna\\.seg$", "", basename(file))
    seg_corrected_call <- seg_data %>%
      dplyr::select(chr, start, end, Corrected_Call = all_of(corrected_call_cols)) %>%
      dplyr::mutate(chr = normalize_chromosome_id(chr))
    colnames(seg_corrected_call)[4] <- sample_name
    seg_corrected_call
  })
  
  # A full_join retains every unique genomic segment present in any sample.
  # Samples without an exactly matching segment receive NA for that genomic
  # interval. Arm-level summaries preserve those NAs; the probe-level and
  # hyperdiploidy sections below currently convert them to binary 0.
  combined_seg_data <- purrr::reduce(seg_data_list, full_join, by = c("chr", "start", "end"))
  
  ## Add the arm info to the directory
  # Cytoband definitions supply hg38 chromosomal band coordinates. The arm label
  # (for example, "1p" or "13q") is assigned by overlap, enabling arm-level
  # summarisation downstream.
  cytoband_txt <- "cytoband.txt"
  if (!file.exists(cytoband_txt)) {
    stop("Missing required hg38 cytoband file: ", cytoband_txt,
         ". This file is needed to map segment-level CNA calls to chromosome arms.")
  }
  cb_frame <- readr::read_tsv(
    cytoband_txt,
    col_names = c("chr", "start", "end", "band", "stain"),
    show_col_types = FALSE
  )
  
  cb_frame_arm <- cb_frame %>%
    mutate(arm = paste0(gsub("chr", "", chr), substr(band, 1, 1))) %>%
    dplyr::select(chr, start, end, arm) %>%
    mutate(chr = normalize_chromosome_id(chr))
  
  check_overlap <- function(chr, start1, end1, cb_frame_arm) {
    cb_matches <- cb_frame_arm %>%
      filter(chr == !!chr & start <= !!end1 & end >= !!start1)
    
    if (nrow(cb_matches) > 0) {
      cb_matches$arm[1]
    } else {
      NA_character_
    }
  }
  
  combined_seg_data <- combined_seg_data %>%
    rowwise() %>%
    mutate(arm = check_overlap(chr, start, end, cb_frame_arm)) %>%
    ungroup() %>%
    relocate(arm, .after = end)
  
  write.csv(combined_seg_data, file = combined_seg_cache_csv, row.names = FALSE)
  saveRDS(combined_seg_data, file = combined_seg_cache_rds)
  message("Support segment cache written: ", combined_seg_cache_rds)
} else {
  # Review bundles and Code Ocean capsules may preserve derived intermediates
  # rather than all raw ichorCNA SEG files. This fallback is explicit: it uses
  # only the combined segment table produced by the raw-data branch above.
  combined_seg_cache_rds <- first_existing_file(
    c(
      combined_seg_cache_rds,
      "Oct_2024_combined_corrected_calls.rds"
    ),
    "combined ichorCNA segment cache"
  )
  message("No raw ichorCNA SEG files found in ", seg_dir,
          "; using preserved support cache: ", combined_seg_cache_rds)
  combined_seg_data <- readRDS(combined_seg_cache_rds)
}

# Apply the same normalization to preserved caches and fail loudly if a future
# input reintroduces mixed chromosome namespaces.
combined_seg_data <- combined_seg_data %>%
  mutate(chr = normalize_chromosome_id(chr))
if (any(is.na(combined_seg_data$chr))) {
  stop("Combined ichorCNA segment data contain missing chromosome identifiers.",
       call. = FALSE)
}
if (any(grepl("^chr", combined_seg_data$chr, ignore.case = TRUE))) {
  stop("Chromosome normalization failed: 'chr'-prefixed identifiers remain.",
       call. = FALSE)
}

# Regression guard for the exact Spring 2026 failure mode: every supplied
# revision sample must retain at least one non-missing corrected call after the
# cross-sample join. A zero count indicates an identifier/join mismatch, not a
# biologically CNA-negative sample.
if (length(spring2026_seg_files) > 0L) {
  spring2026_sample_names <- gsub(
    "\\.cna\\.seg$", "", basename(spring2026_seg_files)
  )
  missing_revision_columns <- setdiff(
    spring2026_sample_names, names(combined_seg_data)
  )
  if (length(missing_revision_columns) > 0L) {
    stop(
      "Revision ichorCNA samples missing after segment join: ",
      paste(missing_revision_columns, collapse = ", "),
      call. = FALSE
    )
  }
  revision_nonmissing_counts <- vapply(
    combined_seg_data[spring2026_sample_names],
    function(x) sum(!is.na(x)),
    integer(1)
  )
  if (any(revision_nonmissing_counts == 0L)) {
    stop(
      "Revision ichorCNA samples have zero retained segment calls after join: ",
      paste(names(revision_nonmissing_counts)[revision_nonmissing_counts == 0L],
            collapse = ", "),
      call. = FALSE
    )
  }
}



## Build ichorCNA calls at clinical FISH probe loci
# =====================================================================
# Purpose: Recompute per-sample CNA calls at FISH probe cytobands
#          using ichorCNA segments, independent of other scripts.
# Output:  probe_calls_bin_cytoband (in memory) + optional exports
# =====================================================================

fish_probe_xlsx    <- "Clinical data/FISH probe locations.xlsx" # must contain Chromosome, Target, Location in hg38
# EITHER provide cb_frame as an RDS with columns chr,start,end,band
cb_frame_rds       <- NULL  # e.g., "cb_frame_hg38.rds"  # set to a path OR leave NULL
# OR provide a UCSC cytoband txt (no header), which we'll parse
cytoband_txt       <- "cytoband.txt"  # fallback if cb_frame_rds is NULL



## ---- Tunables & label maps ----
# pad_bp: ichorCNA segments do not always extend to the exact cytogenomic
#   band boundaries annotated for FISH probes. ±150 kb padding allows
#   segments that end just outside a probe band to still be matched.
# min_bp: a 1 kb minimum overlap prevents false matches from rounding
#   at segment endpoints.
# max_nearest_bp: ichorCNA may leave gaps (e.g. centromere, low-mappability
#   regions). If no segment overlaps a FISH probe window, the nearest
#   segment within 10 Mb is used as a fallback; beyond that distance the
#   raw probe_call is left uncalled (NA). Under the implemented `%in%` binning
#   below, however, its binary is_altered_at_probe flag becomes 0.
pad_bp  <- 150000L  # ±150 kb padding around band windows
min_bp  <- 1000L    # require >= 1 kb overlap to count
max_nearest_bp  <- 10e6L     # allow nearest-segment fallback within 10 Mb 

# FISH feature mapping by arm label in fish table
# These four chromosomal loci are the canonical high-risk FISH panel in MM:
#   del1p (1p32 CDKN2C/FAF1 locus) – adverse prognosis;
#   amp1q (1q21 CKS1B locus)       – adverse prognosis;
#   del17p (17p13 TP53 locus)       – highest-risk feature;
#   del13q (13q14 RB1 locus)        – high-risk, frequently monosomal.
feature_map <- c("1p" = "del1p", "1q" = "amp1q", "17p" = "del17p", "13q" = "del13q")
# Direction map for binning
# amp1q requires a gain event; all others are defined by copy loss.
dir_map     <- c(amp1q = "gain", del1p = "loss", del13q = "loss", del17p = "loss")

# Call label sets
# Gain labels: GAIN (3 copies), AMP (>=4), HLAMP (focal >=6+);
#   all count as amplification for amp1q scoring.
# Loss labels: HETD (1 copy) and HOMD (0 copies);
#   LOSS and CNLOH added for Sequenza-derived calls (see 1_4A).
gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD")

if (!file.exists(fish_probe_xlsx)) stop("Missing: ", fish_probe_xlsx)
fish_probe_locations <- readxl::read_excel(fish_probe_xlsx)

# Cytobands: prefer cb_frame RDS if provided; otherwise parse UCSC txt
if (!is.null(cb_frame_rds)) {
  if (!file.exists(cb_frame_rds)) stop("Missing: ", cb_frame_rds)
  cb_frame <- readRDS(cb_frame_rds)
  cyto <- cb_frame %>%
    transmute(
      chrom = gsub("^chr","", as.character(chr)),
      start = as.integer(start),
      end   = as.integer(end),
      band  = as.character(band)
    ) %>% arrange(chrom, start)
} else {
  if (!file.exists(cytoband_txt)) stop("Missing: ", cytoband_txt)
  cyto <- readr::read_tsv(
    cytoband_txt,
    col_names = c("chrom", "start", "end", "band", "gieStain"),
    col_types = cols(
      chrom = col_character(),
      start = col_double(),
      end   = col_double(),
      band  = col_character(),
      gieStain = col_character()
    ),
    progress = FALSE
  ) %>%
    mutate(
      chrom = gsub("^chr","", chrom),
      start = as.integer(start),
      end   = as.integer(end)
    ) %>%
    arrange(chrom, start)
}
stopifnot(all(c("chrom","start","end","band") %in% names(cyto)))

## ---- Build GRanges from ichorCNA segments (wide) ----
# Expect first 4 columns: chr, start, end, arm (ichorCNA pipeline output)
need_cols <- c("chr","start","end","arm")
if (!all(need_cols %in% names(combined_seg_data))) {
  stop("combined_seg_data must have columns: ", paste(need_cols, collapse = ", "))
}

sample_cols_seg <- setdiff(names(combined_seg_data), need_cols)
if (length(sample_cols_seg) == 0L) stop("No sample columns found in combined_seg_data")

# GRanges of all merged segments
seg_gr <- GRanges(
  seqnames = as.character(combined_seg_data$chr),
  ranges   = IRanges(as.integer(combined_seg_data$start),
                     as.integer(combined_seg_data$end))
)

# Attach calls (one column per sample) to mcols; uppercase for consistency
calls_df <- combined_seg_data[, sample_cols_seg] %>%
  mutate(across(everything(), ~ toupper(as.character(.))))
if (anyDuplicated(names(calls_df))) {
  names(calls_df) <- make.unique(names(calls_df), sep = "_dup")
}
mcols(seg_gr) <- S4Vectors::DataFrame(as.list(calls_df), check.names = FALSE)
samples <- colnames(mcols(seg_gr))

## ---- Helpers: band parsing & expansion ----
.band_span <- function(cyto_chr, arm, band_tag) {
  tag <- paste0(arm, band_tag)  # e.g., "q21", "p13.1"
  rows <- which(startsWith(cyto_chr$band, tag))
  if (!length(rows)) return(NULL)
  tibble(start = min(cyto_chr$start[rows]),
         end   = max(cyto_chr$end[rows]))
}

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

## ---- Turn FISH Targets into genomic windows via cytobands ----
stopifnot(all(c("Chromosome","Target") %in% names(fish_probe_locations)))

probe_windows_cytoband <- fish_probe_locations %>%
  transmute(
    feature = dplyr::recode(Chromosome, !!!feature_map, .default = NA_character_),
    Target  = Target
  ) %>%
  rowwise() %>%
  mutate(.pt = list(.parse_target(Target))) %>%
  ungroup() %>%
  filter(!purrr::map_lgl(.pt, is.null), !is.na(feature)) %>%
  mutate(
    chr  = purrr::map_chr(.pt, ~ .x$chr),
    arm1 = purrr::map_chr(.pt, ~ .x$arm1),
    b1   = purrr::map_chr(.pt, ~ .x$band1),
    arm2 = purrr::map_chr(.pt, ~ .x$arm2),
    b2   = purrr::map_chr(.pt, ~ .x$band2)
  ) %>%
  select(-.pt) %>%
  group_by(feature, chr, arm1, b1, arm2, b2) %>%
  reframe({
    cy_chr <- dplyr::filter(cyto, chrom == chr) %>% arrange(start)
    sspan  <- .band_span(cy_chr, arm1, b1)
    espan  <- .band_span(cy_chr, arm2, b2)
    if (is.null(sspan) || is.null(espan)) tibble(start = NA_integer_, end = NA_integer_) else
      tibble(start = min(sspan$start, espan$start),
             end   = max(sspan$end,   espan$end))
  }) %>%
  ungroup() %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(
    start_pad = pmax(1L, start - pad_bp),
    end_pad   = end + pad_bp
  ) %>%
  distinct(feature, chr, start_pad, end_pad, .keep_all = TRUE)

if (!nrow(probe_windows_cytoband)) stop("No valid probe windows built from Target field")

## ---- Overlap probe windows with segments ----
probe_gr_cyto <- GRanges(
  seqnames = probe_windows_cytoband$chr,
  ranges   = IRanges(probe_windows_cytoband$start_pad, probe_windows_cytoband$end_pad),
  feature  = probe_windows_cytoband$feature
)

hits_cyto <- findOverlaps(probe_gr_cyto, seg_gr, ignore.strand = TRUE)

ovl_bp_cyto <- if (length(hits_cyto)) width(pintersect(
  ranges(probe_gr_cyto)[queryHits(hits_cyto)],
  ranges(seg_gr)[subjectHits(hits_cyto)]
)) else integer(0)

hits_df_cyto <- tibble::tibble(
  probe_idx = as.integer(queryHits(hits_cyto)),
  seg_idx   = as.integer(subjectHits(hits_cyto)),
  ovl       = as.integer(ovl_bp_cyto),
  src       = "overlap"
) %>%
  dplyr::filter(ovl >= min_bp)

# --- Fallback: nearest segment per probe when no overlaps made it, since ichorCNA skips problematic regions of the genome# ichorCNA uses 1 Mb bins; centromeric regions, segmental duplications,
# and low-mappability windows may be entirely absent from the seg file.
# When a FISH probe window has no overlapping segment, distanceToNearest
# finds the closest segment in bp. Only segments within max_nearest_bp
# (10 Mb) are used; probes beyond that threshold remain NA so that
# centromere-proximal probes are not assigned calls from the wrong arm.# distanceToNearest returns the *closest* segment (ties broken arbitrarily), with distance in bp.
nearest_hits <- GenomicRanges::distanceToNearest(probe_gr_cyto, seg_gr, ignore.strand = TRUE)

nearest_df <- tibble::tibble(
  probe_idx = as.integer(S4Vectors::queryHits(nearest_hits)),
  seg_idx   = as.integer(S4Vectors::subjectHits(nearest_hits)),
  dist_bp   = as.integer(S4Vectors::mcols(nearest_hits)$distance)
) %>%
  # keep only probes that have no overlapping rows already
  dplyr::anti_join(hits_df_cyto %>% dplyr::select(probe_idx) %>% dplyr::distinct(),
                   by = "probe_idx") %>%
  # enforce a maximum distance cap to avoid silly jumps
  dplyr::filter(dist_bp <= max_nearest_bp) %>%
  # encode "overlap score" as negative distance so sorting still prefers overlaps
  dplyr::transmute(
    probe_idx,
    seg_idx,
    ovl = -dist_bp,   # negative = farther; less negative = closer
    src = "nearest"
  )

# Combine, order by probe then "best support" (overlap first, then nearest by distance)
hits_df_cyto <- dplyr::bind_rows(hits_df_cyto, nearest_df) %>%
  dplyr::arrange(probe_idx, dplyr::desc(ovl))

## ---- Choose per-sample call by severity first, then overlap ----
probe_ids_cyto <- seq_along(probe_gr_cyto)
final_calls_mat_cyto <- matrix(
  NA_character_, nrow = length(probe_ids_cyto), ncol = length(samples),
  dimnames = list(NULL, samples)
)

# Severity-ranked label vectors: when multiple ichorCNA segments
# overlap a FISH probe window, the call with the highest clinical
# severity is chosen rather than the one with the greatest overlap.
# For gains:  HLAMP > AMP > GAIN (focal/high-level takes precedence).
# For losses: HOMD > HETD > LOSS > CNLOH (homozygous deletion is worst).
# This ensures that focal high-risk events are not obscured by adjacent
# neutral segments that happen to cover more base pairs.
gain_priority <- c("HLAMP","AMP","GAIN")
loss_priority <- c("HOMD","HETD","LOSS","CNLOH")

feature_vec_cyto <- as.character(mcols(probe_gr_cyto)$feature)

for (p in probe_ids_cyto) {
  seg_rows <- hits_df_cyto %>% filter(probe_idx == p)
  segs     <- seg_rows$seg_idx
  if (!length(segs)) next
  
  seg_calls <- as.matrix(mcols(seg_gr)[segs, samples, drop = FALSE])
  ovl_here  <- seg_rows$ovl
  
  feat <- feature_vec_cyto[p]
  dir  <- unname(dir_map[feat])   # "gain" or "loss"
  pri  <- if (identical(dir, "gain")) gain_priority else loss_priority
  
  lab <- toupper(seg_calls)  # same shape as seg_calls
  sev_rank <- array(
    match(lab, pri, nomatch = length(pri) + 1L),  # unknown => worst rank
    dim       = dim(lab),
    dimnames  = dimnames(lab)
  )
  
  for (j in seq_along(samples)) {
    calls_j <- seg_calls[, j]
    valid   <- which(!is.na(calls_j) & nzchar(calls_j))
    if (!length(valid)) next
    
    rnk  <- sev_rank[valid, j]
    best <- which(rnk == min(rnk))
    if (length(best) > 1L) best <- best[which.max(ovl_here[valid][best])]
    final_calls_mat_cyto[p, j] <- calls_j[valid][best]
  }
}

## ---- Tidy & bin by expected direction ----
probe_meta_cyto <- tibble(
  probe_idx = probe_ids_cyto,
  feature   = as.character(mcols(probe_gr_cyto)$feature)
)

probe_calls_long_cyto <- as_tibble(final_calls_mat_cyto, .name_repair = "minimal") %>%
  mutate(probe_idx = probe_ids_cyto) %>%
  pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "probe_call") %>%
  left_join(probe_meta_cyto, by = "probe_idx") %>%
  select(Sample, feature, probe_call)

probe_calls_bin_cytoband <- probe_calls_long_cyto %>%
  mutate(
    direction = unname(dir_map[feature]),
    is_altered_at_probe = if_else(
      direction == "gain",
      as.integer(probe_call %in% gain_labels),
      as.integer(probe_call %in% loss_labels)
    )
  ) %>%
  select(-direction) %>%
  pivot_wider(
    names_from  = feature,
    values_from = c(probe_call, is_altered_at_probe),
    names_sep   = "_"
  )

# View head
print(head(probe_calls_bin_cytoband, 10))

## ---- Export FISH-probe-level ichorCNA helper table ----
# This helper table is used by downstream concordance logic. It is not itself a
# final manuscript table, but it supports Extended Data Figure 2C and
# Supplementary Table 2 through the later 2_3 concordance script.
saveRDS(probe_calls_bin_cytoband, file = file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"))
write.table(probe_calls_bin_cytoband,
            file = file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Active CNA/FISH helper written: ",
        file.path(export_dir, "FISH_probe_calls_bin_cytoband_ichorCNA.rds"))


### Now continue with regular data transformations

# Pivot the data (assuming combined_seg_data is already defined)
long_data <- combined_seg_data %>%
  pivot_longer(cols = -(1:4), names_to = "Sample", values_to = "Value") %>%
  mutate(Value = toupper(Value))


# Define the list of chromosomal arms with corrected chr values (no "chr" prefix)
myeloma_cna_arms <- list(
  del1p = list(chr = "1", arm = "1p"),    # Deletion in 1p
  amp1q = list(chr = "1", arm = "1q"),      # Amplification in 1q
  del13q = list(chr = "13", arm = "13q"),   # Deletion in 13q
  del17p = list(chr = "17", arm = "17p")     # Deletion in 17p
)

# Flatten the list into a data frame
myeloma_arms_df <- data.frame(
  Arm = names(myeloma_cna_arms),
  chr = sapply(myeloma_cna_arms, function(x) x$chr),
  arm = sapply(myeloma_cna_arms, function(x) x$arm),
  stringsAsFactors = FALSE
)

# Define gain and loss labels
gain_labels <- c("GAIN", "AMP", "HLAMP")
loss_labels <- c("HOMD", "HETD")

# Get the list of samples
samples <- unique(long_data$Sample)

# Initialize the results data frame
results <- data.frame(Sample = samples, stringsAsFactors = FALSE)

# Loop over each chromosomal arm
for (i in 1:nrow(myeloma_arms_df)) {
  arm_name <- myeloma_arms_df$Arm[i]
  chr_value <- myeloma_arms_df$chr[i]
  arm_value <- myeloma_arms_df$arm[i]
  
  # Filter the data for the specific chromosome and arm
  data_arm <- long_data %>%
    filter(chr == chr_value, arm == arm_value)
  
  # Debug: print dimensions and unique values for one arm (optional)
  # print(dim(data_arm))
  # print(unique(data_arm$Value))
  
  # Determine if it's an amplification or deletion
  if (startsWith(arm_name, "amp")) {
    data_arm <- data_arm %>%
      mutate(IsAltered = case_when(
        is.na(Value) ~ NA_real_,
        Value %in% gain_labels ~ 1,
        TRUE ~ 0
      ))
  } else if (startsWith(arm_name, "del")) {
    data_arm <- data_arm %>%
      mutate(IsAltered = case_when(
        is.na(Value) ~ NA_real_,
        Value %in% loss_labels ~ 1,
        TRUE ~ 0
      ))
  } else {
    next
  }
  
  # Calculate the proportion of altered segments per sample
  sample_proportions <- data_arm %>%
    group_by(Sample) %>%
    summarise(
      NumSegmentsEvaluated = sum(!is.na(IsAltered)),
      ProportionAltered = if_else(
        NumSegmentsEvaluated > 0,
        mean(IsAltered, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )
  
  # Ensure all samples are included (set missing samples to 0)
  missing_samples <- setdiff(samples, sample_proportions$Sample)
  if (length(missing_samples) > 0) {
    sample_proportions <- bind_rows(
      sample_proportions,
      data.frame(
        Sample = missing_samples,
        NumSegmentsEvaluated = 0L,
        ProportionAltered = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Merge with results and assign 1 if proportion > 1/3
  # The 33% threshold (>1/3 of arm segments altered) was chosen to
  # implement the project-specific manuscript rule for WGS-based arm calls.
  # The denominator is the number of evaluated ichorCNA segments on the arm;
  # it must not be interpreted as genomic coverage or percent abnormal cells.
  results <- left_join(results, sample_proportions, by = "Sample")
  results[[arm_name]] <- ifelse(
    is.na(results$ProportionAltered),
    NA_integer_,
    as.integer(results$ProportionAltered > 1/3)
  )
  results[[paste0(arm_name, "_segments_evaluated")]] <-
    results$NumSegmentsEvaluated
  results[[paste0(arm_name, "_proportion_altered")]] <-
    results$ProportionAltered
  results <- results %>% select(-NumSegmentsEvaluated, -ProportionAltered)
}

# View results to check the output
print(results)

myeloma_CNA_matrix_with_HRD <- results



## Add hyperdiploid info:

## Edit based on Suzanne new criteria 
## Define canonical hyperdiploidy chromosomes
myeloma_cna_HRD <- list(
  hyperdiploid_chr3  = list(chr = "3"),
  hyperdiploid_chr5  = list(chr = "5"),
  hyperdiploid_chr7  = list(chr = "7"),
  hyperdiploid_chr9  = list(chr = "9"),
  hyperdiploid_chr11 = list(chr = "11"),
  hyperdiploid_chr15 = list(chr = "15"),
  hyperdiploid_chr19 = list(chr = "19"),
  hyperdiploid_chr21 = list(chr = "21")
)

# Flatten the hyperdiploid list into a data frame
hyperdiploid_arms_df <- data.frame(
  Chr = names(myeloma_cna_HRD),
  chr = sapply(myeloma_cna_HRD, function(x) x$chr),
  stringsAsFactors = FALSE
)

# Loop over hyperdiploid chromosomes and assess broad gains (threshold >65%)
for (i in 1:nrow(hyperdiploid_arms_df)) {
  chr_value <- hyperdiploid_arms_df$chr[i]
  
  # Filter the data for the specific chromosome
  data_chr <- long_data %>%
    filter(chr == chr_value)
  
  # Mark segments as gained if they match the gain labels. `%in%` returns FALSE
  # for NA here, so cross-sample missing segment calls enter the chromosome
  # denominator as 0 (not gained) rather than remaining unevaluable.
  data_chr <- data_chr %>%
    mutate(IsGain = ifelse(Value %in% gain_labels, 1, 0))
  
  # Calculate the proportion of segments gained per sample for this chromosome
  sample_gains <- data_chr %>%
    group_by(Sample) %>%
    summarise(ProportionGained = mean(IsGain))
  
  # Ensure all samples are included (assign 0 if no data)
  missing_samples <- setdiff(samples, sample_gains$Sample)
  if (length(missing_samples) > 0) {
    sample_gains <- bind_rows(
      sample_gains,
      data.frame(Sample = missing_samples, ProportionGained = 0)
    )
  }
  
  # Merge with the results table and assign gain indicator if >65% of segments are gained
  # The 65% threshold for per-chromosome gain is intentionally permissive:
  # in cfDNA at low tumor fractions ichorCNA may not call every bin as
  # GAIN even on a truly trisomic chromosome. Requiring >65% (rather than
  # >80%) improves sensitivity while a downstream requirement for >=5/8
  # canonical HRD chromosomes provides specificity at the hyperdiploid
  # phenotype level. Updated from 80% threshold per Suzanne's criteria.
  results <- left_join(results, sample_gains, by = "Sample")
  results[[hyperdiploid_arms_df$Chr[i]]] <- ifelse(results$ProportionGained > 0.65, 1, 0)
  results <- results %>% select(-ProportionGained)
}

# Mark a sample as hyperdiploid when >=5 canonical chromosomes are broadly gained
# Hyperdiploid MM (HRD) is defined by gains of odd-numbered chromosomes
# 3,5,7,9,11,15,19,21. Requiring >=5/8 (rather than all 8) reflects that
# WGS may miss gains on individual chromosomes at low tumor fractions,
# with >=5 chosen here as the project-specific specificity threshold.
# Intermediate values (3-4) are treated as not HRD.
results <- results %>%
  mutate(
    hyperdiploid = if_else(
      rowSums(select(., starts_with("hyperdiploid_"))) >= 5,
      1, 0
    )
  )

myeloma_CNA_matrix_with_HRD <- results


# Prepare CNA data
# Final output columns:
#   Sample       – BAM basename (without .bam) matching Bam_clean_tmp in metadata
#   del1p        – 1 if >33% of chr1p arm segments show loss (HETD/HOMD)
#   amp1q        – 1 if >33% of chr1q arm segments show gain (GAIN/AMP/HLAMP)
#   del13q       – 1 if >33% of chr13q arm segments show loss
#   del17p       – 1 if >33% of chr17p arm segments show loss
#   hyperdiploid – 1 if >=5/8 canonical HRD chromosomes show >65% gain
# Consumed by: 1_5_Integrate_WGS_Feature_Data.R (cna_data_ichorCNA.rds)
# We'll select relevant columns and set Sample as row names
cna_data <- myeloma_CNA_matrix_with_HRD %>%
  dplyr::select(Sample, del1p, amp1q, del13q, del17p, hyperdiploid) %>%
  mutate_at(vars(-Sample), as.character)  # Convert numeric to character for consistency

# Keep the full hyperdiploidy chromosome-by-chromosome matrix as support/QC so
# analysts can audit why a sample did or did not meet the >=5/8 HRD rule.
cna_data_backup <- cna_data
cna_data_hyperdiploid_qc <- myeloma_CNA_matrix_with_HRD %>%
  mutate_at(vars(-Sample), as.character)

# Arm-level audit table preserves the evaluated-segment denominator and raw
# altered proportion behind every binary call.
cna_arm_call_qc <- myeloma_CNA_matrix_with_HRD %>%
  dplyr::select(
    Sample,
    matches("^(del1p|amp1q|del13q|del17p)(_segments_evaluated|_proportion_altered)?$")
  )
readr::write_tsv(cna_arm_call_qc, support_file("cna_arm_call_qc.tsv"))

## Export final ichorCNA arm-level feature table
# This is the main output consumed by 1_5_Integrate_WGS_Feature_Data.R.
# cna_data_ichorCNA.rds  – compact binary arm-level call table; read in
#   1_5_Integrate_WGS_Feature_Data.R as the ichorCNA-derived CNA layer
#   before merging with Sequenza calls (1_4A output). The RDS format
#   preserves column types exactly and loads ~10× faster than the TSV.
saveRDS(cna_data, file = file.path(export_dir, "cna_data_ichorCNA.rds"))

write.table(cna_data, file = file.path(export_dir, "cna_data_ichorCNA.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Active CNA feature table written: ", file.path(export_dir, "cna_data_ichorCNA.rds"))

# Support-only exports. These are useful for auditing and historical comparison
# but are not read by downstream manuscript scripts.
write.table(
  cna_data_backup,
  file = support_file("cna_data_compact_backup.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  cna_data_hyperdiploid_qc,
  file = support_file("cna_data_hyperdiploid_chromosome_qc.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message("Support CNA QC tables written under: ", support_table_dir)
