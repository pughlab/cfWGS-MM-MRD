# =============================================================================
# Script:   3_1_Optimize_cfWGS_thresholds.R
# Project:  cfWGS MRD detection in multiple myeloma (MM)
#           Part of Abelman et al. (2025) manuscript
# How to run from the repository root:
#   Rscript 3_1_Optimize_cfWGS_thresholds.R
#
# Current manuscript outputs created/updated by the default preserved-model path:
#   - Figure 3B: BM model confusion matrices.
#   - Figure 4B: blood/cfDNA model confusion matrices.
#   - Extended Data Figure 5B-C: BM training- and test-cohort performance.
#   - Extended Data Figure 7B-C and 7E: blood/cfDNA performance and ROC panels.
#   - Extended Data Figure 9C-F: fragmentomics performance and confusion
#     panels. The test-cohort panel in 9F is calculated here from the current
#     37-sample scored set.
#   - Supplementary Table 5: final all-model metrics workbook. Note that the
#     historical source file is named Supplementary_Table_4_All_Model_Metrics_*
#     because of a manuscript-numbering correction; the artifact map records it
#     as Supplementary Table 5.
#   - Classifier-performance inputs used by
#     3_1C_Expanded_test_clustered_sensitivity.R to assemble Supplementary
#     Table 6. The final Table 6 workbook is written by that later script.
#
# Current grouped-CV outputs made elsewhere:
#   Figure 3A, Figure 4A, Extended Data Figure 5A, Extended Data Figure 7A,
#   and Extended Data Figure 9A-B come from the patient-grouped repeated
#   nested-CV workflow in scripts 6_12, 6_13, 6_14_Generate, and 6_17. Final
#   Supplementary Table 4 is also assembled from that workflow. This script
#   retains the older sample-row CV panels and tables only as analysis history.
#
# Pipeline role and reproducibility note:
#   The default path reads the February 2026 models and thresholds, applies them
#   to the current cohort, and updates the outputs listed above without fitting
#   the models again. Setting CFWGS_RETRAIN_MODELS=1 enters the older development
#   and refitting code below the preserved-model exit. That older branch is not
#   the current patient-grouped repeated nested-CV analysis and should not be
#   used to recreate final Supplementary Table 4 or the six grouped-CV figure
#   panels listed above.
#
# Author:   Dory Abelman
# Date:     May 28, 2025 (updated through Feb 2026)
#
# OVERVIEW:
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# This is the CORE model optimization and evaluation script for the cfWGS MRD detection project.
# It performs the following major tasks:
#
#   1. DATA IMPORT & COHORT DEFINITION
#      Load clinical features, cfWGS metrics, and define patient cohorts
#      (Frontline=training, Non-frontline=validation/hold-out)
#
#   2. GROUND-TRUTH MRD LABELING
#      Create binary MRD labels from clonoSEQ and MFC assays
#      (Hierarchy: clonoSEQ > MFC; missing = NA)
#
#   3. THRESHOLD OPTIMIZATION (Sections 4-5)
#      Evaluate univariate and simple ridge-regression thresholds
#      NOTE: These are in the disabled historical branch. Current model
#      comparison uses patient-grouped repeated nested CV, while current
#      thresholded calls use the preserved February 2026 fitted models.
#
#   4. HISTORICAL NESTED CROSS-VALIDATION (Section 7)
#      Earlier elastic-net analysis using sample-row resampling:
#        • Outer fold: one stratified five-fold split
#        • Inner fold: five-fold CV repeated five times for tuning
#      This does not keep repeated samples from one patient in the same fold.
#      The current patient-grouped repeated analysis is in script 6_12.
#      Three model "families":
#        - BM models: Trained on bone-marrow-derived WGS features
#        - Blood models: Trained on blood cfDNA WGS features  
#        - Fragmentomics models: Trained on 3 different cohorts
#          * Full cohort (all samples with fragmentomics data)
#          * BM-restricted (only samples with BM data, using BM folds)
#          * Blood-restricted (only samples with blood data, using blood folds)
#
#   5. MODEL APPLICATION & SCORING (Section 9)
#      Apply trained models to full cohort with patient-eligibility masking:
#        • BM models → score only patients with BM sample
#        • Blood models → score only patients with blood sample
#        • Fragmentomics → score all patients (no masking)
#
#   6. PERFORMANCE METRICS & VALIDATION (Sections 9B-10)
#      Calculate:
#        • ROC/AUC with confidence intervals (via Delong)
#        • Sensitivity, specificity at fixed thresholds
#        • Calibration, Brier score, partial AUC
#        • Sensitivity @95% specificity, specificity @95% sensitivity
#
#   7. TABLE & FIGURE GENERATION (Sections 10+)
#      Export supplementary tables and create manuscript figures
#
# KEY DESIGN DECISIONS:
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# • HISTORICAL NESTED CV: The outer fold is independent of inner tuning, but
#   splitting is by sample row rather than patient. Use script 6_12 for the
#   current patient-grouped repeated estimate.
# • FOLD REUSE: Fragmentomics models on restricted cohorts use exact same folds as original
#   BM/blood models to ensure fair comparison.
# • MODEL NAMING: Fragmentomics models include cohort suffix (_Full, _BM_restricted, _Blood_restricted)
#   to preserve cohort distinction through entire pipeline.
# • MASKING: Patient eligibility by sample type prevents scoring BM models on patients with only
#   blood samples (and vice versa).
#
# OUTPUTS FROM THE NORMAL PRESERVED-MODEL RUN:
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   • all_patients_with_BM_and_blood_calls_updated6.rds|csv (masked scores)
#   • all_patients_with_BM_and_blood_calls_updated6_full.rds|csv (unmasked scores)
#   • Current test-cohort metric and audit CSV files used by later table scripts
#   • Figure 3B, Figure 4B, and Extended Data performance/confusion panels
#
# The all_nested_cv_metrics, all_cfWGS_models_list, and all_model_thresholds
# files named later in this script belong to the disabled historical refitting
# branch. They are not created by the normal run and are not the 50-repeat
# patient-grouped nested-CV results.
#
# =============================================================================
# Pipeline status:
#   Active but cache-sensitive. Routine runs use the preserved-model path.
#   Deliberate model validation is performed by script 6_12 rather than by the
#   historical refitting branch below.
#

# -----------------------------------------------------------------------------
# 1. Setup & Package Loading
# -----------------------------------------------------------------------------

# -------- 0.  Load packages --------------------------------------------------
# Core data manipulation & statistical packages
library(dplyr)          # data wrangling (filter, select, mutate, join, etc.)
library(tidyr)          # reshape (pivot_wider, pivot_longer, drop_na, etc.)
library(purrr)          # functional programming (map_dfr, imap, walk2, etc.)
library(stringr)        # string manipulation (str_detect, str_remove, str_replace)
library(forcats)        # factor operations (fct_reorder, fct_recode)
library(janitor)        # data cleaning (tabyl, clean_names)
library(glue)           # string interpolation for logging

# Statistical & ML packages
library(caret)          # machine learning framework
                        #   - train(): fit models with CV control
                        #   - trainControl(): configure CV strategy (5x5 CV, metrics, etc.)
                        #   - createFolds(): stratified CV splits
library(glmnet)         # elastic-net regularized regression
                        #   - cv.glmnet(): cross-validated Ridge/Lasso/ElasticNet
library(pROC)           # ROC analysis and AUC computation
                        #   - roc(): build ROC curve object
                        #   - auc(): extract AUC value with CI
                        #   - coords(): get sensitivity/specificity at thresholds
library(Matrix)         # sparse matrix operations (dgCMatrix)
                        #   - as(, "dgCMatrix"): convert to sparse format for large X matrices

# Visualization & output packages
library(viridis)        # color-blind friendly palettes for ggplot
library(patchwork)      # combine multiple ggplot figures into panels
library(scales)         # axis formatting (percent_format, comma, etc.)

# Specialized packages
library(DescTools)      # calibration plots and additional statistical functions
library(PRROC)          # precision-recall curve analysis
library(readr)          # fast CSV reading
library(readxl)         # Excel file reading

# Shared helper for final manuscript-organized outputs.
# The model script retains its historical filenames and cached model artifacts,
# while these helper calls copy the final figures and tables into
# final_manuscript_objects/ with final labels such as Figure_3B, Figure_4B,
# Extended_Data_Figure_5B, and Supplementary_Table_5.
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

# -----------------------------------------------------------------------------
# 2. Data Import & Cohort Definition
# -----------------------------------------------------------------------------

### Set paths 
outdir   <- "Output_tables_2025"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

### Load data 
file <- readRDS("Final_aggregate_table_cfWGS_features_with_clinical_and_demographics_updated9.rds")
cohort_df <- load_final_cohort_assignment()

dat <- file 

# 1.  Join cohort_df and keep frontline only -------------------------------------
dat <- dat %>%                # <‑‑ the master data frame
  left_join(cohort_df, by = "Patient") 


# =============================================================================
# 3. Ground-Truth MRD Label Construction
# =============================================================================
#
# WHAT:   Define the binary outcome variable (MRD_truth) used to train
#         all models. MRD truth comes from two sources:
#         • clonoSEQ: primary molecular MRD reference assay
#         • MFC:     flow cytometry MRD assay (backup)
#
# HIERARCHY: Use clonoSEQ if available; otherwise use MFC; missing=NA
#           This ensures we use the most reliable test for each patient.
#
# WHY:    We need a reliable binary outcome (MRD positive/negative) to:
#         (1) Train supervised learning models
#         (2) Evaluate model performance
#         (3) Optimize decision thresholds
#         
# FILTERING: We then exclude Baseline and Diagnosis timepoints since we
#            want to evaluate MRD detection at post-treatment visits.
#
# =============================================================================

#    Create MRD_truth reference label
# Choose which external test(s) are most reliable:
# "Adaptive_Binary"-> clonoSEQ (high sensitivity molecular assay)
# "Flow_Binary"    -> MFC (flow cytometry assay)
# Prefer clonoSEQ, fall back to MFC if missing
dat <- dat %>% 
  mutate(
    MRD_truth = case_when(
      !is.na(Adaptive_Binary)            ~ Adaptive_Binary,             # use clonoSEQ if available
      is.na(Adaptive_Binary)
      & !is.na(Flow_Binary)              ~ Flow_Binary,                 # else use MFC 
      TRUE                               ~ NA_real_                     # missing if neither assay run
    )
  )

## Keep only MRD timepoints to optimize timepoint 
dat_mrd <- dat %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline"))

USE_PRESERVED_MODELS_ONLY <- Sys.getenv("CFWGS_RETRAIN_MODELS", unset = "0") != "1"

# Rebuild the manuscript ED7E panel from the frozen full-training/refit models.
# This is deterministic and does not refit a classifier or select a threshold.
# Colors match the patient-grouped CV model palette used by script 6_14:
# Sites = blue and Combined = purple.
write_full_training_ed7e_panel <- function(output_png) {
  model_path <- "nested_blood_validation_updated3.rds"
  frame_path <- "nested_blood_validation_updated3_run2.rds"
  missing_inputs <- c(model_path, frame_path)[!file.exists(c(model_path, frame_path))]
  if (length(missing_inputs)) {
    stop(
      "Cannot rebuild ED7E; missing frozen input(s): ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }

  refit_object <- readRDS(model_path)
  frame_object <- readRDS(frame_path)
  evaluation_data <- frame_object$models$Blood_plus_fragment_min$trainingData
  required_models <- c("Blood_zscore_only_sites", "Blood_plus_fragment_min")
  if (is.null(evaluation_data) || nrow(evaluation_data) != 34L ||
      !".outcome" %in% names(evaluation_data)) {
    stop("ED7E requires the preserved 34-row full-training evaluation frame.", call. = FALSE)
  }
  if (!all(required_models %in% names(refit_object$models))) {
    stop(
      "ED7E frozen refit object is missing model(s): ",
      paste(setdiff(required_models, names(refit_object$models)), collapse = ", "),
      call. = FALSE
    )
  }

  roc_df <- purrr::map_dfr(required_models, function(model_name) {
    probability <- predict(
      refit_object$models[[model_name]],
      newdata = evaluation_data,
      type = "prob"
    )[, "pos"]
    roc_object <- pROC::roc(
      response = evaluation_data$.outcome,
      predictor = probability,
      levels = c("neg", "pos"),
      direction = "<",
      quiet = TRUE
    )
    tibble::tibble(
      combo = model_name,
      fpr = 1 - roc_object$specificities,
      tpr = roc_object$sensitivities
    )
  }) %>%
    dplyr::mutate(combo = factor(.data$combo, levels = required_models))

  # These coordinates were calculated from the canonical scored table in the
  # original analysis and are intentionally kept separate from the ROC inputs.
  marker_pts <- tibble::tribble(
    ~combo, ~threshold, ~fpr, ~tpr,
    "Blood_zscore_only_sites", 0.432, 0.500, 0.750,
    "Blood_zscore_only_sites", 0.380, 1.000, 0.938,
    "Blood_plus_fragment_min", 0.435, 0.167, 0.875
  ) %>%
    dplyr::mutate(combo = factor(.data$combo, levels = required_models))

  model_colours <- c(
    Blood_zscore_only_sites = "#2C7FB8",
    Blood_plus_fragment_min = "#6A3D9A"
  )
  model_labels <- c(
    Blood_zscore_only_sites = "Sites model",
    Blood_plus_fragment_min = "Combined model"
  )

  ed7e_plot <- ggplot2::ggplot(
    roc_df,
    ggplot2::aes(x = .data$fpr, y = .data$tpr, colour = .data$combo)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(linetype = 2, colour = "grey60") +
    ggplot2::geom_point(
      data = marker_pts,
      ggplot2::aes(x = .data$fpr, y = .data$tpr),
      fill = "white", colour = "darkgrey", stroke = 0.9, size = 2,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(
      name = "Model",
      values = model_colours,
      breaks = required_models,
      labels = model_labels,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "False-positive rate (1 - specificity)",
      y = "True-positive rate (sensitivity)",
      title = "Full-Cohort ROC (Refit on all Training Samples)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.title.position = "plot"
    )

  dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output_png, plot = ed7e_plot, width = 5.5, height = 4, dpi = 600)
  invisible(output_png)
}

if (USE_PRESERVED_MODELS_ONLY) {
  message("Using preserved cfWGS models and thresholds; no model retraining will be run.")

  selected_models <- readRDS(file.path(outdir, "selected_combo_models_2026-02-16.rds"))
  selected_thr <- readRDS(file.path(outdir, "selected_combo_thresholds_2026-02-16.rds"))

  if (!all(names(selected_models) %in% names(selected_thr))) {
    stop(
      "Preserved selected model names do not all have matching thresholds: ",
      paste(setdiff(names(selected_models), names(selected_thr)), collapse = ", ")
    )
  }

  required_predictors <- unique(unlist(lapply(
    selected_models,
    function(fit) setdiff(names(fit$trainingData), ".outcome")
  )))
  missing_predictors <- setdiff(required_predictors, names(dat))
  if (length(missing_predictors) > 0L) {
    stop(
      "Cannot apply preserved cfWGS models because predictors are missing from the aggregate table: ",
      paste(missing_predictors, collapse = ", ")
    )
  }

  apply_preserved_selected <- function(dat, models, thresholds, positive_class = "pos") {
    out <- dat
    n <- nrow(dat)
    common <- intersect(names(models), names(thresholds))

    for (cmb in common) {
      fit <- models[[cmb]]
      thr <- thresholds[[cmb]]
      preds <- setdiff(names(fit$trainingData), ".outcome")
      can_score <- stats::complete.cases(dat[, preds, drop = FALSE])

      probs <- rep(NA_real_, n)
      calls <- rep(NA_integer_, n)

      if (any(can_score)) {
        probs[can_score] <- predict(
          fit,
          newdata = dat[can_score, preds, drop = FALSE],
          type = "prob"
        )[[positive_class]]
        calls[can_score] <- as.integer(probs[can_score] >= thr)
      }

      out[[paste0(cmb, "_prob")]] <- probs
      out[[paste0(cmb, "_call")]] <- calls
    }

    out
  }

  Good_pts <- read.csv("baseline_high_quality_patients_updated.csv",
                       stringsAsFactors = FALSE)
  bm_good_patients <- Good_pts %>%
    filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
    pull(Patient) %>%
    unique()
  cfDNA_good_patients <- Good_pts %>%
    filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
    pull(Patient) %>%
    unique()

  data_scored <- apply_preserved_selected(
    dat = dat,
    models = selected_models,
    thresholds = selected_thr,
    positive_class = "pos"
  )

  # The February 2026 preserved fragmentomics models carry explicit training-
  # cohort suffixes. Older downstream survival scripts use the historical
  # unsuffixed names for the full-cohort fragmentomics model, so keep aliases.
  full_fragmentomics_aliases <- c(
    "Fragmentomics_full",
    "Fragmentomics_min",
    "Fragmentomics_FS_only",
    "Fragmentomics_mean_coverage_only",
    "Fragmentomics_prop_short_only",
    "Fragmentomics_tumor_fraction_only"
  )
  for (base_name in full_fragmentomics_aliases) {
    for (suffix in c("prob", "call")) {
      full_col <- paste0(base_name, "_Full_", suffix)
      alias_col <- paste0(base_name, "_", suffix)
      if (full_col %in% names(data_scored) && !alias_col %in% names(data_scored)) {
        data_scored[[alias_col]] <- data_scored[[full_col]]
      }
    }
  }

  score_cols <- grep("_(prob|call)$", names(data_scored), value = TRUE, ignore.case = TRUE)
  model_map <- tibble(
    model_col = score_cols,
    model_base = str_remove(score_cols, "_(?i:prob|call)$")
  ) %>%
    mutate(family = case_when(
      str_detect(model_base, regex("^BM_", ignore_case = TRUE)) ~ "BM",
      str_detect(model_base, regex("^Blood_", ignore_case = TRUE)) ~ "BLOOD",
      str_detect(model_base, regex("^Fragmentomics_", ignore_case = TRUE)) ~ "FRAG",
      TRUE ~ "OTHER"
    ))

  typed_na <- function(x) {
    if (is.factor(x)) return(NA)
    if (is.logical(x)) return(NA)
    if (is.integer(x)) return(NA_integer_)
    if (is.double(x)) return(NA_real_)
    if (is.numeric(x)) return(NA_real_)
    if (is.character(x)) return(NA_character_)
    if (inherits(x, "Date")) return(as.Date(NA))
    if (inherits(x, "POSIXt")) return(as.POSIXct(NA))
    NA
  }

  data_scored_masked <- data_scored
  walk2(model_map$model_col, model_map$family, function(col, fam) {
    allowed <- switch(
      fam,
      "BM" = bm_good_patients,
      "BLOOD" = cfDNA_good_patients,
      "FRAG" = unique(data_scored$Patient),
      "OTHER" = unique(data_scored$Patient)
    )
    x <- data_scored[[col]]
    x[!data_scored$Patient %in% allowed] <- typed_na(x)
    data_scored_masked[[col]] <<- x
  })

  saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.rds"))
  write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.csv"))
  saveRDS(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.rds"))
  write_csv(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.csv"))

  preserved_model_sample_group <- function(model_name) {
    dplyr::case_when(
      stringr::str_detect(model_name, "^BM_") ~ "BM_restricted",
      stringr::str_detect(model_name, "^Blood_") ~ "Blood_restricted",
      stringr::str_detect(model_name, "_BM_restricted$") ~ "BM_restricted",
      stringr::str_detect(model_name, "_Blood_restricted$") ~ "Blood_restricted",
      stringr::str_detect(model_name, "_Full$") ~ "Full",
      TRUE ~ "Full"
    )
  }

  safe_div <- function(num, den) {
    ifelse(is.na(den) | den == 0, NA_real_, num / den)
  }

  metric_at_threshold <- function(truth, prob, threshold) {
    pred <- as.integer(prob >= threshold)
    tp <- sum(pred == 1 & truth == 1, na.rm = TRUE)
    tn <- sum(pred == 0 & truth == 0, na.rm = TRUE)
    fp <- sum(pred == 1 & truth == 0, na.rm = TRUE)
    fn <- sum(pred == 0 & truth == 1, na.rm = TRUE)
    sensitivity <- safe_div(tp, tp + fn)
    specificity <- safe_div(tn, tn + fp)
    ppv <- safe_div(tp, tp + fp)
    npv <- safe_div(tn, tn + fn)
    accuracy <- safe_div(tp + tn, tp + tn + fp + fn)
    f1 <- ifelse(is.na(ppv) | is.na(sensitivity) | (ppv + sensitivity) == 0,
                 NA_real_,
                 2 * ppv * sensitivity / (ppv + sensitivity))
    list(
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv,
      accuracy = accuracy,
      bal_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
      f1 = f1
    )
  }

  preserved_test_eval <- data_scored_masked %>%
    filter(!timepoint_info %in% c("Baseline", "Diagnosis"),
           Cohort == "Non-frontline",
           !is.na(MRD_truth))

  preserved_test_metrics <- purrr::map_dfr(names(selected_models), function(model_name) {
    prob_col <- paste0(model_name, "_prob")
    threshold <- selected_thr[[model_name]]
    if (!prob_col %in% names(preserved_test_eval)) {
      return(tibble(
        combo = model_name,
        auc_mean = NA_real_, sens_mean = NA_real_, spec_mean = NA_real_,
        sens95_mean = NA_real_, spec95_mean = NA_real_, balacc_mean = NA_real_,
        ppv_mean = NA_real_, npv_mean = NA_real_, f1_mean = NA_real_,
        brier_mean = NA_real_, pAUC90_mean = NA_real_, acc_mean = NA_real_,
        sample_group = preserved_model_sample_group(model_name),
        eval_cohort = "Testing"
      ))
    }

    eval_df <- preserved_test_eval %>%
      transmute(
        MRD_truth = as.integer(MRD_truth),
        prob = as.numeric(.data[[prob_col]])
      ) %>%
      filter(!is.na(MRD_truth), !is.na(prob))

    if (nrow(eval_df) < 2L || dplyr::n_distinct(eval_df$MRD_truth) < 2L) {
      return(tibble(
        combo = model_name,
        auc_mean = NA_real_, sens_mean = NA_real_, spec_mean = NA_real_,
        sens95_mean = NA_real_, spec95_mean = NA_real_, balacc_mean = NA_real_,
        ppv_mean = NA_real_, npv_mean = NA_real_, f1_mean = NA_real_,
        brier_mean = if (nrow(eval_df)) mean((eval_df$prob - eval_df$MRD_truth)^2, na.rm = TRUE) else NA_real_,
        pAUC90_mean = NA_real_, acc_mean = NA_real_,
        sample_group = preserved_model_sample_group(model_name),
        eval_cohort = "Testing"
      ))
    }

    roc_obj <- pROC::roc(
      response = eval_df$MRD_truth,
      predictor = eval_df$prob,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    fixed_metrics <- metric_at_threshold(eval_df$MRD_truth, eval_df$prob, threshold)
    roc_all <- pROC::coords(
      roc_obj,
      x = "all",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    ) %>%
      as_tibble()

    spec95_row <- roc_all %>%
      filter(specificity >= 0.95) %>%
      arrange(desc(sensitivity), desc(specificity)) %>%
      slice(1)
    if (!nrow(spec95_row)) {
      spec95_row <- roc_all %>%
        arrange(desc(specificity), desc(sensitivity)) %>%
        slice(1)
    }

    tibble(
      combo = model_name,
      auc_mean = as.numeric(pROC::auc(roc_obj)),
      sens_mean = fixed_metrics$sensitivity,
      spec_mean = fixed_metrics$specificity,
      sens95_mean = spec95_row$sensitivity[[1]],
      spec95_mean = spec95_row$specificity[[1]],
      balacc_mean = fixed_metrics$bal_accuracy,
      ppv_mean = fixed_metrics$ppv,
      npv_mean = fixed_metrics$npv,
      f1_mean = fixed_metrics$f1,
      brier_mean = mean((eval_df$prob - eval_df$MRD_truth)^2, na.rm = TRUE),
      pAUC90_mean = as.numeric(pROC::auc(
        roc_obj,
        partial.auc = c(0.9, 1),
        partial.auc.correct = TRUE
      )),
      acc_mean = fixed_metrics$accuracy,
      sample_group = preserved_model_sample_group(model_name),
      eval_cohort = "Testing"
    )
  }) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  preserved_test_metric_audit <- purrr::map_dfr(names(selected_models), function(model_name) {
    prob_col <- paste0(model_name, "_prob")
    if (!prob_col %in% names(preserved_test_eval)) {
      return(tibble(combo = model_name, n_evaluable = 0L, n_mrd_positive = NA_integer_, n_mrd_negative = NA_integer_))
    }
    eval_df <- preserved_test_eval %>%
      transmute(MRD_truth = as.integer(MRD_truth), prob = as.numeric(.data[[prob_col]])) %>%
      filter(!is.na(MRD_truth), !is.na(prob))
    tibble(
      combo = model_name,
      n_evaluable = nrow(eval_df),
      n_mrd_positive = sum(eval_df$MRD_truth == 1, na.rm = TRUE),
      n_mrd_negative = sum(eval_df$MRD_truth == 0, na.rm = TRUE)
    )
  })

  write_csv(
    preserved_test_metrics,
    "Final Tables and Figures/Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv"
  )
  write_csv(
    preserved_test_metrics,
    file.path(outdir, "preserved_model_current_test_cohort_metrics_for_supplementary_table_6.csv")
  )
  write_csv(
    preserved_test_metric_audit,
    file.path(outdir, "preserved_model_current_test_cohort_metrics_audit.csv")
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv",
    artifact_id = "STABLE6",
    role = "table_csv",
    description = "Current testing-cohort model-performance table computed by applying preserved February 2026 models and thresholds to the expanded test cohort.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )

  # Rebuild the Figure 3B BM confusion matrices from the same preserved models
  # applied above. The cache-sensitive branch exits before the historical
  # plotting block, so these panels must be generated here to keep the manuscript
  # object tree synchronized with the expanded test cohort without retraining.
  selected_bm_models <- c(
    "BM_zscore_only_detection_rate",
    "BM_base_zscore"
  )
  selected_bm_labels <- c(
    BM_zscore_only_detection_rate = "cVAF model",
    BM_base_zscore = "Combined model"
  )

  missing_bm_models <- setdiff(selected_bm_models, names(selected_models))
  missing_bm_thresholds <- setdiff(selected_bm_models, names(selected_thr))
  if (length(missing_bm_models) > 0L || length(missing_bm_thresholds) > 0L) {
    stop(
      "Cannot regenerate Figure 3B preserved-model confusion matrices. ",
      "Missing models: ", paste(missing_bm_models, collapse = ", "),
      "; missing thresholds: ", paste(missing_bm_thresholds, collapse = ", "),
      call. = FALSE
    )
  }

  build_preserved_confusion_df <- function(eval_dat, cohort_label, model_names, model_labels) {
    purrr::map_dfr(model_names, function(model_name) {
      prob_col <- paste0(model_name, "_prob")
      if (!prob_col %in% names(eval_dat)) {
        stop("Missing preserved-model probability column: ", prob_col, call. = FALSE)
      }

      model_eval <- eval_dat %>%
        dplyr::transmute(
          MRD_truth = as.integer(.data$MRD_truth),
          prob = as.numeric(.data[[prob_col]])
        ) %>%
        dplyr::filter(!is.na(.data$MRD_truth), !is.na(.data$prob)) %>%
        dplyr::mutate(
          Obs = factor(dplyr::if_else(.data$MRD_truth == 1L, "pos", "neg"),
                       levels = c("neg", "pos")),
          Pred = factor(dplyr::if_else(.data$prob >= selected_thr[[model_name]], "pos", "neg"),
                        levels = c("neg", "pos"))
        )

      if (nrow(model_eval) == 0L) {
        stop(
          "No evaluable rows for Figure 3B ", cohort_label, " confusion matrix: ",
          model_name,
          call. = FALSE
        )
      }

      cm <- caret::confusionMatrix(model_eval$Pred, model_eval$Obs, positive = "pos")
      as_tibble(cm$table) %>%
        dplyr::rename(Pred = Prediction, Obs = Reference, Count = n) %>%
        tidyr::complete(
          Pred = factor(c("neg", "pos"), levels = c("neg", "pos")),
          Obs = factor(c("neg", "pos"), levels = c("neg", "pos")),
          fill = list(Count = 0L)
        ) %>%
        dplyr::mutate(
          model = model_labels[[model_name]],
          model_key = model_name,
          cohort = cohort_label,
          threshold = as.numeric(selected_thr[[model_name]]),
          n_evaluable = nrow(model_eval),
          n_mrd_positive = sum(model_eval$Obs == "pos"),
          n_mrd_negative = sum(model_eval$Obs == "neg"),
          PPV = unname(cm$byClass["Pos Pred Value"]),
          NPV = unname(cm$byClass["Neg Pred Value"])
        )
    }) %>%
      dplyr::mutate(model = factor(.data$model, levels = unname(model_labels[model_names])))
  }

  build_frozen_training_confusion_df <- function(
      models, thresholds, model_names, model_labels, positive_class = "pos") {
    # Training-cohort panels must be reproducible from the locked model-training
    # cohort, not from the current aggregate table.  caret::train objects retain
    # their original trainingData matrix, including the locked .outcome labels.
    # This prevents test-cohort metadata refreshes from changing training-panel
    # denominators while still allowing test panels to be regenerated live.
    purrr::map_dfr(model_names, function(model_name) {
      fit <- models[[model_name]]
      if (is.null(fit$trainingData) || !".outcome" %in% names(fit$trainingData)) {
        stop(
          "Cannot reconstruct locked training confusion matrix for ",
          model_name, ": model object does not contain caret trainingData/.outcome.",
          call. = FALSE
        )
      }
      if (!model_name %in% names(thresholds)) {
        stop("Missing preserved threshold for locked training model: ", model_name, call. = FALSE)
      }

      probs <- predict(fit, newdata = fit$trainingData, type = "prob")[[positive_class]]
      preds <- factor(dplyr::if_else(probs >= thresholds[[model_name]], "pos", "neg"),
                      levels = c("neg", "pos"))
      obs <- factor(as.character(fit$trainingData$.outcome), levels = c("neg", "pos"))
      cm <- caret::confusionMatrix(preds, obs, positive = "pos")

      as_tibble(cm$table) %>%
        dplyr::rename(Pred = Prediction, Obs = Reference, Count = n) %>%
        tidyr::complete(
          Pred = factor(c("neg", "pos"), levels = c("neg", "pos")),
          Obs = factor(c("neg", "pos"), levels = c("neg", "pos")),
          fill = list(Count = 0L)
        ) %>%
        dplyr::mutate(
          model = model_labels[[model_name]],
          model_key = model_name,
          cohort = "Training Cohort",
          threshold = as.numeric(thresholds[[model_name]]),
          n_evaluable = length(obs),
          n_mrd_positive = sum(obs == "pos"),
          n_mrd_negative = sum(obs == "neg"),
          PPV = unname(cm$byClass["Pos Pred Value"]),
          NPV = unname(cm$byClass["Neg Pred Value"])
        )
    }) %>%
      dplyr::mutate(model = factor(.data$model, levels = unname(model_labels[model_names])))
  }

  plot_preserved_confusion_df <- function(cm_df, plot_title) {
    plot_df <- cm_df %>%
      dplyr::group_by(.data$model) %>%
      dplyr::mutate(text_col = dplyr::if_else(
        .data$Count >= 0.7 * max(.data$Count, na.rm = TRUE),
        "white",
        "black"
      )) %>%
      dplyr::ungroup()

    ggplot(plot_df, aes(x = Pred, y = Obs, fill = Count)) +
      geom_tile(color = "white", linewidth = 0) +
      geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
      facet_wrap(~ model) +
      scale_fill_gradient(
        name = "Count",
        low = "#f2f2f2",
        high = "#4a4a4a",
        limits = c(0, max(plot_df$Count, na.rm = TRUE)),
        oob = scales::squish
      ) +
      scale_color_identity() +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = c("pos", "neg")) +
      labs(
        x = "cfWGS-MRD",
        y = "BM-MRD",
        title = plot_title
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, vjust = 0),
        axis.title = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  }

  # -----------------------------------------------------------------------
  # Annotate the sequencing platform used for fragmentomics evaluation
  #
  # Why: the Spring 2026 expansion introduced genuinely new NovaSeq XPlus
  # libraries, whereas the original analysis cohort was generated on NovaSeq
  # 6000. Platform labels let us report performance separately and detect a
  # platform-associated shift without changing the frozen model predictions.
  # -----------------------------------------------------------------------
  FRAGMENTOMICS_PLATFORM_NOVASEQ6000 <- "NovaSeq 6000"
  FRAGMENTOMICS_PLATFORM_XPLUS <- "NovaSeq XPlus"
  FRAGMENTOMICS_PLATFORM_COMBINED <- "NovaSeq 6000 + NovaSeq XPlus"

  # Construct a case- and separator-insensitive patient/sample join key.
  # Sample_Code is preferred; patient-timepoint is a fallback for older rows.
  # Plasma/DNA suffixes are removed because they are inconsistently present
  # across the scored table and the revision metadata.
  normalize_fragmentomics_platform_key <- function(patient, sample_code, timepoint = NULL) {
    key_sample <- dplyr::coalesce(as.character(sample_code), paste(patient, timepoint, sep = "-"))
    key_sample <- stringr::str_replace_all(toupper(key_sample), "_", "-")
    key_sample <- stringr::str_replace_all(key_sample, "\\s+", "")
    key_sample <- stringr::str_remove(key_sample, "-(P|PL|PLASMA|CFDNA|DNA)$")
    paste(stringr::str_replace_all(toupper(as.character(patient)), "_", "-"), key_sample, sep = "||")
  }

  # Add auditable platform columns to a scored sample table.
  #
  # Classification rule: only revision rows explicitly marked as having no
  # pre-existing patient-timepoint-type overlap are new XPlus samples. All
  # original or overlapping rows remain NovaSeq 6000. If optional revision
  # metadata is unavailable, retain the original-cohort default and record the
  # fallback in fragmentomics_platform_reason rather than silently guessing.
  annotate_fragmentomics_platform <- function(df) {
    revision <- load_spring2026_revision_metadata(required = FALSE)
    out <- df %>%
      dplyr::mutate(
        fragmentomics_platform_key = normalize_fragmentomics_platform_key(
          .data$Patient,
          .data$Sample_Code,
          .data$Timepoint
        ),
        fragmentomics_revision_metadata_match = FALSE,
        fragmentomics_revision_new_xplus = FALSE,
        fragmentomics_sequencing_platform = FRAGMENTOMICS_PLATFORM_NOVASEQ6000,
        fragmentomics_platform_reason = "No matching new Spring 2026 revision row"
      )

    if (is.null(revision) || nrow(revision) == 0L) {
      out$fragmentomics_platform_reason <- "Spring 2026 revision metadata unavailable; defaulted to original cohort"
      return(out)
    }

    required_cols <- c("Patient", "Sample_ID", "existing_repo_patient_timepoint_type_overlap")
    missing_cols <- setdiff(required_cols, names(revision))
    if (length(missing_cols) > 0L) {
      stop(
        "Cannot classify fragmentomics sequencing platform because revision metadata is missing columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    # Collapse duplicate revision records to one key. `any(FALSE overlap)` is
    # intentionally conservative: one confirmed new library is sufficient to
    # classify the matched scored sample as XPlus.
    revision_keys <- revision %>%
      dplyr::transmute(
        fragmentomics_platform_key = normalize_fragmentomics_platform_key(
          .data$Patient,
          .data$Sample_ID,
          NA_character_
        ),
        revision_existing_overlap = .data$existing_repo_patient_timepoint_type_overlap,
        revision_bam = if ("Bam" %in% names(revision)) as.character(.data$Bam) else NA_character_
      ) %>%
      dplyr::group_by(.data$fragmentomics_platform_key) %>%
      dplyr::summarise(
        revision_metadata_match = TRUE,
        revision_new_xplus = any(.data$revision_existing_overlap %in% FALSE, na.rm = TRUE),
        revision_overlap_missing = any(is.na(.data$revision_existing_overlap)),
        revision_bam = paste(unique(stats::na.omit(.data$revision_bam)), collapse = ";"),
        .groups = "drop"
      )

    out <- out %>%
      dplyr::left_join(revision_keys, by = "fragmentomics_platform_key") %>%
      dplyr::mutate(
        fragmentomics_revision_metadata_match = dplyr::coalesce(.data$revision_metadata_match, FALSE),
        fragmentomics_revision_new_xplus = dplyr::coalesce(.data$revision_new_xplus, FALSE),
        fragmentomics_sequencing_platform = dplyr::case_when(
          .data$fragmentomics_revision_new_xplus ~ FRAGMENTOMICS_PLATFORM_XPLUS,
          TRUE ~ FRAGMENTOMICS_PLATFORM_NOVASEQ6000
        ),
        fragmentomics_platform_reason = dplyr::case_when(
          .data$fragmentomics_revision_new_xplus ~
            "Matched Spring 2026 revision row with no existing repo patient-timepoint-type overlap",
          .data$fragmentomics_revision_metadata_match ~
            "Matched Spring 2026 revision row that overlaps the original repo cohort",
          TRUE ~
            "No matching new Spring 2026 revision row"
        )
      ) %>%
      dplyr::select(-dplyr::any_of(c("revision_metadata_match", "revision_new_xplus")))

    out
  }

  # Persist the annotated scored objects because all downstream platform-split
  # tables must be traceable to the same sample-level predictions.  Apply the
  # non-model fragmentomics provenance columns to both the full and masked
  # artifacts so audit comparisons can enforce identical schemas; the two
  # tables should differ only in model-score masking, not in metadata presence.
  data_scored <- annotate_fragmentomics_platform(data_scored)
  data_scored_masked <- annotate_fragmentomics_platform(data_scored_masked)
  saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.rds"))
  write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.csv"))
  saveRDS(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.rds"))
  write_csv(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.csv"))

  # Write a row-level audit containing the join outcome, platform rationale,
  # BAM provenance (when available), and whether any full model was evaluable.
  fragmentomics_probability_cols <- grep("^Fragmentomics_.*_Full_prob$", names(data_scored_masked), value = TRUE)
  fragmentomics_platform_audit <- data_scored_masked %>%
    dplyr::mutate(
      fragmentomics_any_full_model_evaluable = if (length(fragmentomics_probability_cols) == 0L) {
        FALSE
      } else {
        rowSums(!is.na(dplyr::pick(dplyr::all_of(fragmentomics_probability_cols)))) > 0
      }
    ) %>%
    dplyr::select(
      dplyr::any_of(c(
        "Patient", "Timepoint", "Sample_Code", "timepoint_info", "Cohort", "Study",
        "MRD_truth", "fragmentomics_sequencing_platform", "fragmentomics_platform_reason",
        "fragmentomics_revision_metadata_match", "fragmentomics_revision_new_xplus",
        "fragmentomics_any_full_model_evaluable", "revision_bam"
      ))
    )
  readr::write_csv(
    fragmentomics_platform_audit,
    file.path(outdir, "fragmentomics_platform_split_audit.csv")
  )

  # -----------------------------------------------------------------------
  # Platform-correction impact and robust sensitivity audit
  #
  # Reapply the same frozen models to two alternate XPlus feature versions:
  #   1) uncorrected raw XPlus values;
  #   2) median/MAD harmonized XPlus values.
  # The primary scored table above remains the prespecified mean/SD mapping.
  # -----------------------------------------------------------------------
  fragmentomics_feature_map <- tibble::tribble(
    ~primary, ~unharmonized, ~robust,
    "FS", "FS_unharmonized", "FS_harmonized_median_mad",
    "Proportion.Short", "Proportion.Short_unharmonized", "Proportion.Short_harmonized_median_mad",
    "Mean.Coverage", "Mean.Coverage_unharmonized", "Mean.Coverage_harmonized_median_mad"
  )
  required_impact_columns <- unique(unlist(fragmentomics_feature_map, use.names = FALSE))
  missing_impact_columns <- setdiff(required_impact_columns, names(dat))
  if (length(missing_impact_columns)) {
    stop(
      "Cannot build fragmentomics platform impact audit; aggregate table is missing: ",
      paste(missing_impact_columns, collapse = ", "),
      call. = FALSE
    )
  }

  xplus_index <- data_scored$fragmentomics_sequencing_platform == FRAGMENTOMICS_PLATFORM_XPLUS
  uncorrected_dat <- dat
  robust_dat <- dat
  for (i in seq_len(nrow(fragmentomics_feature_map))) {
    primary <- fragmentomics_feature_map$primary[[i]]
    unharmonized <- fragmentomics_feature_map$unharmonized[[i]]
    robust <- fragmentomics_feature_map$robust[[i]]
    uncorrected_dat[[primary]][xplus_index] <- uncorrected_dat[[unharmonized]][xplus_index]
    robust_dat[[primary]][xplus_index] <- robust_dat[[robust]][xplus_index]
  }

  uncorrected_scored <- apply_preserved_selected(
    dat = uncorrected_dat,
    models = selected_models,
    thresholds = selected_thr,
    positive_class = "pos"
  )
  robust_scored <- apply_preserved_selected(
    dat = robust_dat,
    models = selected_models,
    thresholds = selected_thr,
    positive_class = "pos"
  )

  affected_models <- names(selected_models)[vapply(
    selected_models,
    function(fit) {
      predictors <- setdiff(names(fit$trainingData), ".outcome")
      length(intersect(predictors, fragmentomics_feature_map$primary)) > 0L
    },
    logical(1)
  )]
  if (!length(affected_models)) {
    stop("No preserved models were identified as using fragmentomics predictors.", call. = FALSE)
  }

  fragmentomics_platform_impact <- purrr::map_dfr(affected_models, function(model_name) {
    prob_col <- paste0(model_name, "_prob")
    call_col <- paste0(model_name, "_call")
    tibble::tibble(
      Patient = data_scored$Patient[xplus_index],
      Timepoint = data_scored$Timepoint[xplus_index],
      Sample_Code = data_scored$Sample_Code[xplus_index],
      Date = data_scored$Date[xplus_index],
      timepoint_info = data_scored$timepoint_info[xplus_index],
      Cohort = data_scored$Cohort[xplus_index],
      fragmentomics_sequencing_platform = data_scored$fragmentomics_sequencing_platform[xplus_index],
      model = model_name,
      FS_uncorrected = dat$FS_unharmonized[xplus_index],
      FS_mean_sd = dat$FS[xplus_index],
      FS_median_mad = dat$FS_harmonized_median_mad[xplus_index],
      Proportion_Short_uncorrected = dat$Proportion.Short_unharmonized[xplus_index],
      Proportion_Short_mean_sd = dat$Proportion.Short[xplus_index],
      Proportion_Short_median_mad = dat$Proportion.Short_harmonized_median_mad[xplus_index],
      Mean_Coverage_uncorrected = dat$Mean.Coverage_unharmonized[xplus_index],
      Mean_Coverage_mean_sd = dat$Mean.Coverage[xplus_index],
      Mean_Coverage_median_mad = dat$Mean.Coverage_harmonized_median_mad[xplus_index],
      probability_uncorrected = uncorrected_scored[[prob_col]][xplus_index],
      probability_mean_sd = data_scored[[prob_col]][xplus_index],
      probability_median_mad = robust_scored[[prob_col]][xplus_index],
      call_uncorrected = uncorrected_scored[[call_col]][xplus_index],
      call_mean_sd = data_scored[[call_col]][xplus_index],
      call_median_mad = robust_scored[[call_col]][xplus_index]
    ) %>%
      dplyr::mutate(
        probability_delta_mean_sd_vs_uncorrected = .data$probability_mean_sd - .data$probability_uncorrected,
        probability_delta_median_mad_vs_uncorrected = .data$probability_median_mad - .data$probability_uncorrected,
        call_flip_mean_sd_vs_uncorrected = dplyr::if_else(
          !is.na(.data$call_uncorrected) & !is.na(.data$call_mean_sd),
          .data$call_uncorrected != .data$call_mean_sd,
          NA
        ),
        call_flip_median_mad_vs_uncorrected = dplyr::if_else(
          !is.na(.data$call_uncorrected) & !is.na(.data$call_median_mad),
          .data$call_uncorrected != .data$call_median_mad,
          NA
        ),
        call_flip_median_mad_vs_mean_sd = dplyr::if_else(
          !is.na(.data$call_mean_sd) & !is.na(.data$call_median_mad),
          .data$call_mean_sd != .data$call_median_mad,
          NA
        )
      )
  })

  fragmentomics_platform_impact_summary <- fragmentomics_platform_impact %>%
    dplyr::group_by(.data$model) %>%
    dplyr::summarise(
      n_xplus_samples = dplyr::n_distinct(.data$Sample_Code),
      n_evaluable_uncorrected = sum(!is.na(.data$probability_uncorrected)),
      n_evaluable_mean_sd = sum(!is.na(.data$probability_mean_sd)),
      n_evaluable_median_mad = sum(!is.na(.data$probability_median_mad)),
      median_probability_delta_mean_sd_vs_uncorrected = median(
        .data$probability_delta_mean_sd_vs_uncorrected,
        na.rm = TRUE
      ),
      median_probability_delta_median_mad_vs_uncorrected = median(
        .data$probability_delta_median_mad_vs_uncorrected,
        na.rm = TRUE
      ),
      n_call_flips_mean_sd_vs_uncorrected = sum(.data$call_flip_mean_sd_vs_uncorrected, na.rm = TRUE),
      n_call_flips_median_mad_vs_uncorrected = sum(.data$call_flip_median_mad_vs_uncorrected, na.rm = TRUE),
      n_call_flips_median_mad_vs_mean_sd = sum(.data$call_flip_median_mad_vs_mean_sd, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(
    fragmentomics_platform_impact,
    file.path(outdir, "fragmentomics_platform_correction_sample_model_impact_audit.csv")
  )
  readr::write_csv(
    fragmentomics_platform_impact_summary,
    file.path(outdir, "fragmentomics_platform_correction_model_summary.csv")
  )

  # Define mutually exclusive evaluation cohorts after excluding baseline and
  # diagnosis rows, where MRD response classification is not meaningful.
  # "Frontline" preserves the development/training cohort; "Non-frontline"
  # is the independent test cohort, including newly incorporated samples.
  preserved_train_eval <- data_scored_masked %>%
    dplyr::filter(
      !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
      .data$Cohort == "Frontline",
      !is.na(.data$MRD_truth)
    )
  preserved_test_eval <- data_scored_masked %>%
    dplyr::filter(
      !.data$timepoint_info %in% c("Baseline", "Diagnosis"),
      .data$Cohort == "Non-frontline",
      !is.na(.data$MRD_truth)
    )

  # Regenerate Figure 3B by applying the frozen February 2026 BM models and
  # thresholds. These calls are evaluation only: no fitting or threshold
  # optimization occurs in preserved-model mode.
  figure3b_training_confusion <- build_frozen_training_confusion_df(
    models = selected_models,
    thresholds = selected_thr,
    model_names = selected_bm_models,
    model_labels = selected_bm_labels
  )
  figure3b_test_confusion <- build_preserved_confusion_df(
    preserved_test_eval,
    cohort_label = "Test Cohort",
    model_names = selected_bm_models,
    model_labels = selected_bm_labels
  )

  readr::write_csv(
    dplyr::bind_rows(figure3b_training_confusion, figure3b_test_confusion),
    file.path(outdir, "preserved_model_current_figure3b_confusion_source_data.csv")
  )
  readr::write_csv(
    figure3b_test_confusion,
    "Final Tables and Figures/Figure_3B_BM_test_confusion_source_data.csv"
  )

  figure3b_training_plot <- plot_preserved_confusion_df(
    figure3b_training_confusion,
    "Confusion Matrix at Youden Index in Training Cohort"
  )
  figure3b_test_plot <- plot_preserved_confusion_df(
    figure3b_test_confusion,
    "Confusion Matrix at Youden Index in Test Cohort"
  )

  ggsave(
    "Final Tables and Figures/Fig4C_confusion_tables_primary_updated5.png",
    plot = figure3b_training_plot,
    width = 5,
    height = 2.75,
    dpi = 600,
    bg = "white"
  )
  ggsave(
    "Final Tables and Figures/Fig4C_confusion_tables_test5.png",
    plot = figure3b_test_plot,
    width = 5,
    height = 2.75,
    dpi = 600,
    bg = "white"
  )

  ms_copy_artifact(
    source_path = "Final Tables and Figures/Fig4C_confusion_tables_primary_updated5.png",
    artifact_id = "FIG3B",
    role = "figure_panel_png_training",
    description = "Locked February 2026 BM model training-cohort confusion matrix.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Fig4C_confusion_tables_test5.png",
    artifact_id = "FIG3B",
    role = "figure_panel_png_test",
    description = "BM model test-cohort confusion matrix regenerated by applying preserved February 2026 models to the expanded test cohort.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = file.path(outdir, "preserved_model_current_figure3b_confusion_source_data.csv"),
    artifact_id = "FIG3B",
    role = "source_data_csv",
    description = "Source confusion counts for Figure 3B. Training-cohort counts are locked February 2026 manuscript values; test-cohort counts are regenerated from the current expanded test cohort using preserved models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Figure_3B_BM_test_confusion_source_data.csv",
    artifact_id = "FIG3B",
    role = "supporting_data_csv_test_confusion",
    description = "Test-cohort source confusion counts for Figure 3B after applying preserved February 2026 BM models to the current cohort.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )

  # Repeat the frozen-model confusion workflow for the blood/cfDNA models used
  # in Figure 4B, failing early if either a model or threshold is absent.
  selected_blood_models <- c("Blood_zscore_only_sites", "Blood_plus_fragment")
  selected_blood_labels <- c(
    Blood_zscore_only_sites = "Sites model",
    Blood_plus_fragment = "Combined model"
  )
  missing_blood_models <- setdiff(selected_blood_models, names(selected_models))
  missing_blood_thresholds <- setdiff(selected_blood_models, names(selected_thr))
  if (length(missing_blood_models) > 0L || length(missing_blood_thresholds) > 0L) {
    stop(
      "Cannot regenerate Figure 4B preserved-model confusion matrices. ",
      "Missing models: ", paste(missing_blood_models, collapse = ", "),
      "; missing thresholds: ", paste(missing_blood_thresholds, collapse = ", "),
      call. = FALSE
    )
  }

  figure4b_training_confusion <- build_frozen_training_confusion_df(
    models = selected_models,
    thresholds = selected_thr,
    model_names = selected_blood_models,
    model_labels = selected_blood_labels
  )
  figure4b_test_confusion <- build_preserved_confusion_df(
    preserved_test_eval,
    cohort_label = "Test Cohort",
    model_names = selected_blood_models,
    model_labels = selected_blood_labels
  )
  readr::write_csv(
    dplyr::bind_rows(figure4b_training_confusion, figure4b_test_confusion),
    file.path(outdir, "preserved_model_current_figure4b_confusion_source_data.csv")
  )
  readr::write_csv(
    figure4b_test_confusion,
    "Final Tables and Figures/Figure_4B_blood_test_confusion_source_data.csv"
  )
  ggsave(
    "Final Tables and Figures/Fig5C_confusion_tables_primary_blood6.png",
    plot = plot_preserved_confusion_df(
      figure4b_training_confusion,
      "Confusion Matrix at Youden Index in Training Cohort"
    ),
    width = 5, height = 2.75, dpi = 600, bg = "white"
  )
  ggsave(
    "Final Tables and Figures/Fig5C_confusion_tables_test_blood6.png",
    plot = plot_preserved_confusion_df(
      figure4b_test_confusion,
      "Confusion Matrix at Youden Index in Test Cohort"
    ),
    width = 5, height = 2.75, dpi = 600, bg = "white"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Fig5C_confusion_tables_primary_blood6.png",
    artifact_id = "FIG4B", role = "figure_panel_png_training",
    description = "Locked February 2026 blood model training-cohort confusion matrix.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Fig5C_confusion_tables_test_blood6.png",
    artifact_id = "FIG4B", role = "figure_panel_png_test",
    description = "Blood model test-cohort confusion matrix regenerated from preserved models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = file.path(outdir, "preserved_model_current_figure4b_confusion_source_data.csv"),
    artifact_id = "FIG4B", role = "source_data_csv",
    description = "Source confusion counts for Figure 4B. Training-cohort counts are locked February 2026 manuscript values; test-cohort counts are regenerated from the current expanded test cohort using preserved models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Figure_4B_blood_test_confusion_source_data.csv",
    artifact_id = "FIG4B", role = "supporting_data_csv_test_confusion",
    description = "Test-cohort source confusion counts for Figure 4B from preserved blood models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )

  # Select the two fragmentomics models used in the ED Figure 9 confusion
  # panels. The broader three-model set below is used for performance summaries.
  selected_fragmentomics_models <- c(
    "Fragmentomics_mean_coverage_only_Full",
    "Fragmentomics_min_Full"
  )
  selected_fragmentomics_labels <- c(
    Fragmentomics_mean_coverage_only_Full = "MM-sites coverage",
    Fragmentomics_min_Full = "Fragment score + MM-sites"
  )
  missing_fragmentomics_models <- setdiff(selected_fragmentomics_models, names(selected_models))
  missing_fragmentomics_thresholds <- setdiff(selected_fragmentomics_models, names(selected_thr))
  if (length(missing_fragmentomics_models) > 0L || length(missing_fragmentomics_thresholds) > 0L) {
    stop(
      "Cannot regenerate Extended Data Figure 9E/F preserved-model confusion matrices. ",
      "Missing models: ", paste(missing_fragmentomics_models, collapse = ", "),
      "; missing thresholds: ", paste(missing_fragmentomics_thresholds, collapse = ", "),
      call. = FALSE
    )
  }

  edfigure9_training_confusion <- build_frozen_training_confusion_df(
    models = selected_models,
    thresholds = selected_thr,
    model_names = selected_fragmentomics_models,
    model_labels = selected_fragmentomics_labels
  )
  edfigure9_test_confusion <- build_preserved_confusion_df(
    preserved_test_eval,
    cohort_label = "Test Cohort",
    model_names = selected_fragmentomics_models,
    model_labels = selected_fragmentomics_labels
  )

  # Models shown in the platform-aware performance comparison: the combined
  # model and its two principal single-feature counterparts.
  fragmentomics_platform_models <- c(
    "Fragmentomics_full_Full",
    "Fragmentomics_mean_coverage_only_Full",
    "Fragmentomics_prop_short_only_Full"
  )
  fragmentomics_platform_labels <- c(
    Fragmentomics_full_Full = "Combined model",
    Fragmentomics_mean_coverage_only_Full = "Coverage model",
    Fragmentomics_prop_short_only_Full = "Prop. short model"
  )
  missing_platform_models <- setdiff(fragmentomics_platform_models, names(selected_models))
  missing_platform_thresholds <- setdiff(fragmentomics_platform_models, names(selected_thr))
  if (length(missing_platform_models) > 0L || length(missing_platform_thresholds) > 0L) {
    stop(
      "Cannot regenerate platform-split fragmentomics performance panels. ",
      "Missing models: ", paste(missing_platform_models, collapse = ", "),
      "; missing thresholds: ", paste(missing_platform_thresholds, collapse = ", "),
      call. = FALSE
    )
  }

  # Calculate threshold-dependent binary performance after complete-case
  # filtering. Undefined denominators remain NA through safe_div(); an empty
  # subset returns a typed one-row result so downstream binding is stable.
  calc_binary_metrics <- function(truth, prob, threshold) {
    keep <- !is.na(truth) & !is.na(prob)
    truth <- as.integer(truth[keep])
    prob <- as.numeric(prob[keep])
    if (length(truth) == 0L) {
      return(tibble(
        threshold = threshold, n_evaluable = 0L, n_mrd_positive = NA_integer_,
        n_mrd_negative = NA_integer_, sensitivity = NA_real_, specificity = NA_real_,
        bal_accuracy = NA_real_, f1 = NA_real_
      ))
    }
    pred <- as.integer(prob >= threshold)
    tp <- sum(pred == 1L & truth == 1L)
    tn <- sum(pred == 0L & truth == 0L)
    fp <- sum(pred == 1L & truth == 0L)
    fn <- sum(pred == 0L & truth == 1L)
    sensitivity <- safe_div(tp, tp + fn)
    specificity <- safe_div(tn, tn + fp)
    precision <- safe_div(tp, tp + fp)
    f1 <- ifelse(
      is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
      NA_real_,
      2 * precision * sensitivity / (precision + sensitivity)
    )
    tibble(
      threshold = threshold,
      n_evaluable = length(truth),
      n_mrd_positive = sum(truth == 1L),
      n_mrd_negative = sum(truth == 0L),
      sensitivity = sensitivity,
      specificity = specificity,
      bal_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
      f1 = f1
    )
  }

  # Derive the fixed high-sensitivity operating point from the preserved
  # training data only. The independent test cohort must never be used to
  # choose its own threshold: doing so would leak outcome information into the
  # evaluation and force the displayed test sensitivity toward the target by
  # construction. The selected point maximizes training specificity among ROC
  # coordinates reaching at least 95% empirical training sensitivity.
  training_95sens_threshold <- function(model_name) {
    fit <- selected_models[[model_name]]
    if (is.null(fit) || is.null(fit$trainingData) || !".outcome" %in% names(fit$trainingData)) {
      stop(
        "Cannot derive training-fixed 95% sensitivity threshold for ",
        model_name,
        ": preserved model does not contain caret trainingData/.outcome.",
        call. = FALSE
      )
    }

    training_eval <- tibble::tibble(
      MRD_truth = dplyr::case_when(
        as.character(fit$trainingData$.outcome) == "pos" ~ 1L,
        as.character(fit$trainingData$.outcome) == "neg" ~ 0L,
        TRUE ~ NA_integer_
      ),
      prob = as.numeric(
        predict(fit, newdata = fit$trainingData, type = "prob")[["pos"]]
      )
    ) %>%
      dplyr::filter(!is.na(.data$MRD_truth), !is.na(.data$prob))

    if (nrow(training_eval) < 2L || dplyr::n_distinct(training_eval$MRD_truth) < 2L) {
      stop(
        "Cannot derive training-fixed 95% sensitivity threshold for ",
        model_name,
        ": preserved training data do not contain both outcome classes.",
        call. = FALSE
      )
    }

    roc_obj <- pROC::roc(
      response = training_eval$MRD_truth,
      predictor = training_eval$prob,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
    roc_df <- pROC::coords(
      roc_obj,
      x = "all",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    ) %>%
      as_tibble()

    candidate <- roc_df %>%
      dplyr::filter(.data$sensitivity >= 0.95) %>%
      dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$sensitivity)) %>%
      dplyr::slice(1)
    if (nrow(candidate) == 0L) {
      stop(
        "Cannot find a training ROC coordinate reaching 95% sensitivity for ",
        model_name,
        ".",
        call. = FALSE
      )
    }
    as.numeric(candidate$threshold[[1]])
  }

  # Generate one current test-cohort performance panel and all of its source
  # data/provenance copies. Bars show preserved Youden-threshold metrics;
  # triangles show test-cohort metrics at a high-sensitivity operating point
  # fixed exclusively in the preserved training data. Legacy copies are
  # refreshed because manuscript assembly still reads those paths, while
  # ms_copy_artifact registers the labeled manuscript copies.
  write_current_test_performance_panel <- function(eval_dat,
                                                   model_names,
                                                   model_labels,
                                                   output_png,
                                                   source_csv,
                                                   artifact_id,
                                                   description,
                                                   legacy_component_png = NULL,
                                                   legacy_component_source_csv = NULL,
                                                   legacy_export_sidecar_csv = NULL) {
    missing_models <- setdiff(model_names, names(selected_models))
    missing_thresholds <- setdiff(model_names, names(selected_thr))
    if (length(missing_models) > 0L || length(missing_thresholds) > 0L) {
      stop(
        "Cannot regenerate ", artifact_id, " current test-cohort performance panel. ",
        "Missing models: ", paste(missing_models, collapse = ", "),
        "; missing thresholds: ", paste(missing_thresholds, collapse = ", "),
        call. = FALSE
      )
    }

    # Keep one row per model and operating point in the source CSV before
    # pivoting to the long plotting representation.
    perf_df <- purrr::map_dfr(model_names, function(model_name) {
      prob_col <- paste0(model_name, "_prob")
      if (!prob_col %in% names(eval_dat)) {
        stop("Missing preserved-model probability column: ", prob_col, call. = FALSE)
      }
      dplyr::bind_rows(
        calc_binary_metrics(
          truth = eval_dat$MRD_truth,
          prob = eval_dat[[prob_col]],
          threshold = as.numeric(selected_thr[[model_name]])
        ) %>% dplyr::mutate(operating_point = "Youden"),
        calc_binary_metrics(
          truth = eval_dat$MRD_truth,
          prob = eval_dat[[prob_col]],
          threshold = training_95sens_threshold(model_name)
        ) %>% dplyr::mutate(operating_point = "Training-derived 95% sensitivity")
      ) %>%
        dplyr::mutate(model = model_name, combo = unname(model_labels[[model_name]]))
    })

    plot_df <- perf_df %>%
      dplyr::select(
        .data$combo, .data$operating_point, .data$sensitivity,
        .data$specificity, .data$bal_accuracy, .data$f1
      ) %>%
      tidyr::pivot_longer(
        cols = c(.data$sensitivity, .data$specificity, .data$bal_accuracy, .data$f1),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        metric = dplyr::recode(
          .data$metric,
          sensitivity = "Sens",
          specificity = "Spec",
          bal_accuracy = "Acc",
          f1 = "F1"
        )
      )

    readr::write_csv(perf_df, source_csv)
    ggsave(
      output_png,
      plot = ggplot(
        plot_df %>% dplyr::filter(.data$operating_point == "Youden"),
        aes(x = .data$combo, y = .data$value, fill = .data$combo)
      ) +
        geom_col(width = 0.6) +
        geom_point(
          data = plot_df %>% dplyr::filter(.data$operating_point == "Training-derived 95% sensitivity"),
          aes(y = .data$value, shape = "95% sensitivity"),
          size = 2,
          colour = "black"
        ) +
        scale_shape_manual(name = NULL, values = c("95% sensitivity" = 17)) +
        facet_wrap(
          ~ metric,
          nrow = 1,
          labeller = labeller(metric = c(
            F1 = "F1",
            Sens = "Sensitivity",
            Spec = "Specificity",
            Acc = "Bal. Accuracy"
          ))
        ) +
        scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
        labs(
          x = NULL,
          y = NULL,
          title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
        ) +
        theme_classic(base_size = 9) +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          strip.text = element_text(face = "bold", size = 8),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          panel.spacing = unit(0.8, "lines"),
          legend.position = "bottom",
          legend.box.margin = margin(t = -10),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
        ),
      width = 5,
      height = 3.5,
      dpi = 600,
      bg = "white"
    )

    # Create parent directories before synchronizing compatibility copies; any
    # failed copy is fatal to prevent a partially updated submission package.
    legacy_paths <- c(legacy_component_png, legacy_component_source_csv, legacy_export_sidecar_csv)
    legacy_paths <- legacy_paths[!is.na(legacy_paths)]
    for (legacy_path in legacy_paths) {
      dir.create(dirname(legacy_path), recursive = TRUE, showWarnings = FALSE)
    }
    if (!is.null(legacy_component_png)) {
      copied_png <- file.copy(output_png, legacy_component_png, overwrite = TRUE)
      if (!isTRUE(copied_png)) {
        stop("Failed to update legacy generated-component PNG for ", artifact_id, ": ", legacy_component_png, call. = FALSE)
      }
    }
    if (!is.null(legacy_component_source_csv)) {
      for (legacy_source_csv in legacy_component_source_csv) {
        copied_source <- file.copy(source_csv, legacy_source_csv, overwrite = TRUE)
        if (!isTRUE(copied_source)) {
          stop("Failed to update legacy generated-component source CSV for ", artifact_id, ": ", legacy_source_csv, call. = FALSE)
        }
      }
    }
    if (!is.null(legacy_export_sidecar_csv)) {
      for (legacy_sidecar_csv in legacy_export_sidecar_csv) {
        copied_sidecar <- file.copy(source_csv, legacy_sidecar_csv, overwrite = TRUE)
        if (!isTRUE(copied_sidecar)) {
          stop("Failed to update legacy manuscript-export sidecar CSV for ", artifact_id, ": ", legacy_sidecar_csv, call. = FALSE)
        }
      }
    }

    ms_copy_artifact(
      source_path = output_png,
      artifact_id = artifact_id,
      role = "figure_panel_png",
      description = description,
      script_name = "3_1_Optimize_cfWGS_thresholds.R"
    )
    ms_copy_artifact(
      source_path = source_csv,
      artifact_id = artifact_id,
      role = "source_data_csv_current_test_cohort",
      description = paste(description, "Source data generated from the current preserved-model test cohort."),
      script_name = "3_1_Optimize_cfWGS_thresholds.R"
    )
    invisible(perf_df)
  }

  # Register training/refit panels that are intentionally not recomputed by
  # the expanded test-cohort workflow. This preserves the scientific boundary
  # between frozen training results and newly evaluated test results.
  sync_existing_panel_copy <- function(source_png,
                                       artifact_id,
                                       description,
                                       legacy_component_png = NULL) {
    if (!file.exists(source_png)) {
      stop("Cannot synchronize existing panel for ", artifact_id, ". Missing file: ", source_png, call. = FALSE)
    }
    if (!is.null(legacy_component_png)) {
      for (legacy_png in legacy_component_png) {
        dir.create(dirname(legacy_png), recursive = TRUE, showWarnings = FALSE)
        copied_png <- file.copy(source_png, legacy_png, overwrite = TRUE)
        if (!isTRUE(copied_png)) {
          stop("Failed to update legacy generated-component PNG for ", artifact_id, ": ", legacy_png, call. = FALSE)
        }
      }
    }
    ms_copy_artifact(
      source_path = source_png,
      artifact_id = artifact_id,
      role = "figure_panel_png",
      description = description,
      script_name = "3_1_Optimize_cfWGS_thresholds.R"
    )
    invisible(source_png)
  }

  sync_existing_panel_copy(
    source_png = "Final Tables and Figures/Supp5A_classifier_performance_bar_updated3.png",
    artifact_id = "EDFIG5B",
    description = "BM training-cohort performance panel copied from the retained manuscript output; not recomputed for test-cohort expansion.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_5/panel_B",
      "Supp5A_classifier_performance_bar_updated3.png"
    )
  )

  sync_existing_panel_copy(
    source_png = "Final Tables and Figures/Supp7A_classifier_performance_bar_updated_blood_muts2.png",
    artifact_id = "EDFIG7B",
    description = "Blood/cfDNA training-cohort performance panel copied from the retained manuscript output; not recomputed for test-cohort expansion.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_7/panel_B",
      "Supp7A_classifier_performance_bar_updated_blood_muts2.png"
    )
  )

  sync_existing_panel_copy(
    source_png = "Final Tables and Figures/Supp9D_classifier_performance_bar_updated_frag2.png",
    artifact_id = "EDFIG9C",
    description = "Fragmentomics training-cohort performance panel copied from the retained manuscript output; not recomputed for test-cohort expansion.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_9/panel_C",
      "Supp9D_classifier_performance_bar_updated_frag2.png"
    )
  )

  ed7e_main_cohort_png <- file.path(
    "Final Tables and Figures",
    "Supp7D_ROC_performance_blood_refit_main_cohort.png"
  )
  write_full_training_ed7e_panel(ed7e_main_cohort_png)
  sync_existing_panel_copy(
    source_png = ed7e_main_cohort_png,
    artifact_id = "EDFIG7E",
    description = "Blood/cfDNA full-cohort refit ROC panel copied from the retained manuscript output; this training/refit panel is not recomputed for test-cohort expansion.",
    legacy_component_png = c(
      file.path(
        "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_7/panel_E",
        "Supp7D_ROC_performance_blood_updated3.png"
      ),
      file.path(
        "Scripts_2025/Final_Scripts/final_manuscript_objects/02_extended_data_figures/Extended_Data_Figure_7/Extended_Data_Figure_7E",
        "ED7E_ROC_refit_main_cohort.png"
      )
    )
  )

  # Recompute only the three independent test-cohort performance panels (BM,
  # blood/cfDNA, and fragmentomics) from the current scored sample table.
  write_current_test_performance_panel(
    eval_dat = preserved_test_eval,
    model_names = c("BM_base_zscore", "BM_zscore_only_detection_rate", "BM_zscore_only_sites"),
    model_labels = c(
      BM_base_zscore = "Combined model",
      BM_zscore_only_detection_rate = "cVAF model",
      BM_zscore_only_sites = "Sites model"
    ),
    output_png = "Final Tables and Figures/Supp5A_classifier_performance_bar_test_cohort_updated3.png",
    source_csv = file.path(outdir, "current_test_cohort_EDFIG5C_BM_performance_source_data.csv"),
    artifact_id = "EDFIG5C",
    description = "BM classifier test-cohort performance summary regenerated from current preserved-model scored calls for Extended Data Figure 5C.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_5/panel_C",
      "Supp5A_classifier_performance_bar_test_cohort_updated3.png"
    ),
    legacy_component_source_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_5/panel_C",
        "CodeOcean_GitHub_Staging/reproducible_workflow/outputs/generated/figure_components/Extended_Data_Figure_5/panel_C"
      ),
      "Extended_Data_Figure_5C_BM_test_performance_source_data.csv"
    ),
    legacy_export_sidecar_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_5/generated_outputs",
        "CodeOcean_GitHub_Staging/Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_5/generated_outputs"
      ),
      "EDFIG5C__C__generated_sidecar_file__Extended_Data_Figure_5C_BM_test_performance_source_data.csv"
    )
  )

  write_current_test_performance_panel(
    eval_dat = preserved_test_eval,
    model_names = c("Blood_plus_fragment", "Blood_rate_only", "Blood_zscore_only_sites"),
    model_labels = c(
      Blood_plus_fragment = "Combined model",
      Blood_rate_only = "cVAF model",
      Blood_zscore_only_sites = "Sites model"
    ),
    output_png = "Final Tables and Figures/Supp_7B_classifier_performance_bar_test_cohort_updated3.png",
    source_csv = file.path(outdir, "current_test_cohort_EDFIG7C_blood_performance_source_data.csv"),
    artifact_id = "EDFIG7C",
    description = "Blood/cfDNA classifier test-cohort performance summary regenerated from current preserved-model scored calls for Extended Data Figure 7C.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_7/panel_C",
      "Supp_7B_classifier_performance_bar_test_cohort_updated3.png"
    ),
    legacy_component_source_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_7/panel_C",
        "CodeOcean_GitHub_Staging/reproducible_workflow/outputs/generated/figure_components/Extended_Data_Figure_7/panel_C"
      ),
      "Extended_Data_Figure_7C_blood_test_performance_source_data.csv"
    ),
    legacy_export_sidecar_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_7/generated_outputs",
        "CodeOcean_GitHub_Staging/Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_7/generated_outputs"
      ),
      "EDFIG7C__C__generated_sidecar_file__Extended_Data_Figure_7C_blood_test_performance_source_data.csv"
    )
  )

  write_current_test_performance_panel(
    eval_dat = preserved_test_eval,
    model_names = c(
      "Fragmentomics_full_Full",
      "Fragmentomics_mean_coverage_only_Full",
      "Fragmentomics_prop_short_only_Full"
    ),
    model_labels = c(
      Fragmentomics_full_Full = "Combined model",
      Fragmentomics_mean_coverage_only_Full = "Coverage model",
      Fragmentomics_prop_short_only_Full = "Prop. short model"
    ),
    output_png = "Final Tables and Figures/Supp_Fig9F_classifier_performance_bar_test_cohort_updated2_frag2.png",
    source_csv = file.path(outdir, "current_test_cohort_EDFIG9D_fragmentomics_performance_source_data.csv"),
    artifact_id = "EDFIG9D",
    description = "Fragmentomics classifier test-cohort performance summary regenerated from current preserved-model scored calls for Extended Data Figure 9D.",
    legacy_component_png = file.path(
      "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_9/panel_D",
      "Supp_Fig9F_classifier_performance_bar_test_cohort_updated2_frag2.png"
    ),
    legacy_component_source_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/final_manuscript_objects/generated/figure_components/Extended_Data_Figure_9/panel_D",
        "CodeOcean_GitHub_Staging/reproducible_workflow/outputs/generated/figure_components/Extended_Data_Figure_9/panel_D"
      ),
      "Extended_Data_Figure_9D_fragmentomics_test_performance_source_data.csv"
    ),
    legacy_export_sidecar_csv = file.path(
      c(
        "Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_9/generated_outputs",
        "CodeOcean_GitHub_Staging/Scripts_2025/Final_Scripts/outputs/manuscript/02_extended_data_figures/Extended_Data_Figure_9/generated_outputs"
      ),
      "EDFIG9D__D__generated_sidecar_file__Extended_Data_Figure_9D_fragmentomics_test_performance_source_data.csv"
    )
  )

  # For platform comparisons, derive the high-sensitivity thresholds once in
  # the training cohort and freeze them before evaluating either platform.
  # This avoids test-set-specific threshold optimization and information leak.
  training_95sens_thresholds <- purrr::map_dbl(
    fragmentomics_platform_models,
    training_95sens_threshold
  )
  names(training_95sens_thresholds) <- fragmentomics_platform_models

  # Apply both fixed operating points to every model in a supplied cohort or
  # platform subset, returning tidy metrics with explicit cohort labels.
  build_fragmentomics_platform_performance_df <- function(eval_dat, cohort_label, platform_label) {
    purrr::map_dfr(fragmentomics_platform_models, function(model_name) {
      prob_col <- paste0(model_name, "_prob")
      if (!prob_col %in% names(eval_dat)) {
        stop("Missing preserved-model probability column: ", prob_col, call. = FALSE)
      }
      fixed_youden <- calc_binary_metrics(
        truth = eval_dat$MRD_truth,
        prob = eval_dat[[prob_col]],
        threshold = as.numeric(selected_thr[[model_name]])
      ) %>%
        dplyr::mutate(operating_point = "Youden")
      fixed_95sens <- calc_binary_metrics(
        truth = eval_dat$MRD_truth,
        prob = eval_dat[[prob_col]],
        threshold = as.numeric(training_95sens_thresholds[[model_name]])
      ) %>%
        dplyr::mutate(operating_point = "Training-derived 95% sensitivity")
      dplyr::bind_rows(fixed_youden, fixed_95sens) %>%
        dplyr::mutate(
          model = model_name,
          combo = fragmentomics_platform_labels[[model_name]],
          cohort = cohort_label,
          platform = platform_label
        )
    })
  }

  # Render the same metric layout for each platform subset so visual
  # comparisons reflect data differences rather than plotting differences.
  plot_fragmentomics_platform_performance <- function(perf_df, plot_title) {
    plot_df <- perf_df %>%
      dplyr::select(
        .data$combo, .data$operating_point, .data$sensitivity,
        .data$specificity, .data$bal_accuracy, .data$f1
      ) %>%
      tidyr::pivot_longer(
        cols = c(.data$sensitivity, .data$specificity, .data$bal_accuracy, .data$f1),
        names_to = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        metric = dplyr::recode(
          .data$metric,
          sensitivity = "Sens",
          specificity = "Spec",
          bal_accuracy = "Acc",
          f1 = "F1"
        ),
        combo = factor(.data$combo, levels = unname(fragmentomics_platform_labels[fragmentomics_platform_models]))
      )

    ggplot(
      plot_df %>% dplyr::filter(.data$operating_point == "Youden"),
      aes(x = .data$combo, y = .data$value, fill = .data$combo)
    ) +
      geom_col(width = 0.6) +
      geom_point(
        data = plot_df %>% dplyr::filter(.data$operating_point == "Training-derived 95% sensitivity"),
        aes(y = .data$value, shape = "95% sensitivity"),
        size = 2,
        colour = "black"
      ) +
      scale_shape_manual(name = NULL, values = c("95% sensitivity" = 17)) +
      facet_wrap(
        ~ metric,
        nrow = 1,
        labeller = labeller(metric = c(
          F1 = "F1",
          Sens = "Sensitivity",
          Spec = "Specificity",
          Acc = "Bal. Accuracy"
        ))
      ) +
      scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      labs(x = NULL, y = NULL, title = plot_title) +
      theme_classic(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.spacing = unit(0.8, "lines"),
        legend.position = "bottom",
        legend.box.margin = margin(t = -10),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
      )
  }

  # Convert human-readable platform labels into stable filename tokens.
  platform_slug <- function(x) {
    dplyr::case_when(
      identical(x, FRAGMENTOMICS_PLATFORM_COMBINED) ~ "combined",
      stringr::str_detect(x, "XPlus") ~ "novaseqxplus",
      TRUE ~ "novaseq6000"
    )
  }

  # Return every observed platform subset and, only when both platforms are
  # present, an additional combined set. This avoids empty or misleading
  # "combined" panels for single-platform cohorts.
  platform_eval_options <- function(cohort_df) {
    observed_platforms <- sort(unique(cohort_df$fragmentomics_sequencing_platform))
    options <- purrr::map(observed_platforms, function(platform_label) {
      list(
        platform_label = platform_label,
        eval_df = cohort_df %>%
          dplyr::filter(.data$fragmentomics_sequencing_platform == platform_label)
      )
    })
    if (length(observed_platforms) > 1L) {
      options <- c(
        options,
        list(list(
          platform_label = FRAGMENTOMICS_PLATFORM_COMBINED,
          eval_df = cohort_df
        ))
      )
    }
    options
  }

  # Evaluate training, test, and pooled descriptive cohorts. The pooled result
  # is for platform context only and is not used to tune thresholds or claim
  # independent validation. Source rows are accumulated before one final write.
  platform_split_performance_source <- list()
  platform_split_confusion_source <- list()
  platform_eval_sets <- list(
    "Training Cohort" = preserved_train_eval,
    "Test Cohort" = preserved_test_eval,
    "Training + Test Cohort" = dplyr::bind_rows(preserved_train_eval, preserved_test_eval)
  )

  for (cohort_label in names(platform_eval_sets)) {
    cohort_df <- platform_eval_sets[[cohort_label]]
    for (platform_option in platform_eval_options(cohort_df)) {
      platform_label <- platform_option$platform_label
      platform_df <- platform_option$eval_df
      if (nrow(platform_df) == 0L) next

      # Skip subsets without an observed truth/probability pair. Retaining the
      # class count in this check also makes future one-class handling visible,
      # although fixed-threshold metrics themselves do not require an ROC fit.
      has_fragmentomics_eval <- platform_df %>%
        dplyr::summarise(
          n = sum(!is.na(.data$MRD_truth) & !is.na(.data$Fragmentomics_full_Full_prob)),
          n_classes = dplyr::n_distinct(.data$MRD_truth[!is.na(.data$MRD_truth) & !is.na(.data$Fragmentomics_full_Full_prob)])
        )
      if (has_fragmentomics_eval$n[[1]] == 0L) next

      slug <- platform_slug(platform_label)
      cohort_slug <- dplyr::case_when(
        cohort_label == "Training Cohort" ~ "train",
        cohort_label == "Test Cohort" ~ "test",
        cohort_label == "Training + Test Cohort" ~ "train_test",
        TRUE ~ stringr::str_replace_all(stringr::str_to_lower(cohort_label), "[^a-z0-9]+", "_")
      )
      # Produce paired metric and confusion outputs from the identical subset
      # so panel values and exported source data share one denominator.
      perf_df <- build_fragmentomics_platform_performance_df(platform_df, cohort_label, platform_label)
      platform_split_performance_source[[paste(cohort_slug, slug, sep = "_")]] <- perf_df

      perf_title <- paste0(
        cohort_label,
        " Fragmentomics Performance\n",
        platform_label
      )
      perf_path <- file.path(
        "Final Tables and Figures",
        paste0("fragmentomics_", cohort_slug, "_performance_bar_", slug, ".png")
      )
      ggsave(
        perf_path,
        plot = plot_fragmentomics_platform_performance(perf_df, perf_title),
        width = 5,
        height = 3.5,
        dpi = 600,
        bg = "white"
      )

      cm_df <- build_preserved_confusion_df(
        platform_df,
        cohort_label = paste(cohort_label, platform_label, sep = " - "),
        model_names = selected_fragmentomics_models,
        model_labels = selected_fragmentomics_labels
      ) %>%
        dplyr::mutate(platform = platform_label)
      platform_split_confusion_source[[paste(cohort_slug, slug, sep = "_")]] <- cm_df

      cm_path <- file.path(
        "Final Tables and Figures",
        paste0("fragmentomics_", cohort_slug, "_confusion_tables_", slug, ".png")
      )
      ggsave(
        cm_path,
        plot = plot_preserved_confusion_df(
          cm_df,
          paste0("Confusion Matrix at Youden Index\n", cohort_label, " - ", platform_label)
        ),
        width = 5,
        height = 2.75,
        dpi = 600,
        bg = "white"
      )
    }
  }

  # Write consolidated source tables after all subsets complete successfully;
  # this prevents partial per-platform CSVs from looking like a full analysis.
  readr::write_csv(
    dplyr::bind_rows(platform_split_performance_source),
    file.path(outdir, "fragmentomics_platform_split_performance_source_data.csv")
  )
  readr::write_csv(
    dplyr::bind_rows(platform_split_confusion_source),
    file.path(outdir, "fragmentomics_platform_split_confusion_source_data.csv")
  )

  readr::write_csv(
    dplyr::bind_rows(edfigure9_training_confusion, edfigure9_test_confusion),
    file.path(outdir, "preserved_model_current_edfigure9_confusion_source_data.csv")
  )
  ggsave(
    "Final Tables and Figures/Supp_Fig9E_confusion_tables_primary_fragmentomics4.png",
    plot = plot_preserved_confusion_df(
      edfigure9_training_confusion,
      "Confusion Matrix at Youden Index in Training Cohort"
    ),
    width = 5, height = 2.75, dpi = 600, bg = "white"
  )
  ggsave(
    "Final Tables and Figures/Supp_Fig9F_confusion_tables_val_fragmentomics3.png",
    plot = plot_preserved_confusion_df(
      edfigure9_test_confusion,
      "Confusion Matrix at Youden Index in Test Cohort"
    ),
    width = 5, height = 2.75, dpi = 600, bg = "white"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Supp_Fig9E_confusion_tables_primary_fragmentomics4.png",
    artifact_id = "EDFIG9E", role = "figure_panel_png",
    description = "Fragmentomics training-cohort confusion matrix regenerated from preserved models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )
  ms_copy_artifact(
    source_path = "Final Tables and Figures/Supp_Fig9F_confusion_tables_val_fragmentomics3.png",
    artifact_id = "EDFIG9F", role = "figure_panel_png",
    description = "Fragmentomics test-cohort confusion matrix regenerated from preserved models.",
    script_name = "3_1_Optimize_cfWGS_thresholds.R"
  )

  message("Preserved-model cfWGS scoring complete. Wrote all_patients_with_BM_and_blood_calls_updated6*.rds/csv.")
  message("Current preserved-model testing metrics complete. Wrote Supplementary Table 6 source and audit outputs without model retraining.")
  quit(save = "no", status = 0)
}


# create train and test df - train on frontline, apply on hold out df
# frontline = training set, everything else = hold‑out
# DESIGN RATIONALE: Frontline patients are used as the TRAINING set because
# they have the largest number of matched BM/cfDNA time-points and represent
# the primary clinical scenario the assay targets (post-ASCT MRD assessment).
# Non-frontline patients are reserved as an independent HOLD-OUT set to test
# whether the thresholds generalise beyond the training cohort; they were NOT
# used to tune any hyperparameter in the model.
# DESIGN RATIONALE: Frontline patients form the TRAINING set because they have
# the largest number of matched BM/cfDNA time-points and represent the primary
# clinical scenario the assay targets (post-ASCT MRD monitoring). Non-frontline
# patients are held out as an independent TEST set to verify that thresholds
# generalise beyond the training cohort; they were never used to tune any model
# hyperparameter.
train_df <- dat_mrd %>% filter(Cohort == "Frontline")
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")


# ──────────────────────────────────────────────────────────────────────────────
# 2.  Helper: build sparse X & y
# DESIGN NOTE: The model matrix is converted to a sparse dgCMatrix before
# passing to glmnet. This matters because:
#   - The feature matrix can have >1000 columns after one-hot encoding
#   - Many patients have NA for some assay features (e.g. not all have
#     fragmentomics), making the effective rows sparse
#   - Sparse format saves memory and speeds up glmnet's coordinate descent
# The intercept column (column 1) is always dropped because glmnet adds its
# own intercept internally.
# ──────────────────────────────────────────────────────────────────────────────
make_sparse_xy <- function(df, predictors) {
  mm <- model.matrix(
    reformulate(predictors, response = "MRD_truth"),
    data = df
  )
  list(
    x = as(mm[, -1], "dgCMatrix"),  # drop intercept
    y = df$MRD_truth
  )
}

# -----------------------------------------------------------------------------
# 4. Exploratory univariate/ridge threshold screening
#
# Manuscript role:
#   Retained for provenance and model-development audit only. These exploratory
#   thresholds are not the manuscript-reported model results; the active
#   manuscript model performance comes from the nested cross-validation and
#   downstream scoring sections later in this script.
# -----------------------------------------------------------------------------

# 3.  Train ridge‐penalized combos ------------------------------------------
# 3a) Classic BM combo
bm_complete    <- train_df %>% drop_na(MRD_truth, zscore_BM, detect_rate_BM)
bm_xy          <- make_sparse_xy(bm_complete, c("zscore_BM", "detect_rate_BM"))
cv_bm          <- cv.glmnet(bm_xy$x, bm_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3b) Extended BM combo (+ FS, coverage, BM_tf)
bm_ext_complete <- train_df %>% 
  drop_na(MRD_truth, zscore_BM, detect_rate_BM, FS, Mean.Coverage, WGS_Tumor_Fraction_Blood_plasma_cfDNA)
bm_ext_xy       <- make_sparse_xy(
  bm_ext_complete,
  c("zscore_BM", "detect_rate_BM", "FS", "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)
cv_bm_ext      <- cv.glmnet(bm_ext_xy$x, bm_ext_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3c) Classic blood combo
bl_complete    <- train_df %>% drop_na(MRD_truth, zscore_blood, detect_rate_blood)
bl_xy          <- make_sparse_xy(bl_complete, c("zscore_blood", "detect_rate_blood"))
cv_bl          <- cv.glmnet(bl_xy$x, bl_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 3d) Extended blood combo (+ FS, coverage, blood_tf)
bl_ext_complete <- train_df %>%
  drop_na(MRD_truth, zscore_blood, detect_rate_blood, FS, Mean.Coverage, WGS_Tumor_Fraction_Blood_plasma_cfDNA)
bl_ext_xy       <- make_sparse_xy(
  bl_ext_complete,
  c("zscore_blood", "detect_rate_blood", "FS",
    "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)
cv_bl_ext      <- cv.glmnet(bl_ext_xy$x, bl_ext_xy$y, family="binomial", alpha=0, nfolds=5, nlambda=30)

# 4.  Add combo probabilities to all patients -------------------------------
add_prob <- function(df, predictors, fit) {
  # logical mask: rows with ALL predictors present
  keep <- stats::complete.cases(df[, predictors])
  
  # initialise output vector
  prob <- rep(NA_real_, nrow(df))
  
  if (any(keep)) {
    mm   <- model.matrix(reformulate(predictors), data = df[keep, ])
    newx <- as(mm[, -1], "dgCMatrix")
    prob[keep] <- as.numeric(
      predict(fit, newx = newx, s = "lambda.min", type = "response")
    )
  }
  prob
}


## Add it to the dataframe
dat <- dat %>%
  mutate(
    combo_BM_prob        = add_prob(., c("zscore_BM", "detect_rate_BM"),          cv_bm),
    combo_BM_ext_prob    = add_prob(., c("zscore_BM", "detect_rate_BM", "FS",
                                         "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"), cv_bm_ext),
    combo_blood_prob     = add_prob(., c("zscore_blood", "detect_rate_blood"),    cv_bl),
    combo_blood_ext_prob = add_prob(., c("zscore_blood", "detect_rate_blood", "FS",
                                         "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"), cv_bl_ext)
  )

# 5.  Define metric sets ----------------------------------------------------
base_metrics     <- c("zscore_BM","detect_rate_BM",
                      "zscore_blood","detect_rate_blood",
                      "combo_BM_prob","combo_blood_prob")
extended_metrics <- c(base_metrics,
                      "FS","Mean.Coverage",
                      "WGS_Tumor_Fraction_Blood_plasma_cfDNA",
                      "combo_BM_ext_prob","combo_blood_ext_prob")

# 6.  Optimiser function ----------------------------------------------------
optimise_metric <- function(df, metric, truth_col="MRD_truth",
                            trial = c("accuracy","sens>=0.9")) {
  trial <- match.arg(trial)
  tmp   <- df %>% select(all_of(c(metric, truth_col))) %>% drop_na()
  roc   <- pROC::roc(tmp[[truth_col]], tmp[[metric]])
  best  <- if (trial=="accuracy") {
    pROC::coords(roc, "best", best.method="youden",
                 ret=c("threshold","sensitivity","specificity","accuracy"),
                 transpose=FALSE)
  } else {
    pROC::coords(roc, "all",
                 ret=c("threshold","sensitivity","specificity","accuracy"),
                 transpose=FALSE) %>%
      filter(sensitivity>=0.9) %>%
      slice_max(specificity, n=1)
  }
  tibble(metric   = metric,
         trial    = trial,
         threshold= best[["threshold"]],
         sens     = best[["sensitivity"]],
         spec     = best[["specificity"]],
         accuracy = best[["accuracy"]],
         auc      = as.numeric(pROC::auc(roc)))
}

# -------------- Re-slice out train & hold‑out ----------------------------
train_df <- dat %>%
  filter(Cohort == "Frontline",
         !timepoint_info %in% c("Diagnosis", "Baseline"))   # keep MRD visits only

hold_df  <- dat %>%
  filter(Cohort == "Non-frontline",
         !timepoint_info %in% c("Diagnosis", "Baseline"))

run_trial <- function(metrics, trial_tag, df_train) {
  map_dfr(metrics, optimise_metric,
          df        = df_train,
          truth_col = "MRD_truth",
          trial     = trial_tag)
}

# 7.  Run trials on FRONTLINE training set ---------------------------------
results <- bind_rows(
  run_trial(base_metrics,     "accuracy",  train_df) %>% mutate(set="base"),
  run_trial(base_metrics,     "sens>=0.9", train_df) %>% mutate(set="base"),
  run_trial(extended_metrics, "accuracy",  train_df) %>% mutate(set="extended"),
  run_trial(extended_metrics, "sens>=0.9", train_df) %>% mutate(set="extended")
)
print(results)



# -----------------------------------------------------------------------------
# 5. Legacy single 5-fold CV model-screening block
#
# Manuscript role:
#   Retained for provenance because it records earlier candidate feature-set
#   screening. It is not the source of the final manuscript performance values;
#   the nested-CV section below is the active model-analysis section.
# -----------------------------------------------------------------------------

#### See if different combination of features would be better
# 1.  Define candidate combos ------------------------------------------------
combos <- list(
  # - BM combos with zscore -
  BM_zscore_only         = c("zscore_BM"),
  BM_base             = c("zscore_BM",    "detect_rate_BM"),
  BM_plus_FS          = c("zscore_BM",    "detect_rate_BM", "FS", "Proportion.Short"),
  BM_plus_MeanCov     = c("zscore_BM",    "detect_rate_BM", "Mean.Coverage"),
  BM_plus_TF          = c("zscore_BM",    "detect_rate_BM", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_all_extras       = c("zscore_BM",    "detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  # - BM combos without zscore -
  BM_rate_base        = c("detect_rate_BM"),
  BM_rate_plus_FS     = c("detect_rate_BM", "FS", "Proportion.Short"),
  BM_rate_plus_MeanCov= c("detect_rate_BM", "Mean.Coverage"),
  BM_rate_plus_TF     = c("detect_rate_BM", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_rate_all_extras  = c("detect_rate_BM", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),

  # - Blood combos with zscore -
  Blood_zscore_only      = c("zscore_blood"),
  Blood_base          = c("zscore_blood", "detect_rate_blood"),
  Blood_plus_FS       = c("zscore_blood", "detect_rate_blood", "FS", "Proportion.Short"),
  Blood_plus_MeanCov  = c("zscore_blood", "detect_rate_blood", "Mean.Coverage"),
  Blood_plus_TF       = c("zscore_blood", "detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_all_extras    = c("zscore_blood", "detect_rate_blood", "FS",
                          "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  # - Blood combos without zscore -
  Blood_rate_base     = c("detect_rate_blood"),
  Blood_rate_plus_FS  = c("detect_rate_blood", "FS"),
  Blood_rate_plus_MeanCov = c("detect_rate_blood", "Mean.Coverage"),
  Blood_rate_plus_TF  = c("detect_rate_blood",
                          "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_rate_all_extras   = c("detect_rate_blood", "FS",
                              "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  
  ### BM and blood together 
  BM_blood_ultimate = c("zscore_BM", "zscore_blood", "detect_rate_blood", "FS",
                        "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"), 
  
  ### Just with fragmentomics 
  Just_fragmentomics_full = c("FS", "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  Just_fragmentomics_small = c("FS", "Mean.Coverage","WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Just_fragmentomics_small2 = c("FS", "Mean.Coverage")
)


# 2.  Utility: evaluate one combo ------------------------------------------
eval_combo <- function(preds, label) {
  # complete-case slice
  df_cc <- train_df %>% drop_na(MRD_truth, all_of(preds))
  
  # skip if too small or only one class
  if (nrow(df_cc) < 20 || n_distinct(df_cc$MRD_truth) < 2)
    return(tibble())
  
  # --------------------------------------------------------------------
  # ➊  FIT MODEL  (glmnet if ≥2 preds, else glm) ------------------------
  # --------------------------------------------------------------------
  if (length(preds) >= 2) {
    # ridge CV
    mm  <- model.matrix(reformulate(preds, response = "MRD_truth"), df_cc)
    X   <- as(mm[, -1, drop = FALSE], "dgCMatrix")
    y   <- df_cc$MRD_truth
    
    fit <- cv.glmnet(X, y,
                     family  = "binomial",
                     alpha   = 0,
                     nfolds  = 5,
                     nlambda = 30)
    
    prob <- as.numeric(
      predict(fit, newx = X, s = "lambda.min", type = "response")
    )
  } else {
    # single predictor → plain logistic regression
    form <- reformulate(preds, response = "MRD_truth")
    fit  <- glm(form, data = df_cc, family = binomial)
    prob <- predict(fit, type = "response")
  }
  
  # --------------------------------------------------------------------
  # ➋  ROC & performance metrics ---------------------------------------
  # --------------------------------------------------------------------
  roc_obj <- pROC::roc(df_cc$MRD_truth, prob)
  
  # Accuracy‑optimal threshold
  acc_best <- pROC::coords(
    roc_obj, "best", best.method = "youden",
    ret = c("threshold","sensitivity","specificity","accuracy"),
    transpose = FALSE
  )
  # Highest specificity with sens ≥ 0.9
  sens90 <- pROC::coords(
    roc_obj, "all",
    ret = c("threshold","sensitivity","specificity","accuracy"),
    transpose = FALSE
  ) %>%
    filter(sensitivity >= 0.9) %>%
    slice_max(specificity, n = 1)
  
  bind_rows(
    tibble(combo = label, trial = "accuracy",  !!!acc_best),
    tibble(combo = label, trial = "sens>=0.9", !!!sens90)
  ) %>%
    mutate(auc = as.numeric(pROC::auc(roc_obj)))
}


# 3.  Run across all combos --------------------------------------------------
set.seed(2025)
combo_results <- imap_dfr(combos, eval_combo)

print(combo_results)

# Save as an RDS:
saveRDS(combo_results,
        file = file.path(outdir, "combo_results_updated.rds"))

# Export as CSV (no row names):
write.csv(combo_results,
          file = file.path(outdir, "combo_results_updated.csv"),
          row.names = FALSE)



#  Save fitted combo models for later scoring ─────────────────────────

# Collect the exact fits from the trained models
model_list <- list(
  BM_zscore_only   = cv_bm,       # classic BM combo
  BM_ext           = cv_bm_ext,   # extended BM combo
  Blood_base       = cv_bl,       # classic blood combo
  Blood_ext        = cv_bl_ext    # extended blood combo
)

# Write them to disk for downstream analysis in dilution-series evaluation
saveRDS(model_list,
        file = file.path(outdir, "models_cfWGS_final.rds"))





# -----------------------------------------------------------------------------
# 6. Repeated-CV model-screening block
#
# Manuscript role:
#   Retained to document the intermediate model-selection path. Final reported
#   values are taken from the later nested-CV/scoring/export sections that write
#   the manuscript figures and Supplementary Tables 4-6.
# -----------------------------------------------------------------------------

###### Now do a repeated 5x5 cross validation to improve generalizability 
###### Final numbers reported in the manuscript

# parallel backend (optional, speeds up train())
library(doParallel)
detected_cores <- parallel::detectCores(logical = TRUE)
parallel_workers <- if (is.na(detected_cores) || detected_cores < 2L) {
  1L
} else {
  max(1L, detected_cores - 1L)
}
registerDoParallel(cores = parallel_workers)
message(sprintf("Registered doParallel backend with %d worker(s)", parallel_workers))

# Training/validation splits
train_df <- dat_mrd %>% filter(Cohort == "Frontline") 
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")

train_df <- train_df %>% filter(!is.na(MRD_truth))
hold_df <- hold_df %>% filter(!is.na(MRD_truth))

# convert binary to factor with levels "neg"/"pos"
train_df <- train_df %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))
hold_df  <- hold_df  %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))

positive_class <- "pos"
# ──────────────────────────────────────────────────────────────────────────────
# 1.  caret::train control (5×5 repeated CV) ------------------------------------
# ──────────────────────────────────────────────────────────────────────────────
ctrl <- trainControl(
  method            = "repeatedcv",
  number            = 5,
  repeats           = 5,
  classProbs        = TRUE,
  summaryFunction   = twoClassSummary,
  savePredictions   = "final",   # store OOF preds
  allowParallel     = TRUE
)

# ──────────────────────────────────────────────────────────────────────────────
# 2.  Fit each combo, get CV + hold‐out metrics -------------------------------
# ──────────────────────────────────────────────────────────────────────────────
## Helper function 
youden_metrics <- function(response, prob, thr) {
  roc_obj <- roc(response, prob, quiet = TRUE, levels = c("neg", "pos"))
  tibble(
    auc = as.numeric(auc(roc_obj)),
    sens = coords(roc_obj, thr, ret = "sensitivity"),
    spec = coords(roc_obj, thr, ret = "specificity"),
    acc  = mean((prob >= thr) == (response == positive_class))
  )
}

## Run on cohort 
results <- imap(combos, function(preds, label) {
  
  # complete cases in train
  df_cc <- train_df %>% drop_na(all_of(preds), MRD_truth)
  if(nrow(df_cc) < 20 || n_distinct(df_cc$MRD_truth) < 2) {
    message("Skipping ", label, " (n=", nrow(df_cc), ")")
    return(NULL)
  }
  
  # fit the model (glm or glmnet)
  form <- reformulate(preds, "MRD_truth")
  fit  <- suppressWarnings(
    if (length(preds) == 1) {
      train(form, df_cc, method="glm", family=binomial,
            metric="ROC", trControl=ctrl, preProc=c("center","scale"))
    } else {
      train(form, df_cc, method="glmnet", metric="ROC",
            tuneLength=10, trControl=ctrl, preProc=c("center","scale"))
    }
  )
  
  # pull out OOF preds
  oof <- if (length(preds)==1) fit$pred else {
    fit$pred %>%
      filter(alpha==fit$bestTune$alpha,
             lambda==fit$bestTune$lambda)
  }
  roc_oof <- roc(oof$obs, oof[[positive_class]], quiet=TRUE,
                 levels=c("neg","pos"))
  
  # helper to pick threshold by criterion
  pick_thr <- function(mode){
    df <- coords(
      roc_oof, "all",
      ret = c("threshold","sensitivity","specificity","accuracy"),
      transpose = FALSE
    ) %>% as_tibble()
    if(mode=="youden"){
      coords(
        roc_oof, "best", best.method="youden",
        ret="threshold", transpose=FALSE
      )[[1]]
    } else if(mode=="sens90"){
      df %>%
        filter(sensitivity >= .9) %>%
        slice_max(specificity, n = 1, with_ties = FALSE) %>%
        pull(threshold)
    } else {
      df %>%
        filter(specificity >= .9) %>%
        slice_max(sensitivity, n = 1, with_ties = FALSE) %>%
        pull(threshold)
    }
  }
  
  # now build a numeric vector of exactly three thresholds:
  thr_list <- map_dbl(c("youden","sens90","spec90"), pick_thr)
  
  
  # build a little function to compute metrics at a given thr
  eval_at <- function(roc_obj, response, prob, thr){
    coords(roc_obj, x=thr,
           ret=c("sensitivity","specificity","accuracy"),
           transpose=FALSE) %>%
      as_tibble() %>%
      mutate(
        auc = as.numeric(auc(roc_obj)),
        thr = thr
      )
  }
  
  # 1) CV metrics for each mode
  cv_metrics <- map2_dfr(
    thr_list,
    c("youden","sens>=0.9","spec>=0.9"),
    ~ eval_at(roc_oof, oof$obs, oof[[positive_class]], .x) %>%
      mutate(
        combo = label,    # <- the feature-combo name
        mode  = .y,       # <- which thresholding mode
        .before = 1
      )
  )
  
  # now do hold-out
  hold_cc    <- hold_df %>% drop_na(all_of(preds), MRD_truth)
  if (nrow(hold_cc) < 2L || n_distinct(hold_cc$MRD_truth) < 2L) {
    message("Skipping hold-out metrics for ", label, " (n=", nrow(hold_cc),
            ", classes=", n_distinct(hold_cc$MRD_truth), ")")
    hold_metrics <- map2_dfr(
      thr_list,
      c("youden","sens>=0.9","spec>=0.9"),
      ~ tibble(
        sensitivity = NA_real_,
        specificity = NA_real_,
        accuracy = NA_real_,
        auc = NA_real_,
        thr = .x,
        combo = label,
        mode = .y,
        .before = 1
      )
    )
  } else {
    probs_hold <- predict(fit, newdata=hold_cc, type="prob")[[positive_class]]
    roc_hold   <- roc(hold_cc$MRD_truth, probs_hold, quiet=TRUE,
                      levels=c("neg","pos"))

    hold_metrics <- map2_dfr(
      thr_list,
      c("youden","sens>=0.9","spec>=0.9"),
      ~ eval_at(roc_hold, hold_cc$MRD_truth, probs_hold, .x) %>%
        mutate(
          combo = label,    # <- again, the feature-combo name
          mode  = .y,
          .before = 1
        )
    )
  }
  
  # return
  list(cv   = cv_metrics,
       hold = hold_metrics,
       fit  = fit                   # <<< return the caret model object
  )
})


# bind into tables
nested_metrics <- map_dfr(results, "cv")
hold_metrics   <- map_dfr(results, "hold")

# Save for manuscript ------------------------------------------------------
write_csv(nested_metrics, "nested_metrics.csv")
write_csv(hold_metrics,   "hold_metrics.csv")


### Now select best one and apply to everything 
# ────────────────────────────────────────────────────────────────────────────
# 4.  SELECT THE BEST PERFORMING MODELS ----------------------------------------
# Helper to pull top‑n accuracy rows for a given combo
pick_rows <- function(cmb, n = 1) {
  nested_metrics %>%
    filter(combo == cmb) %>%
    slice_max(accuracy, n = n, with_ties = FALSE)
}

bm_row        <- pick_rows("BM_zscore_only", n = 1)
blood_rows    <- pick_rows("Blood_rate_all_extras", n = 1)
blood_row_add <- pick_rows("Blood_base", n = 1)
selected_rows <- bind_rows(bm_row, blood_rows, blood_row_add)   # 4 rows in total
print(selected_rows)


# ────────────────────────────────────────────────────────────────────────────
# 4b.  Fit & save only the selected combo models ----------------------------

# Extract the selected combo names
selected_combos <- selected_rows$combo

# Get the models from before
selected_models <- map(results[selected_combos], "fit")

# Save both the models and the thresholds
saveRDS(selected_models,
        file = file.path(outdir, "selected_combo_models_updated.rds"))
saveRDS(selected_rows,
        file = file.path(outdir, "selected_combo_thresholds_updated.rds"))



# ────────────────────────────────────────────────────────────────────────────
# 5. NOW APPLY TO REST OF COHORT  --------------------------------------------

# Define function 
apply_selected <- function(row, idx, dat, combos, models, positive_class){
  cmb   <- row$combo
  thr   <- row$threshold      # ensure selected_rows contains this column
  preds <- combos[[cmb]]
  fit   <- models[[cmb]]
  
  # build a prediction function
  if (inherits(fit, "cv.glmnet")) {
    predict_fun <- function(df_sub) {
      mm   <- model.matrix(reformulate(preds), data = df_sub)
      Xnew <- as(mm[,-1, drop=FALSE], "dgCMatrix")
      as.numeric(predict(fit, newx = Xnew, s="lambda.min", type="response"))
    }
  } else if (inherits(fit, "glm")) {
    predict_fun <- function(df_sub) {
      predict(fit, newdata = df_sub, type="response")
    }
  } else if (inherits(fit, "train")) {
    predict_fun <- function(df_sub) {
      prob_df <- predict(fit, newdata = df_sub, type="prob")
      prob_df[[ positive_class ]]
    }
  } else {
    stop("Unknown model type for combo: ", cmb)
  }
  
  # apply to all rows with non-missing preds
  keep <- complete.cases(dat[, preds, drop=FALSE])
  prob <- rep(NA_real_, nrow(dat))
  prob[keep] <- predict_fun(dat[keep, ])
  
  # optional suffix if multiple instances of the same combo exist
  suffix <- if(cmb=="Blood_all_extras") paste0("_acc", idx) else ""
  pcol   <- paste0(cmb, suffix, "_prob")
  ccol   <- paste0(cmb, suffix, "_call")
  
  dat[[pcol]] <- prob
  dat[[ccol]] <- if_else(prob >= thr, 1L, 0L, NA_integer_)
  
  dat
}

selected_rows <- selected_rows %>%
  mutate(threshold = thr)    # or rename(threshold = threshold_cv)

# Loop through all selected rows
for(i in seq_len(nrow(selected_rows))) {
  dat <- apply_selected(
    selected_rows[i, ], 
    i, 
    dat, 
    combos, 
    selected_models, 
    positive_class
  )
}





### Historical model-development and full-training refit section
# -----------------------------------------------------------------------------
# 7. Elastic-Net Classifier Training & Historical Nested CV
# -----------------------------------------------------------------------------

### One five-fold outer split with five-fold CV repeated five times inside.
### The current 50-repeat patient-grouped analysis is in script 6_12.
positive_class <- "pos"

### Since doing hyperparameter tuning, can use less features 
combos_small <- list(
  # Bone marrow signal comparisons
  BM_zscore_only_sites     = c("zscore_BM"),
  BM_zscore_only_detection_rate     = c("z_score_detection_rate_BM"),
  BM_rate_only = c("detect_rate_BM"),
  BM_base            = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM"),
  BM_base_zscore            = c("zscore_BM", "z_score_detection_rate_BM"),
  BM_plus_fragment   = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM", "FS", "Mean.Coverage", "Proportion.Short", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  BM_plus_fragment_min   = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM", "FS", "Mean.Coverage"),
  
  # Blood signal comparisons
  Blood_zscore_only_sites  = c("zscore_blood"),
  Blood_zscore_only_detection_rate  = c("z_score_detection_rate_blood"),
  Blood_rate_only = c("detect_rate_blood"),
  Blood_base         = c("zscore_blood", "detect_rate_blood", "z_score_detection_rate_blood"),
  Blood_base_zscore         = c("zscore_blood", "z_score_detection_rate_blood"),
  Blood_plus_fragment = c("zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood", "FS", "Mean.Coverage", "Proportion.Short", "WGS_Tumor_Fraction_Blood_plasma_cfDNA"),
  Blood_plus_fragment_min = c("zscore_blood", "z_score_detection_rate_blood", "detect_rate_blood", "FS", "Mean.Coverage"),
  
  # Combined signals
  BM_blood_combo     = c("zscore_BM", "zscore_blood", "detect_rate_blood", "FS", "Mean.Coverage"),
  
  # Fragmentomics-only
  Fragmentomics_full = c("FS", "Mean.Coverage", "WGS_Tumor_Fraction_Blood_plasma_cfDNA", "Proportion.Short"),
  Fragmentomics_min  = c("FS", "Mean.Coverage"),
  Fragmentomics_FS_only = c("FS"), 
  Fragmentomics_mean_coverage_only = c("Mean.Coverage"), 
  Fragmentomics_prop_short_only = c("Proportion.Short"), 
  Fragmentomics_tumor_fraction_only = c("WGS_Tumor_Fraction_Blood_plasma_cfDNA")
)

# 1) Identify which predictors belong to BM vs Blood combos
# BM combos: anything starting “BM_” but drop “BM_blood_combo”
bm_names <- names(combos_small)[
  startsWith(names(combos_small), "BM_") &
    names(combos_small) != "BM_blood_combo"
]
combos_bm    <- combos_small[bm_names]

# Blood combos: anything starting “Blood_” OR containing “Just_fragmentomics”
blood_names <- names(combos_small)[
  startsWith(names(combos_small), "Blood_") |
    grepl("Fragmentomics_", names(combos_small))
]
combos_blood <- combos_small[blood_names]

# Combos just fragmentomics - train seperately to get bigger cohort 
fragmentomics_names <- names(combos_small)[
    grepl("Fragmentomics_", names(combos_small))
]
combos_fragmentomics <- combos_small[fragmentomics_names]


# 2) Build train_bm / train_blood by dropping any row missing *any* BM or Blood predictors
bm_preds    <- unique(unlist(combos_bm))
blood_preds <- unique(unlist(combos_blood))

## Add filter that evidence of disease at baseline 
# 2) Which cohort patients have a “good” baseline BM sample?
Good_pts <- read.csv("baseline_high_quality_patients_updated.csv",
                     stringsAsFactors = FALSE) ## From 2_0 script

bm_good_patients <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_BM_cells == 1) %>%
  select(Patient) %>%         # keep as a tibble with a Patient column
  distinct()                  # remove duplicates

cfDNA_good_patients <- Good_pts %>%
  filter(WGS_Evidence_of_Disease_Blood_plasma_cfDNA_Relaxed == 1) %>%
  select(Patient) %>%         # keep as a tibble with a Patient column
  distinct()                  # remove duplicates

# train_bm_original <- train_bm
# train_blood_original <- train_blood
# hold_bm_original <- hold_bm
# hold_blood_original <- hold_blood
# Training/validation splits
train_df <- dat_mrd %>% filter(Cohort == "Frontline") 
hold_df  <- dat_mrd %>% filter(Cohort == "Non-frontline")

train_df <- train_df %>% filter(!is.na(MRD_truth))
hold_df <- hold_df %>% filter(!is.na(MRD_truth))

# convert binary to factor with levels "neg"/"pos"
train_df <- train_df %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))
hold_df  <- hold_df  %>% mutate(MRD_truth = factor(MRD_truth, levels = c(0,1), labels=c("neg","pos")))

positive_class <- "pos"
train_bm    <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds))) %>% filter(Patient %in% bm_good_patients$Patient)
train_blood <- train_df %>% drop_na(all_of(c("MRD_truth", blood_preds))) %>% filter(Patient %in% cfDNA_good_patients$Patient)

hold_bm    <- hold_df %>% drop_na(all_of(c("MRD_truth", bm_preds)))  %>% filter(Patient %in% bm_good_patients$Patient)
hold_blood <- hold_df %>% drop_na(all_of(c("MRD_truth", blood_preds))) %>% filter(Patient %in% cfDNA_good_patients$Patient)


## Get more info on dates for manuscript 
Baseline_dates <- read_csv(file = "Final Tables and Figures/Baseline dates for samples.csv")

make_sentence <- function(df, cohort_label) {
  df %>%
    left_join(Baseline_dates, by = c("Patient" = "patient")) %>%
    mutate(months_from_dx = interval(start, Date) / months(1)) %>%
    summarise(
      n_samples  = n(),
      n_patients = n_distinct(Patient),
      min_m      = round(min(months_from_dx), 1),
      max_m      = round(max(months_from_dx), 1),
      median_m   = round(median(months_from_dx), 1)
    ) %>%
    mutate(
      text = paste0(
        n_samples, " samples from ", n_patients,
        " patients were available with matched clinical MRD in the ",
        cohort_label, " cohort, collected over ",
        min_m, "-", max_m, " months following diagnosis (median ",
        median_m, ")."
      )
    ) %>%
    pull(text)
}

cat(make_sentence(train_bm,  "training (bone marrow)"), "\n")
cat(make_sentence(hold_bm,   "test (bone marrow)"),    "\n")
cat(make_sentence(train_blood,"training (blood)"),     "\n")
cat(make_sentence(hold_blood, "test (blood)"),         "\n")

### see what has blood but not BM 
# Restrict to frontline samples in each training set
fb <- train_blood %>% filter(Cohort == "Frontline")
fbm <- train_bm    %>% filter(Cohort == "Frontline")

# 1) Identify blood samples with no BM call
no_bm <- fb %>%
  # drop any that have a matching (Patient, Timepoint) in the BM table
  anti_join(fbm, by = c("Patient", "Timepoint")) %>%
  # keep only non‐baseline/diagnosis timepoints
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 


n_samples_no_bm  <- nrow(no_bm)
n_patients_no_bm <- no_bm %>% distinct(Patient) %>% nrow()

# Summarize
tibble(
  samples_without_bm   = n_samples_no_bm,
  patients_without_bm  = n_patients_no_bm
)

## Get reason
# `1_6_Identify_High_Quality_Patient_Pairs.R` writes the Figure 1B source
# table to Final Tables and Figures/. Keep a legacy root-file fallback for older
# local runs, but the organized path is the pipeline source going forward.
fig1b_source_table <- file.path("Final Tables and Figures", "Table for creating sample flowchart updated3.csv")
if (!file.exists(fig1b_source_table)) {
  fig1b_source_table <- "Table for creating sample flowchart updated3.csv"
}
failures <- read_csv(fig1b_source_table)
no_bm <- left_join(no_bm)

patients_no_bm <- no_bm %>%
  select(Patient) %>%
  distinct()

failures_joined <- patients_no_bm %>%
  left_join(failures, by = "Patient")

failures_joined

## Get total number of samples
dat_mrd %>%
  filter(Cohort == "Frontline") %>%
  drop_na(all_of(c(blood_preds))) %>%
  filter(Patient %in% cfDNA_good_patients$Patient) %>%
  summarise(
    total_patients = n_distinct(Patient),
    total_samples  = dplyr::n()
  )


# 3) Build train/hold sets for fragmentomics-only models
frag_preds       <- unique(unlist(combos_fragmentomics))

train_fragmentomics <- train_df %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

hold_fragmentomics  <- hold_df %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

# =============================================================================
# SECTION 7: HISTORICAL SAMPLE-ROW NESTED CROSS-VALIDATION
# =============================================================================
# 
# WHAT:    Earlier nested-CV function that trains elastic-net models and
#          evaluates them on held-out sample rows.
#
# HOW IT WORKS:
#   1. Outer loop: Split training rows once into five folds, stratified by
#      MRD_truth. Repeated samples from one patient are not grouped here.
#   2. Inner loop: For each outer fold, use caret::train() with 5x5 repeated CV
#      to optimize elastic-net hyperparameters (alpha, lambda).
#   3. Final model: Train on entire train_data with optimized parameters,
#      then score on both validation data (held-out cohort) and outer-fold preds.
#
# OUTPUTS: List containing:
#   - nested_metrics: Per-fold CV performance (5 rows x ~10 metrics)
#   - models: Fitted caret::train objects (one per combo, trained on full data)
#   - validation_metrics: Test-set performance (1 row per combo)
#   - outer_predictions: OOF predictions for outer folds (unbiased for AUC/calibration)
#   - thresholds: Youden-optimal threshold per combo. NOTE (2026-07-29 audit):
#       these are computed on the FULL-TRAINING REFIT, not from the OOF
#       predictions (see the `roc_full`/`youden_full` block below). They are
#       therefore apparent (resubstitution) operating points derived from
#       training data only. This is safe with respect to the test cohort - no
#       test data is involved - but the values are optimistic as descriptions of
#       training-set performance, and the Methods now say so. The comment
#       previously claimed these came from OOF predictions.
#   - fold_indices: The outer fold structure (useful for reusing folds across models)
#
# CURRENT MANUSCRIPT ANALYSIS:
#   Script 6_12 repeats the outer five-fold split 50 times and groups both outer
#   and inner folds by patient. Scripts 6_13 through 6_20 assemble its 32-model
#   table and current figure panels. The code below remains for the original
#   model-development record and preserved full-training model objects.
#
# =============================================================================

## Define the main nested CV function
# Outer 5-fold nested CV
run_nested_with_validation <- function(train_data,
                                       valid_data,           # ◀ NEW: hold-out validation cohort
                                       combo_list,
                                       positive_class = "pos",
                                       fold_indices = NULL) {  # ◀ OPTIONAL: pass pre-defined fold indices for exact reproducibility
  # 0) Ensure `ctrl` (trainControl) is defined in the parent environment
  #    exactly as before (5×5 repeated CV, classProbs=TRUE, summaryFunction=twoClassSummary, etc.)
  
  # 1) Outer 5‐fold indices (stratified on MRD_truth)
  # If fold_indices provided, use those; otherwise create fresh ones
  if (is.null(fold_indices)) {
    outer_folds <- createFolds(train_data$MRD_truth,
                               k = 5,
                               returnTrain = TRUE)
  } else {
    outer_folds <- fold_indices
    message("Using pre-defined fold indices for exact reproducibility")
  }
  
  outer_preds <- list()    # one element per combo
  
  # 2) Loop over each feature-combo
  results <- purrr::imap(combo_list, function(preds, label) {
    
    combo_fold_preds <- vector("list", length(outer_folds))  # collect per fold
    
    # 2a) Collect per‐fold metrics
    combo_metrics <- purrr::map_dfr(seq_along(outer_folds), function(i) {
      # split out vs. in
      train_out <- train_data[ outer_folds[[i]], ]
      test_out  <- train_data[-outer_folds[[i]], ]
      
      # complete‐case
      train_cc <- train_out %>% drop_na(MRD_truth, all_of(preds))
      test_cc  <- test_out  %>% drop_na(MRD_truth, all_of(preds))
      
      # Export diagnostics
      message(sprintf("[%s | fold %d] n_train=%d  n_test=%d  classes(train)=%d  classes(test)=%d",
                      label, i,
                      nrow(train_cc), nrow(test_cc),
                      n_distinct(train_cc$MRD_truth),
                      n_distinct(test_cc$MRD_truth)))
      
      # SKIP if either side is too small or single‐class
    #  if ( # nrow(train_cc) < 20 || removed this criteria
    #      n_distinct(train_cc$MRD_truth) < 2 ||
    #      n_distinct(test_cc$MRD_truth)  < 2) {
    #    return(NULL)
    #  }
      
      if (n_distinct(train_cc$MRD_truth) < 2 ||
          n_distinct(test_cc$MRD_truth)  < 2) {
        message("‑‑ fold dropped: missing class after drop_na()")
        return(NULL)
      }
      
      message("‑‑ fold kept: fitting model now")
      
      # build formula
      f <- reformulate(preds, response = "MRD_truth")
      
      # 2b) inner CV via caret::train
      fit <- if (length(preds) == 1) {
        train(f, data    = train_cc,
              method  = "glm",
              family  = binomial,
              metric  = "ROC",
              trControl = ctrl,
              preProc = c("center","scale"))
      } else {
        train(f, data    = train_cc,
              method    = "glmnet",
              metric    = "ROC",
              # ◀ EXPANDED GRID over alpha ∈ {0,0.25,0.5,0.75,1} and 50 log‐lambda values
              tuneGrid  = expand.grid(
                alpha  = seq(0,1, by=0.25),
                lambda = 10^seq(-3, 1, length.out = 50)
              ),
              trControl = ctrl,
              preProc   = c("center","scale"))
      }
      
      # 2c) predict on outer test
      probs    <- predict(fit, newdata = test_cc, type = "prob")[[positive_class]]
      message("    --> length(probs) = ", length(probs))
      
      roc_o    <- pROC::roc(test_cc$MRD_truth, probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
      youden_t <- unlist(
        pROC::coords(roc_o, "best",
                     best.method = "youden",
                     ret         = "threshold")
      )
      youden_t <- as.numeric(youden_t[1])
      
      ## --- store the raw predictions for later plotting ----
      combo_fold_preds[[i]] <<- tibble(
        combo = label,
        fold  = i,
        truth = test_cc$MRD_truth,
        prob  = probs
      )
      
      # pull out single‐value sens/spec
      sens_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "sensitivity",
                             transpose = FALSE)[[1]]
      spec_o <- pROC::coords(roc_o, x = youden_t,
                             ret       = "specificity",
                             transpose = FALSE)[[1]]
      
      # sensitivity at 95% specificity
      sens95_o <- pROC::coords(roc_o,
                               x      = 0.95,
                               input  = "specificity",
                               ret    = "sensitivity",
                               transpose = FALSE)[[1]]
      
      # specificity at 95% sensitivity
      spec95_o <- pROC::coords(roc_o,
                               x      = 0.95,
                               input  = "sensitivity",
                               ret    = "specificity",
                               transpose = FALSE)[[1]]
      
      # PPV / NPV at Youden threshold
      ppv_o    <- pROC::coords(roc_o,
                               x      = youden_t,
                               input  = "threshold",
                               ret    = "ppv",
                               transpose = FALSE)[[1]]
      npv_o    <- pROC::coords(roc_o,
                               x      = youden_t,
                               input  = "threshold",
                               ret    = "npv",
                               transpose = FALSE)[[1]]
      
      # F1 and balanced accuracy
      f1_o     <- 2 * ppv_o * sens_o / (ppv_o + sens_o)
      balacc_o <- (sens_o + spec_o) / 2
      
      # Brier score
      brier_o  <- mean((probs - as.numeric(test_cc$MRD_truth == positive_class))^2)
      
      # partial AUC over specificity ≥ 0.9
      pauc90_o <- as.numeric(
        pROC::auc(roc_o,
                  partial.auc = c(0.9, 1),
                  partial.auc.correct = TRUE)
      )
      
      # fold sizes
      n_train_i <- nrow(train_cc)
      n_test_i  <- nrow(test_cc)
      
      tibble(
        combo       = label,
        auc         = as.numeric( pROC::auc(roc_o) ),
        sensitivity = sens_o,
        specificity = spec_o,
        sens_at_95_spec       = sens95_o,
        spec_at_95_sen       = spec95_o,
        bal_accuracy  = balacc_o,
        ppv           = ppv_o,
        npv           = npv_o,
        f1            = f1_o,
        brier         = brier_o,
        pAUC90        = pauc90_o,
        n_train       = n_train_i,
        n_test        = n_test_i,
        accuracy    = mean((probs >= youden_t) ==
                             (test_cc$MRD_truth == positive_class))
      )
    })
    
    # --- NEW: bind the fold‑level predictions for this combo
    preds_df <- bind_rows(combo_fold_preds)
    message("   --> preds_df rows for ", label, " = ", nrow(preds_df))
    
    # 2d) summarize across *kept* folds
    nested_summary <- combo_metrics %>%
      summarise(
        combo            = first(combo),
        n_folds          = dplyr::n(),
        auc_mean         = mean(auc,           na.rm = TRUE),
        auc_sd           = sd(auc,             na.rm = TRUE),
        sens_mean        = mean(sensitivity,   na.rm = TRUE),
        sens_sd          = sd(sensitivity,     na.rm = TRUE),
        spec_mean        = mean(specificity,   na.rm = TRUE),
        spec_sd          = sd(specificity,     na.rm = TRUE),
        balacc_mean      = mean(bal_accuracy,   na.rm = TRUE),
        balacc_sd        = sd(bal_accuracy,     na.rm = TRUE),
        ppv_mean         = mean(ppv,            na.rm = TRUE),
        npv_mean         = mean(npv,            na.rm = TRUE),
        f1_mean          = mean(f1,             na.rm = TRUE),
        brier_mean       = mean(brier,          na.rm = TRUE),
        pAUC90_mean      = mean(pAUC90,         na.rm = TRUE),
        sens95_mean      = mean(sens_at_95_spec,         na.rm = TRUE),
        spec95_mean      = mean(spec_at_95_sen,         na.rm = TRUE),
        mean_n_train     = mean(n_train,        na.rm = TRUE),
        mean_n_test      = mean(n_test,         na.rm = TRUE),
        acc_mean         = mean(accuracy,       na.rm = TRUE),
        .groups          = "drop"
      )
    
    # 3) Re‐fit on the *entire* train_data
    full_cc   <- train_data %>% drop_na(MRD_truth, all_of(preds))
    f_full    <- reformulate(preds, "MRD_truth")
    final_fit <- if (length(preds) == 1) {
      train(f_full, data    = full_cc,
            method  = "glm", family = binomial,
            metric  = "ROC", trControl = ctrl,
            preProc = c("center","scale"))
    } else {
      train(f_full, data    = full_cc,
            method    = "glmnet",
            metric    = "ROC",
            tuneGrid  = expand.grid(
              alpha  = seq(0,1, by=0.25),
              lambda = 10^seq(-3, 1, length.out = 10)
            ),
            trControl = ctrl,
            preProc   = c("center","scale"))
    }
    
    # 4) Compute Youden threshold on full‐train ROC
    full_probs <- predict(final_fit, newdata = full_cc, type="prob")[[positive_class]]
    roc_full   <- pROC::roc(full_cc$MRD_truth, full_probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
    youden_full <- unlist(
      pROC::coords(roc_full, "best",
                   best.method = "youden",
                   ret         = "threshold")
    )
    youden_full <- as.numeric(youden_full[1])
    
    # keep that threshold so we can reuse it later 
    threshold_full <- youden_full
    
    # 5) EVALUATE on the *external* validation set
    valid_cc   <- valid_data %>% drop_na(MRD_truth, all_of(preds))
    valid_probs<- predict(final_fit, newdata = valid_cc, type="prob")[[positive_class]]
    roc_val    <- pROC::roc(valid_cc$MRD_truth, valid_probs,
                            quiet  = TRUE,
                            levels = c("neg","pos"))
    sens_val   <- pROC::coords(roc_val, x = youden_full,
                               ret = "sensitivity")[[1]]
    spec_val   <- pROC::coords(roc_val, x = youden_full,
                               ret = "specificity")[[1]]
    acc_val    <- mean((valid_probs >= youden_full) ==
                         (valid_cc$MRD_truth == positive_class))
    
    # sensitivity/spec at 95% on the hold‐out
    sens95_val <- pROC::coords(roc_val,
                               x      = 0.95,
                               input  = "specificity",
                               ret    = "sensitivity",
                               transpose = FALSE)[1]
    
    spec95_val <- pROC::coords(roc_val,
                               x      = 0.95,
                               input  = "sensitivity",
                               ret    = "specificity",
                               transpose = FALSE)[1]
    
    # PPV / NPV at the full-train Youden cutoff
    ppv_val    <- pROC::coords(roc_val,
                               x      = youden_full,
                               input  = "threshold",
                               ret    = "ppv",
                               transpose = FALSE)[1]
    npv_val    <- pROC::coords(roc_val,
                               x      = youden_full,
                               input  = "threshold",
                               ret    = "npv",
                               transpose = FALSE)[1]
    
    # F1 and balanced accuracy
    f1_val     <- 2 * ppv_val * sens_val / (ppv_val + sens_val)
    balacc_val <- (sens_val + spec_val) / 2
    
    # Brier score
    brier_val  <- mean((valid_probs - as.numeric(valid_cc$MRD_truth == positive_class))^2)
    
    # Partial AUC over specificity ≥ 0.9
    pauc90_val <- as.numeric(
      pROC::auc(roc_val,
                partial.auc           = c(0.9, 1),
                partial.auc.correct   = TRUE)
    )
    
    
    validation_summary <- tibble(
      combo       = label,
      auc_valid   = as.numeric( pROC::auc(roc_val) ),
      sens_valid  = sens_val,
      spec_valid  = spec_val,
      sens_at_95_spec_valid  = sens95_val,
      spec_at_95_sens_valid  = spec95_val,
      bal_accuracy   = balacc_val,
      ppv            = ppv_val,
      npv            = npv_val,
      f1             = f1_val,
      brier          = brier_val,
      pAUC90         = pauc90_val,
      accuracy   = acc_val
    )
    
    # return all pieces
    list(
      nested     = nested_summary,
      final_fit  = final_fit,
      bestTune   = final_fit$bestTune,
      validation = validation_summary,
      threshold  = threshold_full,
      preds_df   = preds_df           # <‑‑  NEW: stash the preds in the result
    )
  })
  
  # pull them back out
  nested_metrics     <- map_dfr(results, "nested")
  models             <- map(results, "final_fit")
  best_tunes         <- map(results, "bestTune")
  validation_metrics <- map_dfr(results, "validation")
  thresholds         <- map_dbl(results, "threshold") # collect thresholds
  outer_predictions  <- map_dfr(results, "preds_df") # ---  NEW: bind all combos’ outer‑fold predictions
  
  message("Rows in outer_predictions = ", nrow(outer_predictions))
  
  list(
    nested_metrics     = nested_metrics,
    models             = models,
    best_tunes         = best_tunes,
    validation_metrics = validation_metrics,
    thresholds         = thresholds,
    outer_predictions  = outer_predictions,   ### NEW -- outer predictions are here
    fold_indices       = outer_folds          ### NEW -- return the fold indices for reuse
  )
}


# 4) Run nested CV separately
# OPTION 1: Run BM and blood models as originally planned (no changes here)
set.seed(2025)
nested_bm_validation_updated3 <- run_nested_with_validation(
  train_data = train_bm,
  valid_data = hold_bm,
  combo_list = combos_bm,
  positive_class = "pos"
)

set.seed(2025)
nested_blood_validation_updated3 <- run_nested_with_validation(
  train_data = train_blood,
  valid_data = hold_blood,
  combo_list = combos_blood,
  positive_class = "pos"
)

# OPTION 1 ORIGINAL: Keep the original fragmentomics run on full cohort (unchanged)
set.seed(2025)
nested_fragmentomics_validation_updated3 <- run_nested_with_validation(
  train_data    = train_fragmentomics,
  valid_data    = hold_fragmentomics,
  combo_list    = combos_fragmentomics,
  positive_class = "pos"
)

# OPTION 1 NEW: ALSO run FRAGMENTOMICS on the SAME restricted cohorts with EXACT SAME FOLDS
# This allows comparison of fragmentomics on the restricted evidence-of-disease cohorts
# Extract the fold indices from the original runs
frag_preds <- unique(unlist(combos_fragmentomics))

# Filter BM and blood datasets to keep only rows with fragmentomics data
train_bm_with_frag <- train_bm %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

hold_bm_with_frag <- hold_bm %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

train_blood_with_frag <- train_blood %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

hold_blood_with_frag <- hold_blood %>%
  drop_na(all_of(c("MRD_truth", frag_preds)))

# Now run fragmentomics on BM cohort using the EXACT SAME FOLDS as the original BM run
set.seed(2025)
nested_fragmentomics_bm_validation_updated3 <- run_nested_with_validation(
  train_data = train_bm_with_frag,
  valid_data = hold_bm_with_frag,
  combo_list = combos_fragmentomics,
  positive_class = "pos",
  fold_indices = nested_bm_validation_updated3$fold_indices  # ◀ USE EXACT SAME FOLDS
)

# And run fragmentomics on blood cohort using the EXACT SAME FOLDS as the original blood run
set.seed(2025)
nested_fragmentomics_blood_validation_updated3 <- run_nested_with_validation(
  train_data = train_blood_with_frag,
  valid_data = hold_blood_with_frag,
  combo_list = combos_fragmentomics,
  positive_class = "pos",
  fold_indices = nested_blood_validation_updated3$fold_indices  # ◀ USE EXACT SAME FOLDS
)

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTANT: Add fold_indices to existing BM/blood results from previous runs
# ─────────────────────────────────────────────────────────────────────────────
# Since the function was updated to return fold_indices, previously captured
# results need fold indices re-captured from the nested CV process
# We'll re-run BM and blood with seed(2025) - since randomness is seeded identically,
# results will be IDENTICAL but now we'll have fold_indices captured
# This is safe: same seed + same data = identical results + new fold_indices

message("Capturing fold_indices from BM and blood models...")

# Re-run BM with same seed and parameters (will produce identical results)
set.seed(2025)
nested_bm_validation_updated3_with_folds <- run_nested_with_validation(
  train_data = train_bm,
  valid_data = hold_bm,
  combo_list = combos_bm,
  positive_class = "pos"
)

# Re-run blood with same seed and parameters (will produce identical results)  
set.seed(2025)
nested_blood_validation_updated3_with_folds <- run_nested_with_validation(
  train_data = train_blood,
  valid_data = hold_blood,
  combo_list = combos_blood,
  positive_class = "pos"
)

# Verify results are identical (everything except fold_indices should match)
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("VERIFICATION: Comparing original vs re-run results (should be IDENTICAL)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

cat("\n── BM VALIDATION ──\n")
cat("  Nested metrics rows:", nrow(nested_bm_validation_updated3$nested_metrics), 
    "→", nrow(nested_bm_validation_updated3_with_folds$nested_metrics), "\n")
cat("  AUC match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$nested_metrics$auc_mean,
                     nested_bm_validation_updated3_with_folds$nested_metrics$auc_mean)), "\n")
cat("  Sensitivity match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$nested_metrics$sens_mean,
                     nested_bm_validation_updated3_with_folds$nested_metrics$sens_mean)), "\n")
cat("  Specificity match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$nested_metrics$spec_mean,
                     nested_bm_validation_updated3_with_folds$nested_metrics$spec_mean)), "\n")
cat("  Balanced accuracy match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$nested_metrics$balacc_mean,
                     nested_bm_validation_updated3_with_folds$nested_metrics$balacc_mean)), "\n")
cat("  Thresholds match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$thresholds,
                     nested_bm_validation_updated3_with_folds$thresholds)), "\n")

cat("\n── BLOOD VALIDATION ──\n")
cat("  Nested metrics rows:", nrow(nested_blood_validation_updated3$nested_metrics), 
    "→", nrow(nested_blood_validation_updated3_with_folds$nested_metrics), "\n")
cat("  AUC match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$nested_metrics$auc_mean,
                     nested_blood_validation_updated3_with_folds$nested_metrics$auc_mean)), "\n")
cat("  Sensitivity match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$nested_metrics$sens_mean,
                     nested_blood_validation_updated3_with_folds$nested_metrics$sens_mean)), "\n")
cat("  Specificity match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$nested_metrics$spec_mean,
                     nested_blood_validation_updated3_with_folds$nested_metrics$spec_mean)), "\n")
cat("  Balanced accuracy match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$nested_metrics$balacc_mean,
                     nested_blood_validation_updated3_with_folds$nested_metrics$balacc_mean)), "\n")
cat("  Thresholds match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$thresholds,
                     nested_blood_validation_updated3_with_folds$thresholds)), "\n")

cat("\n── VALIDATION METRICS (hold-out cohort) ──\n")
cat("  BM validation metrics match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$validation_metrics,
                     nested_bm_validation_updated3_with_folds$validation_metrics)), "\n")
cat("  Blood validation metrics match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$validation_metrics,
                     nested_blood_validation_updated3_with_folds$validation_metrics)), "\n")

cat("\n── OUTER FOLD PREDICTIONS ──\n")
cat("  BM predictions match:", 
    isTRUE(all.equal(nested_bm_validation_updated3$outer_predictions,
                     nested_bm_validation_updated3_with_folds$outer_predictions)), "\n")
cat("  Blood predictions match:", 
    isTRUE(all.equal(nested_blood_validation_updated3$outer_predictions,
                     nested_blood_validation_updated3_with_folds$outer_predictions)), "\n")

cat("\n✓ If all above are TRUE, results are IDENTICAL between runs\n")
cat("  (Only difference: fold_indices now captured for fragmentomics)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

# If verification passed, replace with versions that have fold_indices
nested_bm_validation_updated3 <- nested_bm_validation_updated3_with_folds
nested_blood_validation_updated3 <- nested_blood_validation_updated3_with_folds

rm(nested_bm_validation_updated3_with_folds, nested_blood_validation_updated3_with_folds)

message("✓ Fold indices successfully captured for BM and blood models")

bm_preds <- nested_bm_validation_updated3$outer_predictions


# -------------------------------------------------------------------------
# Optional development-only quick-iteration example (disabled)
#
# Purpose:
#   This commented block was used during method development to run the nested-CV
#   function on a single BM feature combination before launching the full model
#   grid. It is intentionally not active in the manuscript pipeline and should
#   not be uncommented for routine manuscript regeneration.
#
# Reproducibility note:
#   The manuscript model artifacts are generated by the full nested-CV/refit
#   sections below. This small example is retained only to document how the
#   model function can be smoke-tested during future code maintenance.
# -------------------------------------------------------------------------
# # 1) Minimal combos list with only BM_base
# combos_small_test <- list(
#   BM_base = c("zscore_BM", "detect_rate_BM", "z_score_detection_rate_BM")
# )
# 
# # 2) For clarity we’ll treat all entries in combos_small_test as BM combos
# combos_bm_test <- combos_small_test
# 
# # 3) Pull out the predictor names and build train / hold sets
# bm_preds_test   <- unique(unlist(combos_bm_test))
# train_bm_test   <- train_df %>% drop_na(all_of(c("MRD_truth", bm_preds_test)))
# hold_bm_test    <- hold_df  %>% drop_na(all_of(c("MRD_truth", bm_preds_test)))
# 
# # 4) Run nested CV just on BM_base
# test_out <- run_nested_with_validation(
#   train_data     = train_bm_test,
#   valid_data     = hold_bm_test,
#   combo_list     = combos_bm_test,
#   positive_class = positive_class
# )

# 5) Export
## Export the models that are updated
# saveRDS(nested_bm_validation_updated2, file = "nested_bm_validation_updated2.rds")
# saveRDS(nested_blood_validation_updated2, file = "nested_blood_validation_updated2.rds")
# saveRDS(nested_fragmentomics_validation_updated2, file = "nested_fragmentomics_validation_updated2.rds")

## Export the models that are updated
saveRDS(nested_bm_validation_updated3, file = "nested_bm_validation_updated5.rds")
saveRDS(nested_blood_validation_updated3, file = "nested_blood_validation_updated5.rds")
# Original fragmentomics (on full cohort) - unchanged
saveRDS(nested_fragmentomics_validation_updated3, file = "nested_fragmentomics_validation_updated3_original.rds")
# New fragmentomics (on restricted BM/blood cohorts with their folds)
saveRDS(nested_fragmentomics_bm_validation_updated3, file = "nested_fragmentomics_bm_validation_updated5.rds")
saveRDS(nested_fragmentomics_blood_validation_updated3, file = "nested_fragmentomics_blood_validation_updated5.rds")

## Load back in (optional)
# nested_bm_validation_updated3 <- readRDS("nested_bm_validation_updated5.rds")
# nested_blood_validation_updated3 <- readRDS("nested_blood_validation_updated5.rds")
# nested_fragmentomics_validation_updated3 <- readRDS("nested_fragmentomics_validation_updated3_original.rds")
# nested_fragmentomics_bm_validation_updated3 <- readRDS("nested_fragmentomics_bm_validation_updated3.rds")
# nested_fragmentomics_blood_validation_updated3 <- readRDS("nested_fragmentomics_blood_validation_updated3.rds")


## For consistency 
nested_bm_validation_updated2 <- nested_bm_validation_updated3
nested_blood_validation_updated2 <- nested_blood_validation_updated3

# -----------------------------------------------------------------------------
# 8.Get Model Metrics Together and Exported
# -----------------------------------------------------------------------------

# All validation metrics in one tibble
all_val <- bind_rows(
  nested_bm_validation_updated2$validation_metrics,
  nested_blood_validation_updated2$validation_metrics
)

## All primary
all_primary <- bind_rows(
  nested_bm_validation_updated2$nested_metrics,
  nested_blood_validation_updated2$nested_metrics
)

# ---- Step 3: Combine all fitted models into one named list ----
# Each element is a caret::train object (the elastic-net model trained on full cohort)
# Models are named by combo (e.g., "BM_zscore_only", "Blood_base", etc.)
# These trained models are what we apply to new patients for prediction
all_models <- c(
  nested_bm_validation_updated2$models,      # BM family models
  nested_blood_validation_updated2$models    # Blood family models
)

# ---- Step 4: Combine Youden-optimal thresholds ----
# For each model, this vector holds the threshold that maximizes (sens + spec - 1)
# These thresholds will be used to convert predicted probabilities to binary calls
all_thresholds <- c(
  nested_bm_validation_updated2$thresholds,   # BM family thresholds
  nested_blood_validation_updated2$thresholds # Blood family thresholds
)


### CRITICAL STEP: Handle Fragmentomics Cohort Variants ===================================
# BACKGROUND: We trained fragmentomics models in THREE ways:
#   1. "Full":         On ALL patients with fragmentomics data (largest cohort)
#   2. "BM_restricted": On ONLY patients with BM data (same subset as BM models)
#   3. "Blood_restricted": On ONLY patients with blood data (same subset as blood models)
#
# WHY THREE VARIANTS?
# - Full model: Shows what we can do if we use all available data
# - Restricted models: Fair comparison with BM/blood; same patient set, same fold structure
#   * Uses IDENTICAL folds from the BM/blood nested CV runs (prevents confounding due to different fold assignments)
#   * Isolates the effect of fragmentomics features vs. other sample types
#
# THE PROBLEM: Model names collide!
# If we just concatenate:  all_models <- c(all_models, frag_full, frag_bm, frag_blood)
# We get THREE models named "Fragmentomics_tumor_fraction_only" (and others)
# The named list only keeps the LAST one, silently losing the other two!
#
# THE SOLUTION: Add cohort suffix to fragmentomics model names
# "Fragmentomics_tumor_fraction_only_Full"
# "Fragmentomics_tumor_fraction_only_BM_restricted"
# "Fragmentomics_tumor_fraction_only_Blood_restricted"
# This preserves all three variants and makes the cohort explicit in the name.

# Remove old fragmentomics entries (if they exist) from BM/blood results
frag_names <- names(combos_fragmentomics)

all_val_clean <- all_val %>%
  filter(!combo %in% frag_names) %>%
  mutate(
    cohort = case_when(
      startsWith(combo, "BM_") ~ "BM_restricted",
      startsWith(combo, "Blood_") ~ "Blood_restricted",
      TRUE ~ NA_character_
    )
  )

all_primary_clean <- all_primary %>%
  filter(!combo %in% frag_names) %>%
  mutate(
    cohort = case_when(
      startsWith(combo, "BM_") ~ "BM_restricted",
      startsWith(combo, "Blood_") ~ "Blood_restricted",
      TRUE ~ NA_character_
    )
  )

# Now REBIND with all THREE fragmentomics variants
# Notice the mutate(cohort = ...) which labels which training cohort was used
# CRITICAL: Also update combo names to include cohort suffix so they match model names in all_models
all_val <- bind_rows(
  all_val_clean,
  # Original fragmentomics on full cohort
  nested_fragmentomics_validation_updated3$validation_metrics %>%
    mutate(combo = paste0(combo, "_Full"), cohort = "Full"),
  # Fragmentomics trained on restricted BM cohort with BM folds
  nested_fragmentomics_bm_validation_updated3$validation_metrics %>%
    mutate(combo = paste0(combo, "_BM_restricted"), cohort = "BM_restricted"),
  # Fragmentomics trained on restricted blood cohort with blood folds
  nested_fragmentomics_blood_validation_updated3$validation_metrics %>%
    mutate(combo = paste0(combo, "_Blood_restricted"), cohort = "Blood_restricted")
)

all_primary <- bind_rows(
  all_primary_clean,
  # Original fragmentomics on full cohort
  nested_fragmentomics_validation_updated3$nested_metrics %>%
    mutate(combo = paste0(combo, "_Full"), cohort = "Full"),
  # Fragmentomics trained on restricted BM cohort with BM folds
  nested_fragmentomics_bm_validation_updated3$nested_metrics %>%
    mutate(combo = paste0(combo, "_BM_restricted"), cohort = "BM_restricted"),
  # Fragmentomics trained on restricted blood cohort with blood folds
  nested_fragmentomics_blood_validation_updated3$nested_metrics %>%
    mutate(combo = paste0(combo, "_Blood_restricted"), cohort = "Blood_restricted")
)

# ---- Rebuild all_models with unique fragmentomics names ----
all_models_clean <- all_models[ ! names(all_models) %in% frag_names ]

# Create unique names for fragmentomics models by cohort
frag_models_full <- nested_fragmentomics_validation_updated3$models
names(frag_models_full) <- paste0(names(frag_models_full), "_Full")

frag_models_bm <- nested_fragmentomics_bm_validation_updated3$models
names(frag_models_bm) <- paste0(names(frag_models_bm), "_BM_restricted")

frag_models_blood <- nested_fragmentomics_blood_validation_updated3$models
names(frag_models_blood) <- paste0(names(frag_models_blood), "_Blood_restricted")

all_models <- c(
  all_models_clean,
  # Original fragmentomics models (on full cohort) - suffixed with _Full
  frag_models_full,
  # Fragmentomics models trained on restricted BM cohort - suffixed with _BM_restricted
  frag_models_bm,
  # Fragmentomics models trained on restricted blood cohort - suffixed with _Blood_restricted
  frag_models_blood
)

# CRITICAL: Update all_thresholds to include fragmentomics thresholds with cohort suffixes
# Remove old fragmentomics thresholds (if any existed)
all_thresholds_clean <- all_thresholds[ ! names(all_thresholds) %in% frag_names ]

# Create unique names for fragmentomics thresholds by cohort
frag_thresholds_full <- nested_fragmentomics_validation_updated3$thresholds
names(frag_thresholds_full) <- paste0(names(frag_thresholds_full), "_Full")

frag_thresholds_bm <- nested_fragmentomics_bm_validation_updated3$thresholds
names(frag_thresholds_bm) <- paste0(names(frag_thresholds_bm), "_BM_restricted")

frag_thresholds_blood <- nested_fragmentomics_blood_validation_updated3$thresholds
names(frag_thresholds_blood) <- paste0(names(frag_thresholds_blood), "_Blood_restricted")

all_thresholds <- c(
  all_thresholds_clean,
  # Original fragmentomics thresholds (on full cohort) - suffixed with _Full
  frag_thresholds_full,
  # Fragmentomics thresholds trained on restricted BM cohort - suffixed with _BM_restricted
  frag_thresholds_bm,
  # Fragmentomics thresholds trained on restricted blood cohort - suffixed with _Blood_restricted
  frag_thresholds_blood
)

## Add additional metrics

# ====== SECTION: POOL OUTER-FOLD PREDICTIONS FOR ROC ANALYSIS =========================================
# WHAT IS HAPPENING:
# The nested CV produces "outer-fold predictions" - these are predictions on the held-out outer folds
# that are COMPLETELY INDEPENDENT of the inner hyperparameter tuning. This is why they're so valuable:
# they give an UNBIASED estimate of model performance.
#
# We combine all outer-fold predictions from all three model families (BM, blood, fragmentomics)
# into one big data frame. Then we compute pooled ROC curves and AUC values across all folds.
#
# WHY POOL PREDICTIONS?
# - Individual fold estimates are noisy (small sample size per fold)
# - Pooling across folds gives a robust, single AUC estimate with confidence interval
# - The pooled AUC is what goes into the manuscript tables
#
# KEY POINT: Fragmentomics have THREE variants (Full, BM_restricted, Blood_restricted)
# They should NOT be pooled together! We group by (combo, cohort) to keep them separate.

# ---- 1. Gather all outer-fold predictions from all families ---
outer_preds <- bind_rows(
  # BM predictions (drop fragmentomics)
  nested_bm_validation_updated2$outer_predictions %>%
    filter(! combo %in% frag_names) %>%
    mutate(cohort = "BM_restricted"),
  # Blood predictions (drop fragmentomics)
  nested_blood_validation_updated2$outer_predictions %>%
    filter(! combo %in% frag_names) %>%
    mutate(cohort = "Blood_restricted"),
  # Original fragmentomics (on full cohort)
  nested_fragmentomics_validation_updated3$outer_predictions %>%
    mutate(combo = paste0(combo, "_Full"), cohort = "Full"),
  # Fragmentomics trained on BM cohort using BM folds
  nested_fragmentomics_bm_validation_updated3$outer_predictions %>%
    mutate(combo = paste0(combo, "_BM_restricted"), cohort = "BM_restricted"),
  # Fragmentomics trained on blood cohort using blood folds
  nested_fragmentomics_blood_validation_updated3$outer_predictions %>%
    mutate(combo = paste0(combo, "_Blood_restricted"), cohort = "Blood_restricted")
) %>%
  mutate(truth = factor(truth, levels = c("neg","pos"))) # make sure the response is a factor with neg first, pos second

# ---- 2. Compute POOLED AUC & 95% CI per combo and cohort ----
# WHAT IS HAPPENING:
# We take ALL the outer-fold predictions for a given model combo,
# pool them together, and compute a SINGLE AUC with 95% CI.
#
# WHY THIS APPROACH?
# - We have outer-fold predictions from 5 outer folds
# - Pooling them gives us ~50-100 test samples per combo (if training on 200-300 samples)
# - This is much more stable than computing AUC separately per fold
# - The DeLong method computes confidence intervals that account for variability
#
# CRITICAL: group_by(combo, cohort)
# For fragmentomics, this ensures we don't pool Full + BM_restricted + Blood_restricted together.
# They're different training cohorts and should have separate AUC values.
#
# OUTPUT: pooled_auc_tbl
# - One row per combo (and cohort if fragmentomics)
# - Columns: combo, cohort, auc_pooled, ci_low, ci_high, auc_pooled_ci (formatted)
#
# EXAMPLE:
#   combo = \"BM_zscore_only\", cohort = NA, auc_pooled = 0.75, ci = (0.68-0.82)
#   combo = \"Fragmentomics_tumor_fraction_only\", cohort = \"Full\", auc_pooled = 0.82, ci = (0.76-0.88)
#   combo = \"Fragmentomics_tumor_fraction_only\", cohort = \"BM_restricted\", auc_pooled = 0.80, ci = (0.74-0.86)
pooled_auc_tbl <- outer_preds %>%
  group_by(combo, cohort) %>%
  summarise(
    roc_obj     = list(roc(truth, prob,
                           levels    = c("neg","pos"),
                           direction = "<",
                           quiet     = TRUE)),
    auc_pooled  = as.numeric(auc(roc_obj[[1]])),
    ci_low      = ci.auc(roc_obj[[1]], method = "delong")[1],
    ci_high     = ci.auc(roc_obj[[1]], method = "delong")[3],
    .groups     = "drop"
  ) %>%
  mutate(
    auc_pooled_ci = sprintf("%.2f (%.2f-%.2f)", auc_pooled, ci_low, ci_high)
  ) %>%
  select(-roc_obj)   # drop the heavy objects

# ---- 3. Add to primary metrics table ----------------------------------
# For left_join, we need to match on both combo and cohort
all_primary <- all_primary %>%
  left_join(pooled_auc_tbl, by = c("combo", "cohort"))

# ═══════════════════════════════════════════════════════════════════════════════════
# SECTION: DETECT & FLAG INVERTED PREDICTIONS
# ═══════════════════════════════════════════════════════════════════════════════════
# WHAT DOES THIS DO?
# Sometimes outer-fold predictions are completely inverted/flipped. This manifests as:
#   - Nested CV AUC:  0.68 (legitimately high, within-fold performance)
#   - Pooled AUC:     0.33 (inverted, worse than random, outer-fold predictions backwards)
#   - Difference:     0.35 (large discrepancy indicates instability)
#   - Pooled < 0.50:  Confirms inversion (AUC should be > 0.50 always)
#
# WHY DOES THIS HAPPEN?
# Small sample size (n=42) in restricted cohorts causes fold-to-fold instability.
# Model finds real signal within folds, but outer-fold predictions systematically flip.
#
# HOW DO WE DETECT IT?
# Threshold: abs(auc_mean - auc_pooled) > 0.20 AND auc_pooled < 0.50
#   - 0.20 point difference = clearly meaningful discrepancy
#   - AUC < 0.50 = predictions are worse than random (inverted direction)
#
# KEY INSIGHT:
# Nested CV AUC is reliable (proper cross-validation).
# Pooled AUC is unreliable for these models (use nested CV instead).

# Step 1: Compute difference and flag inverted models
all_primary_tmp <- all_primary %>%
  mutate(
    # Difference between nested CV and pooled outer-fold AUC
    auc_discrepancy = abs(auc_mean - auc_pooled),
    # Flag as inverted if large discrepancy + poor pooled AUC
    metric_status = case_when(
      is.na(auc_pooled) ~ "no_pooled",
      auc_discrepancy > 0.20 & auc_pooled < 0.50 ~ "INVERTED_FLAG",
      auc_discrepancy > 0.15 & auc_pooled < 0.50 ~ "caution_inverted",
      TRUE ~ "check_ok"
    ),
    # Recommendation for metric use
    metric_recommendation = case_when(
      metric_status == "INVERTED_FLAG" ~ "Use NESTED CV only; pooled AUC unreliable",
      metric_status == "caution_inverted" ~ "Use NESTED CV with caution; pooled AUC unstable",
      TRUE ~ "Both metrics reliable"
    )
  )

# Step 2: Create diagnostic summary table
inverted_summary <- all_primary_tmp %>%
  filter(metric_status %in% c("INVERTED_FLAG", "caution_inverted")) %>%
  select(combo, cohort, auc_mean, auc_pooled, auc_discrepancy, metric_status, metric_recommendation)

# Step 3: Print diagnostic report
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║         METRIC QUALITY ASSESSMENT - INVERTED PREDICTION DETECTION             ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

if (nrow(inverted_summary) > 0) {
  cat("\n⚠ WARNING: INVERTED METRICS DETECTED\n")
  cat("─────────────────────────────────────────────────────────────────────────────────\n")
  cat("The following models show discrepancy between nested CV and pooled AUC:\n")
  cat("PROBABLE CAUSE: Small sample size in restricted cohort causing fold instability\n")
  cat("IMPACT: Outer-fold predictions are INVERTED (negatives → higher scores than positives)\n")
  cat("SOLUTION: Use nested CV AUC only; disregard pooled AUC for these models\n\n")
  
  # Print detailed table
  print_tbl <- inverted_summary %>%
    mutate(
      auc_mean = sprintf("%.3f", auc_mean),
      auc_pooled = sprintf("%.3f", auc_pooled),
      auc_discrepancy = sprintf("%.3f", auc_discrepancy)
    )
  
  print(print_tbl %>% as.data.frame(), row.names = FALSE)
  
  cat("\n✓ Flagged models identified and marked in all_primary table\n")
  cat("  Use 'metric_status' and 'metric_recommendation' columns for guidance\n\n")
  
} else {
  cat("\n✓ GOOD NEWS: No inverted metrics detected\n")
  cat("  All models show consistent nested CV and pooled AUC performance\n")
  cat("  (auc_discrepancy ≤ 0.20 OR auc_pooled ≥ 0.50)\n\n")
}

# Summary statistics
n_inverted_flag <- sum(all_primary_tmp$metric_status == "INVERTED_FLAG", na.rm = TRUE)
n_caution_inverted <- sum(all_primary_tmp$metric_status == "caution_inverted", na.rm = TRUE)
n_ok <- sum(all_primary_tmp$metric_status == "check_ok", na.rm = TRUE)

cat("SUMMARY:\n")
cat(sprintf("  Models with INVERTED FLAG:           %2d\n", n_inverted_flag))
cat(sprintf("  Models with CAUTION (borderline):    %2d\n", n_caution_inverted))
cat(sprintf("  Models with OK metrics:              %2d\n", n_ok))
cat(sprintf("  Models with no pooled AUC:          %2d\n", sum(all_primary_tmp$metric_status == "no_pooled", na.rm = TRUE)))
cat("\n")

# Step 4: Show which specific models are problematic
if (n_inverted_flag > 0) {
  cat("INVERTED MODELS (pooled AUC < 0.50 + large discrepancy):\n")
  problematic <- all_primary_tmp %>%
    filter(metric_status == "INVERTED_FLAG") %>%
    pull(combo) %>%
    unique()
  for (m in problematic) {
    cat(sprintf("  • %s\n", m))
  }
  cat("\n")
}

cat("═════════════════════════════════════════════════════════════════════════════════\n\n")



## Pick winners after reviewin the tables
wanted <- c(
  # BM-only models
  "BM_zscore_only_sites",              # best BM by AUC
  "BM_zscore_only_detection_rate",     # best BM by sensitivity
  "BM_base_zscore",     # best BM by sensitivity
  
  
  # Blood-only models
  "Blood_base",   #  best blood by AUC
  "Blood_rate_only",                   # second best blood by AUC
  "Blood_plus_fragment_min",           # best blood by accuracy
  "Blood_zscore_only_detection_rate",  # best blood by sensitivity
  "Blood_zscore_only_sites",  # best blood in validation cohort
  
  # Fragmentomics-only models
  "Fragmentomics_tumor_fraction_only_Full",             # best fragmentomics in validation cohort and good sens
  "Fragmentomics_mean_coverage_only_Full"   # best auc in primary and accuracy
)

# Then pull just those rows
selected_primary <- all_primary # %>%
#  filter(combo %in% wanted)

## Save
selected_combos   <- selected_primary$combo
selected_models   <- all_models[selected_combos]
selected_thr      <- all_thresholds[selected_combos]

# saveRDS(selected_models,
#         file = file.path(outdir, "selected_combo_models_2025-07-25.rds"))
# saveRDS(selected_thr,
#         file = file.path(outdir, "selected_combo_thresholds_2025-07-25.rds"))
# 
# saveRDS(selected_models,
#         file = file.path(outdir, "selected_combo_models_2025-07-31.rds"))
# saveRDS(selected_thr,
#         file = file.path(outdir, "selected_combo_thresholds_2025-07-31.rds"))


saveRDS(selected_models,
        file = file.path(outdir, "selected_combo_models_2026-02-16.rds"))
saveRDS(selected_thr,
        file = file.path(outdir, "selected_combo_thresholds_2026-02-16.rds"))


# -----------------------------------------------------------------------------
# 8b. Export Combined Metrics & Models
# -----------------------------------------------------------------------------
# 1) Validation‐set metrics (hold‐out cohort)
saveRDS(all_val,
        file = file.path(outdir, "all_validation_metrics_v2_with_fragmentomics_restricted_cohorts.rds"))
write.csv(all_val,
          file = file.path(outdir, "all_validation_metrics_v2_with_fragmentomics_restricted_cohorts.csv"),
          row.names = FALSE)

# 2) Nested‐CV metrics (primary/frontline cohort)
saveRDS(all_primary,
        file = file.path(outdir, "all_nested_cv_metrics_v2_with_fragmentomics_restricted_cohorts.rds"))
write.csv(all_primary,
          file = file.path(outdir, "all_nested_cv_metrics_v2_with_fragmentomics_restricted_cohorts.csv"),
          row.names = FALSE)

# 3) Full list of fitted models
saveRDS(all_models,
        file = file.path(outdir, "all_cfWGS_models_list_v2_with_fragmentomics_restricted_cohorts.rds"))

# (Optional) Export the thresholds vector:
saveRDS(all_thresholds,
        file = file.path(outdir, "all_model_thresholds_v2_with_fragmentomics_restricted_cohorts.rds"))
write.csv(
  tibble(
    combo     = names(all_thresholds),
    threshold = unname(all_thresholds)
  ),
  file = file.path(outdir, "all_model_thresholds_v2_with_fragmentomics_restricted_cohorts.csv"),
  row.names = FALSE
)

message("Exported validation metrics, nested‐CV metrics, models list, and thresholds to ", outdir)


### Export cleaned version for supplementary tables
# 1) Pull only the non-NA columns from the training results
primary_clean <- all_primary %>%
  select(where(~ !all(is.na(.)))) %>%   # drop any column that is entirely NA
  rename(sample_group = cohort) %>%      # cohort column shows fragmentomics training cohort (Full/BM_restricted/Blood_restricted)
  mutate(eval_cohort = "Training")       # eval_cohort distinguishes Training vs Testing

# 2) do the same for the validation set, and rename its cols to match the *_mean convention
val_clean <- all_val %>%
  select(where(~ !all(is.na(.)))) %>%
  rename(
    auc_mean       = auc_valid,
    sens_mean      = sens_valid,
    spec_mean      = spec_valid,
    balacc_mean    = bal_accuracy,
    ppv_mean       = ppv,
    npv_mean       = npv,
    f1_mean        = f1,
    brier_mean     = brier,
    pAUC90_mean    = pAUC90,
    acc_mean       = accuracy,
    sens95_mean    = sens_at_95_spec_valid,
    spec95_mean    = spec_at_95_sens_valid,
    sample_group   = cohort              # cohort column shows fragmentomics training cohort (Full/BM_restricted/Blood_restricted)
  ) %>%
  mutate(eval_cohort = "Testing")        # eval_cohort distinguishes Training vs Testing



# 1. Helper to pull the first (and only) column out of any data.frame columns
flatten_df_cols <- function(df) {
  is_df <- vapply(df, is.data.frame, logical(1))
  df[is_df] <- lapply(df[is_df], function(x) x[[1]])
  df
}

# 2. Helper to clean up the sens95/spec95 names
fix_names <- function(df) {
  names(df) <- sub("sens95_mean\\$sensitivity", "sens95_mean", names(df))
  names(df) <- sub("spec95_mean\\$specificity", "spec95_mean", names(df))
  df
}

# Apply to both data sets
primary_clean <- primary_clean %>%
  flatten_df_cols() %>%
  fix_names()

val_clean <- val_clean %>%
  flatten_df_cols() %>%
  fix_names()

# 3) stack them into one combined table
combined_metrics <- bind_rows(primary_clean, val_clean) %>%
  select(
    combo,
    sample_group,      # ← Distinguishes Full vs BM_restricted vs Blood_restricted (fragmentomics training cohorts)
    eval_cohort,       # ← Distinguishes Training vs Testing
    auc_mean,
    sens_mean,
    spec_mean,
    balacc_mean,
    brier_mean,
    pAUC90_mean,
    acc_mean,
    everything()       # keep all other columns
  ) %>%
  select(
    -sens95_mean,
    -spec95_mean,
    -ppv_mean,
    -npv_mean,
    -f1_mean
  )

# 4) Export only the Fragmentomics models
frag_tbl <- combined_metrics %>%
  filter(startsWith(combo, "Fragmentomics_"))

write_csv(
  frag_tbl,
  "Output_tables_2025/Supplementary_Table_4_Fragmentomics_Performance_v3_Feb2026_with_restricted_cohorts.csv"
)

# 2) Export all models (includes sample_group and eval_cohort to track fragmentomics training cohorts)
# Columns:
#   - combo: Model name
#   - sample_group: For fragmentomics, shows training cohort (Full/BM_restricted/Blood_restricted)
#                   For BM/blood combos, this will be NA (not applicable)
#   - eval_cohort: Training or Testing (evaluation dataset)
#   - auc_mean, sens_mean, etc.: Performance metrics

# NOTE: This file exports as "Supplementary_Table_3_All_Model_performance_nested_CV_v3_Feb2026_with_restricted_fragmentomics.csv"
#       But in the final manuscript, this corresponds to: Supplementary_Table_4_All_Model_performance_nested_CV.csv
#       (Table numbering adjusts in final manuscript due to rearrangement of supplementary tables)

write_csv(
  primary_clean %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  "Final Tables and Figures/Supplementary_Table_3_All_Model_performance_nested_CV_v3_Feb2026_with_restricted_fragmentomics.csv"
)

# -------------------------------------------------------------------------
# Historical nested-CV table (superseded as Supplementary Table 4)
#
# What this is:
#   The earlier sample-row nested-CV model-performance table across BM,
#   blood/cfDNA, and fragmentomics model families.
#
# Why it is here:
#   This export was previously registered as Supplementary Table 4 and is kept
#   so the original analysis can be traced. The current final Supplementary
#   Table 4 is the 32-model, 50-repeat patient-grouped result assembled from
#   script 6_12 outputs. Do not use this older file in place of that table.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supplementary_Table_3_All_Model_performance_nested_CV_v3_Feb2026_with_restricted_fragmentomics.csv",
  artifact_id = "STABLE4",
  role = "table_csv",
  description = "Nested cross-validation model-performance table used as Supplementary Table 4.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)

# NOTE: This file exports as "Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv"
#       But in the final manuscript, this corresponds to: Supplementary_Table_6_All_Model_performance_nested_CV_on_test_cohort.csv
#       (Table numbering adjusts in final manuscript due to rearrangement of supplementary tables)

write_csv(
  val_clean %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  "Final Tables and Figures/Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv"
)

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 6
#
# What this is:
#   Test/non-frontline cohort model-performance table using the preserved
#   training-derived models and thresholds.
#
# Why it is here:
#   The historical filename says Supplementary Table 5, but the audited
#   manuscript map identifies this export as final Supplementary Table 6.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supplementary_Table_5_All_Model_performance_testing_cohort_v3_Feb2026_with_restricted_fragmentomics.csv",
  artifact_id = "STABLE6",
  role = "table_csv",
  description = "Testing-cohort model-performance table used as Supplementary Table 6.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)




# Training-cohort operating-point diagnostics for all model probability columns.
# These summaries feed the full-refit metrics table copied as Supplementary
# Table 5 below.
# all_models is a named list of caret::train objects
sens_at_95spec <- map_dbl(all_models, function(model) {
  # 1) Grab the best‐tuning parameters
  best <- model$bestTune
  
  # 2) Subset the out‐of‐fold predictions to only those rows matching bestTune
  preds_best <- inner_join(
    model$pred,
    best,
    by = names(best)
  )
  
  # 3) Build ROC on those held‐out probabilities
  roc_obj <- roc(
    response  = preds_best$obs,    # observed neg/pos
    predictor = preds_best$pos,    # predicted probability of “pos”
    levels    = c("neg", "pos"),
    direction = "<"
  )
  
  # 4) Extract the sensitivity at 95% specificity
  sens <- coords(
    roc_obj,
    x         = 0.95,
    input     = "specificity",
    ret       = "sensitivity",
    transpose = FALSE
  )
  
  as.numeric(sens)
})

# Combine into a dataframe for reporting
tibble(
  model            = names(sens_at_95spec),
  sens_at_95_spec  = sens_at_95spec
)


# Compute specificity at 95% sensitivity for each model
spec_at_95sen <- map_dbl(all_models, function(model) {
  # 1) Optimal tuning parameters
  best <- model$bestTune
  
  # 2) Subset to out-of-fold predictions at bestTune
  preds_best <- inner_join(
    model$pred,
    best,
    by = names(best)
  )
  
  # 3) Build ROC on those held-out predictions
  roc_obj <- roc(
    response  = preds_best$obs,   # true class: “neg”/“pos”
    predictor = preds_best$pos,   # predicted probability for “pos”
    levels    = c("neg","pos"),
    direction = "<"
  )
  
  # 4) Extract specificity when sensitivity = 0.95
  spec <- coords(
    roc_obj,
    x         = 0.95,
    input     = "sensitivity",
    ret       = "specificity",
    transpose = FALSE
  )
  
  as.numeric(spec)
})

# Create a summary table
table <- tibble(
  model                  = names(spec_at_95sen),
  spec_at_95_sensitivity = spec_at_95sen
)


# ===================================================================================================
# SECTION 9A: APPLY TRAINED MODELS TO FULL COHORT & SCORE ALL PATIENTS
# ===================================================================================================
# WHAT HAPPENS IN THIS SECTION:
# We take the trained elastic-net models and apply them to EVERY PATIENT in the dataset,
# generating predicted probability of MRD+ (score) and binary calls (MRD pos/neg).
#
# WHY IS THIS IMPORTANT?
# - Models were trained on Frontline cohort only (to avoid data leakage)
# - But we want predictions for all patients, including Non-frontline (validation) and future patients
# - This section shows how to operationalize the models for deployment
#
# WORKFLOW:
# 1. apply_selected():     Helper function to score multiple models at once
# 2. Patient masking:      Set predictions to NA for ineligible patients
#    - BM models: only score patients WITH bone marrow samples
#    - Blood models: only score patients WITH blood samples
#    - Fragmentomics: score everyone (universal applicability)
# 3. Export two versions:
#    - data_scored_masked:  Clinical use (ineligible patients = NA, can't use their predictions)
#    - data_scored_full:    Research use (all patients scored, documents model capability)
#
# KEY INSIGHT:
# Not all patients have all sample types! A patient without BM data can't be scored by a BM model
# (it would use missing feature values = garbage output). Masking prevents inappropriate predictions.

# ---- Helper Function: apply_selected() ----
# PURPOSE: Apply a batch of trained models to new data
# INPUTS:
#   - dat:              New data frame (with same column structure as training data)
#   - models:           Named list of fitted caret::train objects
#   - thresholds:       Named numeric vector of decision thresholds
#   - positive_class:   Which class label means "MRD positive" (usually "pos")
# OUTPUTS:
#   - Same data frame WITH NEW COLUMNS appended:
#     * model_name_prob:  Predicted probability of MRD+ (numeric 0-1)
#     * model_name_call:  Binary prediction (0=MRD-, 1=MRD+) based on threshold
#
# HOW IT WORKS:
# For each model combo in the list:
#   1. Extract the set of predictors (features) the model needs
#   2. Find complete-case rows (no missing values in required features)
#   3. Predict probabilities ONLY for those rows (others remain NA)
#   4. Apply threshold to convert probabilities to binary calls
#   5. Add two new columns to output: _prob and _call



# -----------------------------------------------------------------------------
# 9A. Full-Cohort & Prep for Dilution-Series Scoring
# -----------------------------------------------------------------------------
### Apply on entire dataframe 
## Define function
apply_selected <- function(dat, models, thresholds, positive_class = "pos") {
  out <- dat            # start with the full data frame
  n   <- nrow(dat)
  
  # Check that models and thresholds align
  if (!all(names(models) %in% names(thresholds))) {
    stop("Model names and threshold names do not match!")
  }
  common <- intersect(names(models), names(thresholds))
  
  for (cmb in common) {
    fit <- models[[cmb]]
    thr <- thresholds[[cmb]]
    
    # which predictors the model actually needs
    preds <- setdiff(names(fit$trainingData), ".outcome")
    
    # figure out which rows we can actually score
    can_score <- complete.cases(dat[preds])
    
    # pre‐allocate
    probs <- rep(NA_real_,    n)
    calls <- rep(NA_integer_, n)
    
    # only predict where allowed
    if (any(can_score)) {
      probs[can_score] <- predict(fit,
                                  newdata = dat[can_score, preds, drop = FALSE],
                                  type    = "prob")[[positive_class]]
      calls[can_score] <- as.integer(probs[can_score] >= thr)
    }
    
    # stick them back on
    out[[paste0(cmb, "_prob")]] <- probs
    out[[paste0(cmb, "_call")]] <- calls
  }
  
  out
}

## Run function
### Primary cohort
data_scored <- apply_selected(
  dat           = dat,            # the master data frame with all features
  models        = selected_models,
  thresholds    = selected_thr,
  positive_class= "pos"
)


# ===================================================================================================
# SECTION 9A-CONTINUED: PATIENT ELIGIBILITY MASKING
# ===================================================================================================
# CRITICAL STEP: After generating predictions for all patients, we need to "mask" (set to NA)
# predictions that are inappropriate for certain patients.
#
# THE PROBLEM:
# - BM models use bone marrow-specific features (e.g., BM z-scores, BM detection rates)
# - Blood models use cfDNA-specific features (e.g., Blood z-scores)
# - A patient WITHOUT bone marrow data cannot be scored by a BM model
#   (the model would get missing feature values = garbage output = meaningless probability)
# - Similarly for blood models
# - Fragmentomics models use universal features (fragment score, tumor fraction) applicable to anyone
#
# THE SOLUTION:
# 1. Identify which patients have complete data for each sample type
#    - good_bm_patients: Patients with all BM features available
#    - good_blood_patients: Patients with all blood features available
#
# 2. Create a "model_map" that assigns each model to a family (BM/BLOOD/FRAG)
#
# 3. For each patient/model combination:
#    - IF patient is NOT eligible for that model family
#    - THEN set their prediction to NA (using typed_na to preserve column type)
#
# RESULT:
# - data_scored_masked: Safe for clinical use (guarantees only valid predictions are present)
# - data_scored_full:   Research version (unmasked, shows all predictions)

### Now restrict to eligible patients 
## Patient eligibility vectors from earlier 
good_blood_patients <- cfDNA_good_patients$Patient
good_bm_patients    <- bm_good_patients$Patient

## All probability and score columns produced by apply_selected
# These are columns like "BM_zscore_only_prob", "Blood_base_call", "Fragmentomics_tumor_fraction_only_prob"
score_cols <- grep("_(prob|call)$", names(data_scored), value = TRUE, ignore.case = TRUE)

## Map prob column -> model base name (strip the trailing "_prob" or "_call")
# Example: "BM_zscore_only_prob" -> "BM_zscore_only"
model_base <- tibble(
  model_col  = score_cols,
  model_base = str_remove(score_cols, "_(?i:prob|call)$")   # strip _prob or _call
)

## Assign each model to its family (BM / BLOOD / FRAG)
# This is a heuristic: models starting with "BM_" belong to BM family, etc.
# This determines which patients are eligible for which models
model_map <- model_base %>%
  mutate(family = case_when(
    str_detect(model_base, regex("^BM_", ignore_case = TRUE))            ~ "BM",
    str_detect(model_base, regex("^Blood_", ignore_case = TRUE))         ~ "BLOOD",
    str_detect(model_base, regex("^Fragmentomics_", ignore_case = TRUE)) ~ "FRAG",
    TRUE                                                                 ~ "OTHER"
  ))

# --- 4) Mask disallowed patients for BOTH prob+call columns --------
data_scored_masked <- data_scored

typed_na <- function(x) {
  if (is.factor(x))         return(NA)                 # factor NA
  if (is.logical(x))        return(NA)                 # logical NA
  if (is.integer(x))        return(NA_integer_)        # integer NA
  if (is.double(x))         return(NA_real_)           # numeric NA
  if (is.numeric(x))        return(NA_real_)
  if (is.character(x))      return(NA_character_)      # character NA
  if (inherits(x, "Date"))  return(as.Date(NA))        # Date NA
  if (inherits(x, "POSIXt"))return(as.POSIXct(NA))     # POSIXct NA
  NA
}

walk2(model_map$model_col, model_map$family, function(col, fam) {
  # Step 2a: Determine eligible patients for this model family
  allowed <- switch(
    fam,
    "BM"    = good_bm_patients,        # Only BM patients (have complete BM features)
    "BLOOD" = good_blood_patients,     # Only blood patients (have complete cfDNA features)
    "FRAG"  = unique(data_scored$Patient),  # Everyone (fragmentomics universal)
    "OTHER" = unique(data_scored$Patient)   # Default: everyone
  )
  
  # Step 2b: Create mask (TRUE = eligible, FALSE = ineligible)
  x <- data_scored[[col]]
  mask <- data_scored$Patient %in% allowed
  
  # Step 2c: Set ineligible patient predictions to typed NA
  x[!mask] <- typed_na(x)        # preserve the column's original type
  data_scored_masked[[col]] <<- x   # <<- assigns to parent environment
})

# ---- Step 3: Diagnostic summary: N before/after masking ----
# These messages show how many patients have valid predictions for each model
# BEFORE masking: includes invalid predictions (garbage from missing features)
# AFTER masking:  only includes valid predictions (from eligible patients)
#
# If a model shows large N changes, it means many patients lack required sample types

message("=== Effective N BEFORE masking (all predicted, some invalid) ===")
for (m in score_cols) {
  n_eff <- data_scored %>% tidyr::drop_na(all_of(m)) %>% nrow()
  message(sprintf("%-40s N = %d", m, n_eff))
}

message("=== Effective N AFTER masking (only eligible patients) ===")
for (m in score_cols) {
  n_eff <- data_scored_masked %>% tidyr::drop_na(all_of(m)) %>% nrow()
  message(sprintf("%-40s N = %d", m, n_eff))
}


### Save this 
# saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated3.rds"))
# write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated3.csv"))

## Another export
saveRDS(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.rds"))
write_csv(data_scored, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6_full.csv"))

saveRDS(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.rds"))
write_csv(data_scored_masked, file = file.path(outdir, "all_patients_with_BM_and_blood_calls_updated6.csv"))





# -----------------------------------------------------------------------------
# 9B. Now re-get model metrics on the whole cohort 
# -----------------------------------------------------------------------------
### See if can optimize 
# --- Build a "whole cohort" evaluation set (exclude Dx/Baseline timepoints) ---
whole <- data_scored_masked %>%
  filter(!timepoint_info %in% c("Baseline","Diagnosis")) %>%
  # keep rows that have either label or any model prob; we will drop NAs per-model anyway
  mutate(MRD_truth = as.integer(MRD_truth)) %>% filter(Cohort == "Frontline")

prev_whole <- mean(whole$MRD_truth == 1, na.rm = TRUE)  # prevalence matters for accuracy/PPV/NPV

# Identify all probability columns generated by model predictions
prob_cols <- grep("_prob$", names(whole), value = TRUE)

# Utility to compute derived metrics from sens/spec plus prevalence
.add_metrics <- function(df, prev) {
  df %>%
    mutate(
      accuracy      = sensitivity * prev + specificity * (1 - prev),
      bal_accuracy  = (sensitivity + specificity) / 2,
      precision     = ifelse((sensitivity*prev + (1 - specificity)*(1 - prev)) == 0,
                             NA_real_,
                             (sensitivity*prev) / (sensitivity*prev + (1 - specificity)*(1 - prev))),
      f1            = ifelse((precision + sensitivity) == 0,
                             NA_real_,
                             2 * precision * sensitivity / (precision + sensitivity))
    )
}

# ===================================================================================================
# THRESHOLD TUNING FUNCTIONS
# ===================================================================================================
# THE PROBLEM WITH THRESHOLDS:
# - A model predicts a probability P(MRD+) for each patient
# - But we need a binary decision: MRD+ or MRD-
# - We do this by choosing a threshold T: if P >= T, predict MRD+, else predict MRD-
# - DIFFERENT thresholds give DIFFERENT sensitivity/specificity tradeoffs:
#   * Threshold = 0.1: Very low bar to call MRD+, high sensitivity, low specificity
#   * Threshold = 0.5: Middle ground
#   * Threshold = 0.9: Very high bar, low sensitivity, high specificity
#
# THE SOLUTION: Define threshold objectives
# We evaluate thresholds at several OBJECTIVES and report all of them:
# 1. Youden-optimal:     Maximizes (sensitivity + specificity - 1)
#    - Interpretation: Balanced, reasonable default
# 2. Max accuracy:       Maximizes overall % correct
#    - Interpretation: Best if costs of FP and FN are equal, accounting for prevalence
# 3. Max balanced acc.:  Maximizes mean(sensitivity, specificity)
#    - Interpretation: Balanced across both classes, ignores prevalence
# 4. Weighted score:     Maximizes (0.7*sensitivity + 0.3*specificity)
#    - Interpretation: Prioritizes sensitivity (catching positives) 70% of the time
# 5. Constrained:        Given a constraint (sen >= 90%, spec >= 90%), optimize the other
#    - Interpretation: Clinical use cases with hard requirements

# ---- Function: tune_thresholds_one() ----
# PURPOSE: For ONE model, evaluate it at multiple thresholds and find optimal operating points
# INPUTS:
#   - prob_col: Name of the probability column (e.g., "BM_zscore_only_prob")
#   - dat:      Data frame containing the labels and probabilities
#   - y_col:    Name of the label column (default: "MRD_truth")
#   - prev:     Prevalence of positives (needed for accuracy calculation)
#   - min_spec / target_sens: Constraints for constrained optimization
#   - sens_weight: Weight for sens vs spec in weighted score (0.7 = 70% sens, 30% spec)
# OUTPUTS:
#   - Tibble with ONE ROW PER OBJECTIVE, showing:
#     * Model, Objective, Threshold (the value)
#     * Sensitivity, Specificity, Accuracy, BalancedAccuracy, Precision, F1
#
# HOW IT WORKS:
# 1. Filter to complete cases (no NA in label or probability)
# 2. Build ROC curve on these cases
# 3. Extract all (threshold, sensitivity, specificity) pairs
# 4. Compute derived metrics (accuracy, F1, etc.) for each threshold
# 5. For each objective, find the best threshold
# 6. Return all 5 objectives as separate rows

# Core function: sweep thresholds and choose by several objectives
tune_thresholds_one <- function(prob_col,
                                dat,
                                y_col = "MRD_truth",
                                prev  = prev_whole,
                                min_spec      = 0.95,   # for "sens↑ @ spec≥min_spec"
                                target_sens   = 0.95,   # for "acc↑ @ sens≥target_sens"
                                sens_weight   = 0.7,    # α in weighted score α*sens + (1-α)*spec
                                direction     = "<") {
  
  # Drop rows where either label or predictor is NA for this model
  dd <- dat %>%
    select(all_of(c(y_col, prob_col))) %>%
    rename(y = !!y_col, p = !!prob_col) %>%
    filter(!is.na(y) & !is.na(p))
  
  if (nrow(dd) == 0) {
    return(tibble(
      Model     = prob_col, Objective = c("max_accuracy","max_balacc","sens_weighted",
                                          paste0("sens_at_spec≥", min_spec),
                                          paste0("acc_at_sens≥", target_sens)),
      Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
      Accuracy  = NA_real_, BalancedAccuracy = NA_real_, Precision = NA_real_, F1 = NA_real_
    ))
  }
  
  roc_obj <- roc(response = dd$y, predictor = dd$p, levels = c(0,1), direction = direction)
  
  roc_df <- coords(roc_obj, x = "all",
                   ret = c("threshold","sensitivity","specificity"),
                   transpose = FALSE) %>%
    as_tibble() %>%
    .add_metrics(prev) %>%
    # convenience columns for ranking
    mutate(weighted_score = sens_weight * sensitivity + (1 - sens_weight) * specificity)
  
  # 1) Global maximizers
  best_acc     <- roc_df %>% slice_max(order_by = accuracy, n = 1, with_ties = FALSE)
  best_balacc  <- roc_df %>% slice_max(order_by = bal_accuracy, n = 1, with_ties = FALSE)
  best_weight  <- roc_df %>% slice_max(order_by = weighted_score, n = 1, with_ties = FALSE)
  
  # 2) Constrained: maximize sensitivity with specificity ≥ min_spec
  cand_spec <- roc_df %>% filter(specificity >= min_spec)
  best_sens_at_spec <- if (nrow(cand_spec)) {
    cand_spec %>% arrange(desc(sensitivity), desc(specificity)) %>% slice(1)
  } else {
    tibble(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_,
           accuracy = NA_real_, bal_accuracy = NA_real_, precision = NA_real_, f1 = NA_real_)
  }
  
  # 3) Constrained: maximize accuracy with sensitivity ≥ target_sens
  cand_sens <- roc_df %>% filter(sensitivity >= target_sens)
  best_acc_at_sens <- if (nrow(cand_sens)) {
    cand_sens %>% arrange(desc(accuracy), desc(specificity)) %>% slice(1)
  } else {
    tibble(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_,
           accuracy = NA_real_, bal_accuracy = NA_real_, precision = NA_real_, f1 = NA_real_)
  }
  
  bind_rows(
    best_acc     %>% mutate(Objective = "max_accuracy"),
    best_balacc  %>% mutate(Objective = "max_balacc"),
    best_weight  %>% mutate(Objective = paste0("sens_weighted_α=", sens_weight)),
    best_sens_at_spec %>% mutate(Objective = paste0("sens_at_spec≥", min_spec)),
    best_acc_at_sens  %>% mutate(Objective = paste0("acc_at_sens≥", target_sens))
  ) %>%
    transmute(
      Model            = prob_col,
      Objective,
      Threshold        = threshold,
      Sensitivity      = sensitivity,
      Specificity      = specificity,
      Accuracy         = accuracy,
      BalancedAccuracy = bal_accuracy,
      Precision        = precision,
      F1               = f1
    )
}

# Run for all models on the WHOLE cohort
tuned_whole <- map_dfr(
  prob_cols,
  ~tune_thresholds_one(.x,
                       dat         = whole,
                       prev        = prev_whole,
                       min_spec    = 0.8,   # adjust these to taste
                       target_sens = 0.8,
                       sens_weight = 0.7)
)

# Round for display and write alongside existing outputs
tuned_whole_out <- tuned_whole %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(tuned_whole_out)

# If all_perf_metrics is in memory, write a single workbook with both
write_xlsx(
  list(
    `All Performance Metrics (old)` = all_perf_metrics,
    `Whole Cohort Tuned Thresholds` = tuned_whole_out
  ),
  path = "Final Tables and Figures/Supplementary_Table_All_Model_Metrics_with_Tuned_Whole_updated_Feb2026.xlsx"
)


### See fixed threshold 
# ------------------------------------------------------------
# Helper: metrics from a prob column, fixed threshold
# ------------------------------------------------------------

# ============================================================================================
# SIMPLIFIED METRICS FUNCTION
# ============================================================================================
# WHAT THIS DOES:
# Sometimes a single specific threshold is desired (not a range to search)
# This function computes performance metrics at that SINGLE fixed threshold
#
# USE CASE:
# - A threshold T = 0.38 may be chosen based on clinical reasoning or prior validation
# - The goal is to know: sensitivity, specificity, accuracy, F1, etc. AT THAT THRESHOLD
# - Evaluation can occur across multiple cohorts (training, validation, whole, etc.)
#
# KEY METRICS COMPUTED:
# - Confusion matrix: TP, FP, TN, FN (how many of each classification type)
# - Sensitivity = TP/(TP+FN)     (recall, % of positives detected)
# - Specificity = TN/(TN+FP)     (% of negatives correctly identified)
# - Precision = TP/(TP+FP)       (PPV, if we predict positive, % correct)
# - NPV = TN/(TN+FN)             (if we predict negative, % correct)
# - Accuracy = (TP+TN)/N         (overall % correct)
# - Balanced Accuracy = mean(sensitivity, specificity)
# - F1 = harmonic mean of precision and sensitivity
#
# ---- Function: metrics_at_threshold() ----
# PURPOSE: Compute confusion matrix and derived metrics for ONE threshold
# INPUTS:
#   - dat:       Data frame with labels and probabilities
#   - prob_col:  Name of probability column
#   - thr:       Fixed threshold (e.g., 0.38)
#   - label_col: Name of label column (default: "MRD_truth")
# OUTPUTS:
#   - Tibble with ONE ROW showing all confusion matrix metrics

metrics_at_threshold <- function(dat, prob_col, thr, label_col = "MRD_truth") {
  stopifnot(prob_col %in% names(dat), label_col %in% names(dat))
  
  df <- dat %>%
    select(all_of(c(label_col, prob_col))) %>%
    rename(y = !!label_col, p = !!prob_col) %>%
    filter(!is.na(y) & !is.na(p)) %>%
    mutate(pred = as.integer(p >= thr))
  
  if (nrow(df) == 0) {
    return(tibble(
      Model = prob_col, Threshold = thr, N = 0,
      TP = NA_integer_, FP = NA_integer_, TN = NA_integer_, FN = NA_integer_,
      Prevalence = NA_real_,
      Sensitivity = NA_real_, Specificity = NA_real_,
      Precision = NA_real_, NPV = NA_real_,
      Accuracy = NA_real_, BalancedAccuracy = NA_real_, F1 = NA_real_
    ))
  }
  
  TP <- sum(df$pred == 1 & df$y == 1)
  FP <- sum(df$pred == 1 & df$y == 0)
  TN <- sum(df$pred == 0 & df$y == 0)
  FN <- sum(df$pred == 0 & df$y == 1)
  N  <- TP + FP + TN + FN
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  acc  <- (TP + TN) / N
  bal  <- mean(c(sens, spec), na.rm = TRUE)
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) {
    2 * prec * sens / (prec + sens)
  } else NA_real_
  
  tibble(
    Model = prob_col,
    Threshold = thr,
    N = N,
    TP = TP, FP = FP, TN = TN, FN = FN,
    Prevalence = mean(df$y == 1),
    Sensitivity = sens,
    Specificity = spec,
    Precision = prec,
    NPV = npv,
    Accuracy = acc,
    BalancedAccuracy = bal,
    F1 = f1
  )
}

# ------------------------------------------------------------
# Cohorts: reuse available objects; build a "whole, post-Dx" slice
# ------------------------------------------------------------
whole_postdx <- data_scored_masked %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis"))

# Target model and threshold
mdl <- "BM_base_zscore_prob"
thr <- 0.38

# Compute for each cohort
res_whole     <- metrics_at_threshold(whole_postdx, mdl, thr)
res_frontline <- metrics_at_threshold(frontline,      mdl, thr)
res_testing   <- metrics_at_threshold(test_cohort,    mdl, thr)

metrics_fixed_038 <- bind_rows(
  res_whole     %>% mutate(Cohort = "Whole (post-Dx)"),
  res_frontline %>% mutate(Cohort = "Frontline"),
  res_testing   %>% mutate(Cohort = "Non-frontline (test)")
) %>%
  relocate(Cohort, .before = Model) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(metrics_fixed_038)


# ============================================================================================
# DETAILED THRESHOLD EVALUATION - ALL THRESHOLDS, ALL METRICS
# ============================================================================================
# The sections below use a more comprehensive approach:
# 1. For EACH model, we evaluate at EVERY possible threshold
# 2. We extract sensitivity, specificity, and other metrics at each threshold
# 3. We then identify special thresholds (Youden-optimal, 95% specificity, etc.)
# 4. We create comprehensive tables showing model performance across thresholds
#
# This allows us to:
# - Plot ROC curves (sensitivity vs 1-specificity)
# - Select thresholds for specific operating points
# - Compare models at the same operating point
# - Create manuscript tables with multiple evaluation criteria

# Training-cohort threshold grid.
# This active manuscript section evaluates every observed probability threshold
# for each model, then extracts the Youden, >=95% specificity, and >=95%
# sensitivity operating points used in the metrics exports.

# 1) Filter to frontline cohort and compute prevalence
frontline <- data_scored_masked %>% 
  filter(Cohort == "Frontline") %>% 
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 

test_cohort <- data_scored_masked %>% 
  filter(Cohort == "Non-frontline") %>%
  filter(! timepoint_info %in% c("Baseline", "Diagnosis")) 

prev <- mean(frontline$MRD_truth == 1, na.rm = TRUE)

# 2) Identify all “_prob” columns (i.e. the different models)
prob_cols <- grep("_prob$", colnames(data_scored), value = TRUE)

### 3A) First get metrics for entire cohort 
# 3) Loop over models, build ROC & compute metrics at every threshold
all_metrics_rescored_primary <- map_dfr(prob_cols, function(prob_col) {
  # ---- Build ROC curve for this model ----
  # ROC curve shows: X-axis = 1 - specificity (false positive rate)
  #                  Y-axis = sensitivity (true positive rate)
  # As we vary the threshold from 1 (almost all negatives) to 0 (almost all positives)
  # sensitivity increases and specificity decreases
  # This creates the characteristic ROC curve shape
  roc_obj <- roc(
    response  = frontline$MRD_truth,    # TRUE labels (0 or 1)
    predictor = frontline[[prob_col]],  # Predicted probabilities
    levels    = c(0, 1),                # Class labels
    direction = "<"                     # Higher probability = positive class
  )
  
  # ---- Extract metrics at EVERY possible threshold ----
  # pROC's coords() with x="all" gives us all unique thresholds observed in the data
  # and sensitivity, specificity, PPV, NPV at each threshold
  # For ~100 samples, we might get 20-50 unique thresholds
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  ) %>% 
    as_tibble()
  
  # ---- Compute derived metrics from sensitivity/specificity ----
  # These are metrics that can be derived once we know sens and spec
  # accuracy depends on prevalence, balanced_accuracy doesn't
  # f1 is harmonic mean of precision (ppv) and sensitivity
  roc_df %>% 
    mutate(
      model        = prob_col,
      accuracy     = sensitivity * prev + specificity * (1 - prev),
      bal_accuracy = (sensitivity + specificity) / 2,
      f1           = 2 * sensitivity * ppv / (sensitivity + ppv)
    ) %>%
    select(model, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
})

# Get Youden cutoff
threshold_df <- tibble(
  model = paste0(names(all_thresholds), "_prob"),
  youden = unname(all_thresholds)
)

# Filter all_metrics_rescored_primary to the Youden-optimal row for each model
metrics_youden <- all_metrics_rescored_primary %>%
  inner_join(threshold_df, by = "model") %>%
  group_by(model) %>%
  # pick the row whose threshold minimizes |threshold - youden|
  slice_min(order_by = abs(threshold - youden), n = 1) %>%
  ungroup() %>%
  select(-youden)

print(metrics_youden)

### 3B) Now for each model, find the threshold with ≥95% specificity that maximizes sensitivity,
#    then compute sensitivity, specificity and accuracy at that threshold
metrics_at_95spec <- map_dfr(prob_cols, function(prob_col) {
  # build ROC
  roc_obj <- roc(
    response  = frontline$MRD_truth,
    predictor = frontline[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # get all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x        = "all",
    ret      = c("threshold","sensitivity","specificity"),
    transpose= FALSE
  ) %>% as_tibble()
  
  # filter to spec ≥ 0.95
  cand <- roc_df %>%
    filter(specificity >= 0.95)
  
  # if no threshold meets the bar, return NA
  if (nrow(cand) == 0) {
    return(tibble(
      model        = prob_col,
      threshold    = NA_real_,
      sensitivity  = NA_real_,
      specificity  = NA_real_,
      accuracy     = NA_real_
    ))
  }
  
  # pick the one with highest sensitivity
  best <- cand %>%
    arrange(desc(sensitivity)) %>%
    slice(1)
  
  # compute accuracy = sens*prev + spec*(1-prev)
  best <- best %>%
    mutate(
      accuracy = sensitivity * prev + specificity * (1 - prev),
      model    = prob_col
    ) %>%
    select(model, threshold, sensitivity, specificity, accuracy)
  
  return(best)
})

print(metrics_at_95spec)

## Now for sens


### Get more metrics 
# Prevalence: proportion of positives in data
prev <- mean(frontline$MRD_truth == 1, na.rm = TRUE)

metrics_at_95sens <- map_dfr(prob_cols, function(prob_col) {
  # 1) Build ROC
  roc_obj <- roc(
    response  = frontline$MRD_truth,
    predictor = frontline[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # 2) Extract all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold","sensitivity","specificity"),
    transpose = FALSE
  ) %>% as_tibble()
  
  # 3) Filter to sens ≥ 0.95
  #    2026-07-29 audit fix: the filter read `>= 0.94` while the section, the
  #    object name, and the downstream panel titles all say 95%.
  cand <- roc_df %>%
    filter(sensitivity >= 0.95)
  
  # 4) If none meet the bar, fall back to highest observed sens
  if (nrow(cand) == 0) {
    best <- roc_df %>% slice_max(sensitivity, n = 1)
  } else {
    # 5) Otherwise pick the one with highest specificity
    best <- cand %>% arrange(desc(specificity)) %>% slice(1)
  }
  
  # 6) Compute precision, F1, and balanced accuracy at that threshold
  best <- best %>%
    mutate(
      model      = prob_col,
      accuracy   = sensitivity * prev + specificity * (1 - prev),
      precision  = ifelse(
        (sensitivity * prev + (1 - specificity) * (1 - prev)) == 0,
        NA,
        (sensitivity * prev) / (sensitivity * prev + (1 - specificity) * (1 - prev))
      ),
      f1         = ifelse(
        (precision + sensitivity) == 0,
        NA,
        2 * (precision * sensitivity) / (precision + sensitivity)
      ),
      bal_accuracy = (sensitivity + specificity) / 2
    )
  
  best %>% select(model, threshold, sensitivity, specificity, accuracy, precision, f1, bal_accuracy)
})

print(metrics_at_95sens)


### Save this 
# Save metrics evaluated at Youden threshold
write_rds(metrics_youden, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold3_Feb2026.rds"))
write_csv(metrics_youden, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold3_Feb2026.csv"))

# Save metrics evaluated at fixed 95% specificity (if applicable)
write_rds(metrics_at_95spec, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec3_Feb2026.rds"))
write_csv(metrics_at_95spec, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec3_Feb2026.csv"))

write_rds(metrics_at_95sens, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens3_Feb2026.rds"))
write_csv(metrics_at_95sens, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens3_Feb2026.csv"))


## Tidy into one thing 
# 1. Tidy up and rename columns on each
youden_tbl <- metrics_youden %>%
  rename(
    Model                 = model,
    Threshold_Youden      = threshold,
    Sensitivity_Youden    = sensitivity,
    Specificity_Youden    = specificity,
    PPV_Youden            = ppv,
    NPV_Youden            = npv,
    Accuracy_Youden       = accuracy,
    BalancedAcc_Youden    = bal_accuracy,
    F1_Youden             = f1
  )

spec95_tbl <- metrics_at_95spec %>%
  rename(
    Model                   = model,
    Threshold_95Spec        = threshold,
    Sensitivity_95Spec      = sensitivity,
    Specificity_95Spec      = specificity,
    Accuracy_95Spec         = accuracy
  )

sens95_tbl <- metrics_at_95sens %>%
  rename(
    Model                   = model,
    Threshold_95Sens        = threshold,
    Sensitivity_95Sens      = sensitivity,
    Specificity_95Sens      = specificity,
    Accuracy_95Sens         = accuracy,
    Precision_95Sens        = precision,
    F1_95Sens               = f1,
    BalancedAccuracy_95Sens = bal_accuracy
  )


# 2. Join them all together
all_perf_metrics <- youden_tbl %>%
  left_join(spec95_tbl, by = "Model") %>%
  left_join(sens95_tbl, by = "Model")

all_perf_metrics <- all_perf_metrics %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# 3. Export as a single supplementary Excel (or CSV)
library(writexl)    # or use write.csv as preferred output format
write_xlsx(
  list(`All Performance Metrics` = all_perf_metrics),
  path = "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.xlsx"
)

# - or, if prefer CSV:
write.csv(all_perf_metrics, 
          "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.csv", 
          row.names = FALSE)

# -------------------------------------------------------------------------
# Manuscript output: Supplementary Table 5
#
# What this is:
#   Full refit/training-cohort model metrics table at Youden, fixed 95%
#   specificity, and fixed 95% sensitivity operating points.
#
# Why it is here:
#   This is the corrected source for final Supplementary Table 5. The old
#   filename says Supplementary Table 4; that naming drift is preserved because
#   it is documented in docs/supplementary_table_5_name_comparison.tsv.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.xlsx",
  artifact_id = "STABLE5",
  role = "workbook_xlsx",
  description = "Correct full-refit model metrics workbook used as Supplementary Table 5.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supplementary_Table_4_All_Model_Metrics_Refit6_Feb2026.csv",
  artifact_id = "STABLE5",
  role = "table_csv_companion",
  description = "CSV companion for the corrected Supplementary Table 5 full-refit model metrics.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)




# Non-frontline validation cohort threshold grid.
# This mirrors the training-cohort grid above, but evaluates model probabilities
# on the independent non-frontline/test cohort. The final exported table is
# copied as Supplementary Table 6.
### 3A) First get metrics for the non-frontline validation cohort
# 3) Loop over models, build ROC & compute metrics at every threshold
all_metrics_rescored_testing <- map_dfr(prob_cols, function(prob_col) {
  # build ROC object
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0, 1),
    direction = "<"
  )
  
  # extract thresholds with sens, spec, ppv, npv
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  ) %>% 
    as_tibble()
  
  # compute additional metrics
  roc_df %>% 
    mutate(
      model        = prob_col,
      accuracy     = sensitivity * prev + specificity * (1 - prev),
      bal_accuracy = (sensitivity + specificity) / 2,
      f1           = 2 * sensitivity * ppv / (sensitivity + ppv)
    ) %>%
    select(model, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
})

# Filter all_metrics_rescored_testing to the Youden-optimal row for each test model
metrics_youden_testing <- all_metrics_rescored_testing %>%
  inner_join(threshold_df, by = "model") %>%
  group_by(model) %>%
  # pick the row whose threshold minimizes |threshold - youden|
  slice_min(order_by = abs(threshold - youden), n = 1) %>%
  ungroup() %>%
  select(-youden)


# ============================================================================ #
# LEGACY / DIAGNOSTIC ONLY - DO NOT USE FOR MANUSCRIPT VALUES
# ---------------------------------------------------------------------------- #
# 2026-07-29 audit note.
#
# The two blocks below (`metrics_at_95spec_test`, `metrics_at_95sens_test`)
# search the *test-cohort* ROC for an operating point. That is threshold
# selection on the independent test set, and any sensitivity/specificity read
# off it is optimistically biased. It must never reach the manuscript.
#
# Retained only because the derived tables are useful as a diagnostic on how
# separable the test cohort is at all. They are now written to an explicitly
# quarantined output directory, and the panels that consumed them no longer
# overwrite manuscript figures. Previously the bar panel built from these
# tables wrote to
#   Final Tables and Figures/Supp_7B_classifier_performance_bar_test_cohort_updated3.png
# which is the same path the frozen-threshold branch writes at line ~1456, so
# running with CFWGS_RETRAIN_MODELS=1 silently replaced a manuscript figure
# with a test-set-optimised one under an identical filename.
#
# Manuscript operating points come exclusively from
#   selected_combo_thresholds_2026-02-16.rds  (Youden, full-training refit)
#   training_95sens_threshold()               (fixed sensitivity, training only)
# ============================================================================ #

quarantine_dir <- file.path(
  "Output_tables_2025", "LEAKY_test_threshold_diagnostics_DO_NOT_PUBLISH"
)
dir.create(quarantine_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(
  c(
    "Files in this directory were produced by searching the independent",
    "test-cohort ROC for an operating point. They are optimistically biased",
    "and must not be used for any reported value, figure, or table.",
    "See the audit note in 3_1_Optimize_cfWGS_thresholds.R."
  ),
  file.path(quarantine_dir, "README_DO_NOT_PUBLISH.txt")
)

# Test-cohort prevalence, defined before first use (see audit note below).
prev <- mean(test_cohort$MRD_truth == 1, na.rm = TRUE)

### 3B) Now for each model, find the threshold with ≥95% specificity that maximizes sensitivity,
#    then compute sensitivity, specificity and accuracy at that threshold
#    LEAKY: threshold chosen on the test cohort. Diagnostic use only.
metrics_at_95spec_test <- map_dfr(prob_cols, function(prob_col) {
  # build ROC
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # get all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x        = "all",
    ret      = c("threshold","sensitivity","specificity"),
    transpose= FALSE
  ) %>% as_tibble()
  
  # filter to spec ≥ 0.95
  cand <- roc_df %>%
    filter(specificity >= 0.95)
  
  # if no threshold meets the bar, return NA
  if (nrow(cand) == 0) {
    return(tibble(
      model        = prob_col,
      threshold    = NA_real_,
      sensitivity  = NA_real_,
      specificity  = NA_real_,
      accuracy     = NA_real_
    ))
  }
  
  # pick the one with highest sensitivity
  best <- cand %>%
    arrange(desc(sensitivity)) %>%
    slice(1)
  
  # compute accuracy = sens*prev + spec*(1-prev)
  best <- best %>%
    mutate(
      accuracy = sensitivity * prev + specificity * (1 - prev),
      model    = prob_col
    ) %>%
    select(model, threshold, sensitivity, specificity, accuracy)
  
  return(best)
})


## Now for sens

### Get more metrics
# Prevalence: proportion of positives in the TEST cohort.
#
# 2026-07-29 audit fix: this assignment previously sat *after* the
# `metrics_at_95spec_test` block that uses `prev`, so that block silently used
# the frontline/training prevalence still in scope from line ~4554. Moved above
# both blocks so each uses the test-cohort prevalence its label implies.
prev <- mean(test_cohort$MRD_truth == 1, na.rm = TRUE)

# LEAKY: threshold chosen on the test cohort. Diagnostic use only.
metrics_at_95sens_test <- map_dfr(prob_cols, function(prob_col) {
  # 1) Build ROC
  roc_obj <- roc(
    response  = test_cohort$MRD_truth,
    predictor = test_cohort[[prob_col]],
    levels    = c(0,1),
    direction = "<"
  )
  
  # 2) Extract all sens/spec pairs
  roc_df <- coords(
    roc_obj,
    x         = "all",
    ret       = c("threshold","sensitivity","specificity"),
    transpose = FALSE
  ) %>% as_tibble()
  
  # 3) Filter to sens ≥ 0.95
  #    2026-07-29 audit fix: the filter read `>= 0.94` while the section, the
  #    object name, and the downstream panel titles all say 95%.
  cand <- roc_df %>%
    filter(sensitivity >= 0.95)
  
  # 4) If none meet the bar, fall back to highest observed sens
  if (nrow(cand) == 0) {
    best <- roc_df %>% slice_max(sensitivity, n = 1)
  } else {
    # 5) Otherwise pick the one with highest specificity
    best <- cand %>% arrange(desc(specificity)) %>% slice(1)
  }
  
  # 6) Compute precision, F1, and balanced accuracy at that threshold
  best <- best %>%
    mutate(
      model      = prob_col,
      accuracy   = sensitivity * prev + specificity * (1 - prev),
      precision  = ifelse(
        (sensitivity * prev + (1 - specificity) * (1 - prev)) == 0,
        NA,
        (sensitivity * prev) / (sensitivity * prev + (1 - specificity) * (1 - prev))
      ),
      f1         = ifelse(
        (precision + sensitivity) == 0,
        NA,
        2 * (precision * sensitivity) / (precision + sensitivity)
      ),
      bal_accuracy = (sensitivity + specificity) / 2
    )
  
  best %>% select(model, threshold, sensitivity, specificity, accuracy, precision, f1, bal_accuracy)
})

print(metrics_at_95sens_test)


### Save this 
# Save metrics evaluated at Youden threshold
write_rds(metrics_youden_testing, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold_test_cohort3_Feb2026.rds"))
write_csv(metrics_youden_testing, file = file.path(outdir, "cfWGS_model_metrics_youden_threshold_test_cohort3_Feb2026.csv"))

# Save metrics evaluated at fixed 95% specificity (if applicable)
write_rds(metrics_at_95spec_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec_test_cohort3_Feb2026.rds"))
write_csv(metrics_at_95spec_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95spec_test_cohort3_Feb2026.csv"))

write_rds(metrics_at_95sens_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens_test_cohort3_Feb2026.rds"))
write_csv(metrics_at_95sens_test, file = file.path(outdir, "cfWGS_model_metrics_fixed_95sens_test_cohort3_Feb2026.csv"))


## Tidy into one thing 
# 1. Tidy up and rename columns on each
youden_tbl <- metrics_youden_testing %>%
  rename(
    Model                 = model,
    Threshold_Youden      = threshold,
    Sensitivity_Youden    = sensitivity,
    Specificity_Youden    = specificity,
    PPV_Youden            = ppv,
    NPV_Youden            = npv,
    Accuracy_Youden       = accuracy,
    BalancedAcc_Youden    = bal_accuracy,
    F1_Youden             = f1
  )

spec95_tbl <- metrics_at_95spec_test %>%
  rename(
    Model                   = model,
    Threshold_95Spec        = threshold,
    Sensitivity_95Spec      = sensitivity,
    Specificity_95Spec      = specificity,
    Accuracy_95Spec         = accuracy
  )

sens95_tbl <- metrics_at_95sens_test %>%
  rename(
    Model                   = model,
    Threshold_95Sens        = threshold,
    Sensitivity_95Sens      = sensitivity,
    Specificity_95Sens      = specificity,
    Accuracy_95Sens         = accuracy,
    Precision_95Sens        = precision,
    F1_95Sens               = f1,
    BalancedAccuracy_95Sens = bal_accuracy
  )


# 2. Join them all together
all_perf_metrics <- youden_tbl %>%
  left_join(spec95_tbl, by = "Model") %>%
  left_join(sens95_tbl, by = "Model")

# 3. Export as a single supplementary Excel (or CSV)
library(writexl)    # or use write.csv as preferred output format

all_perf_metrics <- all_perf_metrics %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_xlsx(
  list(`All Performance Metrics` = all_perf_metrics),
  path = "Final Tables and Figures/Supplementary_Table_5_All_Model_Metrics_Refit_Test_Cohort6_Feb2026.xlsx"
)

# - or, if prefer CSV:
write.csv(all_perf_metrics, 
          "Final Tables and Figures/Supplementary_Table_5_All_Model_Metrics_Refit_Test_Cohort6_Feb2026.csv", 
          row.names = FALSE)








### Now get whole cohort performance on rescored dataframe - skipped
# Define the specific model and a second fixed threshold
# model_name <- "BM_zscore_only_sites_prob"
# fixed_thr  <- 0.2211
# 
# # Pull out the Youden threshold for this model
# youden_thr <- threshold_df %>%
#   filter(model == model_name) %>%
#   pull(youden)
# 
# # Subset the full metrics table to just this model
# metrics_model <- all_metrics_rescored_primary %>%
#   filter(model == model_name)
# 
# # Find the row nearest the Youden cutoff
# metrics_youden <- metrics_model %>%
#   slice_min(order_by = abs(threshold - youden_thr), n = 1) %>%
#   mutate(type = "youden")
# 
# # Find the row nearest the fixed cutoff
# metrics_fixed <- metrics_model %>%
#   slice_min(order_by = abs(threshold - fixed_thr), n = 1) %>%
#   mutate(type = "fixed")
# 
# # Combine and print
# metrics_two <- bind_rows(metrics_youden, metrics_fixed) %>%
#   select(model, type, threshold, sensitivity, specificity, ppv, npv, accuracy, bal_accuracy, f1)
# 
# print(metrics_two)
# 
# write_rds(metrics_two, file = file.path(outdir, "cfWGS_model_metrics_fixed_z_score_model_deeper_LOD.rds"))
# write_csv(metrics_two, file = file.path(outdir, "cfWGS_model_metrics_fixed_z_score_model_deeper_LOD.csv"))
# write.csv(all_metrics_rescored_primary, file = "Final Tables and Figures/Suplementary_Table_3_All_BM_Zscore_Model_Performance_Cutoffs.csv")


# -----------------------------------------------------------------------------
# 10A. Performance Summaries & Exports - BM 
# -----------------------------------------------------------------------------

### First get combined color mapping 
# Combined colour mapping: Blood_* and BM_* share identical colours
col_map_all <- c(
  Blood_base                       = "#000000",  # black
  BM_base                          = "#000000",  # black
  
  Blood_base_zscore                = "#E69F00",  # orange
  BM_base_zscore                   = "#E69F00",  # orange
  
  Blood_plus_fragment              = "#56B4E9",  # sky blue
  BM_plus_fragment                 = "#56B4E9",  # sky blue
  
  Blood_plus_fragment_min          = "#009E73",  # bluish green
  BM_plus_fragment_min             = "#009E73",  # bluish green
  
  Blood_rate_only                  = "#CC79A7",  # magenta 
  BM_rate_only                     = "#CC79A7",  # magenta
  
  Blood_zscore_only_detection_rate = "#0072B2",  # blue
  BM_zscore_only_detection_rate    = "#0072B2",  # blue
  
  Blood_zscore_only_sites          = "#D55E00",  # vermillion
  BM_zscore_only_sites             = "#D55E00",   # vermillion
  
  # ---- Fragmentomics (re-use existing colours) ----
  Fragmentomics_full               = "#000000",  # black
  Fragmentomics_min                = "#E69F00",  # orange
  Fragmentomics_FS_only            = "#56B4E9",  # sky blue
  Fragmentomics_mean_coverage_only = "#009E73",  # bluish green
  Fragmentomics_prop_short_only    = "#0072B2",  # blue
  Fragmentomics_tumor_fraction_only= "#D55E00"   # vermillion

)

# 2) nice human-readable labels 
pretty_combo_names <- c(BM_base = "All Mut. Features", 
                        BM_base_zscore = "BM Sites + cVAF Z-score", 
                        BM_plus_fragment = "BM Mut. + Fragmentomics", 
                        BM_plus_fragment_min = "BM Mut. + Fragments (min)", 
                        BM_rate_only = "BM cVAF", 
                        BM_zscore_only_detection_rate = "BM cVAF Z-score", 
                        BM_zscore_only_sites = "BM Sites Z-score",
                        Blood_base                          = "All Mut. Features",
                        Blood_base_zscore                   = "Blood Sites + cVAF Z-score",
                        Blood_plus_fragment                 = "Blood Mut. + Fragmentomics",
                        Blood_plus_fragment_min             = "Blood Mut. + Fragments (min)",
                        Blood_rate_only                     = "Blood cVAF",
                        Blood_zscore_only_detection_rate    = "Blood cVAF Z-score",
                        Blood_zscore_only_sites             = "Blood Sites Z-score",
                        Fragmentomics_full                  = "FS + MeanCov + TF + PropShort",
                        Fragmentomics_min                   = "FS + MeanCov",
                        Fragmentomics_FS_only               = "Fragment Size Only",
                        Fragmentomics_mean_coverage_only    = "Mean Coverage Only",
                        Fragmentomics_prop_short_only       = "Prop. Short Fragments Only",
                        Fragmentomics_tumor_fraction_only   = "Tumor Fraction Only")

#### Plot to see what the model looks like
# For convenience
bm_obj      <- nested_bm_validation_updated2
models_list <- bm_obj$models            # named list of caret models
valid_df    <- train_bm                  # must contain MRD_truth + all predictors
bm_preds <- nested_bm_validation_updated2$outer_predictions



# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

# 1) Compute one ROC curve per combo from the outer‑fold preds to get the nested CV version
roc_dfs <- bm_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()


roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# ── Make one master level order and one named color map ───────────────────────
combo_levels <- levels(roc_df$combo)  # already AUC-ordered
combo_levels_chr <- as.character(combo_levels)
combo_levels_chr_BM <- combo_levels_chr #Reuse for test cohort(?)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr
)

# Single helper for a shared color scale
scale_cols_shared <- function(labels_vec, title = NULL) {
  scale_colour_manual(
    name   = title,
    values = col_map_all,             # <- fixed colours in AUC order
    limits = combo_levels_chr,      # <- force same legend/data order
    breaks = combo_levels_chr,
    labels = labels_vec,
    drop   = FALSE
  )
}
  
# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve on\nBM-Derived Features (nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_cols_shared(labels_vec = legend_labels)



# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- bm_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity for\nBM-Derived Features (mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  # scale_colour_manual(
  #   values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  # ) +
  scale_cols_shared(labels_vec = legend_labels) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

## Add legend
# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr
)


perf_plot_with_legend <- perf_plot + theme(    legend.position    = c(0.35, 0.2),  # inside, bottom-right
                                               legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
                                               legend.key.size    = unit(0.8, "lines"),
                                               legend.text        = element_text(size = 10)) +
  scale_cols_shared(labels_vec = labels_perf, title = "Model")


# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_updated3.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
    filename = "Final Tables and Figures/Performance_nested_folds_bm_updated2_only_performance.png",
  plot     = perf_plot_with_legend,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 3A
#
# What this is:
#   BM-derived nested-CV sensitivity/specificity performance panel.
#
# Why it is here:
#   This PNG is the script-generated component used for final Main Figure 3A.
#   It is generated from cached nested-CV metrics to avoid stochastic model
#   refitting drift.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Performance_nested_folds_bm_updated2_only_performance.png",
  artifact_id = "FIG3A",
  role = "figure_panel_png",
  description = "BM-derived nested-CV performance panel used as Main Figure 3A.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)

ggsave(
  filename = "Final Tables and Figures/ROC_plot_folds_bm_updated2.png",
  plot     = roc_plot,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 5A
#
# What this is:
#   Standalone pooled outer-fold ROC panel for BM-derived models.
#
# Why it is here:
#   The final Extended Data Figure 5 uses this as panel A.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/ROC_plot_folds_bm_updated2.png",
  artifact_id = "EDFIG5A",
  role = "figure_panel_png",
  description = "BM-derived pooled outer-fold ROC panel used as Extended Data Figure 5A.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## Above is figure 4A and B 

### Now additional stats 

# ──────────────────────────────────────────────────────────────────────────────
# 1. Whole cohort summary: pick just the two BM z-score models, relabel, and rename cols
# supp figure 5A
cv_tbl <- metrics_youden %>%
  filter(model %in% c("BM_zscore_only_sites_prob", 
                      "BM_zscore_only_detection_rate_prob", 
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob          = "Sites model",
                   BM_zscore_only_detection_rate_prob = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
  #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("BM_zscore_only_sites_prob",
                      "BM_zscore_only_detection_rate_prob",
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob           = "Sites model",
                   BM_zscore_only_detection_rate_prob  = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
   # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)

  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp5A_classifier_performance_bar_updated3.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 5B
#
# What this is:
#   BM classifier performance summary in the training/primary cohort.
#
# Why it is here:
#   The final Extended Data Figure 5 uses this as panel B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp5A_classifier_performance_bar_updated3.png",
  artifact_id = "EDFIG5B",
  role = "figure_panel_png",
  description = "BM classifier performance summary used as Extended Data Figure 5B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


 ## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("BM_zscore_only_sites_prob", 
                      "BM_zscore_only_detection_rate_prob", 
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob          = "Sites model",
                   BM_zscore_only_detection_rate_prob = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("BM_zscore_only_sites_prob",
                      "BM_zscore_only_detection_rate_prob",
                      "BM_base_zscore_prob")) %>%
  mutate(
    combo = recode(model,
                   BM_zscore_only_sites_prob           = "Sites model",
                   BM_zscore_only_detection_rate_prob  = "cVAF model",
                   BM_base_zscore_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp5A_classifier_performance_bar_test_cohort_updated3.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 5C
#
# What this is:
#   BM classifier performance summary in the test/non-frontline cohort.
#
# Why it is here:
#   The final Extended Data Figure 5 uses this as panel C.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp5A_classifier_performance_bar_test_cohort_updated3.png",
  artifact_id = "EDFIG5C",
  role = "figure_panel_png",
  description = "BM classifier test-cohort performance summary used as Extended Data Figure 5C.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


### Add points on ROC where models selected - skipped
# # ────────────────────────────────────────────────────────────────
# # A.  Build a tibble of (model, threshold) pairs to mark
# # ────────────────────────────────────────────────────────────────
# mark_tbl <- tribble(
#   ~combo,                          ~threshold,
#   "BM_zscore_only_sites",           0.501,
#   "BM_zscore_only_sites",           0.221,
#   "BM_zscore_only_detection_rate",  0.436
# )
# 
# # ────────────────────────────────────────────────────────────────
# # B.  For each row in mark_tbl, compute FPR & TPR
# # ────────────────────────────────────────────────────────────────
# marker_pts <- pmap_dfr(mark_tbl, function(combo, threshold) {
#   
#   # retrieve fitted caret model & its predicted probs
#   prob <- predict(models_list[[combo]], newdata = valid_df, type = "prob")[ , positive_class]
#   
#   # truth vector (pos/neg as in earlier code)
#   truth <- valid_df$MRD_truth
#   
#   # binary calls at this threshold
#   pred_pos <- prob >= threshold
#   truth_pos <- truth == "pos"
#   
#   # confusion parts
#   tp <- sum(pred_pos & truth_pos)
#   fn <- sum(!pred_pos & truth_pos)
#   fp <- sum(pred_pos & !truth_pos)
#   tn <- sum(!pred_pos & !truth_pos)
#   
#   tibble(
#     combo      = combo,
#     threshold  = threshold,
#     fpr        = fp / (fp + tn),
#     tpr        = tp / (tp + fn)
#   )
# })
# 
# # ensure same factor order as curves
# marker_pts$combo <- factor(marker_pts$combo, levels = levels(roc_df$combo))
# 
# # ────────────────────────────────────────────────────────────────
# # C.  Add to roc_plot
# # ────────────────────────────────────────────────────────────────
# selected_models <- c("BM_zscore_only_sites", "BM_zscore_only_detection_rate")
# 
# roc_plot_tmp <- roc_df %>%
#   filter(combo %in% selected_models) %>%               # keep only those two
#   ggplot(aes(x = fpr, y = tpr, colour = combo)) +
#   geom_line(size = 1) +
#   geom_abline(lty = 2, colour = "grey60") +
#   labs(
#     x      = "False-positive rate (1 − specificity)",
#     y      = "True-positive rate (sensitivity)",
#     colour = NULL,
#     title  = "Cross-validated ROC curves (nested 5×5 folds)"
#   ) +
#   theme_bw(12) +
#   theme(
#     legend.position    = c(0.63, 0.15),
#     legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
#     legend.key.size    = unit(0.8, "lines"),
#     legend.text        = element_text(size = 10),
#     panel.grid         = element_blank()
#   ) +
#   scale_colour_manual(
#     # pick the same palette entries but only for the two models
#     values = okabe_ito8[ match(selected_models, levels(roc_df$combo)) ],
#     labels = legend_labels[selected_models]
#   )
# roc_plot_final <- roc_plot_tmp +
#   geom_point(
#     data    = marker_pts,
#     aes(x = fpr, y = tpr, shape = combo),
#     fill    = "white",
#     colour  = "black",
#     stroke  = 0.9,
#     size    = 3,
#     inherit.aes = FALSE
#   ) +
#   scale_shape_manual(
#     name   = "Marked thresholds",
#     values = c(
#       BM_zscore_only_sites           = 21,  # circle
#       BM_zscore_only_detection_rate  = 24   # triangle
#     ),
#     labels = c(
#       BM_zscore_only_sites          = "Sites model",
#       BM_zscore_only_detection_rate = "cVAF model"
#     )
#   ) +
#   guides(
#     # a) the shape legend for marked thresholds
#     shape = guide_legend(
#       title       = "Marked thresholds",
#       order       = 1,       # draw this first
#       ncol        = 1,
#       byrow       = FALSE,
#       keywidth    = unit(1.2, "lines"),
#       keyheight   = unit(1.2, "lines"),
#       override.aes = list(
#         fill   = "white",
#         colour = "black",
#         size   = 4,
#         stroke = 0.9
#       )
#     ),
#     # b) the colour legend for ROC curves
#     colour = guide_legend(
#       title       = "Model (AUC)",
#       order       = 2,       # draw this second
#       ncol        = 1,
#       byrow       = FALSE,
#       keywidth    = unit(2,   "lines"),
#       keyheight   = unit(0.5, "lines"),
#       override.aes = list(
#         size = 1.5           # slightly thicker lines in the key
#       )
#     )
#   ) +
#   theme(
#     # remove legend background box
#     legend.background = element_blank(),
#     legend.position       = c(0.72, 0.24),   # x=0.85 (near right), y=0.25 (higher)
#     # shrink the space between items
#     legend.spacing.y  = unit(0.2, "cm"),
#     legend.key        = element_blank(),
#     # style legend titles & text
#     legend.title      = element_text(face = "bold", size = 10),
#     legend.text       = element_text(size = 9), 
#     # title bold 
#     plot.title   = element_text(face = "bold", size = 13)
#   )
# 
# 
# ggsave(
#   file.path("Final Tables and Figures/Supp5A_ROC_performance.png"),
#   plot   = roc_plot_final,
#   width  = 5,
#   height = 4.25,
#   dpi    = 600
# )
# 
# 





## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
bm_obj <- nested_bm_validation_updated2
thresh <- bm_obj$thresholds
mods   <- bm_obj$models[c("BM_zscore_only_detection_rate",
                          "BM_base_zscore")]

# nice labels
model_labs <- c(
  "cVAF model" = "BM_zscore_only_detection_rate",
  "Combined model"  = "BM_base_zscore"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_bm, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_bm$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs)) %>%
  mutate(model = fct_relevel(model, "cVAF model", "Combined model"))

## make grey 
col_low  <- "#f2f2f2"
col_high <- "#4a4a4a"

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  # scale_fill_viridis_c(
  #   option = "D",
  #   name   = "Count",
  #   begin  = 0.3,      # shift palette toward its lighter end
  #   end    = 0.9       # avoid the very darkest purples
  # ) +  
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = c("pos", "neg")) +
  labs(
    x = "cfWGS-MRD",
    y = "BM-MRD",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig4C_confusion_tables_primary_updated5.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 3B
#
# What this is:
#   BM model confusion-matrix panel for the training/primary cohort.
#
# Why it is here:
#   This PNG is one of the two confusion-matrix components used for final Main
#   Figure 3B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4C_confusion_tables_primary_updated5.png",
  artifact_id = "FIG3B",
  role = "figure_panel_png_training",
  description = "BM model training-cohort confusion matrix used as part of Main Figure 3B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_bm, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_bm$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs)) %>%
  mutate(model = fct_relevel(model, "cVAF model", "Combined model"))

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  # geom_text(aes(label = Count), size = 5) +
  # facet_wrap(~ model) +
  # scale_fill_viridis_c(
  #   option = "D",
  #   name   = "Count",
  #   begin  = 0.3,      # shift palette toward its lighter end
  #   end    = 0.9       # avoid the very darkest purples
  # ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "cfWGS-MRD",
    y = "BM-MRD",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  scale_y_discrete(limits = c("pos", "neg")) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig4C_confusion_tables_test5.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 3B
#
# What this is:
#   BM model confusion-matrix panel for the test/non-frontline cohort.
#
# Why it is here:
#   This PNG is the second confusion-matrix component used for final Main
#   Figure 3B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig4C_confusion_tables_test5.png",
  artifact_id = "FIG3B",
  role = "figure_panel_png_test",
  description = "BM model test-cohort confusion matrix used as part of Main Figure 3B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)




#### Now do the performance plot on validation cohort
# For convenience
valid_df    <- hold_bm                  # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)

# DO NOT reorder by validation AUC - lock to training order
roc_df$combo <- factor(roc_df$combo, levels = combo_levels_chr_BM)

auc_tbl <- roc_df %>% distinct(combo, auc)

# # 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr_BM], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr_BM, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr_BM
)


# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Cohort ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_cols_shared(labels_vec = legend_labels)



# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- bm_obj$validation_metrics %>%
  select(combo, 
         sens_mean = sens_valid,
         spec_mean = spec_valid) %>%
  # Add placeholder SD columns (set to 0 so error bars are invisible)
  mutate(
    sens_sd = 0,
    spec_sd = 0
  ) %>%
  # Keep only combos that made it into the ROC dataframe
  filter(combo %in% levels(roc_df$combo)) %>%
  # Ensure identical factor ordering for plotting consistency
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr
)
    
# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_cols_shared(labels_vec = labels_perf) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. Specificity of cfWGS\nModels in Test Cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 16)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_validation_updated5.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/4E_performance_nested_folds_bm_validation_updated5.png",
  plot     = perf_plot,
  width    = 5,
  height   = 4,
  dpi      = 500
)






#### Now re-do for blood derived muts 
## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
blood_obj      <- nested_blood_validation_updated2

thresh <- blood_obj$thresholds
mods   <- blood_obj$models[c("Blood_plus_fragment",
                          "Blood_zscore_only_sites")]

# nice labels
model_labs <- c(
  "Combined model" = "Blood_plus_fragment",
  "Sites model"  = "Blood_zscore_only_sites"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_blood, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_blood$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs)) %>%
  mutate(model = fct_relevel(model, "Sites model", "Combined model"))

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  # scale_fill_viridis_c(
  #   option = "D",
  #   name   = "Count",
  #   begin  = 0.3,      # shift palette toward its lighter end
  #   end    = 0.9       # avoid the very darkest purples
  # ) +  
  scale_x_discrete(position = "top") +
  labs(
    x = "cfWGS-MRD",
    y = "BM-MRD",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig5C_confusion_tables_primary_blood6.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 4B
#
# What this is:
#   Blood/cfDNA model confusion-matrix panel for the training/primary cohort.
#
# Why it is here:
#   This PNG is one of the two confusion-matrix components used for final Main
#   Figure 4B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5C_confusion_tables_primary_blood6.png",
  artifact_id = "FIG4B",
  role = "figure_panel_png_training",
  description = "Blood/cfDNA model training-cohort confusion matrix used as part of Main Figure 4B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_blood, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_blood$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs)) %>%
  mutate(model = fct_relevel(model, "Sites model", "Combined model"))

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  # scale_fill_viridis_c(
  #   option = "D",
  #   name   = "Count",
  #   begin  = 0.3,      # shift palette toward its lighter end
  #   end    = 0.9       # avoid the very darkest purples
  # ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "cfWGS-MRD",
    y = "BM-MRD",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Fig5C_confusion_tables_test_blood6.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 4B
#
# What this is:
#   Blood/cfDNA model confusion-matrix panel for the test/non-frontline cohort.
#
# Why it is here:
#   This PNG is the second confusion-matrix component used for final Main
#   Figure 4B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Fig5C_confusion_tables_test_blood6.png",
  artifact_id = "FIG4B",
  role = "figure_panel_png_test",
  description = "Blood/cfDNA model test-cohort confusion matrix used as part of Main Figure 4B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)



# Blood/cfDNA sensitivity/specificity bar plots for Figure 4B.
# This section formats the selected blood model metrics at the manuscript
# operating points and writes the panel component copied below.
cv_tbl <- metrics_youden %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_prob = "Combined model")  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp7A_classifier_performance_bar_updated_blood_muts2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 7B
#
# What this is:
#   Blood/cfDNA classifier performance summary in the training/primary cohort.
#
# Why it is here:
#   The final Extended Data Figure 7 uses this as panel B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp7A_classifier_performance_bar_updated_blood_muts2.png",
  artifact_id = "EDFIG7B",
  role = "figure_panel_png",
  description = "Blood/cfDNA classifier training-cohort performance summary used as Extended Data Figure 7B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the two models, relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("Blood_zscore_only_sites_prob", 
                      "Blood_rate_only_prob",
                      "Blood_plus_fragment_prob")) %>%
  mutate(
    combo = recode(model,
                   Blood_zscore_only_sites_prob           = "Sites model",
                   Blood_rate_only_prob  = "cVAF model",
                   Blood_plus_fragment_prob = "Combined model")
    ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc")) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
#
# 2026-07-29 audit fix: this panel is built from `metrics_at_95spec_test` /
# `metrics_at_95sens_test`, whose operating points were chosen on the test
# cohort. It previously wrote to
#   Final Tables and Figures/Supp_7B_classifier_performance_bar_test_cohort_updated3.png
# which is the identical path used by the frozen-threshold branch (line ~1456),
# so a run with CFWGS_RETRAIN_MODELS=1 silently substituted a test-optimised
# panel for the manuscript one. Redirected to the quarantine directory.
ggsave(
  file.path(quarantine_dir,
            "LEAKY_Supp_7B_classifier_performance_bar_test_cohort_DO_NOT_PUBLISH.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 7C
#
# What this is:
#   Blood/cfDNA classifier performance summary in the test/non-frontline cohort.
#
# Why it is here:
#   The final Extended Data Figure 7 uses this as panel C.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_7B_classifier_performance_bar_test_cohort_updated3.png",
  artifact_id = "EDFIG7C",
  role = "figure_panel_png",
  description = "Blood/cfDNA classifier test-cohort performance summary used as Extended Data Figure 7C.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)




### Make barplots for fragmentomic features (ORIGINAL full-cohort models only)
cv_tbl <- metrics_youden %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_Full_prob", 
                      "Fragmentomics_prop_short_only_Full_prob", 
                      "Fragmentomics_full_Full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_Full_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_Full_prob = "Prop. short model",
                   Fragmentomics_full_Full_prob = "Combined model")
  ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only the fragmentomics models (FULL COHORT), relabel, and rename
fix95_tbl <- metrics_at_95sens %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_Full_prob", 
                      "Fragmentomics_prop_short_only_Full_prob", 
                      "Fragmentomics_full_Full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_Full_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_Full_prob = "Prop. short model",
                   Fragmentomics_full_Full_prob = "Combined model")
  ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Train Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )
# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp9D_classifier_performance_bar_updated_frag2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9C
#
# What this is:
#   Fragmentomics classifier performance summary in the training/primary cohort.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel C.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp9D_classifier_performance_bar_updated_frag2.png",
  artifact_id = "EDFIG9C",
  role = "figure_panel_png",
  description = "Fragmentomics classifier training-cohort performance summary used as Extended Data Figure 9C.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## do for validation as well 
# ──────────────────────────────────────────────────────────────────────────────
cv_tbl <- metrics_youden_testing %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_Full_prob", 
                      "Fragmentomics_prop_short_only_Full_prob", 
                      "Fragmentomics_full_Full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_Full_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_Full_prob = "Prop. short model",
                   Fragmentomics_full_Full_prob = "Combined model")
    ) %>%
  # rename so we have *_mean and *_sd
  rename(
    #  AUC_mean  = auc_mean,   AUC_sd  = auc_sd,
    Sens_mean = sensitivity, # Sens_sd = sens_sd,
    Spec_mean = specificity, # Spec_sd = spec_sd,
    Acc_mean  = bal_accuracy,    # no sd for accuracy? if there is, rename _sd too,
    F1_mean = f1
  ) %>%
  select(combo, ends_with("_mean"), ends_with("_sd"))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Fixed-95% sensitivity metrics: pull only fragmentomics (FULL COHORT), relabel, and rename
fix95_tbl <- metrics_at_95sens_test %>%
  filter(model %in% c("Fragmentomics_mean_coverage_only_Full_prob", 
                      "Fragmentomics_prop_short_only_Full_prob", 
                      "Fragmentomics_full_Full_prob")) %>%
  mutate(
    combo = recode(model,
                   Fragmentomics_mean_coverage_only_Full_prob          = "Coverage model",
                   Fragmentomics_prop_short_only_Full_prob = "Prop. short model",
                   Fragmentomics_full_Full_prob = "Combined model")
    ) %>%
  select(combo, sensitivity, specificity, accuracy, bal_accuracy, f1) %>%
  rename(
    Sens_95 = sensitivity,
    Spec_95 = specificity,
    Acc_95  = bal_accuracy, 
    F1_95 = f1
  )

# ──────────────────────────────────────────────────────────────────────────────
# 3. Pivot CV summary to long form (one row per combo × metric)
cv_long <- cv_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = c("Metric","Stat"),
    names_sep = "_",        # splits e.g. "Sens_mean" → Metric="Sens", Stat="mean"
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = Stat,
    values_from = value
  )
# cv_long now has columns: combo, Metric, mean, sd

# ──────────────────────────────────────────────────────────────────────────────
# 4. Pivot fixed-95 table to long form
fix95_long <- fix95_tbl %>%
  pivot_longer(
    cols      = -combo,
    names_to  = "Metric95",
    values_to = "fixed"
  ) %>%
  mutate(
    Metric = case_when(
      Metric95 == "Sens_95" ~ "Sens",
      Metric95 == "Spec_95" ~ "Spec",
      Metric95 == "Acc_95"  ~ "Acc",
      Metric95 == "F1_95" ~ "F1",
      TRUE                  ~ NA_character_
    )
  ) %>% 
  select(combo, Metric, fixed)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Join CV + fixed-95
plot_df <- left_join(cv_long, fix95_long, by = c("combo","Metric"))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Restrict to only the metrics you want to plot:
#    here: F1, Sens, Spec, Acc
plot_df <- plot_df %>%
  filter(Metric %in% c("F1","Sens","Spec","Acc"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Give nicer facet labels
metric_labs <- c(
  AUC  = "F1",
  Sens = "Sensitivity",
  Spec = "Specificity",
  Acc  = "Bal. Accuracy"
)

# ──────────────────────────────────────────────────────────────────────────────
# 8. Build the plot
p_perf <- ggplot(plot_df, aes(x = combo, y = mean, fill = combo)) +
  # CV bar + error
  geom_col(width = 0.6) +
  # fixed‐95% dots, now mapped to a shape legend
  geom_point(aes(y = fixed, shape = "95% sensitivity"),
             size = 2, colour = "black") +
  
  # Tell ggplot how to draw that shape
  scale_shape_manual(
    name   = NULL,         # no legend title
    values = c("95% sensitivity" = 17)
  ) +
  
  facet_wrap(~ Metric, nrow = 1, labeller = labeller(Metric = metric_labs)) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.8, guide = "none") +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Test Cohort Performance at Youden and 95%\nSensitivity Thresholds"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10), # not as far down,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20)
    
  )

# ──────────────────────────────────────────────────────────────────────────────
# 9. Save
ggsave(
  file.path("Final Tables and Figures/Supp_Fig9F_classifier_performance_bar_test_cohort_updated2_frag2.png"),
  plot   = p_perf,
  width  = 5,
  height = 3.5,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9D
#
# What this is:
#   Fragmentomics classifier performance summary in the test/non-frontline
#   cohort.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel D.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Fig9F_classifier_performance_bar_test_cohort_updated2_frag2.png",
  artifact_id = "EDFIG9D",
  role = "figure_panel_png",
  description = "Fragmentomics classifier test-cohort performance summary used as Extended Data Figure 9D.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)




##### Now make the contingency tables
## add contingency table 
# ──────────────────────────────────────────────────────────────────────────────
# 2) pull thresholds and models
fragmentomics_obj      <- nested_fragmentomics_validation_updated3
models_list <- fragmentomics_obj$models            # named list of caret models
valid_df    <- train_fragmentomics                  # must contain MRD_truth + all predictors
fragmentomics_preds <- nested_fragmentomics_validation_updated3$outer_predictions

mods   <- fragmentomics_obj$models[c("Fragmentomics_mean_coverage_only",                             "Fragmentomics_min")]

# nice labels
model_labs <- c(
  "MM-sites coverage" = "Fragmentomics_mean_coverage_only",
  "Fragment score + MM-sites"  = "Fragmentomics_min"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = train_fragmentomics, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, train_fragmentomics$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Training Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Supp_Fig9E_confusion_tables_primary_fragmentomics4.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9E
#
# What this is:
#   Fragmentomics model confusion-matrix panel for the training/primary cohort.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel E.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Fig9E_confusion_tables_primary_fragmentomics4.png",
  artifact_id = "EDFIG9E",
  role = "figure_panel_png",
  description = "Fragmentomics model training-cohort confusion matrix used as Extended Data Figure 9E.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


## Now validation
# ──────────────────────────────────────────────────────────────────────────────
# 3) predict + build confusion‐table tibbles
cm_list <- imap(mods, function(mod, nm){
  # get threshold
  th <- thresh[[nm]]
  # predict probabilities
  probs <- predict(mod, newdata = hold_fragmentomics, type = "prob")[[ positive_class ]]
  # call class by threshold
  preds <- factor(if_else(probs >= th, "pos","neg"), levels = c("neg","pos"))
  # confusionMatrix
  cm <- confusionMatrix(preds, hold_fragmentomics$MRD_truth, positive = "pos")
  # turn table to tibble
  as_tibble(cm$table) %>%
    rename(Obs = Reference, Pred = Prediction, Count = n) %>%
    mutate(
      model = nm,
      PPV   = cm$byClass["Pos Pred Value"],
      NPV   = cm$byClass["Neg Pred Value"]
    )
})

cm_df <- bind_rows(cm_list) %>%
  mutate(model = fct_recode(model, !!!model_labs))

# make text switch to white on darker tiles for legibility (by facet)
cm_df <- cm_df %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(text_col = if_else(Count >= 0.7 * max(Count), "white", "black")) %>%
  dplyr::ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4) plot 2×2 tiles + add PPV/NPV text
p_tables <- ggplot(cm_df, aes(x = Pred, y = Obs, fill = Count)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(aes(label = Count, color = text_col), size = 5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_fill_gradient(
    name   = "Count",
    low    = col_low,
    high   = col_high,
    limits = c(0, max(cm_df$Count)),   # consistent scaling across facets
    oob    = squish
  )  +
  scale_color_identity() +
  # scale_fill_viridis_c(
  #   option = "D",
  #   name   = "Count",
  #   begin  = 0.3,      # shift palette toward its lighter end
  #   end    = 0.9       # avoid the very darkest purples
  # ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Predicted MRD status",
    y = "Observed MRD status",
    title = "Confusion Matrix at Youden Index in Test Cohort"
  ) +
  theme_minimal(base_size = 10) +
  # reverse the Obs axis so neg is at the top
  scale_y_discrete(limits = c("pos", "neg")) +
  theme(
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(size = 9, vjust = 0),
    axis.title        = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", hjust = 0.5)
  )

# ──────────────────────────────────────────────────────────────────────────────
# 5) save
ggsave(
  "Final Tables and Figures/Supp_Fig9F_confusion_tables_val_fragmentomics3.png",
  plot   = p_tables,
  width  = 5,
  height = 2.75,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9F
#
# What this is:
#   Fragmentomics model confusion-matrix panel for the test/non-frontline cohort.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel F.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Supp_Fig9F_confusion_tables_val_fragmentomics3.png",
  artifact_id = "EDFIG9F",
  role = "figure_panel_png",
  description = "Fragmentomics model test-cohort confusion matrix used as Extended Data Figure 9F.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)









## DO ROC curve on primary blood muts
# -----------------------------------------------------------------------------
# 10B. Performance Summaries & Exports - Blood 
# -----------------------------------------------------------------------------
#### Plot to see what the model looks like
# For convenience
models_list <- blood_obj$models            # named list of caret models
valid_df    <- train_blood                  # must contain MRD_truth + all predictors
blood_preds <- blood_obj$outer_predictions

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

roc_dfs <- blood_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

# Then recompute AUC table
auc_tbl <- roc_df %>%
  distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# ── Make one master level order and one named color map ───────────────────────
combo_levels <- levels(roc_df$combo)  # already AUC-ordered
combo_levels_chr <- as.character(combo_levels)

## Save for hold out test set 
combo_levels_blood <- combo_levels
combo_levels_chr_blood <- combo_levels_chr

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr
)


# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve on\ncfDNA-Derived Features (nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  ) +
  # scale_colour_manual(
  #   values = okabe_ito8[1:length(levels(roc_df$combo))],
  #   labels = legend_labels
  # )
  scale_cols_shared(labels_vec = legend_labels)



# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- blood_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # keep only combos that actually made it into the roc_df
  filter(combo %in% levels(roc_df$combo)) %>%
  # now drop any leftover unused levels
  mutate(combo = factor(combo, levels = levels(roc_df$combo))) 

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity for\ncfDNA-Derived Features (mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  # scale_colour_manual(
  #   values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  # ) +
  scale_cols_shared(labels_vec = legend_labels) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

## Add legend
# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr
)


perf_plot_with_legend <- perf_plot + theme(    legend.position    = c(0.37, 0.2),  # inside, bottom-right
                                               legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
                                               legend.key.size    = unit(0.8, "lines"),
                                               legend.text        = element_text(size = 10)) +
  scale_cols_shared(labels_vec = labels_perf, title = "Model")

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Fig_5A_updated_combined_ROC_and_performance_nested_folds_blood_features_updated3.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/Performance_nested_folds_blood_updated2_only_performance2.png",
  plot     = perf_plot_with_legend,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Main Figure 4A
#
# What this is:
#   Blood/cfDNA-derived nested-CV sensitivity/specificity performance panel.
#
# Why it is here:
#   This PNG is the script-generated component used for final Main Figure 4A.
#   It is generated from cached nested-CV metrics to avoid stochastic model
#   refitting drift.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Performance_nested_folds_blood_updated2_only_performance2.png",
  artifact_id = "FIG4A",
  role = "figure_panel_png",
  description = "Blood/cfDNA-derived nested-CV performance panel used as Main Figure 4A.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)

ggsave(
  filename = "Final Tables and Figures/ROC_plot_folds_blood_updated3.png",
  plot     = roc_plot,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 7A
#
# What this is:
#   Standalone pooled outer-fold ROC panel for blood/cfDNA-derived models.
#
# Why it is here:
#   The final Extended Data Figure 7 uses this as panel A.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/ROC_plot_folds_blood_updated3.png",
  artifact_id = "EDFIG7A",
  role = "figure_panel_png",
  description = "Blood/cfDNA-derived pooled outer-fold ROC panel used as Extended Data Figure 7A.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


### Add the metrics from dilution series
### Add points on ROC where models selected 
analysis_df <- data_scored_masked %>% filter(Sample_Code %in% train_blood$Sample_Code) %>%
  filter(!timepoint_info %in% c("Baseline", "Diagnosis"))

# Identify probability columns present in analysis_df
prob_cols <- grep("_prob$", names(analysis_df), value = TRUE)

# Build long-form predictions: one row per sample × model
# Map "Blood_zscore_only_sites_prob" -> combo "Blood_zscore_only_sites", etc.
pred_long <- analysis_df %>%
  select(MRD_truth, all_of(prob_cols)) %>%
  pivot_longer(cols = all_of(prob_cols),
               names_to = "model_prob",
               values_to = "prob") %>%
  mutate(
    combo = sub("_prob$", "", model_prob),
    # robust truth: accept 1/"pos" as positive
    truth_pos = MRD_truth %in% c(1, "pos"),
    truth_neg = MRD_truth %in% c(0, "neg")
  ) %>%
  select(combo, prob, truth_pos, truth_neg)

# Do ROC df on whole cohort refit
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# Keep factor levels consistent with the ROC plot (if roc_df exists already)
if (exists("roc_df")) {
  pred_long$combo <- factor(pred_long$combo, levels = levels(roc_df$combo))
}

# ------------------------------------------------------------------
# Thresholds to mark (same names as in the ROC 'combo' column)
# ------------------------------------------------------------------
mark_tbl <- tribble(
  ~combo,                      ~threshold,
  "Blood_zscore_only_sites",      0.5166693,
  "Blood_zscore_only_sites",      0.380,
  "Blood_plus_fragment",      0.4716880
)

# ------------------------------------------------------------------
# Compute sensitivity/specificity at each (combo, threshold)
# ------------------------------------------------------------------
marker_pts <- pmap_dfr(mark_tbl, function(combo, threshold) {
  df <- pred_long %>% filter(combo == !!combo)
  
  # Confusion counts using >= threshold rule
  pred_pos <- df$prob >= threshold
  tp <- sum(pred_pos & df$truth_pos, na.rm = TRUE)
  fn <- sum(!pred_pos & df$truth_pos, na.rm = TRUE)
  fp <- sum(pred_pos & df$truth_neg, na.rm = TRUE)
  tn <- sum(!pred_pos & df$truth_neg, na.rm = TRUE)
  
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  
  tibble(
    combo       = combo,
    threshold   = threshold,
    sensitivity = sens,
    specificity = spec,
    tpr         = sens,
    fpr         = 1 - spec
  )
})

# Align factor levels for plotting aesthetics
if (exists("roc_df")) {
  marker_pts$combo <- factor(marker_pts$combo, levels = levels(roc_df$combo))
}

# Optional: round for reporting
marker_pts <- marker_pts %>%
  mutate(across(c(sensitivity, specificity, tpr, fpr), ~ round(.x, 3)))

print(marker_pts)

# 2026-07-29 audit fix: this line previously reassigned `selected_models`,
# clobbering the named list of fitted caret objects held in that variable
# earlier in the script with a plain character vector. Renamed so the model list
# survives; any later use of `selected_models` as a model list would otherwise
# fail or misbehave silently.
ed7e_plot_combos <- c("Blood_zscore_only_sites", "Blood_plus_fragment")

roc_plot_tmp <- roc_df %>%
  filter(combo %in% ed7e_plot_combos) %>%              # keep only those two
  ggplot(aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,
    title  = "Full-Cohort ROC (Refit on all Training Cohort Samples)"
  ) +
  theme_bw(12) +
  theme(
    legend.position    = c(0.63, 0.15),
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid         = element_blank()
  ) +
  scale_colour_manual(
    # pick the same palette entries but only for the two models
    values = col_map_all,
    labels = c(
      Blood_zscore_only_sites          = "Sites model",
      Blood_plus_fragment = "Combined model"
    )
  )
roc_plot_final <- roc_plot_tmp +
  geom_point(
    data    = marker_pts,
    aes(x = fpr, y = tpr),
    fill    = "white",
    colour  = "darkgrey",
    stroke  = 0.9,
    size    = 2,
    inherit.aes = FALSE
  ) +
  guides(
    # a) the shape legend for marked thresholds
    shape = guide_legend(
      title       = "Marked thresholds",
      order       = 1,       # draw this first
      ncol        = 1,
      byrow       = FALSE,
      keywidth    = unit(1.2, "lines"),
      keyheight   = unit(1.2, "lines"),
      override.aes = list(
        fill   = "white",
        colour = "black",
        size   = 4,
        stroke = 0.9
      )
    ),
    # b) the colour legend for ROC curves
    colour = guide_legend(
      title       = "Model",
      order       = 2,       # draw this second
      ncol        = 1,
      byrow       = FALSE,
      keywidth    = unit(2,   "lines"),
      keyheight   = unit(0.5, "lines"),
      override.aes = list(
        size = 1.5           # slightly thicker lines in the key
      )
    )
  ) +
  theme(
    # remove legend background box
    legend.background = element_blank(),
    #  legend.position       = c(0.72, 0.24),   # x=0.85 (near right), y=0.25 (higher)
    legend.position       = "right",
    # shrink the space between items
    #    legend.spacing.y  = unit(0.2, "cm"),
    legend.key        = element_blank(),
    # style legend titles & text
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9), 
    # title bold 
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.title.position = "plot"
  )


ggsave(
  file.path("Final Tables and Figures/Supp7D_ROC_performance_blood_updated4.png"),
  plot   = roc_plot_final,
  width  = 5.5,
  height = 4,
  dpi    = 600
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 7E
#
# What this is:
#   Blood/cfDNA ROC performance panel with the fixed confirmatory marker used
#   in the final manuscript.
#
# Why it is here:
#   The final Extended Data Figure 7 uses this as panel E.
# -------------------------------------------------------------------------
ms_copy_artifact(
  # Stage a deterministic redraw from the frozen 34-row full-training/refit
  # inputs. The historical `updated4` file is not a safe manuscript source
  # because `valid_df` is mutable in this development branch.
  source_path = write_full_training_ed7e_panel(file.path(
    "Final Tables and Figures",
    "Supp7D_ROC_performance_blood_refit_main_cohort.png"
  )),
  artifact_id = "EDFIG7E",
  role = "figure_panel_png",
  description = "Blood/cfDNA full-training/refit ROC performance panel used as Extended Data Figure 7E.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


### Do ROC curve on validation for blood derived muts
#### Now do the performance plot on validation cohort
# For convenience
valid_df    <- hold_blood                # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)

# Rebuild roc_df without any Fragmentomics_ combos as re-trained later
roc_df <- bind_rows(roc_dfs) %>%
  filter(!grepl("^Fragmentomics_", combo))

# DO NOT reorder by validation AUC - lock to training order
roc_df$combo <- factor(roc_df$combo, levels = combo_levels_chr_blood)

auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr_blood], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr_blood, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr_blood
)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Cohort ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_cols_shared(labels_vec = legend_labels)
  # scale_colour_manual(
  #   values = okabe_ito8[1:length(levels(roc_df$combo))],
  #   labels = legend_labels
  # )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- blood_obj$validation_metrics %>%
  select(combo, 
         sens_mean = sens_valid,
         spec_mean = spec_valid) %>%
  # Add placeholder SD columns (set to 0 so error bars are invisible)
  mutate(
    sens_sd = 0,
    spec_sd = 0
  ) %>%
  # Keep only combos that made it into the ROC dataframe
  filter(combo %in% levels(roc_df$combo)) %>%
  # Ensure identical factor ordering for plotting consistency
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))


# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr_blood], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr_blood, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr_blood, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr_blood
)

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_cols_shared(labels_vec = labels_perf) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. Specificity of cfWGS\nModels in Test Cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 16)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))


# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_Fig_7C_combined_ROC_and_performance_nested_folds_blood_validation4.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/Supp_Fig_7C_performance_nested_folds_blood_validation_updated4.png",
  plot     = perf_plot,
  width    = 5,
  height   = 4,
  dpi      = 500
)


















# #### Now do ROC curve on validation cohort for BM muts - old 
# # For convenience
# valid_df    <- hold_bm                  # must contain MRD_truth + all predictors
# 
# # ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })
# 
# roc_df  <- bind_rows(roc_dfs)
# auc_tbl <- roc_df %>% distinct(combo, auc)
# 
# # 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)
# 
# # 4) make the plot, using a qualitative brewer palette
# roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
#   geom_line(size = 1) +
#   geom_abline(lty = 2, colour = "grey60") +
#   labs(
#     x      = "False-positive rate (1 − specificity)",
#     y      = "True-positive rate (sensitivity)",
#     colour = NULL,              # no title above the legend
#     title  = "Test Samples ROC Curves"
#   ) +
#   theme_bw(14) +
#   theme(
#     legend.position    = c(0.63, 0.15),  # inside, bottom-right
#     legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
#     legend.key.size    = unit(0.8, "lines"),
#     legend.text        = element_text(size = 10),
#     panel.grid = element_blank()
#   ) +
#   scale_colour_manual(
#     values = okabe_ito8[1:length(levels(roc_df$combo))],
#     labels = legend_labels
#   )
# 
# 
# # ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
# perf_df <- bm_obj$nested_metrics %>%
#   select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
#   # force the same ordering of combos
#   mutate(combo = factor(combo, levels = levels(roc_df$combo)))
# 
# # Make a named vector of labels
# labels_perf <- perf_df %>%
#   transmute(
#     combo,
#     pretty = pretty_combo_names[combo],
#     lbl = paste0(
#       pretty,
#       " (", percent(sens_mean, 1),
#       " sens, ", percent(spec_mean, 1),
#       " spec)"
#     )
#   ) %>%
#   { setNames(.$lbl, .$combo) }
# 
# # 1) Build performance plot
# perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
#   geom_point(size=3) +
#   geom_errorbarh(aes(
#     xmin = pmax(0, sens_mean - sens_sd),
#     xmax = pmin(1, sens_mean + sens_sd)
#   ), height=0.015) +
#   geom_errorbar(aes(
#     ymin = pmax(0, spec_mean - spec_sd),
#     ymax = pmin(1, spec_mean + spec_sd)
#   ), width=0.015) +
#   geom_vline(xintercept=0.5, lty=2, colour="grey80") +
#   geom_hline(yintercept=0.5, lty=2, colour="grey80") +
#   scale_colour_manual(
#     values = okabe_ito8[1:length(levels(perf_df$combo))],
#     labels = labels_perf
#   ) +
#   scale_x_continuous(limits=c(0,1)) +
#   scale_y_continuous(limits=c(0,1)) +
#   labs(
#     x     = "Mean sensitivity",
#     y     = "Mean specificity",
#     title = "Sensitivity vs. specificity of cfWGS models\nin the test cohort",
#     colour = NULL
#   ) +
#   theme_bw(12) +
#   theme(
#     panel.grid      = element_blank(),
#     legend.position = c(0.05, 0.05),
#     legend.justification = c(0,0),
#     legend.background = element_rect(fill = alpha("white",0.7), colour=NA),
#     plot.title      = element_text(hjust=0.5, face = "bold", size = 14)
#   )
# 
# # ── 3) Combine with roc_plot ────────────────────────────────────────────────
# combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))
# 
# # ── 4) Export ───────────────────────────────────────────────────────────────
# ggsave(
#   filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_bm_validation.png",
#   plot     = combined_plot,
#   width    = 12,
#   height   = 6,
#   dpi      = 500
# )
# 
# ggsave(
#   filename = "Final Tables and Figures/4E_performance_nested_folds_bm_validation_updated.png",
#   plot     = perf_plot,
#   width    = 6,
#   height   = 4,
#   dpi      = 500
# )
# 
# 
# 











## Add the other metrics
#### Now do on validation for blood
# # For convenience
# valid_df    <- hold_blood                  # must contain MRD_truth + all predictors
# 
# # ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })
# 
# # Rebuild roc_df without any Fragmentomics_ combos as re-trained later
# roc_df <- bind_rows(roc_dfs) %>%
#   filter(!grepl("^Fragmentomics_", combo))
# 
# auc_tbl <- roc_df %>% distinct(combo, auc)
# 
# # 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)
# 
# # 4) make the plot, using a qualitative brewer palette
# roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
#   geom_line(size = 1) +
#   geom_abline(lty = 2, colour = "grey60") +
#   labs(
#     x      = "False-positive rate (1 − specificity)",
#     y      = "True-positive rate (sensitivity)",
#     colour = NULL,              # no title above the legend
#     title  = "Hold-out samples ROC curves"
#   ) +
#   theme_bw(14) +
#   theme(
#     legend.position    = c(0.63, 0.15),  # inside, bottom-right
#     legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
#     legend.key.size    = unit(0.8, "lines"),
#     legend.text        = element_text(size = 10),
#     panel.grid = element_blank()
#   ) +
#   scale_colour_manual(
#     values = okabe_ito8[1:length(levels(roc_df$combo))],
#     labels = legend_labels
#   )
# 
# 
# # ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
# perf_df <- bm_obj$nested_metrics %>%
#   select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
#   # force the same ordering of combos
#   mutate(combo = factor(combo, levels = levels(roc_df$combo)))
# 
# # 1) Build performance plot
# perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(
#     xmin = pmax(0, sens_mean - sens_sd),
#     xmax = pmin(1, sens_mean + sens_sd)
#   ), height = 0.015) +
#   geom_errorbar(aes(
#     ymin = pmax(0, spec_mean - spec_sd),
#     ymax = pmin(1, spec_mean + spec_sd)
#   ), width = 0.015) +
#   geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
#   geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
#   scale_x_continuous(limits = c(0,1)) +
#   scale_y_continuous(limits = c(0,1)) +
#   labs(
#     x     = "Mean sensitivity",
#     y     = "Mean specificity",
#     title = "Hold-out samples performance"
#   ) +
#   # reuse the Okabe-Ito palette
#   scale_colour_manual(
#     values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
#   ) +
#   theme_bw(base_size = 14) +
#   theme(
#     panel.grid      = element_blank(),
#     legend.position = "none",
#     plot.title      = element_text(hjust = 0.5),
#     axis.title      = element_text(size = 13),
#     axis.text       = element_text(size = 11)
#   )
# 
# # ── 3) Combine with roc_plot ────────────────────────────────────────────────
# combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))
# 
# # ── 4) Export ───────────────────────────────────────────────────────────────
# ggsave(
#   filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_blood_features_validation.png",
#   plot     = combined_plot,
#   width    = 12,
#   height   = 6,
#   dpi      = 500
# )






# -----------------------------------------------------------------------------
# 10C. Performance Summaries & Exports - Fragmentomics
# -----------------------------------------------------------------------------
#### Plot to see what the model looks like
# For convenience
fragmentomics_obj      <- nested_fragmentomics_validation_updated3
models_list <- fragmentomics_obj$models            # named list of caret models
valid_df    <- train_fragmentomics                  # must contain MRD_truth + all predictors
fragmentomics_preds <- nested_fragmentomics_validation_updated3$outer_predictions

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
# roc_dfs <- imap(models_list,
#                 function(fit, label) {
#                   # caret::predict returns a data-frame; pull out the column for the + class
#                   prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
#                   
#                   roc_obj <- roc(response   = valid_df$MRD_truth,
#                                  predictor  = prob,
#                                  levels     = c("neg", "pos"),   # neg first, pos second
#                                  direction  = "<",
#                                  quiet      = TRUE)
#                   
#                   tibble(
#                     combo = label,
#                     fpr   = 1 - roc_obj$specificities,
#                     tpr   = roc_obj$sensitivities,
#                     auc   = as.numeric(auc(roc_obj))
#                   )
#                 })

## On the pooled outer fold predictions 
# 1) Compute one ROC curve per combo from the outer‑fold preds
# 1) Compute one ROC curve per combo from the outer‑fold preds to get the nested CV version
roc_dfs <- fragmentomics_preds %>% 
  group_by(combo) %>% 
  group_map(~{
    roc_obj <- roc(.x$truth, .x$prob,
                   levels    = c("neg","pos"),
                   direction = "<",
                   quiet     = TRUE)
    tibble(
      combo = .y$combo,
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      auc   = as.numeric(auc(roc_obj))
    )
  }) %>% 
  bind_rows()


roc_df  <- bind_rows(roc_dfs)
auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
auc_tbl <- auc_tbl %>% arrange(desc(auc))
roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# ── Make one master level order and one named color map ───────────────────────
combo_levels <- levels(roc_df$combo)  # already AUC-ordered
combo_levels_chr <- as.character(combo_levels)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr
)


# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Pooled Outer‑Fold ROC Curve\n(nested 5×5 CV)"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  ) +
  # scale_colour_manual(
  #   values = okabe_ito8[1:length(levels(roc_df$combo))],
  #   labels = legend_labels
  # )
  scale_cols_shared(labels_vec = legend_labels)


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- fragmentomics_obj$nested_metrics %>%
  select(combo, sens_mean, sens_sd, spec_mean, spec_sd) %>%
  # force the same ordering of combos
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(x = sens_mean, y = spec_mean, colour = combo)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height = 0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width = 0.015) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey80") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey80") +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x     = "Mean sensitivity (CV)",
    y     = "Mean specificity (CV)",
    title = "Fold‑Wise Sensitivity & Specificity\n(mean ± SD)"
  ) +
  # reuse the Okabe-Ito palette
  # scale_colour_manual(
  #   values = okabe_ito8[ seq_along(levels(perf_df$combo)) ]
  # ) +
  scale_cols_shared(labels_vec = legend_labels) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.title      = element_text(size = 13),
    axis.text       = element_text(size = 11)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))

## Add legend
# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr
)

perf_plot_with_legend <- perf_plot + theme(    legend.position    = c(0.4, 0.15),  # inside, bottom-right
                                               legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
                                               legend.key.size    = unit(0.8, "lines"),
                                               legend.text        = element_text(size = 10)) +
  scale_cols_shared(labels_vec = labels_perf) 


# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_fragmentomics_updated6.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/Performance_nested_folds_fragmentomics_only_performance4.png",
  plot     = perf_plot_with_legend,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9B
#
# What this is:
#   Fragmentomics nested-CV sensitivity/specificity performance diagnostic.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel B.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/Performance_nested_folds_fragmentomics_only_performance4.png",
  artifact_id = "EDFIG9B",
  role = "figure_panel_png",
  description = "Fragmentomics nested-CV performance panel used as Extended Data Figure 9B.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)

ggsave(
  filename = "Final Tables and Figures/ROC_plot_folds_fragmentomics_train3.png",
  plot     = roc_plot,
  width    = 6,
  height   = 6,
  dpi      = 500
)

# -------------------------------------------------------------------------
# Manuscript output: Extended Data Figure 9A
#
# What this is:
#   Standalone pooled outer-fold ROC panel for fragmentomics models.
#
# Why it is here:
#   The final Extended Data Figure 9 uses this as panel A.
# -------------------------------------------------------------------------
ms_copy_artifact(
  source_path = "Final Tables and Figures/ROC_plot_folds_fragmentomics_train3.png",
  artifact_id = "EDFIG9A",
  role = "figure_panel_png",
  description = "Fragmentomics pooled outer-fold ROC panel used as Extended Data Figure 9A.",
  script_name = "3_1_Optimize_cfWGS_thresholds.R"
)


combined_plot <- roc_plot + perf_plot_with_legend + plot_layout(ncol = 2, widths = c(1,1))

ggsave(
  filename = "Final Tables and Figures/combined_ROC_and_performance_nested_folds_fragmentomics_updated7_label.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)


### Redo for validation set 
#### Now do the performance plot on validation cohort
# For convenience
valid_df    <- hold_fragmentomics                # must contain MRD_truth + all predictors

# ── 1. ROC curves on the hold-out cohort ─────────────────────────────────
roc_dfs <- imap(models_list,
                function(fit, label) {
                  # caret::predict returns a data-frame; pull out the column for the + class
                  prob <- predict(fit, newdata = valid_df, type = "prob")[ , positive_class]
                  
                  roc_obj <- roc(response   = valid_df$MRD_truth,
                                 predictor  = prob,
                                 levels     = c("neg", "pos"),   # neg first, pos second
                                 direction  = "<",
                                 quiet      = TRUE)
                  
                  tibble(
                    combo = label,
                    fpr   = 1 - roc_obj$specificities,
                    tpr   = roc_obj$sensitivities,
                    auc   = as.numeric(auc(roc_obj))
                  )
                })

roc_df  <- bind_rows(roc_dfs)

# DO NOT reorder by validation AUC - lock to training order
roc_df$combo <- factor(roc_df$combo, levels = combo_levels_chr)

auc_tbl <- roc_df %>% distinct(combo, auc)

# 1) reorder combos by AUC
# auc_tbl <- auc_tbl %>% arrange(desc(auc))
# roc_df$combo <- factor(roc_df$combo, levels = auc_tbl$combo)

# 3) build a vector of final legend labels: “Pretty Name, AUC = 0.83”
legend_labels <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], ", AUC = ",
    formatC(auc_tbl$auc[match(combo_levels_chr, as.character(auc_tbl$combo))],
            digits = 2, format = "f")
  ),
  combo_levels_chr
)

# 4) make the plot, using a qualitative brewer palette
roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, colour = "grey60") +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = NULL,              # no title above the legend
    title  = "Test Cohort ROC Curves"
  ) +
  theme_bw(14) +
  theme(
    legend.position    = c(0.63, 0.15),  # inside, bottom-right
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size    = unit(0.8, "lines"),
    legend.text        = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  scale_cols_shared(labels_vec = legend_labels)
# scale_colour_manual(
#   values = okabe_ito8[1:length(levels(roc_df$combo))],
#   labels = legend_labels
# )


# ── 2) Prepare perf_df with the same factor‐ordering as roc_df ───────────────
perf_df <- fragmentomics_obj$validation_metrics %>%
  select(combo, 
         sens_mean = sens_valid,
         spec_mean = spec_valid) %>%
  # Add placeholder SD columns (set to 0 so error bars are invisible)
  mutate(
    sens_sd = 0,
    spec_sd = 0
  ) %>%
  # Keep only combos that made it into the ROC dataframe
  filter(combo %in% levels(roc_df$combo)) %>%
  # Ensure identical factor ordering for plotting consistency
  mutate(combo = factor(combo, levels = levels(roc_df$combo)))

# Make a named vector of labels
labels_perf <- setNames(
  paste0(
    pretty_combo_names[combo_levels_chr], " (",
    scales::percent(perf_df$sens_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " sens, ",
    scales::percent(perf_df$spec_mean[match(combo_levels_chr, as.character(perf_df$combo))], 1),
    " spec)"
  ),
  combo_levels_chr
)

# 1) Build performance plot
perf_plot <- ggplot(perf_df, aes(sens_mean, spec_mean, colour = combo)) +
  geom_point(size=3) +
  geom_errorbarh(aes(
    xmin = pmax(0, sens_mean - sens_sd),
    xmax = pmin(1, sens_mean + sens_sd)
  ), height=0.015) +
  geom_errorbar(aes(
    ymin = pmax(0, spec_mean - spec_sd),
    ymax = pmin(1, spec_mean + spec_sd)
  ), width=0.015) +
  geom_vline(xintercept=0.5, lty=2, colour="grey80") +
  geom_hline(yintercept=0.5, lty=2, colour="grey80") +
  scale_cols_shared(labels_vec = labels_perf) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    x     = "Mean sensitivity",
    y     = "Mean specificity",
    title = "Sensitivity vs. Specificity of cfWGS\nModels in Test Cohort",
    colour = NULL
  ) +
  theme_bw(12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0,0),
    legend.background = element_rect(
      fill = alpha("white", 0.7),  # lower alpha → more transparent
      colour = NA
    ),
    plot.title      = element_text(hjust=0.5, face = "bold", size = 16)
  )

# ── 3) Combine with roc_plot ────────────────────────────────────────────────
combined_plot <- roc_plot + perf_plot + plot_layout(ncol = 2, widths = c(1,1))


# ── 4) Export ───────────────────────────────────────────────────────────────
ggsave(
  filename = "Final Tables and Figures/Supp_Fig_9C_combined_ROC_and_performance_nested_folds_fragmentomics_validation3.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  dpi      = 500
)

ggsave(
  filename = "Final Tables and Figures/Supp_Fig_9C_performance_nested_folds_fragmentomics_validation_updated3.png",
  plot     = perf_plot,
  width    = 5,
  height   = 4,
  dpi      = 500
)




















# -----------------------------------------------------------------------------
# 11. Additional Diagnostic Plots
# -----------------------------------------------------------------------------

fig_dir <- file.path(outdir, "Diagnostic_Figures_2025")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Helper: pretty label for each combo
pretty_combo_names <- c(
  BM_zscore_only_sites            = "BM Sites Z-score",
  BM_zscore_only_detection_rate   = "BM cVAF Z-score",
  Blood_rate_only                 = "Blood cVAF",
  Blood_plus_fragment_min         = "Blood + Fragments (min)",
  Blood_zscore_only_detection_rate= "Blood cVAF Z-score",
  Fragmentomics_min               = "Fragmentomics (FS + Cov.)",
  Fragmentomics_mean_coverage_only= "Fragmentomics Coverage"
)

# choose the same frontline MRD set you trained on
valid_df <- train_df          # contains MRD_truth + ALL predictors
positive_class <- "pos"

# loop over every wanted combo
for (cmb in wanted) {
  
  fit <- selected_models[[cmb]]
  thr <- selected_thr[[cmb]]
  
  # figure out which predictors this model requires
  preds <- setdiff(colnames(fit$trainingData), ".outcome")
  
  # restrict to rows that have all predictors
  keep <- complete.cases(valid_df[, preds])
  if (!any(keep)) {
    message("Skipping ", cmb, " (no rows with complete predictors).")
    next
  }
  
  probs <- predict(fit,
                   newdata = valid_df[keep, preds, drop = FALSE],
                   type = "prob")[[positive_class]]
  
  truth_fac <- valid_df$MRD_truth[keep]          # factor "neg"/"pos"
  truth     <- ifelse(truth_fac == "pos", 1, 0)  # numeric 0/1
  
  # ── 11.1 Calibration plot (deciles) ──────────────────────────────────────
  cal_df <- tibble(truth = truth, pred = probs) %>%
    mutate(dec = ntile(pred, 10)) %>%
    group_by(dec) %>%
    summarise(mean_pred = mean(pred), obs_rate = mean(truth), .groups = "drop")
  
  cal_plot <- ggplot(cal_df, aes(mean_pred, obs_rate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_abline(lty = 2, col = "grey60") +
    scale_x_continuous(labels = percent_format(1)) +
    scale_y_continuous(labels = percent_format(1)) +
    labs(
      x = "Mean predicted probability",
      y = "Observed MRD-positive rate",
      title = paste0("Calibration: ", pretty_combo_names[cmb])
    ) +
    theme_bw(14) + theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("calibration_", cmb, ".png")),
    plot = cal_plot, width = 6, height = 5, dpi = 300
  )
  
  # ── 11.2 Precision–Recall curve ─────────────────────────────────────────
  pr <- pr.curve(scores.class0 = probs[truth == 1],
                 scores.class1 = probs[truth == 0],
                 curve = TRUE)
  pr_df <- data.frame(recall = pr$curve[, 1],
                      precision = pr$curve[, 2])
  
  pr_plot <- ggplot(pr_df, aes(recall, precision)) +
    geom_line(size = 1, col = "#E69F00") +                 # orange
    labs(
      x = "Recall (sensitivity)",
      y = "Precision (PPV)",
      title = paste0("PR curve: ", pretty_combo_names[cmb]),
      subtitle = sprintf("AUPRC = %.2f", pr$auc.integral)
    ) +
    theme_bw(14) + theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("pr_curve_", cmb, ".png")),
    plot = pr_plot, width = 6, height = 5, dpi = 300
  )
  
  # ── 11.3 Sensitivity & specificity vs cutoff ────────────────────────────
  roc_obj    <- roc(response = truth_fac, predictor = probs,
                    quiet = TRUE, levels = c("neg", "pos"))
  thresholds <- seq(0, 1, by = 0.01)
  
  # coords() returns a matrix with columns "sensitivity" & "specificity"
  metric_df  <- as_tibble(
    coords(
      roc_obj,
      x         = thresholds,
      input     = "threshold",
      ret       = c("sensitivity", "specificity"),
      transpose = FALSE
    )
  )
  
  # bind the numeric threshold vector back on
  coords_df <- tibble(threshold = thresholds) %>%
    bind_cols(metric_df)
  
  thr_plot <- ggplot(coords_df, aes(x = threshold)) +
    geom_line(aes(y = sensitivity,   col = "Sensitivity"),   size = 1) +
    geom_line(aes(y = specificity,   col = "Specificity"),   size = 1) +
    scale_color_manual(
      NULL,
      values = c("Sensitivity" = "#009E73", "Specificity" = "#0072B2")
    ) +
    scale_x_continuous(labels = percent_format(0.1)) +
    scale_y_continuous(labels = percent_format(1)) +
    labs(
      x     = "Probability cutoff",
      y     = "Value",
      title = paste0("Sens/Spec vs cutoff: ", pretty_combo_names[cmb])
    ) +
    theme_bw(14) +
    theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(fig_dir,
                         paste0("sens_spec_vs_cutoff_", cmb, ".png")),
    plot   = thr_plot,
    width  = 6,
    height = 5,
    dpi    = 300
  )
  
  message("Exported diagnostics for ", cmb)
}














#### Optional exploratory diagnostics
#### Below here is diagnostic ROC smoothing/sample-count code. It is excluded
#### from active manuscript figures, tables, thresholds, and reported model
#### metrics. Keep disabled during command-line manuscript regeneration unless
#### an analyst is intentionally auditing exploratory diagnostics.

run_non_manuscript_diagnostics <- FALSE

if (isTRUE(run_non_manuscript_diagnostics)) {

### Now try smoothed 
# ── 1. Compute binormal-smoothed ROC curves ────────────────────────────────
smoothed_rocs <- imap_dfr(models_list, function(fit, label) {
  # 1a) get predicted probabilities on the primary cohort
  prob <- predict(fit, newdata = valid_df, type = "prob")[, positive_class]
  
  # 1b) raw ROC object
  roc_obj <- roc(
    response  = valid_df$MRD_truth,
    predictor = prob,
    levels    = c("neg","pos"),
    direction = "<",
    quiet     = TRUE
  )
  
  # 1c) binormal smoothing
  roc_s <- smooth(roc_obj, method = "binormal")
  
  # 1d) pack into a tibble
  tibble(
    combo = label,
    fpr   = 1 - roc_s$specificities,
    tpr   =     roc_s$sensitivities,
    auc   = as.numeric(auc(roc_obj))
  )
})

# ── 2. Extract one AUC per combo for legend labels ─────────────────────────
auc_tbl <- smoothed_rocs %>%
  distinct(combo, auc) %>%
  arrange(combo)

# ── 3. Plot ────────────────────────────────────────────────────────────────
ggplot(smoothed_rocs, aes(x = fpr, y = tpr, color = combo)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, color = "grey60") +
  scale_color_manual(
    values = scales::hue_pal()(nrow(auc_tbl)),
    labels = setNames(
      sprintf("%s | AUC = %.2f", auc_tbl$combo, auc_tbl$auc),
      auc_tbl$combo
    )
  ) +
  labs(
    x      = "False-positive rate (1 − specificity)",
    y      = "True-positive rate (sensitivity)",
    colour = "Combo",
    title  = "ROC curves on primary cohort\n(binormal smoothing)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )





### Get more info on total number of samples
### Now check sample counts 

# Define export directory
export_dir <- "Output_tables_2025"

## Edit to ones we have cfDNA analyzed on in primary cohorts
dat_valid <- dat %>%
  filter(Patient %in% cohort_df$Patient) %>% filter(!is.na(WGS_Tumor_Fraction_Blood_plasma_cfDNA)) %>% filter(!is.na(FS))

## 2. Per‑patient summaries
patient_summ <- dat_valid %>%
  group_by(Patient) %>%
  summarise(
    n_samples   = dplyr::n(),                                   # how many cfDNA samples
    first_date  = min(Date, na.rm = TRUE),
    last_date   = max(Date, na.rm = TRUE),
    follow_days = as.numeric(difftime(last_date, first_date, units = "days")),
    .groups = "drop"
  )

## 4. Study‑level numbers
n_patients        <- n_distinct(dat_valid$Patient)        # should be 51 (check)
total_samples     <- nrow(dat_valid)
median_per_pt     <- median(patient_summ$n_samples)
range_per_pt      <- range(patient_summ$n_samples)
median_follow     <- median(patient_summ$follow_days)
range_follow      <- range(patient_summ$follow_days)

## 5. Auto‑generate the sentence
glue(
  "This study comprises cfWGS analysis of {total_samples} blood cfDNA samples ",
  "from {n_patients} MM patients (median {median_per_pt} samples per patient, ",
  "range {range_per_pt[1]}–{range_per_pt[2]}), collected over a median follow‑up ",
  "of {median_follow} days (range {range_follow[1]}–{range_follow[2]})."
)


### See how many additional patients had blood samples but not BM 
dat %>%
  filter(!is.na(zscore_blood), is.na(zscore_BM)) %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  distinct(Patient, Cohort) %>%
  count(Cohort, name = "n_patients_with_blood_only")

dat %>%
  filter(!is.na(zscore_blood), is.na(zscore_BM)) %>%
  filter(!timepoint_info %in% c("Diagnosis", "Baseline")) %>%
  filter(!is.na(MRD_truth)) %>%
  distinct(Sample_Code, Cohort) %>%
  count(Cohort, name = "n_samples_with_blood_only")
}
