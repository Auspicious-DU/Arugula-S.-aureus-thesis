# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# ============================================================
# 01_PRE_DT_12h_phenotype_harmonisation.R
#
# Purpose:
# 1. Combine archived 2rep and 3rep threshold-summary tables
# 2. Recover the original same-plate MeOH threshold references
# 3. Harmonise Delta-T to the predefined 0–12 h assay window
# 4. Preserve strong-delay / no-crossing individuals
# 5. Match only the 589 individuals used in the final
#    LC-MS–phenotype dataset
# 6. Generate a clean Delta-T phenotype table for formal
#    downstream analysis
#
# IMPORTANT:
# - Original phenotype columns are NOT modified.
# - Original batch-processing scripts are NOT rerun.
# - The 12-h transformation is applied at the Individual level.
# ============================================================


# ============================================================
# 0. Setup
# ============================================================

source(file.path(repo_dir, "scripts", "phenotype", "00_DT_setup.R"))


# ============================================================
# 1. Output folders
# ============================================================

DT_preprocessed_dir <- file.path(DT_output_dir, "00_preprocessed_data")

DT_pre_QC_dir <- file.path(DT_preprocessed_dir, "QC")

make_DT_dir(
  DT_preprocessed_dir
)

make_DT_dir(
  DT_pre_QC_dir
)


# ============================================================
# 2. Constants
# ============================================================

DT_ANALYSIS_WINDOW_H <- 12


# ============================================================
# 3. Helper functions
# ============================================================

clean_sample_id_DT <- function(x) {
  
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all(
      "：",
      ":"
    ) %>%
    stringr::str_replace_all(
      "[:/\\\\]",
      "."
    ) %>%
    stringr::str_replace_all(
      "\\s+",
      ""
    ) %>%
    stringr::str_replace_all(
      "\\.+",
      "."
    ) %>%
    stringr::str_replace_all(
      "^\\.|\\.$",
      ""
    )
}


normalise_batch_id_DT <- function(x) {
  
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    
    # Normalise different dash characters
    stringr::str_replace_all(
      "[\u2013\u2014\u2212]",
      "-"
    ) %>%
    
    stringr::str_replace_all(
      "\\s+",
      ""
    ) %>%
    
    stringr::str_to_lower()
}


# ============================================================
# 4. Read archived threshold-summary tables
# ============================================================

threshold_2rep_DT <- read_required_csv_DT(
  DT_threshold_2rep_file,
  "2rep threshold-summary table"
) %>%
  dplyr::mutate(
    replicate_type_source = "2rep"
  )


threshold_3rep_DT <- read_required_csv_DT(
  DT_threshold_3rep_file,
  "3rep threshold-summary table"
) %>%
  dplyr::mutate(
    replicate_type_source = "3rep"
  )


dim(
  threshold_2rep_DT
)

dim(
  threshold_3rep_DT
)


# ============================================================
# 5. Check required columns
# ============================================================

required_threshold_cols_DT <- c(
  "FileName",
  "TestDate",
  "AssayBatch",
  "BatchID",
  "threshold",
  "threshold_name",
  "SampleID",
  "Group",
  "T_mean",
  "T_sd",
  "MeOH_T_mean",
  "dT"
)


missing_2rep_cols_DT <- setdiff(
  required_threshold_cols_DT,
  names(
    threshold_2rep_DT
  )
)

missing_3rep_cols_DT <- setdiff(
  required_threshold_cols_DT,
  names(
    threshold_3rep_DT
  )
)


if (
  length(
    missing_2rep_cols_DT
  ) > 0
) {
  
  stop(
    "Missing columns in 2rep threshold table: ",
    paste(
      missing_2rep_cols_DT,
      collapse = ", "
    )
  )
  
}


if (
  length(
    missing_3rep_cols_DT
  ) > 0
) {
  
  stop(
    "Missing columns in 3rep threshold table: ",
    paste(
      missing_3rep_cols_DT,
      collapse = ", "
    )
  )
  
}


# ============================================================
# 6. Combine 2rep and 3rep threshold data
# ============================================================

threshold_all_DT <- dplyr::bind_rows(
  threshold_2rep_DT,
  threshold_3rep_DT
) %>%
  dplyr::mutate(
    
    sample_id_clean =
      clean_sample_id_DT(
        SampleID
      ),
    
    batch_id_clean =
      normalise_batch_id_DT(
        BatchID
      ),
    
    threshold_code =
      dplyr::case_when(
        
        threshold == 0.05 ~
          "005",
        
        threshold == 0.10 ~
          "010",
        
        threshold == 0.15 ~
          "015",
        
        TRUE ~
          NA_character_
      )
  )


# Check threshold definitions
threshold_definition_check_DT <-
  threshold_all_DT %>%
  dplyr::distinct(
    threshold,
    threshold_name,
    threshold_code
  ) %>%
  dplyr::arrange(
    threshold
  )

threshold_definition_check_DT


if (
  any(
    is.na(
      threshold_all_DT$threshold_code
    )
  )
) {
  
  stop(
    "Unexpected threshold values were found."
  )
  
}


# ============================================================
# 7. Basic combined-data QC
# ============================================================

threshold_combined_summary_DT <-
  threshold_all_DT %>%
  dplyr::summarise(
    
    n_rows =
      dplyr::n(),
    
    n_files =
      dplyr::n_distinct(
        FileName
      ),
    
    n_batches =
      dplyr::n_distinct(
        batch_id_clean
      ),
    
    n_sample_ids =
      dplyr::n_distinct(
        sample_id_clean[
          Group == "Sample"
        ]
      )
  )

threshold_combined_summary_DT


threshold_source_summary_DT <-
  threshold_all_DT %>%
  dplyr::group_by(
    replicate_type_source
  ) %>%
  dplyr::summarise(
    
    n_rows =
      dplyr::n(),
    
    n_files =
      dplyr::n_distinct(
        FileName
      ),
    
    n_batches =
      dplyr::n_distinct(
        batch_id_clean
      ),
    
    .groups = "drop"
  )

threshold_source_summary_DT


# ============================================================
# 8. Check uniqueness of batch/sample/threshold records
# ============================================================

threshold_duplicate_check_DT <-
  threshold_all_DT %>%
  dplyr::count(
    batch_id_clean,
    sample_id_clean,
    threshold,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

threshold_duplicate_check_DT


if (
  nrow(
    threshold_duplicate_check_DT
  ) > 0
) {
  
  stop(
    "Duplicated BatchID + SampleID + threshold records found."
  )
  
}


# ============================================================
# 9. Explicit MeOH-control QC
# ============================================================

meoh_threshold_DT <- threshold_all_DT %>%
  dplyr::filter(
    SampleID == "MeOH",
    Group == "Control"
  )


meoh_control_count_DT <- meoh_threshold_DT %>%
  dplyr::count(
    batch_id_clean,
    threshold,
    name = "n_control_rows"
  )


if (
  any(
    meoh_control_count_DT$n_control_rows != 1
  )
) {
  
  stop(
    "Each batch × threshold must contain exactly one MeOH control row."
  )
  
}


meoh_control_QC_DT <- meoh_threshold_DT %>%
  dplyr::transmute(
    FileName,
    TestDate,
    AssayBatch,
    BatchID,
    batch_id_clean,
    replicate_type_source,
    threshold,
    threshold_name,
    
    MeOH_T_mean =
      T_mean,
    
    MeOH_T_sd =
      T_sd,
    
    stored_MeOH_T_mean =
      .data$MeOH_T_mean,
    
    dT,
    
    meoh_missing =
      is.na(
        T_mean
      ),
    
    meoh_above_12h =
      !is.na(
        T_mean
      ) &
      T_mean >
      DT_ANALYSIS_WINDOW_H,
    
    stored_reference_matches =
      dplyr::near(
        T_mean,
        .data$MeOH_T_mean
      ),
    
    dt_is_zero =
      dplyr::near(
        dT,
        0
      )
  )


meoh_control_QC_DT


if (
  any(
    meoh_control_QC_DT$meoh_missing
  )
) {
  
  stop(
    "At least one MeOH threshold time is missing."
  )
  
}


if (
  any(
    meoh_control_QC_DT$meoh_above_12h
  )
) {
  
  stop(
    "At least one MeOH control crosses after 12 h."
  )
  
}


if (
  any(
    !meoh_control_QC_DT$stored_reference_matches
  )
) {
  
  stop(
    "Stored MeOH_T_mean does not match the explicit MeOH control row."
  )
  
}


if (
  any(
    !meoh_control_QC_DT$dt_is_zero
  )
) {
  
  stop(
    "MeOH dT is not zero for at least one control row."
  )
  
}


# ============================================================
# 10. Check that every sample row carries one consistent
#     same-plate MeOH reference
# ============================================================

sample_threshold_DT <- threshold_all_DT %>%
  dplyr::filter(
    Group == "Sample"
  )


meoh_reference_consistency_DT <-
  sample_threshold_DT %>%
  dplyr::group_by(
    batch_id_clean,
    threshold
  ) %>%
  dplyr::summarise(
    
    n_reference_values =
      dplyr::n_distinct(
        MeOH_T_mean
      ),
    
    min_reference =
      min(
        MeOH_T_mean,
        na.rm = TRUE
      ),
    
    max_reference =
      max(
        MeOH_T_mean,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


meoh_reference_consistency_DT


if (
  any(
    meoh_reference_consistency_DT$n_reference_values != 1
  )
) {
  
  stop(
    "Inconsistent MeOH_T_mean values found within a batch × threshold."
  )
  
}


# ============================================================
# 11. Create individual-level 12-h-capped threshold phenotype
#
# IMPORTANT:
# This transformation is applied to the Individual-level T_mean
# stored in the archived threshold-summary tables.
#
# It is NOT a replicate-level recalculation.
#
# Rules:
# - T_mean <= 12 h: retain observed T_mean
# - T_mean > 12 h: cap at 12 h
# - T_mean missing because threshold was not reached:
#   cap at 12 h for the 0–12 h phenotype score
#
# This preserves strong-delay individuals instead of deleting
# the strongest end of the phenotype distribution.
# ============================================================

sample_threshold_12h_DT <- sample_threshold_DT %>%
  dplyr::mutate(
    
    threshold_12h_status =
      dplyr::case_when(
        
        is.na(
          T_mean
        ) ~
          "No crossing; capped at 12 h",
        
        T_mean >
          DT_ANALYSIS_WINDOW_H ~
          "Crossing after 12 h; capped at 12 h",
        
        TRUE ~
          "Observed within 12 h"
      ),
    
    
    T_mean_12h =
      dplyr::case_when(
        
        is.na(
          T_mean
        ) ~
          DT_ANALYSIS_WINDOW_H,
        
        T_mean >
          DT_ANALYSIS_WINDOW_H ~
          DT_ANALYSIS_WINDOW_H,
        
        TRUE ~
          T_mean
      ),
    
    
    dT_12h =
      T_mean_12h -
      MeOH_T_mean
  )


# ============================================================
# 12. Summarise 12-h transformation before final matching
# ============================================================

threshold_12h_status_summary_all_DT <-
  sample_threshold_12h_DT %>%
  dplyr::count(
    threshold,
    threshold_name,
    threshold_code,
    threshold_12h_status,
    name = "n"
  ) %>%
  dplyr::arrange(
    threshold,
    threshold_12h_status
  )

threshold_12h_status_summary_all_DT


# ============================================================
# 13. Create one wide threshold record per
#     SampleID × BatchID
# ============================================================

threshold_12h_wide_DT <-
  sample_threshold_12h_DT %>%
  
  dplyr::select(
    sample_id_clean,
    batch_id_clean,
    replicate_type_source,
    threshold_code,
    
    MeOH_T_mean,
    T_mean_12h,
    dT_12h,
    threshold_12h_status
  ) %>%
  
  tidyr::pivot_wider(
    names_from =
      threshold_code,
    
    values_from = c(
      MeOH_T_mean,
      T_mean_12h,
      dT_12h,
      threshold_12h_status
    ),
    
    names_glue =
      "{.value}_{threshold_code}"
  )


# Rename the downstream phenotype columns clearly
threshold_12h_wide_DT <-
  threshold_12h_wide_DT %>%
  dplyr::rename(
    
    meoh_t0_05 =
      MeOH_T_mean_005,
    
    meoh_t0_10 =
      MeOH_T_mean_010,
    
    meoh_t0_15 =
      MeOH_T_mean_015,
    
    
    t0_05_t_mean_12h =
      T_mean_12h_005,
    
    t0_10_t_mean_12h =
      T_mean_12h_010,
    
    t0_15_t_mean_12h =
      T_mean_12h_015,
    
    
    dt0_05_12h =
      dT_12h_005,
    
    dt0_10_12h =
      dT_12h_010,
    
    dt0_15_12h =
      dT_12h_015,
    
    
    dt0_05_12h_status =
      threshold_12h_status_005,
    
    dt0_10_12h_status =
      threshold_12h_status_010,
    
    dt0_15_12h_status =
      threshold_12h_status_015
  )


# ============================================================
# 14. Read the final 589-individual phenotype as the anchor
# ============================================================

pheno_anchor_DT <- read_required_csv_DT(
  DT_pheno_file,
  "final phenotype anchor"
)


required_anchor_cols_DT <- c(
  "sample_id_clean",
  "batch_id",
  "replicate_type"
)


missing_anchor_cols_DT <- setdiff(
  required_anchor_cols_DT,
  names(
    pheno_anchor_DT
  )
)


if (
  length(
    missing_anchor_cols_DT
  ) > 0
) {
  
  stop(
    "Final phenotype anchor is missing columns: ",
    paste(
      missing_anchor_cols_DT,
      collapse = ", "
    )
  )
  
}


pheno_anchor_DT <- pheno_anchor_DT %>%
  dplyr::mutate(
    
    sample_id_clean =
      clean_sample_id_DT(
        sample_id_clean
      ),
    
    batch_id_clean =
      normalise_batch_id_DT(
        batch_id
      )
  )


# ============================================================
# 15. Check uniqueness of the final phenotype anchor
# ============================================================

anchor_duplicate_check_DT <-
  pheno_anchor_DT %>%
  dplyr::count(
    sample_id_clean,
    batch_id_clean,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )


anchor_duplicate_check_DT


if (
  nrow(
    anchor_duplicate_check_DT
  ) > 0
) {
  
  stop(
    "Duplicated SampleID + BatchID found in final phenotype anchor."
  )
  
}


# ============================================================
# 16. Match archived threshold records to the 589-individual
#     final phenotype
# ============================================================

anchor_keys_DT <- pheno_anchor_DT %>%
  dplyr::distinct(
    sample_id_clean,
    batch_id_clean
  )


threshold_keys_DT <- threshold_12h_wide_DT %>%
  dplyr::distinct(
    sample_id_clean,
    batch_id_clean
  )


anchor_without_threshold_DT <-
  anchor_keys_DT %>%
  dplyr::anti_join(
    threshold_keys_DT,
    by = c(
      "sample_id_clean",
      "batch_id_clean"
    )
  )


threshold_not_in_anchor_DT <-
  threshold_keys_DT %>%
  dplyr::anti_join(
    anchor_keys_DT,
    by = c(
      "sample_id_clean",
      "batch_id_clean"
    )
  )


anchor_without_threshold_DT
threshold_not_in_anchor_DT


if (
  nrow(
    anchor_without_threshold_DT
  ) > 0
) {
  
  stop(
    "One or more final phenotype individuals could not be matched ",
    "to the archived threshold-summary tables. ",
    "Inspect anchor_without_threshold_DT."
  )
  
}


# ============================================================
# 17. Join 12-h phenotype to the final analysis anchor
# ============================================================

DT_phenotype_12h_clean <-
  pheno_anchor_DT %>%
  dplyr::left_join(
    threshold_12h_wide_DT,
    by = c(
      "sample_id_clean",
      "batch_id_clean"
    )
  )


dim(
  DT_phenotype_12h_clean
)


if (
  nrow(
    DT_phenotype_12h_clean
  ) !=
  nrow(
    pheno_anchor_DT
  )
) {
  
  stop(
    "Row count changed after joining the 12-h phenotype."
  )
  
}


# ============================================================
# 18. Check 2rep / 3rep source consistency
# ============================================================

replicate_type_mismatch_DT <-
  DT_phenotype_12h_clean %>%
  dplyr::filter(
    !is.na(
      replicate_type
    ),
    !is.na(
      replicate_type_source
    ),
    replicate_type !=
      replicate_type_source
  ) %>%
  dplyr::select(
    sample_id_clean,
    batch_id,
    replicate_type,
    replicate_type_source
  )


replicate_type_mismatch_DT


if (
  nrow(
    replicate_type_mismatch_DT
  ) > 0
) {
  
  stop(
    "Replicate-type mismatch found after threshold-table matching."
  )
  
}


# ============================================================
# 19. Final 12-h phenotype completeness
# ============================================================

DT_12h_final_completeness <- tibble::tibble(
  
  metric = c(
    "Delta-T0.05",
    "Delta-T0.10",
    "Delta-T0.15"
  ),
  
  n_total =
    nrow(
      DT_phenotype_12h_clean
    ),
  
  n_nonmissing = c(
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_05_12h
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_10_12h
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_15_12h
      )
    )
  )
) %>%
  dplyr::mutate(
    
    n_missing =
      n_total -
      n_nonmissing
  )


DT_12h_final_completeness


# ============================================================
# 20. Final status counts after restricting to the 589
#     individuals
# ============================================================

DT_12h_final_status_summary <-
  dplyr::bind_rows(
    
    DT_phenotype_12h_clean %>%
      dplyr::count(
        status =
          dt0_05_12h_status,
        name = "n"
      ) %>%
      dplyr::mutate(
        metric =
          "Delta-T0.05"
      ),
    
    DT_phenotype_12h_clean %>%
      dplyr::count(
        status =
          dt0_10_12h_status,
        name = "n"
      ) %>%
      dplyr::mutate(
        metric =
          "Delta-T0.10"
      ),
    
    DT_phenotype_12h_clean %>%
      dplyr::count(
        status =
          dt0_15_12h_status,
        name = "n"
      ) %>%
      dplyr::mutate(
        metric =
          "Delta-T0.15"
      )
  ) %>%
  dplyr::select(
    metric,
    status,
    n
  ) %>%
  dplyr::arrange(
    metric,
    status
  )


DT_12h_final_status_summary


# ============================================================
# 21. Compare unchanged <=12-h values with the original
#     phenotype table
#
# This is QC only.
# The original columns remain untouched.
# ============================================================

DT_original_vs_12h_concordance <- tibble::tibble(
  
  metric = c(
    "Delta-T0.05",
    "Delta-T0.10",
    "Delta-T0.15"
  ),
  
  n_original_nonmissing = c(
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$t0_05_d_t
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$t0_1_d_t
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$t0_15_d_t
      )
    )
  ),
  
  n_12h_nonmissing = c(
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_05_12h
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_10_12h
      )
    ),
    
    sum(
      !is.na(
        DT_phenotype_12h_clean$dt0_15_12h
      )
    )
  )
)


DT_original_vs_12h_concordance


# ============================================================
# 22. Save clean analysis phenotype
# ============================================================

readr::write_csv(
  DT_phenotype_12h_clean,
  DT_pheno_12h_clean_file
)


# ============================================================
# 23. Save QC / provenance tables
# ============================================================

readr::write_csv(
  threshold_definition_check_DT,
  file.path(DT_pre_QC_dir, "DT_threshold_definition_check.csv")
)

readr::write_csv(
  threshold_combined_summary_DT,
  file.path(DT_pre_QC_dir, "DT_threshold_combined_summary.csv")
)

readr::write_csv(
  threshold_source_summary_DT,
  file.path(DT_pre_QC_dir, "DT_threshold_source_summary.csv")
)

readr::write_csv(
  meoh_control_QC_DT,
  file.path(DT_pre_QC_dir, "DT_MeOH_threshold_control_QC.csv")
)

readr::write_csv(
  meoh_reference_consistency_DT,
  file.path(DT_pre_QC_dir, "DT_MeOH_reference_consistency.csv")
)

readr::write_csv(
  threshold_12h_status_summary_all_DT,
  file.path(DT_pre_QC_dir, "DT_12h_status_summary_all_archived_records.csv")
)

readr::write_csv(
  anchor_without_threshold_DT,
  file.path(DT_pre_QC_dir, "DT_final_anchor_without_threshold_record.csv")
)

readr::write_csv(
  threshold_not_in_anchor_DT,
  file.path(DT_pre_QC_dir, "DT_threshold_records_not_in_final_anchor.csv")
)

readr::write_csv(
  replicate_type_mismatch_DT,
  file.path(DT_pre_QC_dir, "DT_replicate_type_mismatch.csv")
)

readr::write_csv(
  DT_12h_final_completeness,
  file.path(DT_pre_QC_dir, "DT_12h_final_completeness.csv")
)

readr::write_csv(
  DT_12h_final_status_summary,
  file.path(DT_pre_QC_dir, "DT_12h_final_status_summary.csv")
)

readr::write_csv(
  DT_original_vs_12h_concordance,
  file.path(DT_pre_QC_dir, "DT_original_vs_12h_completeness_comparison.csv")
)


# ============================================================
# 24. Final messages
# ============================================================

message(
  "01_PRE completed."
)

message(
  "Clean 12-h Delta-T phenotype saved to:"
)

message(
  DT_pheno_12h_clean_file
)

message(
  "QC outputs saved to:"
)

message(
  DT_pre_QC_dir
)