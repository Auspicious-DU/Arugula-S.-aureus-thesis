# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 00_setup.R
# Project setup: paths, packages, and reusable functions
# =========================

# =========================
# 0. Packages
# =========================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(janitor)
  library(forcats)
  library(scales)
})


source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))

# =========================
# Plot settings for RE analysis
# =========================

RE_group_cols <- c(
  "Low_AUC_inhibition_delay" = "#6FA8DC",
  "Intermediate" = "#B4A7D6",
  "High_AUC_growth_enhancing" = "#F4A6B7"
)

RE_group_labels <- c(
  "Low_AUC_inhibition_delay" = "Low AUC\ninhibition/delay",
  "Intermediate" = "Intermediate",
  "High_AUC_growth_enhancing" = "High AUC\ngrowth-enhancing"
)

RE_batch_levels <- paste0("bat", 1:11)
# =========================
# 1. Project folders
# =========================

project_dir <- repo_dir

input_dir <- file.path(repo_dir, "input", "metabolite_association")
output_dir <- file.path(repo_dir, "output", "metabolite_association")
script_dir <- file.path(repo_dir, "scripts", "metabolite_association")

# =========================
# 2. Input file paths
# =========================

metab_file <- file.path(
  input_dir,
  "01_metabolite_matrix_659Arugula_275CorrectedPeak.xlsx"
)

anno_file <- file.path(
  input_dir,
  "02_peak_annotation_275Metabolite_Arugula.xlsx"
)

pheno_file <- file.path(
  input_dir,
  "03_Rank1_phenotype_summary_final_wide_clean.xlsx"
)

# =========================
# 3. Reusable functions
# =========================

# Clean sample IDs into a shared format:
# A10:17:2:1 / A10/17/2/1 / A10.17.2.1 -> a10.17.2.1
clean_sample_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("：", ":") %>%
    stringr::str_replace_all("[:/\\\\]", ".") %>%
    stringr::str_replace_all("\\s+", "") %>%
    stringr::str_replace_all("\\.+", ".") %>%
    stringr::str_replace_all("^\\.|\\.$", "")
}

# Identify metabolite peak columns
get_peak_cols <- function(df) {
  names(df) %>%
    stringr::str_subset("^peak\\d+$")
}

# Create directory only if it does not exist
make_output_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Safe -log10 transform for p-values
safe_neg_log10 <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
}

# =========================
# 4. Common output folders
# =========================

output_01_sample_matching <- file.path(output_dir, "01_sample_matching")
output_02_matched_data <- file.path(output_dir, "02_matched_data")
output_03_peak_qc <- file.path(output_dir, "03_peak_qc")
output_04_preprocessed_data <- file.path(output_dir, "04_preprocessed_data")
output_05_phenotype_qc <- file.path(output_dir, "05_phenotype_qc")
output_05B_metabolite_pca_qc <- file.path(output_dir, "05B_metabolite_pca_qc")
output_06_correlation_DAM <- file.path(output_dir, "06_correlation_DAM")
output_06B_batch_adjusted_sensitivity <- file.path(output_dir, "06B_batch_adjusted_sensitivity")
output_07_PLS_VIP_candidate <- file.path(output_dir, "07_PLS_VIP_candidate")
output_08_final_candidate_tables <- file.path(output_dir, "08_final_candidate_tables")

# =========================
# 5. Biological-threshold re-analysis folders
# =========================

reanalysis_dir <- file.path(
  output_dir,
  "biological_threshold_reanalysis"
)

re_output_03_peak_qc <- file.path(
  reanalysis_dir,
  "03_peak_qc"
)

re_output_04_preprocessed_data <- file.path(
  reanalysis_dir,
  "04_preprocessed_data"
)

re_output_05_group_qc <- file.path(
  reanalysis_dir,
  "05_biological_group_qc"
)

re_output_06_inhibition_analysis <- file.path(
  reanalysis_dir,
  "06_inhibition_delay_analysis"
)

re_output_07_promotion_analysis <- file.path(
  reanalysis_dir,
  "07_growth_enhancing_analysis"
)

make_output_dir(reanalysis_dir)