# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 02_metabolite_preprocessing.R
# Metabolite QC and preprocessing
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))

# =========================
# 1. Read matched data
# =========================

matched_data <- read_csv(
  file.path(output_dir, "02_matched_data", "matched_pheno_metab_data.csv"),
  show_col_types = FALSE
)

anno_clean <- read_csv(
  file.path(output_dir, "02_matched_data", "peak_annotation_clean.csv"),
  show_col_types = FALSE
)

dim(matched_data)
dim(anno_clean)
names(matched_data)

# =========================
# 2. Identify metabolite peak columns
# =========================

peak_cols <- names(matched_data) %>%
  keep(~ str_detect(.x, "^peak\\d+$"))

length(peak_cols)

head(peak_cols)
tail(peak_cols)

anno_peak_ids <- anno_clean %>%
  pull(peak_id_clean)

tibble(
  item = c(
    "Peak columns in matched data",
    "Peak IDs in annotation",
    "Matched peaks"
  ),
  n = c(
    length(peak_cols),
    length(anno_peak_ids),
    length(intersect(peak_cols, anno_peak_ids))
  )
)

# =========================
# 3. Peak missing / zero QC
# =========================

peak_qc_summary <- map_dfr(peak_cols, function(peak) {
  
  x <- matched_data[[peak]]
  
  tibble(
    peak_id = peak,
    n_total = length(x),
    n_missing = sum(is.na(x)),
    missing_rate = mean(is.na(x)),
    n_zero = sum(x == 0, na.rm = TRUE),
    zero_rate = mean(x == 0, na.rm = TRUE),
    min_value = min(x, na.rm = TRUE),
    max_value = max(x, na.rm = TRUE),
    sd_value = sd(x, na.rm = TRUE)
  )
}) %>%
  arrange(desc(missing_rate), desc(zero_rate))

peak_qc_overall <- tibble(
  item = c(
    "Total peaks",
    "Peaks with any missing values",
    "Peaks with missing rate > 30%",
    "Peaks with zero rate > 30%",
    "Peaks with zero rate > 50%",
    "Peaks with zero SD"
  ),
  n = c(
    length(peak_cols),
    sum(peak_qc_summary$n_missing > 0),
    sum(peak_qc_summary$missing_rate > 0.30),
    sum(peak_qc_summary$zero_rate > 0.30),
    sum(peak_qc_summary$zero_rate > 0.50),
    sum(peak_qc_summary$sd_value == 0, na.rm = TRUE)
  )
)

peak_qc_overall

peak_qc_summary %>%
  slice(1:20)

# =========================
# 4. Select peaks for downstream analysis
# =========================

peak_cols_keep <- peak_qc_summary %>%
  filter(
    missing_rate <= 0.30,
    zero_rate <= 0.50,
    !is.na(sd_value),
    sd_value > 0
  ) %>%
  pull(peak_id)

length(peak_cols_keep)

# =========================
# 5. Log2 transform metabolite peaks
# =========================

metab_log <- matched_data %>%
  dplyr::select(sample_id_clean, all_of(peak_cols_keep)) %>%
  mutate(
    across(
      all_of(peak_cols_keep),
      ~ log2(.x + 1)
    )
  )

dim(metab_log)

# =========================
# 6. Z-score scaling
# =========================

metab_scaled <- metab_log %>%
  mutate(
    across(
      all_of(peak_cols_keep),
      ~ as.numeric(scale(.x))
    )
  )

dim(metab_scaled)

# =========================
# 7. Check scaled data
# =========================

scaled_check <- map_dfr(peak_cols_keep, function(peak) {
  
  x <- metab_scaled[[peak]]
  
  tibble(
    peak_id = peak,
    mean_after_scaling = mean(x, na.rm = TRUE),
    sd_after_scaling = sd(x, na.rm = TRUE),
    min_scaled = min(x, na.rm = TRUE),
    max_scaled = max(x, na.rm = TRUE)
  )
})

scaled_check %>%
  slice(1:10)

# =========================
# 8. Create phenotype table for analysis
# =========================

pheno_cols_keep <- c(
  "sample_id_clean",
  "pheno_sample_id_original",
  "accession",
  "assay_batch",
  "batch_id",
  "auc_mean",
  "auc_sd",
  "auc_ratio_mean",
  "auc_ratio_sd",
  "t0_05_t_mean",
  "t0_1_t_mean",
  "t0_15_t_mean",
  "t0_05_t_sd",
  "t0_1_t_sd",
  "t0_15_t_sd",
  "t0_05_d_t",
  "t0_1_d_t",
  "t0_15_d_t",
  "replicate_type",
  "accession_rank1",
  "tracker_status",
  "tracker_test_date",
  "tracker_assay_batch",
  "in_rank1",
  "tracker_tested"
)

pheno_for_analysis <- matched_data %>%
  dplyr::select(any_of(pheno_cols_keep))

dim(pheno_for_analysis)
names(pheno_for_analysis)

# =========================
# 9. Save preprocessing outputs
# =========================

dir.create(file.path(output_dir, "03_peak_qc"), showWarnings = FALSE)
dir.create(file.path(output_dir, "04_preprocessed_data"), showWarnings = FALSE)

write_csv(
  peak_qc_summary,
  file.path(output_dir, "03_peak_qc", "peak_missing_zero_qc_summary.csv")
)

write_csv(
  peak_qc_overall,
  file.path(output_dir, "03_peak_qc", "peak_qc_overall_summary.csv")
)

write_csv(
  metab_log,
  file.path(output_dir, "04_preprocessed_data", "metabolite_log2_matrix.csv")
)

write_csv(
  metab_scaled,
  file.path(output_dir, "04_preprocessed_data", "metabolite_log2_zscaled_matrix.csv")
)

write_csv(
  pheno_for_analysis,
  file.path(output_dir, "04_preprocessed_data", "phenotype_for_analysis.csv")
)

write_csv(
  tibble(peak_id = peak_cols_keep),
  file.path(output_dir, "04_preprocessed_data", "peaks_retained_for_analysis.csv")
)

write_csv(
  scaled_check,
  file.path(output_dir, "04_preprocessed_data", "scaled_data_check.csv")
)
