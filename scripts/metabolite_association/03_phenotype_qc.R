# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 03_phenotype_qc.R
# Phenotype QC for matched samples
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))
source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))
# =========================
# 1. Read phenotype table
# =========================

pheno <- read_csv(
  file.path(output_dir, "04_preprocessed_data", "phenotype_for_analysis.csv"),
  show_col_types = FALSE
)

dim(pheno)
names(pheno)

# =========================
# 2. Check primary phenotype
# =========================

primary_pheno <- "auc_ratio_mean"

pheno %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(.data[[primary_pheno]])),
    min_value = min(.data[[primary_pheno]], na.rm = TRUE),
    q1 = quantile(.data[[primary_pheno]], 0.25, na.rm = TRUE),
    median = median(.data[[primary_pheno]], na.rm = TRUE),
    mean = mean(.data[[primary_pheno]], na.rm = TRUE),
    q3 = quantile(.data[[primary_pheno]], 0.75, na.rm = TRUE),
    max_value = max(.data[[primary_pheno]], na.rm = TRUE),
    sd = sd(.data[[primary_pheno]], na.rm = TRUE)
  )

# =========================
# 3. Summary of phenotype variables
# =========================

pheno_vars <- c(
  "auc_mean",
  "auc_ratio_mean",
  "t0_05_t_mean",
  "t0_1_t_mean",
  "t0_15_t_mean",
  "t0_05_d_t",
  "t0_1_d_t",
  "t0_15_d_t"
)

pheno_summary <- pheno %>%
  dplyr::select(any_of(pheno_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "phenotype",
    values_to = "value"
  ) %>%
  group_by(phenotype) %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(value)),
    min_value = min(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

pheno_summary

# =========================
# 4. AUC-ratio distribution
# =========================

p_auc_hist <- ggplot(pheno, aes(x = auc_ratio_mean)) +
  geom_histogram(bins = 40, color = "black", fill = "grey80") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_classic(base_size = 12) +
  labs(
    title = "Distribution of S. aureus AUC-ratio in matched LC-MS samples",
    x = "AUC-ratio mean",
    y = "Number of samples"
  )

p_auc_hist

# =========================
# 5. Save phenotype QC outputs
# =========================

dir.create(file.path(output_dir, "05_phenotype_qc"), showWarnings = FALSE)

write_csv(
  pheno_summary,
  file.path(output_dir, "05_phenotype_qc", "phenotype_summary_matched_samples.csv")
)

ggsave(
  filename = file.path(output_dir, "05_phenotype_qc", "auc_ratio_distribution_histogram.png"),
  plot = p_auc_hist,
  width = 6,
  height = 4,
  dpi = 300
)

# =========================
# 6. Add phenotype classification
# =========================

pheno_qc <- pheno %>%
  mutate(
    auc_class_0_7_1_5 = case_when(
      auc_ratio_mean < 0.7 ~ "Inhibitor",
      auc_ratio_mean > 1.5 ~ "High AUC / weak inhibition",
      TRUE ~ "Intermediate"
    )
  )

pheno_qc %>%
  count(auc_class_0_7_1_5)

# =========================
# 7. Add quartile-based DAM group
# =========================

auc_q1 <- quantile(pheno_qc$auc_ratio_mean, 0.25, na.rm = TRUE)
auc_q3 <- quantile(pheno_qc$auc_ratio_mean, 0.75, na.rm = TRUE)

auc_q1
auc_q3

pheno_qc <- pheno_qc %>%
  mutate(
    auc_dam_group_q25_q75 = case_when(
      auc_ratio_mean <= auc_q1 ~ "Low_AUC_strong_inhibition",
      auc_ratio_mean >= auc_q3 ~ "High_AUC_weak_inhibition",
      TRUE ~ "Middle"
    )
  )

pheno_qc %>%
  count(auc_dam_group_q25_q75)

# =========================
# 8. Group summaries
# =========================

auc_class_summary <- pheno_qc %>%
  count(auc_class_0_7_1_5) %>%
  arrange(desc(n))

auc_dam_group_summary <- pheno_qc %>%
  count(auc_dam_group_q25_q75) %>%
  arrange(desc(n))

auc_class_summary
auc_dam_group_summary

# =========================
# 9. AUC-ratio distribution with thresholds
# =========================

p_auc_hist_thresholds <- ggplot(pheno_qc, aes(x = auc_ratio_mean)) +
  geom_histogram(
    bins = 40,
    color = "black",
    fill = alpha("#6FA8DC", 0.55),
    linewidth = 0.4
  ) +
  geom_vline(xintercept = 0.7, linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 0.8) +
  theme_prism_user() +
  labs(
    title = "Distribution of S. aureus AUC-ratio in matched LC-MS samples",
    x = "AUC-ratio mean",
    y = "Number of samples"
  )

p_auc_hist_thresholds

# =========================
# 10. Save phenotype QC outputs
# =========================

dir.create(file.path(output_dir, "05_phenotype_qc"), showWarnings = FALSE)

write_csv(
  pheno_summary,
  file.path(output_dir, "05_phenotype_qc", "phenotype_summary_matched_samples.csv")
)

write_csv(
  pheno_qc,
  file.path(output_dir, "05_phenotype_qc", "phenotype_qc_with_auc_groups.csv")
)

write_csv(
  auc_class_summary,
  file.path(output_dir, "05_phenotype_qc", "auc_class_summary_0_7_1_5.csv")
)

write_csv(
  auc_dam_group_summary,
  file.path(output_dir, "05_phenotype_qc", "auc_dam_group_summary_q25_q75.csv")
)

ggsave(
  filename = file.path(output_dir, "05_phenotype_qc", "auc_ratio_distribution_histogram.png"),
  plot = p_auc_hist,
  width = 6,
  height = 4,
  dpi = 300
)

save_prism_plot(
  p_auc_hist_thresholds,
  file.path(output_dir, "05_phenotype_qc", "auc_ratio_distribution_prism_thresholds.png"),
  width = 6,
  height = 4,
  dpi = 600
)