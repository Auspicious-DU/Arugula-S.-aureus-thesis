# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 04B_batch_adjusted_sensitivity.R
# Batch-adjusted sensitivity analysis
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))
source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))

# =========================
# 1. Read required data
# =========================

pheno_qc <- read_csv(
  file.path(output_dir, "05_phenotype_qc", "phenotype_qc_with_auc_groups.csv"),
  show_col_types = FALSE
)

metab_scaled <- read_csv(
  file.path(output_dir, "04_preprocessed_data", "metabolite_log2_zscaled_matrix.csv"),
  show_col_types = FALSE
)

metab_log <- read_csv(
  file.path(output_dir, "04_preprocessed_data", "metabolite_log2_matrix.csv"),
  show_col_types = FALSE
)

matched_data <- read_csv(
  file.path(output_dir, "02_matched_data", "matched_pheno_metab_data.csv"),
  show_col_types = FALSE
)

integrated_results_anno <- read_csv(
  file.path(output_dir, "06_correlation_DAM", "integrated_pearson_DAM_results_with_annotation.csv"),
  show_col_types = FALSE
)

dim(pheno_qc)
dim(metab_scaled)
dim(metab_log)
dim(matched_data)
dim(integrated_results_anno)

# =========================
# 2. Extract LC-MS batch from meta_id
# =========================

sample_batch <- matched_data %>%
  dplyr::select(sample_id_clean, meta_id) %>%
  mutate(
    lcms_batch = str_extract(str_to_lower(meta_id), "^bat\\d+")
  ) %>%
  distinct(sample_id_clean, .keep_all = TRUE)

sample_batch %>%
  count(lcms_batch, sort = TRUE)

sample_batch %>%
  summarise(
    n_total = n(),
    n_missing_lcms_batch = sum(is.na(lcms_batch))
  )

# =========================
# 3. Identify peak columns
# =========================

peak_cols <- names(metab_scaled) %>%
  keep(~ str_detect(.x, "^peak\\d+$"))

length(peak_cols)
head(peak_cols)
tail(peak_cols)

# =========================
# 4. Prepare data for batch-adjusted Pearson-like model
# =========================

primary_pheno <- "auc_ratio_mean"

pearson_adj_data <- pheno_qc %>%
  dplyr::select(sample_id_clean, auc_ratio_mean) %>%
  inner_join(sample_batch, by = "sample_id_clean") %>%
  inner_join(metab_scaled, by = "sample_id_clean") %>%
  mutate(
    lcms_batch = factor(lcms_batch)
  )

dim(pearson_adj_data)

pearson_adj_data %>%
  summarise(
    n_total = n(),
    n_missing_auc = sum(is.na(auc_ratio_mean)),
    n_missing_lcms_batch = sum(is.na(lcms_batch))
  )

# =========================
# 5. Batch-adjusted Pearson-like model for all peaks
# auc_ratio_mean ~ metabolite + lcms_batch
# =========================

batch_adj_pearson_results <- map_dfr(peak_cols, function(peak) {
  
  df <- pearson_adj_data %>%
    dplyr::select(auc_ratio_mean, lcms_batch, all_of(peak)) %>%
    rename(metab_value = all_of(peak)) %>%
    drop_na()
  
  if (
    nrow(df) < 10 ||
    sd(df$metab_value, na.rm = TRUE) == 0 ||
    sd(df$auc_ratio_mean, na.rm = TRUE) == 0
  ) {
    return(
      tibble(
        peak_id = peak,
        n = nrow(df),
        adj_pearson_beta = NA_real_,
        adj_pearson_p = NA_real_,
        model_status = "Skipped"
      )
    )
  }
  
  fit_result <- tryCatch({
    
    fit <- lm(
      auc_ratio_mean ~ metab_value + lcms_batch,
      data = df
    )
    
    coef_df <- summary(fit)$coefficients
    
    tibble(
      peak_id = peak,
      n = nrow(df),
      adj_pearson_beta = coef_df["metab_value", "Estimate"],
      adj_pearson_p = coef_df["metab_value", "Pr(>|t|)"],
      model_status = "OK"
    )
    
  }, error = function(e) {
    
    tibble(
      peak_id = peak,
      n = nrow(df),
      adj_pearson_beta = NA_real_,
      adj_pearson_p = NA_real_,
      model_status = paste0("Error: ", e$message)
    )
    
  })
  
  fit_result
}) %>%
  mutate(
    adj_pearson_fdr = p.adjust(adj_pearson_p, method = "BH"),
    adj_pearson_direction = case_when(
      adj_pearson_beta < 0 ~ "Negative with AUC-ratio after LC-MS batch adjustment",
      adj_pearson_beta > 0 ~ "Positive with AUC-ratio after LC-MS batch adjustment",
      TRUE ~ "No direction"
    )
  ) %>%
  arrange(adj_pearson_p)

batch_adj_pearson_results %>%
  slice(1:20)

batch_adj_pearson_results %>%
  summarise(
    total_peaks = n(),
    model_ok = sum(model_status == "OK"),
    adj_pearson_p_lt_0_05 = sum(adj_pearson_p < 0.05, na.rm = TRUE),
    adj_pearson_fdr_lt_0_05 = sum(adj_pearson_fdr < 0.05, na.rm = TRUE),
    adj_negative_fdr_lt_0_05 = sum(adj_pearson_beta < 0 & adj_pearson_fdr < 0.05, na.rm = TRUE),
    adj_positive_fdr_lt_0_05 = sum(adj_pearson_beta > 0 & adj_pearson_fdr < 0.05, na.rm = TRUE)
  )

# =========================
# 6. Prepare data for batch-adjusted DAM-like model
# =========================

dam_adj_data <- pheno_qc %>%
  dplyr::select(sample_id_clean, auc_ratio_mean, auc_dam_group_q25_q75) %>%
  inner_join(sample_batch, by = "sample_id_clean") %>%
  inner_join(metab_log, by = "sample_id_clean") %>%
  filter(
    auc_dam_group_q25_q75 %in% c(
      "Low_AUC_strong_inhibition",
      "High_AUC_weak_inhibition"
    )
  ) %>%
  mutate(
    lcms_batch = factor(lcms_batch),
    auc_group = factor(
      auc_dam_group_q25_q75,
      levels = c(
        "High_AUC_weak_inhibition",
        "Low_AUC_strong_inhibition"
      )
    )
  )

dim(dam_adj_data)

dam_adj_data %>%
  count(auc_group)

dam_adj_data %>%
  count(lcms_batch, auc_group)

# =========================
# 7. Batch-adjusted DAM-like model for all peaks
# metabolite ~ auc_group + lcms_batch
# =========================

batch_adj_dam_results <- map_dfr(peak_cols, function(peak) {
  
  df <- dam_adj_data %>%
    dplyr::select(auc_group, lcms_batch, all_of(peak)) %>%
    rename(metab_value = all_of(peak)) %>%
    drop_na()
  
  if (
    nrow(df) < 10 ||
    n_distinct(df$auc_group) < 2 ||
    sd(df$metab_value, na.rm = TRUE) == 0
  ) {
    return(
      tibble(
        peak_id = peak,
        n = nrow(df),
        adj_log2FC_low_vs_high = NA_real_,
        adj_dam_p = NA_real_,
        model_status = "Skipped"
      )
    )
  }
  
  fit_result <- tryCatch({
    
    fit <- lm(
      metab_value ~ auc_group + lcms_batch,
      data = df
    )
    
    coef_df <- summary(fit)$coefficients
    
    term_name <- "auc_groupLow_AUC_strong_inhibition"
    
    if (!term_name %in% rownames(coef_df)) {
      return(
        tibble(
          peak_id = peak,
          n = nrow(df),
          adj_log2FC_low_vs_high = NA_real_,
          adj_dam_p = NA_real_,
          model_status = "Term not found"
        )
      )
    }
    
    tibble(
      peak_id = peak,
      n = nrow(df),
      adj_log2FC_low_vs_high = coef_df[term_name, "Estimate"],
      adj_dam_p = coef_df[term_name, "Pr(>|t|)"],
      model_status = "OK"
    )
    
  }, error = function(e) {
    
    tibble(
      peak_id = peak,
      n = nrow(df),
      adj_log2FC_low_vs_high = NA_real_,
      adj_dam_p = NA_real_,
      model_status = paste0("Error: ", e$message)
    )
    
  })
  
  fit_result
}) %>%
  mutate(
    adj_dam_fdr = p.adjust(adj_dam_p, method = "BH"),
    adj_dam_direction = case_when(
      adj_log2FC_low_vs_high > 0 ~ "Higher in Low AUC after LC-MS batch adjustment",
      adj_log2FC_low_vs_high < 0 ~ "Higher in High AUC after LC-MS batch adjustment",
      TRUE ~ "No direction"
    ),
    adj_dam_candidate = if_else(
      abs(adj_log2FC_low_vs_high) > 0.58 & adj_dam_p < 0.05,
      "Adjusted DAM candidate",
      "Other"
    )
  ) %>%
  arrange(adj_dam_p)

batch_adj_dam_results %>%
  slice(1:20)

batch_adj_dam_results %>%
  summarise(
    total_peaks = n(),
    model_ok = sum(model_status == "OK"),
    adj_dam_p_lt_0_05 = sum(adj_dam_p < 0.05, na.rm = TRUE),
    adj_dam_fdr_lt_0_05 = sum(adj_dam_fdr < 0.05, na.rm = TRUE),
    adj_log2fc_abs_gt_0_58 = sum(abs(adj_log2FC_low_vs_high) > 0.58, na.rm = TRUE),
    adj_dam_candidate_nominal = sum(
      abs(adj_log2FC_low_vs_high) > 0.58 & adj_dam_p < 0.05,
      na.rm = TRUE
    ),
    adj_higher_in_low_fdr_lt_0_05 = sum(
      adj_log2FC_low_vs_high > 0 & adj_dam_fdr < 0.05,
      na.rm = TRUE
    ),
    adj_higher_in_high_fdr_lt_0_05 = sum(
      adj_log2FC_low_vs_high < 0 & adj_dam_fdr < 0.05,
      na.rm = TRUE
    )
  )

# =========================
# 8. Integrate adjusted results with original integrated table
# =========================

batch_adjusted_integrated <- integrated_results_anno %>%
  left_join(
    batch_adj_pearson_results %>%
      dplyr::select(
        peak_id,
        adj_pearson_beta,
        adj_pearson_p,
        adj_pearson_fdr,
        adj_pearson_direction,
        adj_pearson_model_status = model_status
      ),
    by = "peak_id"
  ) %>%
  left_join(
    batch_adj_dam_results %>%
      dplyr::select(
        peak_id,
        adj_log2FC_low_vs_high,
        adj_dam_p,
        adj_dam_fdr,
        adj_dam_direction,
        adj_dam_candidate,
        adj_dam_model_status = model_status
      ),
    by = "peak_id"
  ) %>%
  mutate(
    batch_adjusted_support = case_when(
      
      priority == "High priority inhibition-associated" &
        adj_pearson_beta < 0 &
        adj_pearson_fdr < 0.05 &
        adj_log2FC_low_vs_high > 0 &
        adj_dam_fdr < 0.05 &
        abs(adj_log2FC_low_vs_high) > 0.58 ~
        "Robust inhibition-associated after LC-MS batch adjustment",
      
      priority == "High priority inhibition-associated" &
        adj_pearson_beta < 0 &
        adj_pearson_p < 0.05 &
        adj_log2FC_low_vs_high > 0 &
        adj_dam_p < 0.05 &
        abs(adj_log2FC_low_vs_high) > 0.58 ~
        "Nominally supported inhibition-associated after LC-MS batch adjustment",
      
      priority == "High priority high-AUC-associated" &
        adj_pearson_beta > 0 &
        adj_pearson_fdr < 0.05 &
        adj_log2FC_low_vs_high < 0 &
        adj_dam_fdr < 0.05 &
        abs(adj_log2FC_low_vs_high) > 0.58 ~
        "Robust high-AUC-associated after LC-MS batch adjustment",
      
      priority == "High priority high-AUC-associated" &
        adj_pearson_beta > 0 &
        adj_pearson_p < 0.05 &
        adj_log2FC_low_vs_high < 0 &
        adj_dam_p < 0.05 &
        abs(adj_log2FC_low_vs_high) > 0.58 ~
        "Nominally supported high-AUC-associated after LC-MS batch adjustment",
      
      priority %in% c(
        "High priority inhibition-associated",
        "High priority high-AUC-associated"
      ) ~
        "Original high-priority signal weakened after LC-MS batch adjustment",
      
      TRUE ~ "Not original high-priority candidate"
    )
  ) %>%
  arrange(
    factor(
      batch_adjusted_support,
      levels = c(
        "Robust inhibition-associated after LC-MS batch adjustment",
        "Nominally supported inhibition-associated after LC-MS batch adjustment",
        "Robust high-AUC-associated after LC-MS batch adjustment",
        "Nominally supported high-AUC-associated after LC-MS batch adjustment",
        "Original high-priority signal weakened after LC-MS batch adjustment",
        "Not original high-priority candidate"
      )
    ),
    pearson_p
  )

batch_adjusted_integrated %>%
  dplyr::select(
    peak_id,
    priority,
    batch_adjusted_support,
    pearson_r,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_fdr,
    adj_pearson_beta,
    adj_pearson_fdr,
    adj_log2FC_low_vs_high,
    adj_dam_fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:30)

# =========================
# 9. Batch-adjusted summary
# =========================

batch_adjusted_support_summary <- batch_adjusted_integrated %>%
  count(batch_adjusted_support) %>%
  arrange(desc(n))

batch_adjusted_support_summary

high_priority_inhibition_batch_summary <- batch_adjusted_integrated %>%
  filter(priority == "High priority inhibition-associated") %>%
  count(batch_adjusted_support) %>%
  arrange(desc(n))

high_priority_high_auc_batch_summary <- batch_adjusted_integrated %>%
  filter(priority == "High priority high-AUC-associated") %>%
  count(batch_adjusted_support) %>%
  arrange(desc(n))

high_priority_inhibition_batch_summary
high_priority_high_auc_batch_summary

# =========================
# 10. Extract robust candidate tables
# =========================

robust_inhibition_candidates <- batch_adjusted_integrated %>%
  filter(
    batch_adjusted_support ==
      "Robust inhibition-associated after LC-MS batch adjustment"
  ) %>%
  arrange(pearson_p)

nominal_supported_inhibition_candidates <- batch_adjusted_integrated %>%
  filter(
    batch_adjusted_support %in% c(
      "Robust inhibition-associated after LC-MS batch adjustment",
      "Nominally supported inhibition-associated after LC-MS batch adjustment"
    )
  ) %>%
  arrange(pearson_p)

robust_high_auc_candidates <- batch_adjusted_integrated %>%
  filter(
    batch_adjusted_support ==
      "Robust high-AUC-associated after LC-MS batch adjustment"
  ) %>%
  arrange(pearson_p)

robust_inhibition_candidates %>%
  dplyr::select(
    peak_id,
    pearson_r,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_fdr,
    adj_pearson_beta,
    adj_pearson_fdr,
    adj_log2FC_low_vs_high,
    adj_dam_fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:30)

# =========================
# 11. Plot original Pearson r vs adjusted beta
# =========================

p_original_vs_adj_pearson <- ggplot(
  batch_adjusted_integrated,
  aes(x = pearson_r, y = adj_pearson_beta)
) +
  geom_point(
    aes(color = priority),
    size = 2.1,
    alpha = 0.75
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7) +
  theme_prism_user() +
  labs(
    title = "Original vs LC-MS batch-adjusted metabolite association",
    x = "Original Pearson r",
    y = "Adjusted beta from linear model",
    color = NULL
  )

p_original_vs_adj_pearson

# =========================
# 12. Plot original DAM log2FC vs adjusted log2FC
# =========================

p_original_vs_adj_dam <- ggplot(
  batch_adjusted_integrated,
  aes(x = log2FC_low_vs_high, y = adj_log2FC_low_vs_high)
) +
  geom_point(
    aes(color = priority),
    size = 2.1,
    alpha = 0.75
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7) +
  theme_prism_user() +
  labs(
    title = "Original vs LC-MS batch-adjusted DAM effect",
    x = "Original log2FC: Low AUC vs High AUC",
    y = "Adjusted log2FC: Low AUC vs High AUC",
    color = NULL
  )

p_original_vs_adj_dam

# =========================
# 13. Save outputs
# =========================

dir.create(file.path(output_dir, "06B_batch_adjusted_sensitivity"), showWarnings = FALSE)

write_csv(
  batch_adj_pearson_results,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_pearson_like_results.csv")
)

write_csv(
  batch_adj_dam_results,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_DAM_like_results.csv")
)

write_csv(
  batch_adjusted_integrated,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_integrated_results.csv")
)

write_csv(
  batch_adjusted_support_summary,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_support_summary.csv")
)

write_csv(
  high_priority_inhibition_batch_summary,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "high_priority_inhibition_batch_summary.csv")
)

write_csv(
  high_priority_high_auc_batch_summary,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "high_priority_high_AUC_batch_summary.csv")
)

write_csv(
  robust_inhibition_candidates,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "robust_inhibition_candidates_batch_adjusted.csv")
)

write_csv(
  nominal_supported_inhibition_candidates,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "nominal_supported_inhibition_candidates_batch_adjusted.csv")
)

write_csv(
  robust_high_auc_candidates,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "robust_high_AUC_candidates_batch_adjusted.csv")
)

save_prism_plot(
  p_original_vs_adj_pearson,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "original_vs_batch_adjusted_pearson.png"),
  width = 6,
  height = 4.5,
  dpi = 600
)

save_prism_plot(
  p_original_vs_adj_dam,
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "original_vs_batch_adjusted_DAM.png"),
  width = 6,
  height = 4.5,
  dpi = 600
)