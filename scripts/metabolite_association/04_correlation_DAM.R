# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 04_correlation_DAM.R
# Correlation and DAM analysis
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))
source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))

# =========================
# 1. Read preprocessed data
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

anno_clean <- read_csv(
  file.path(output_dir, "02_matched_data", "peak_annotation_clean.csv"),
  show_col_types = FALSE
)

dim(pheno_qc)
dim(metab_scaled)
dim(metab_log)
dim(anno_clean)

# =========================
# 2. Identify peak columns
# =========================

peak_cols <- names(metab_scaled) %>%
  keep(~ str_detect(.x, "^peak\\d+$"))

length(peak_cols)
head(peak_cols)
tail(peak_cols)

# =========================
# 3. Prepare data for Pearson correlation
# =========================

primary_pheno <- "auc_ratio_mean"

cor_data <- pheno_qc %>%
  dplyr::select(
    sample_id_clean,
    auc_ratio_mean,
    auc_class_0_7_1_5,
    auc_dam_group_q25_q75
  ) %>%
  inner_join(metab_scaled, by = "sample_id_clean")

dim(cor_data)

cor_data %>%
  summarise(
    n_total = n(),
    n_missing_pheno = sum(is.na(.data[[primary_pheno]]))
  )

# =========================
# 4. Pearson correlation
# AUC-ratio mean vs each metabolite peak
# =========================

pearson_results <- map_dfr(peak_cols, function(peak) {
  
  df <- cor_data %>%
    dplyr::select(all_of(primary_pheno), all_of(peak)) %>%
    drop_na()
  
  x <- df[[peak]]
  y <- df[[primary_pheno]]
  
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(
      tibble(
        peak_id = peak,
        n = length(x),
        pearson_r = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- cor.test(x, y, method = "pearson")
  
  tibble(
    peak_id = peak,
    n = length(x),
    pearson_r = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_value),
    direction = case_when(
      pearson_r < 0 ~ "Negative with AUC-ratio",
      pearson_r > 0 ~ "Positive with AUC-ratio",
      TRUE ~ "No direction"
    ),
    nominal_significant = if_else(p_value < 0.05, "p < 0.05", "NS"),
    fdr_significant = if_else(fdr < 0.05, "FDR < 0.05", "NS")
  ) %>%
  arrange(p_value)

pearson_results <- pearson_results %>%
  mutate(
    pearson_plot_group = case_when(
      pearson_r < 0 & fdr < 0.05 ~ "Negative with AUC\n(FDR < 0.05)",
      pearson_r > 0 & fdr < 0.05 ~ "Positive with AUC\n(FDR < 0.05)",
      p_value < 0.05 ~ "p < 0.05 only",
      TRUE ~ "Other"
    ),
    pearson_plot_group = factor(
      pearson_plot_group,
      levels = c(
        "Negative with AUC\n(FDR < 0.05)",
        "Positive with AUC\n(FDR < 0.05)",
        "p < 0.05 only",
        "Other"
      )
    )
  )

pearson_results %>%
  slice(1:20)

#pearson_r < 0 = metabolite abundance higher → AUC-ratio lower → stronger inhibition / stronger delay

# =========================
# 5. Add peak annotation
# =========================

pearson_results_anno <- pearson_results %>%
  left_join(
    anno_clean,
    by = c("peak_id" = "peak_id_clean")
  )

pearson_results_anno %>%
  dplyr::select(
    peak_id,
    pearson_r,
    p_value,
    fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:20)

# =========================
# 6. Pearson correlation plot
# =========================
p_pearson_volcano <- ggplot(
  pearson_results,
  aes(x = pearson_r, y = neg_log10_p)
) +
  geom_point(
    aes(color = pearson_plot_group),
    size = 2.3,
    alpha = 0.80
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Negative with AUC\n(FDR < 0.05)" = "#6FA8DC",
      "Positive with AUC\n(FDR < 0.05)" = "#F4A6B7",
      "p < 0.05 only" = "#B4A7D6",
      "Other" = "grey75"
    )
  ) +
  theme_prism_user() +
  labs(
    title = "Metabolite correlation with S. aureus AUC-ratio",
    x = "Pearson r with AUC-ratio",
    y = "-log10(p-value)",
    color = NULL
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1))
  )

p_pearson_volcano

# =========================
# 7. Save Pearson outputs
# =========================

dir.create(file.path(output_dir, "06_correlation_DAM"), showWarnings = FALSE)

write_csv(
  pearson_results,
  file.path(output_dir, "06_correlation_DAM", "pearson_auc_ratio_metabolite_results.csv")
)

write_csv(
  pearson_results_anno,
  file.path(output_dir, "06_correlation_DAM", "pearson_auc_ratio_metabolite_results_with_annotation.csv")
)

save_prism_plot(
  p_pearson_volcano,
  file.path(output_dir, "06_correlation_DAM", "pearson_auc_ratio_volcano.png"),
  width = 6,
  height = 4,
  dpi = 600
)

# =========================
# 8. Prepare data for DAM analysis
# =========================

dam_data <- pheno_qc %>%
  dplyr::select(
    sample_id_clean,
    auc_ratio_mean,
    auc_dam_group_q25_q75
  ) %>%
  inner_join(metab_log, by = "sample_id_clean") %>%
  filter(
    auc_dam_group_q25_q75 %in% c(
      "Low_AUC_strong_inhibition",
      "High_AUC_weak_inhibition"
    )
  )

dim(dam_data)

dam_data %>%
  count(auc_dam_group_q25_q75)

# =========================
# 9. DAM analysis
# Low AUC strong inhibition vs High AUC weak inhibition
# =========================

dam_results <- map_dfr(peak_cols, function(peak) {
  
  df <- dam_data %>%
    dplyr::select(auc_dam_group_q25_q75, all_of(peak)) %>%
    drop_na()
  
  low_values <- df %>%
    filter(auc_dam_group_q25_q75 == "Low_AUC_strong_inhibition") %>%
    pull(all_of(peak))
  
  high_values <- df %>%
    filter(auc_dam_group_q25_q75 == "High_AUC_weak_inhibition") %>%
    pull(all_of(peak))
  
  if (
    length(low_values) < 3 ||
    length(high_values) < 3 ||
    sd(low_values) == 0 ||
    sd(high_values) == 0
  ) {
    return(
      tibble(
        peak_id = peak,
        n_low = length(low_values),
        n_high = length(high_values),
        mean_low_auc = mean(low_values, na.rm = TRUE),
        mean_high_auc = mean(high_values, na.rm = TRUE),
        log2FC_low_vs_high = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- t.test(low_values, high_values)
  
  tibble(
    peak_id = peak,
    n_low = length(low_values),
    n_high = length(high_values),
    mean_low_auc = mean(low_values, na.rm = TRUE),
    mean_high_auc = mean(high_values, na.rm = TRUE),
    log2FC_low_vs_high = mean(low_values, na.rm = TRUE) - mean(high_values, na.rm = TRUE),
    p_value = test$p.value
  )
}) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_value),
    direction = case_when(
      log2FC_low_vs_high > 0 ~ "Higher in Low AUC / strong inhibition",
      log2FC_low_vs_high < 0 ~ "Higher in High AUC / weak inhibition",
      TRUE ~ "No direction"
    ),
    nominal_significant = if_else(p_value < 0.05, "p < 0.05", "NS"),
    fdr_significant = if_else(fdr < 0.05, "FDR < 0.05", "NS"),
    dam_candidate_nominal = if_else(
      abs(log2FC_low_vs_high) > 0.58 & p_value < 0.05, # log2(1.5)=0.58，这是 metabolite abundance 的 effect size threshold，代谢物丰度差 1.5 倍以上
      "DAM candidate",
      "Other"
    )
  ) %>%
  arrange(p_value)

dam_results <- dam_results %>%
  mutate(
    dam_plot_group = case_when(
      log2FC_low_vs_high > 0 &
        abs(log2FC_low_vs_high) > 0.58 &
        p_value < 0.05 ~ "Higher in Low AUC\n(candidate)",
      
      log2FC_low_vs_high < 0 &
        abs(log2FC_low_vs_high) > 0.58 &
        p_value < 0.05 ~ "Higher in High AUC\n(candidate)",
      
      p_value < 0.05 ~ "p < 0.05 only",
      
      TRUE ~ "Other"
    )
  )

dam_results %>%
  slice(1:20)

# =========================
# 10. DAM summary
# =========================

dam_summary <- dam_results %>%
  summarise(
    total_peaks = n(),
    nominal_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
    fdr_lt_0_05 = sum(fdr < 0.05, na.rm = TRUE),
    log2fc_abs_gt_0_58 = sum(abs(log2FC_low_vs_high) > 0.58, na.rm = TRUE),
    dam_candidate_nominal = sum(
      abs(log2FC_low_vs_high) > 0.58 & p_value < 0.05,
      na.rm = TRUE
    ),
    higher_in_low_auc_p_lt_0_05 = sum(
      log2FC_low_vs_high > 0 & p_value < 0.05,
      na.rm = TRUE
    ),
    higher_in_high_auc_p_lt_0_05 = sum(
      log2FC_low_vs_high < 0 & p_value < 0.05,
      na.rm = TRUE
    )
  )

dam_summary

# =========================
# 11. Add annotation to DAM results
# =========================

dam_results_anno <- dam_results %>%
  left_join(
    anno_clean,
    by = c("peak_id" = "peak_id_clean")
  )

dam_results_anno %>%
  dplyr::select(
    peak_id,
    log2FC_low_vs_high,
    p_value,
    fdr,
    direction,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:20)

# =========================
# 12. DAM volcano plot
# =========================
p_dam_volcano <- ggplot(
  dam_results,
  aes(x = log2FC_low_vs_high, y = neg_log10_p)
) +
  geom_point(
    aes(color = dam_plot_group),
    size = 2.3,
    alpha = 0.80
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Higher in Low AUC\n(candidate)" = "#6FA8DC",   # 蓝色
      "Higher in High AUC\n(candidate)" = "#F4A6B7",  # 粉色
      "p < 0.05 only" = "#B4A7D6",                    # 淡紫
      "Other" = "grey75"
    )
  ) +
  theme_prism_user() +
  labs(
    title = "Differential metabolite abundance between AUC groups",
    x = "log2FC: Low AUC vs High AUC",
    y = "-log10(p-value)",
    color = NULL
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1))
  )

p_dam_volcano
# =========================
# 13. Save DAM outputs
# =========================

write_csv(
  dam_results,
  file.path(output_dir, "06_correlation_DAM", "DAM_lowAUC_vs_highAUC_results.csv")
)

write_csv(
  dam_results_anno,
  file.path(output_dir, "06_correlation_DAM", "DAM_lowAUC_vs_highAUC_results_with_annotation.csv")
)

write_csv(
  dam_summary,
  file.path(output_dir, "06_correlation_DAM", "DAM_summary.csv")
)

save_prism_plot(
  p_dam_volcano,
  file.path(output_dir, "06_correlation_DAM", "DAM_lowAUC_vs_highAUC_volcano.png"),
  width = 6,
  height = 4,
  dpi = 600
)

# =========================
# 14. Integrate Pearson and DAM results
# =========================

integrated_results <- pearson_results %>%
  dplyr::select(
    peak_id,
    pearson_n = n,
    pearson_r,
    pearson_p = p_value,
    pearson_fdr = fdr,
    pearson_direction = direction
  ) %>%
  left_join(
    dam_results %>%
      dplyr::select(
        peak_id,
        n_low,
        n_high,
        mean_low_auc,
        mean_high_auc,
        log2FC_low_vs_high,
        dam_p = p_value,
        dam_fdr = fdr,
        dam_direction = direction
      ),
    by = "peak_id"
  ) %>%
  mutate(
    association_pattern = case_when(
      pearson_r < 0 & log2FC_low_vs_high > 0 ~
        "Consistent: inhibition-associated",
      
      pearson_r > 0 & log2FC_low_vs_high < 0 ~
        "Consistent: high-AUC-associated",
      
      pearson_r < 0 & log2FC_low_vs_high < 0 ~
        "Discordant",
      
      pearson_r > 0 & log2FC_low_vs_high > 0 ~
        "Discordant",
      
      TRUE ~ "Unclear"
    ),
    
    priority = case_when(
      pearson_r < 0 &
        log2FC_low_vs_high > 0 &
        pearson_fdr < 0.05 &
        dam_fdr < 0.05 &
        abs(log2FC_low_vs_high) > 0.58 ~
        "High priority inhibition-associated",
      
      pearson_r < 0 &
        log2FC_low_vs_high > 0 &
        pearson_p < 0.05 &
        dam_p < 0.05 &
        abs(log2FC_low_vs_high) > 0.58 ~
        "Medium priority inhibition-associated",
      
      pearson_r > 0 &
        log2FC_low_vs_high < 0 &
        pearson_fdr < 0.05 &
        dam_fdr < 0.05 &
        abs(log2FC_low_vs_high) > 0.58 ~
        "High priority high-AUC-associated",
      
      pearson_p < 0.05 | dam_p < 0.05 ~
        "Statistical signal only",
      
      TRUE ~ "Other"
    ),
    
    combined_score = -log10(pearson_p) +
      -log10(dam_p) +
      abs(pearson_r) * 10 +
      abs(log2FC_low_vs_high)
  ) %>%
  arrange(
    factor(
      priority,
      levels = c(
        "High priority inhibition-associated",
        "Medium priority inhibition-associated",
        "High priority high-AUC-associated",
        "Statistical signal only",
        "Other"
      )
    ),
    desc(combined_score)
  )

integrated_results %>%
  slice(1:20)

# =========================
# 15. Add annotation to integrated candidate table
# =========================

integrated_results_anno <- integrated_results %>%
  left_join(
    anno_clean,
    by = c("peak_id" = "peak_id_clean")
  )

integrated_results_anno %>%
  dplyr::select(
    peak_id,
    priority,
    association_pattern,
    pearson_r,
    pearson_p,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_p,
    dam_fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:25)

# =========================
# 16. Integrated result summary
# =========================

integrated_priority_summary <- integrated_results %>%
  count(priority) %>%
  arrange(desc(n))

integrated_pattern_summary <- integrated_results %>%
  count(association_pattern) %>%
  arrange(desc(n))

integrated_priority_summary
integrated_pattern_summary

# =========================
# 17. Extract top candidate tables
# =========================

top_inhibition_candidates <- integrated_results_anno %>%
  filter(priority == "High priority inhibition-associated") %>%
  arrange(desc(combined_score))

top_high_auc_candidates <- integrated_results_anno %>%
  filter(priority == "High priority high-AUC-associated") %>%
  arrange(desc(combined_score))

top_inhibition_candidates %>%
  dplyr::select(
    peak_id,
    pearson_r,
    pearson_p,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_p,
    dam_fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:20)

top_high_auc_candidates %>%
  dplyr::select(
    peak_id,
    pearson_r,
    pearson_p,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_p,
    dam_fdr,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:20)

# =========================
# 18. Save integrated outputs
# =========================

write_csv(
  integrated_results,
  file.path(output_dir, "06_correlation_DAM", "integrated_pearson_DAM_results.csv")
)

write_csv(
  integrated_results_anno,
  file.path(output_dir, "06_correlation_DAM", "integrated_pearson_DAM_results_with_annotation.csv")
)

write_csv(
  integrated_priority_summary,
  file.path(output_dir, "06_correlation_DAM", "integrated_priority_summary.csv")
)

write_csv(
  integrated_pattern_summary,
  file.path(output_dir, "06_correlation_DAM", "integrated_pattern_summary.csv")
)

write_csv(
  top_inhibition_candidates,
  file.path(output_dir, "06_correlation_DAM", "top_inhibition_associated_candidates.csv")
)

write_csv(
  top_high_auc_candidates,
  file.path(output_dir, "06_correlation_DAM", "top_high_AUC_associated_candidates.csv")
)