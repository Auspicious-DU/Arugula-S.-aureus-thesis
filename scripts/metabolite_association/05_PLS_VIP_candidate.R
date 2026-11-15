# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 05_PLS_VIP_candidate.R
# PLS-VIP and candidate integration
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))
source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))

library(mixOmics)

# Reload tidyverse components after mixOmics/MASS
# to reduce select()/filter() conflicts
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)
library(readr)
library(stringr)

# =========================
# 1. Read data
# =========================

pheno_qc <- read_csv(
  file.path(output_dir, "05_phenotype_qc", "phenotype_qc_with_auc_groups.csv"),
  show_col_types = FALSE
)

metab_scaled <- read_csv(
  file.path(output_dir, "04_preprocessed_data", "metabolite_log2_zscaled_matrix.csv"),
  show_col_types = FALSE
)

anno_clean <- read_csv(
  file.path(output_dir, "02_matched_data", "peak_annotation_clean.csv"),
  show_col_types = FALSE
)

batch_adjusted_integrated <- read_csv(
  file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_integrated_results.csv"),
  show_col_types = FALSE
)

dim(pheno_qc)
dim(metab_scaled)
dim(anno_clean)
dim(batch_adjusted_integrated)

# =========================
# 2. Prepare X and Y for PLS
# =========================

peak_cols <- names(metab_scaled) %>%
  keep(~ str_detect(.x, "^peak\\d+$"))

length(peak_cols)
head(peak_cols)
tail(peak_cols)

pls_data <- pheno_qc %>%
  dplyr::select(
    sample_id_clean,
    auc_ratio_mean,
    auc_dam_group_q25_q75,
    auc_class_0_7_1_5
  ) %>%
  dplyr::inner_join(metab_scaled, by = "sample_id_clean") %>%
  tidyr::drop_na(auc_ratio_mean)

dim(pls_data)

X <- pls_data %>%
  dplyr::select(all_of(peak_cols)) %>%
  as.matrix()

Y <- pls_data$auc_ratio_mean

dim(X)
length(Y)
summary(Y)

# =========================
# 3. Run PLS regression
# =========================

set.seed(123)

pls_ncomp <- 2

pls_model <- mixOmics::pls(
  X = X,
  Y = Y,
  ncomp = pls_ncomp,
  mode = "regression"
)

pls_model

# =========================
# 4. Extract PLS scores
# =========================

pls_score_mat <- pls_model$variates$X[, 1:pls_ncomp, drop = FALSE]

pls_scores <- as_tibble(pls_score_mat)
names(pls_scores) <- paste0("PLS_comp", seq_len(ncol(pls_scores)))

pls_scores <- pls_scores %>%
  mutate(
    sample_id_clean = pls_data$sample_id_clean,
    auc_ratio_mean = pls_data$auc_ratio_mean,
    auc_dam_group_q25_q75 = pls_data$auc_dam_group_q25_q75,
    auc_class_0_7_1_5 = pls_data$auc_class_0_7_1_5
  )

head(pls_scores)

# =========================
# 5. PLS score plot by AUC group
# =========================

p_pls_score_auc <- ggplot(
  pls_scores,
  aes(x = PLS_comp1, y = PLS_comp2, color = auc_dam_group_q25_q75)
) +
  geom_point(size = 2.2, alpha = 0.75) +
  scale_color_manual(
    values = c(
      "Low_AUC_strong_inhibition" = "#6FA8DC",
      "High_AUC_weak_inhibition" = "#F4A6B7",
      "Middle" = "grey75"
    )
  ) +
  theme_prism_user() +
  labs(
    title = "PLS score plot by AUC group",
    x = "PLS component 1",
    y = "PLS component 2",
    color = NULL
  )

p_pls_score_auc

# =========================
# 6. Calculate VIP scores
# =========================

vip_mat <- mixOmics::vip(pls_model)

vip_df <- as_tibble(vip_mat, rownames = "peak_id")

# Force clean VIP column names
vip_value_cols <- names(vip_df)[names(vip_df) != "peak_id"]
names(vip_df)[names(vip_df) != "peak_id"] <- paste0("VIP_comp", seq_along(vip_value_cols))

vip_cols <- names(vip_df) %>%
  keep(~ str_detect(.x, "^VIP_comp"))

vip_df <- vip_df %>%
  mutate(
    vip_max = do.call(pmax, c(dplyr::select(., all_of(vip_cols)), na.rm = TRUE)),
    vip_mean = rowMeans(dplyr::select(., all_of(vip_cols)), na.rm = TRUE),
    vip_tier = case_when(
      vip_max >= 2 ~ "VIP >= 2",
      vip_max >= 1.5 ~ "VIP 1.5–2",
      vip_max >= 1 ~ "VIP 1–1.5",
      TRUE ~ "VIP < 1"
    )
  ) %>%
  arrange(desc(vip_max))

vip_df %>% slice(1:20)

vip_df %>%
  count(vip_tier)

# =========================
# 7. Add annotation to VIP table
# =========================

vip_anno <- vip_df %>%
  left_join(
    anno_clean,
    by = c("peak_id" = "peak_id_clean")
  )

vip_anno %>%
  dplyr::select(
    peak_id,
    vip_max,
    vip_tier,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:25)

# =========================
# 8. Top VIP bar plot
# =========================

top_vip_n <- 25

top_vip_plot_df <- vip_anno %>%
  slice_max(vip_max, n = top_vip_n) %>%
  mutate(
    peak_label = if_else(
      is.na(analog_compound_name) | analog_compound_name == "",
      peak_id,
      paste0(peak_id, "\n", analog_compound_name)
    ),
    peak_label = fct_reorder(peak_id, vip_max)
  )

p_top_vip <- ggplot(
  top_vip_plot_df,
  aes(x = peak_label, y = vip_max, fill = vip_tier)
) +
  geom_col(color = "black", linewidth = 0.25, alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "VIP >= 2" = "#6FA8DC",
      "VIP 1.5–2" = "#B4A7D6",
      "VIP 1–1.5" = "grey70",
      "VIP < 1" = "grey85"
    )
  ) +
  theme_prism_user() +
  labs(
    title = "Top VIP metabolite features from PLS regression",
    x = NULL,
    y = "VIP score",
    fill = NULL
  )

p_top_vip

# =========================
# 9. Integrate VIP with batch-adjusted candidate table
# =========================

candidate_vip_integrated <- batch_adjusted_integrated %>%
  left_join(
    vip_df %>%
      dplyr::select(
        peak_id,
        starts_with("VIP_comp"),
        vip_max,
        vip_mean,
        vip_tier
      ),
    by = "peak_id"
  ) %>%
  mutate(
    final_candidate_tier = case_when(
      batch_adjusted_support ==
        "Robust inhibition-associated after LC-MS batch adjustment" &
        vip_max >= 2 ~
        "Tier 1: robust inhibition + VIP >= 2",
      
      batch_adjusted_support ==
        "Robust inhibition-associated after LC-MS batch adjustment" &
        vip_max >= 1.5 ~
        "Tier 2: robust inhibition + VIP 1.5–2",
      
      batch_adjusted_support ==
        "Nominally supported inhibition-associated after LC-MS batch adjustment" &
        vip_max >= 2 ~
        "Tier 3: nominal inhibition + VIP >= 2",
      
      batch_adjusted_support ==
        "Robust high-AUC-associated after LC-MS batch adjustment" &
        vip_max >= 2 ~
        "High-AUC-associated + VIP >= 2",
      
      priority == "High priority inhibition-associated" &
        vip_max >= 2 ~
        "Original inhibition candidate + VIP >= 2",
      
      TRUE ~
        "Other"
    )
  ) %>%
  arrange(
    factor(
      final_candidate_tier,
      levels = c(
        "Tier 1: robust inhibition + VIP >= 2",
        "Tier 2: robust inhibition + VIP 1.5–2",
        "Tier 3: nominal inhibition + VIP >= 2",
        "High-AUC-associated + VIP >= 2",
        "Original inhibition candidate + VIP >= 2",
        "Other"
      )
    ),
    desc(vip_max)
  )

candidate_vip_integrated %>%
  dplyr::select(
    peak_id,
    final_candidate_tier,
    batch_adjusted_support,
    priority,
    pearson_r,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_fdr,
    adj_pearson_beta,
    adj_pearson_fdr,
    adj_log2FC_low_vs_high,
    adj_dam_fdr,
    vip_max,
    vip_tier,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  slice(1:30)

# =========================
# 10. Final candidate summary
# =========================

final_candidate_summary <- candidate_vip_integrated %>%
  count(final_candidate_tier) %>%
  arrange(desc(n))

final_candidate_summary

# =========================
# 11. Extract final inhibition shortlist
# =========================

final_inhibition_shortlist <- candidate_vip_integrated %>%
  filter(
    final_candidate_tier %in% c(
      "Tier 1: robust inhibition + VIP >= 2",
      "Tier 2: robust inhibition + VIP 1.5–2",
      "Tier 3: nominal inhibition + VIP >= 2"
    )
  ) %>%
  arrange(
    factor(
      final_candidate_tier,
      levels = c(
        "Tier 1: robust inhibition + VIP >= 2",
        "Tier 2: robust inhibition + VIP 1.5–2",
        "Tier 3: nominal inhibition + VIP >= 2"
      )
    ),
    desc(vip_max)
  )

final_inhibition_shortlist %>%
  dplyr::select(
    peak_id,
    final_candidate_tier,
    pearson_r,
    pearson_fdr,
    log2FC_low_vs_high,
    dam_fdr,
    adj_pearson_beta,
    adj_pearson_fdr,
    adj_log2FC_low_vs_high,
    adj_dam_fdr,
    vip_max,
    analog_compound_name,
    cf_class,
    cf_subclass,
    npc_pathway_results
  ) %>%
  print(n = 50)

# =========================
# 12. Integrated candidate plot
# Pearson direction + DAM direction + VIP
# =========================

candidate_plot_df <- candidate_vip_integrated %>%
  mutate(
    final_plot_group = case_when(
      final_candidate_tier == "Tier 1: robust inhibition + VIP >= 2" ~
        "Tier 1 robust inhibition",
      final_candidate_tier == "Tier 2: robust inhibition + VIP 1.5–2" ~
        "Tier 2 robust inhibition",
      final_candidate_tier == "Tier 3: nominal inhibition + VIP >= 2" ~
        "Tier 3 nominal inhibition",
      final_candidate_tier == "High-AUC-associated + VIP >= 2" ~
        "High-AUC + VIP >= 2",
      TRUE ~ "Other"
    ),
    final_plot_group = factor(
      final_plot_group,
      levels = c(
        "Tier 1 robust inhibition",
        "Tier 2 robust inhibition",
        "Tier 3 nominal inhibition",
        "High-AUC + VIP >= 2",
        "Other"
      )
    )
  )

p_candidate_vip_scatter <- ggplot(
  candidate_plot_df,
  aes(
    x = pearson_r,
    y = log2FC_low_vs_high,
    color = final_plot_group,
    size = vip_max
  )
) +
  geom_point(alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Tier 1 robust inhibition" = "#6FA8DC",
      "Tier 2 robust inhibition" = "#B4A7D6",
      "Tier 3 nominal inhibition" = "#8DD3C7",
      "High-AUC + VIP >= 2" = "#F4A6B7",
      "Other" = "grey75"
    )
  ) +
  scale_size_continuous(range = c(1.5, 5)) +
  theme_prism_user() +
  labs(
    title = "Integrated candidate metabolite features",
    x = "Pearson r with AUC-ratio",
    y = "DAM log2FC: Low AUC vs High AUC",
    color = NULL,
    size = "VIP"
  )

p_candidate_vip_scatter

# =========================
# 13. Final shortlist VIP bar plot
# =========================

final_shortlist_plot_df <- final_inhibition_shortlist %>%
  slice_max(vip_max, n = 25) %>%
  mutate(
    peak_label = if_else(
      is.na(analog_compound_name) | analog_compound_name == "",
      peak_id,
      paste0(peak_id, "\n", analog_compound_name)
    ),
    peak_label = fct_reorder(peak_id, vip_max)
  )

p_final_shortlist_vip <- ggplot(
  final_shortlist_plot_df,
  aes(x = peak_label, y = vip_max, fill = final_candidate_tier)
) +
  geom_col(color = "black", linewidth = 0.25, alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Tier 1: robust inhibition + VIP >= 2" = "#6FA8DC",
      "Tier 2: robust inhibition + VIP 1.5–2" = "#B4A7D6",
      "Tier 3: nominal inhibition + VIP >= 2" = "#8DD3C7"
    )
  ) +
  theme_prism_user() +
  labs(
    title = "Final inhibition-associated candidate features",
    x = NULL,
    y = "VIP score",
    fill = NULL
  )

p_final_shortlist_vip

# =========================
# 14. Save outputs
# =========================

dir.create(file.path(output_dir, "07_PLS_VIP_candidate"), showWarnings = FALSE)

write_csv(
  pls_scores,
  file.path(output_dir, "07_PLS_VIP_candidate", "PLS_scores.csv")
)

write_csv(
  vip_df,
  file.path(output_dir, "07_PLS_VIP_candidate", "PLS_VIP_scores.csv")
)

write_csv(
  vip_anno,
  file.path(output_dir, "07_PLS_VIP_candidate", "PLS_VIP_scores_with_annotation.csv")
)

write_csv(
  candidate_vip_integrated,
  file.path(output_dir, "07_PLS_VIP_candidate", "candidate_integrated_with_VIP.csv")
)

write_csv(
  final_candidate_summary,
  file.path(output_dir, "07_PLS_VIP_candidate", "final_candidate_summary.csv")
)

write_csv(
  final_inhibition_shortlist,
  file.path(output_dir, "07_PLS_VIP_candidate", "final_inhibition_candidate_shortlist.csv")
)

save_prism_plot(
  p_pls_score_auc,
  file.path(output_dir, "07_PLS_VIP_candidate", "PLS_score_plot_by_AUC_group.png"),
  width = 6,
  height = 4.5,
  dpi = 600
)

save_prism_plot(
  p_top_vip,
  file.path(output_dir, "07_PLS_VIP_candidate", "Top_25_VIP_features.png"),
  width = 7,
  height = 6,
  dpi = 600
)

save_prism_plot(
  p_candidate_vip_scatter,
  file.path(output_dir, "07_PLS_VIP_candidate", "Integrated_candidate_scatter_Pearson_DAM_VIP.png"),
  width = 7,
  height = 5,
  dpi = 600
)

save_prism_plot(
  p_final_shortlist_vip,
  file.path(output_dir, "07_PLS_VIP_candidate", "Final_inhibition_candidate_VIP_barplot.png"),
  width = 7,
  height = 6,
  dpi = 600
)