# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 06_final_candidate_clean_table.R
# Clean final candidate tables for slide/report
# =========================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))

# =========================
# 1. Read integrated candidate table
# =========================

candidate_vip_integrated <- read_csv(
  file.path(output_07_PLS_VIP_candidate, "candidate_integrated_with_VIP.csv"),
  show_col_types = FALSE
)

dim(candidate_vip_integrated)
names(candidate_vip_integrated)

# =========================
# 2. Define final candidate tiers
# =========================

final_candidate_levels <- c(
  "Tier 1: robust inhibition + VIP >= 2",
  "Tier 2: robust inhibition + VIP 1.5–2",
  "Tier 3: nominal inhibition + VIP >= 2",
  "High-AUC-associated + VIP >= 2",
  "Original inhibition candidate + VIP >= 2",
  "Other"
)

# =========================
# 3. Clean all-candidate evidence table
# =========================

all_candidate_evidence_clean <- candidate_vip_integrated %>%
  mutate(
    final_candidate_tier = factor(
      final_candidate_tier,
      levels = final_candidate_levels
    ),
    
    final_tier_short = case_when(
      final_candidate_tier == "Tier 1: robust inhibition + VIP >= 2" ~ "Tier 1",
      final_candidate_tier == "Tier 2: robust inhibition + VIP 1.5–2" ~ "Tier 2",
      final_candidate_tier == "Tier 3: nominal inhibition + VIP >= 2" ~ "Tier 3",
      final_candidate_tier == "High-AUC-associated + VIP >= 2" ~ "High-AUC VIP",
      final_candidate_tier == "Original inhibition candidate + VIP >= 2" ~ "Original only",
      TRUE ~ "Other"
    ),
    
    association_direction = case_when(
      pearson_r < 0 & log2FC_low_vs_high > 0 ~
        "Inhibition-associated",
      pearson_r > 0 & log2FC_low_vs_high < 0 ~
        "High-AUC-associated",
      TRUE ~ "Other / discordant"
    ),
    
    putative_annotation = case_when(
      !is.na(analog_compound_name) & analog_compound_name != "" ~ analog_compound_name,
      TRUE ~ "Unknown / not assigned"
    ),
    
    chemical_class = case_when(
      !is.na(cf_class) & cf_class != "" ~ cf_class,
      TRUE ~ "Unknown"
    ),
    
    chemical_subclass = case_when(
      !is.na(cf_subclass) & cf_subclass != "" ~ cf_subclass,
      TRUE ~ "Unknown"
    ),
    
    npc_pathway = case_when(
      !is.na(npc_pathway_results) & npc_pathway_results != "" ~ npc_pathway_results,
      TRUE ~ "Unknown"
    )
  ) %>%
  transmute(
    peak_id,
    final_tier_short,
    final_candidate_tier = as.character(final_candidate_tier),
    association_direction,
    
    # Pearson evidence
    pearson_r = round(pearson_r, 3),
    pearson_fdr = signif(pearson_fdr, 3),
    
    # DAM evidence
    log2FC_low_vs_high = round(log2FC_low_vs_high, 3),
    dam_fdr = signif(dam_fdr, 3),
    
    # Batch-adjusted evidence
    adj_pearson_beta = round(adj_pearson_beta, 3),
    adj_pearson_fdr = signif(adj_pearson_fdr, 3),
    adj_log2FC_low_vs_high = round(adj_log2FC_low_vs_high, 3),
    adj_dam_fdr = signif(adj_dam_fdr, 3),
    batch_adjusted_support,
    
    # PLS-VIP evidence
    vip_max = round(vip_max, 3),
    vip_mean = round(vip_mean, 3),
    vip_tier,
    
    # Annotation
    putative_annotation,
    chemical_class,
    chemical_subclass,
    npc_pathway,
    
    # Original metadata from analysis
    priority,
    association_pattern
  ) %>%
  arrange(
    final_candidate_tier,
    desc(vip_max),
    pearson_fdr,
    dam_fdr
  )

all_candidate_evidence_clean %>%
  slice(1:30)

# =========================
# 4. Final inhibition candidate clean table
# =========================

final_inhibition_candidate_clean <- all_candidate_evidence_clean %>%
  filter(
    final_candidate_tier %in% c(
      "Tier 1: robust inhibition + VIP >= 2",
      "Tier 2: robust inhibition + VIP 1.5–2",
      "Tier 3: nominal inhibition + VIP >= 2"
    )
  ) %>%
  arrange(
    final_candidate_tier,
    desc(vip_max),
    pearson_fdr,
    dam_fdr
  )

final_inhibition_candidate_clean %>%
  print(n = 50)

# =========================
# 5. Slide-friendly top candidate table
# =========================

slide_top_n <- 15

final_inhibition_candidate_slide_top <- final_inhibition_candidate_clean %>%
  slice_head(n = slide_top_n) %>%
  transmute(
    Peak = peak_id,
    Tier = final_tier_short,
    `Pearson r` = pearson_r,
    `DAM log2FC` = log2FC_low_vs_high,
    VIP = vip_max,
    `Putative annotation` = putative_annotation,
    Class = chemical_class,
    Pathway = npc_pathway
  )

final_inhibition_candidate_slide_top

final_inhibition_candidate_slide_minimal <- final_inhibition_candidate_clean %>%
  slice_head(n = slide_top_n) %>%
  transmute(
    Peak = peak_id,
    Tier = final_tier_short,
    `r` = pearson_r,
    `log2FC` = log2FC_low_vs_high,
    VIP = vip_max,
    Annotation = putative_annotation
  )

final_inhibition_candidate_slide_minimal
# =========================
# 6. Summary tables
# =========================

final_candidate_summary_clean <- all_candidate_evidence_clean %>%
  count(final_candidate_tier, final_tier_short, name = "n") %>%
  arrange(final_candidate_tier)

final_candidate_summary_clean

final_inhibition_annotation_class_summary <- final_inhibition_candidate_clean %>%
  count(chemical_class, name = "n") %>%
  arrange(desc(n))

final_inhibition_annotation_class_summary

final_inhibition_pathway_summary <- final_inhibition_candidate_clean %>%
  count(npc_pathway, name = "n") %>%
  arrange(desc(n))

final_inhibition_pathway_summary

# =========================
# 7. Save outputs
# =========================

make_output_dir(output_08_final_candidate_tables)

write_csv(
  all_candidate_evidence_clean,
  file.path(output_08_final_candidate_tables, "all_candidate_evidence_clean_table.csv")
)

write_csv(
  final_inhibition_candidate_clean,
  file.path(output_08_final_candidate_tables, "final_inhibition_candidate_clean_table.csv")
)

write_csv(
  final_inhibition_candidate_slide_top,
  file.path(output_08_final_candidate_tables, "final_inhibition_candidate_slide_top15.csv")
)

write_csv(
  final_inhibition_candidate_slide_minimal,
  file.path(output_08_final_candidate_tables, "final_inhibition_candidate_slide_top15_minimal.csv")
)

write_csv(
  final_candidate_summary_clean,
  file.path(output_08_final_candidate_tables, "final_candidate_summary_clean.csv")
)

write_csv(
  final_inhibition_annotation_class_summary,
  file.path(output_08_final_candidate_tables, "final_inhibition_annotation_class_summary.csv")
)

write_csv(
  final_inhibition_pathway_summary,
  file.path(output_08_final_candidate_tables, "final_inhibition_pathway_summary.csv")
)