# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# ============================================================
# THESIS_metabolite_association_figures.R
# Thesis-specific figures for Section 3.4:
# 1) Pearson correlation overview
# 2) DAM overview
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

source(file.path(repo_dir, "scripts", "metabolite_association", "00_setup.R"))

source(file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(patchwork)
})

# ============================================================
# 1. Output folders
# ============================================================

thesis_figure_dir <- file.path(repo_dir, "output", "thesis_figures")

thesis_plot_dir <- file.path(thesis_figure_dir, "plots")

thesis_table_dir <- file.path(thesis_figure_dir, "tables")

dir.create(
  thesis_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  thesis_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Input files
# ============================================================

pearson_file <- file.path(output_dir, "06_correlation_DAM", "pearson_auc_ratio_metabolite_results.csv")

dam_file <- file.path(output_dir, "06_correlation_DAM", "DAM_lowAUC_vs_highAUC_results.csv")

integrated_file <- file.path(output_dir, "06_correlation_DAM", "integrated_pearson_DAM_results.csv")

batch_file <- file.path(output_dir, "06B_batch_adjusted_sensitivity", "batch_adjusted_integrated_results.csv")

# ============================================================
# 3. Read data
# ============================================================

pearson_df <- read_csv(
  pearson_file,
  show_col_types = FALSE
)

dam_df <- read_csv(
  dam_file,
  show_col_types = FALSE
)

integrated_df <- read_csv(
  integrated_file,
  show_col_types = FALSE
)

batch_df <- read_csv(
  batch_file,
  show_col_types = FALSE
)

# ============================================================
# 4. Quick summaries for writing
# ============================================================

pearson_summary <- pearson_df %>%
  count(pearson_plot_group, name = "n_features")

dam_summary <- dam_df %>%
  count(dam_plot_group, name = "n_features")

batch_summary <- batch_df %>%
  count(batch_adjusted_support, name = "n_features")

write_csv(
  pearson_summary,
  file.path(thesis_table_dir, "Figure_3_3_pearson_summary.csv")
)

write_csv(
  dam_summary,
  file.path(thesis_table_dir, "Figure_3_3_DAM_summary.csv")
)

write_csv(
  batch_summary,
  file.path(thesis_table_dir, "Figure_3_3_batch_adjusted_summary.csv")
)

# Optional overall summary table
overall_summary <- tibble(
  analysis = c(
    "Pearson nominal p < 0.05",
    "Pearson FDR < 0.05",
    "Pearson negative with AUC (FDR < 0.05)",
    "Pearson positive with AUC (FDR < 0.05)",
    "DAM nominal p < 0.05",
    "DAM FDR < 0.05",
    "DAM candidate",
    "DAM higher in Low AUC (candidate)",
    "DAM higher in High AUC (candidate)"
  ),
  n_features = c(
    sum(pearson_df$p_value < 0.05, na.rm = TRUE),
    sum(pearson_df$fdr < 0.05, na.rm = TRUE),
    sum(pearson_df$pearson_plot_group == "Negative with AUC\n(FDR < 0.05)", na.rm = TRUE),
    sum(pearson_df$pearson_plot_group == "Positive with AUC\n(FDR < 0.05)", na.rm = TRUE),
    sum(dam_df$p_value < 0.05, na.rm = TRUE),
    sum(dam_df$fdr < 0.05, na.rm = TRUE),
    sum(dam_df$dam_candidate_nominal == "DAM candidate", na.rm = TRUE),
    sum(dam_df$dam_plot_group == "Higher in Low AUC\n(candidate)", na.rm = TRUE),
    sum(dam_df$dam_plot_group == "Higher in High AUC\n(candidate)", na.rm = TRUE)
  )
)

write_csv(
  overall_summary,
  file.path(thesis_table_dir, "Figure_3_3_overall_summary.csv")
)

# ============================================================
# 5. Thresholds and labels
# ============================================================

P_REF <- -log10(0.05)
LOG2FC_REF <- log2(1.5)

pearson_n_label <- paste0(
  "n = ",
  nrow(pearson_df)
)

dam_n_label <- paste0(
  "n = ",
  nrow(dam_df)
)

dam_group_label <- paste0(
  "Low AUC: n = ",
  unique(dam_df$n_low),
  "\nHigh AUC: n = ",
  unique(dam_df$n_high)
)


# ============================================================
# 6. Shared plotting categories
# ============================================================

# Common legend labels across Pearson and DAM panels.
# Statistical criteria remain specific to each analysis.

plot_levels <- c(
  "Low-AUC direction",
  "High-AUC direction",
  "Nominal only",
  "Other"
)

pearson_df <- pearson_df %>%
  mutate(
    plot_class = case_when(
      pearson_plot_group ==
        "Negative with AUC\n(FDR < 0.05)" ~
        "Low-AUC direction",
      
      pearson_plot_group ==
        "Positive with AUC\n(FDR < 0.05)" ~
        "High-AUC direction",
      
      pearson_plot_group == "p < 0.05 only" ~
        "Nominal only",
      
      TRUE ~ "Other"
    ),
    plot_class = factor(plot_class, levels = plot_levels)
  ) %>%
  arrange(
    factor(
      plot_class,
      levels = c(
        "Other", "Nominal only",
        "Low-AUC direction", "High-AUC direction"
      )
    )
  )

dam_df <- dam_df %>%
  mutate(
    plot_class = case_when(
      dam_plot_group ==
        "Higher in Low AUC\n(candidate)" ~
        "Low-AUC direction",
      
      dam_plot_group ==
        "Higher in High AUC\n(candidate)" ~
        "High-AUC direction",
      
      dam_plot_group == "p < 0.05 only" ~
        "Nominal only",
      
      TRUE ~ "Other"
    ),
    plot_class = factor(plot_class, levels = plot_levels)
  ) %>%
  arrange(
    factor(
      plot_class,
      levels = c(
        "Other", "Nominal only",
        "Low-AUC direction", "High-AUC direction"
      )
    )
  )

# ============================================================
# 7. Shared Prism colour scheme
# ============================================================

shared_fill_cols <- c(
  "Low-AUC direction" =
    prism_dataset_colour(2, "main"),
  "High-AUC direction" =
    prism_dataset_colour(5, "main"),
  "Nominal only" =
    prism_dataset_colour(1, "main"),
  "Other" = SCI_COLOUR_OTHER_FILL
)

shared_edge_cols <- c(
  "Low-AUC direction" =
    prism_dataset_colour(2, "edge"),
  "High-AUC direction" =
    prism_dataset_colour(5, "edge"),
  "Nominal only" =
    prism_dataset_colour(1, "edge"),
  "Other" = SCI_COLOUR_OTHER_EDGE
)

# ============================================================
# 8. Panel A — Pearson association overview
# ============================================================

p_pearson <- ggplot(
  pearson_df,
  aes(
    x = pearson_r,
    y = neg_log10_p,
    fill = plot_class,
    colour = plot_class
  )
) +
  geom_hline(
    yintercept = P_REF,
    linetype = "dashed",
    linewidth = SCI_AXIS_LINEWIDTH,
    colour = SCI_COLOUR_GUIDE_LINE
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = SCI_AXIS_LINEWIDTH,
    colour = SCI_COLOUR_ZERO_LINE
  ) +
  geom_point(
    shape = 21,
    size = SCI_POINT_SIZE * 1.7,
    stroke = SCI_SYMBOL_BORDER,
    alpha = SCI_ALPHA
  ) +
  scale_fill_manual(
    values = shared_fill_cols,
    breaks = plot_levels,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.5,
        alpha = 1
      )
    )
  ) +
  scale_colour_manual(
    values = shared_edge_cols,
    guide = "none"
  ) +
  labs(
    x = "Pearson correlation coefficient (r)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_sci_prism(
    legend_position = "bottom"
  )

# ============================================================
# 9. Panel B — DAM overview
# ============================================================

p_dam <- ggplot(
  dam_df,
  aes(
    x = log2FC_low_vs_high,
    y = neg_log10_p,
    fill = plot_class,
    colour = plot_class
  )
) +
  geom_hline(
    yintercept = P_REF,
    linetype = "dashed",
    linewidth = SCI_AXIS_LINEWIDTH,
    colour = SCI_COLOUR_GUIDE_LINE
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = SCI_AXIS_LINEWIDTH,
    colour = SCI_COLOUR_ZERO_LINE
  ) +
  geom_vline(
    xintercept = c(-LOG2FC_REF, LOG2FC_REF),
    linetype = "dotted",
    linewidth = SCI_AXIS_LINEWIDTH,
    colour = SCI_COLOUR_GUIDE_LINE
  ) +
  geom_point(
    shape = 21,
    size = SCI_POINT_SIZE * 1.7,
    stroke = SCI_SYMBOL_BORDER,
    alpha = SCI_ALPHA
  ) +
  scale_fill_manual(
    values = shared_fill_cols,
    breaks = plot_levels,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.5,
        alpha = 1
      )
    )
  ) +
  scale_colour_manual(
    values = shared_edge_cols,
    guide = "none"
  ) +
  labs(
    x = expression(
      log[2]*" fold change (Low AUC / High AUC)"
    ),
    y = expression(-log[10](italic(p)))
  ) +
  theme_sci_prism(
    legend_position = "bottom"
  )

# ============================================================
# 10. Combined Figure 3.3
# ============================================================

p_figure_3_3 <- (
  p_pearson + p_dam
) +
  patchwork::plot_layout(
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# ============================================================
# 11. Save plots
# ============================================================

save_sci_plot(
  plot_object = p_pearson,
  output_dir = thesis_plot_dir,
  file_stem = "Figure_3_3A_Pearson_association_overview",
  width_mm = 90,
  height_mm = 80
)

save_sci_plot(
  plot_object = p_dam,
  output_dir = thesis_plot_dir,
  file_stem = "Figure_3_3B_DAM_association_overview",
  width_mm = 90,
  height_mm = 80
)

save_sci_plot(
  plot_object = p_figure_3_3,
  output_dir = thesis_plot_dir,
  file_stem = "Figure_3_3_combined_metabolite_association_overview",
  width_mm = 180,
  height_mm = 100
)

# ============================================================
# 12. Console messages
# ============================================================

message("THESIS metabolite association figures completed.")
message("Outputs saved under: ", thesis_figure_dir)


# ============================================================
# 13. Figure 3.4 — VIP ranking + NPC pathway annotation
# This section uses the FINAL 29 inhibition-associated features.
# It does not use the older 19-feature biological-threshold set.
# Existing Figure 3.3 code above remains unchanged.
# ============================================================

vip_scores_file <- file.path(
  output_07_PLS_VIP_candidate,
  "PLS_VIP_scores.csv"
)

final_candidates_file <- file.path(
  output_08_final_candidate_tables,
  "final_inhibition_candidate_clean_table.csv"
)

if (!file.exists(vip_scores_file)) {
  stop("Missing PLS VIP file: ", vip_scores_file)
}
if (!file.exists(final_candidates_file)) {
  stop("Missing final candidate file: ", final_candidates_file)
}

vip_all_34 <- readr::read_csv(
  vip_scores_file, show_col_types = FALSE
)
candidates_34 <- readr::read_csv(
  final_candidates_file, show_col_types = FALSE
)

required_candidate_cols <- c(
  "peak_id", "final_tier_short", "vip_max", "npc_pathway"
)
if (!all(required_candidate_cols %in% names(candidates_34))) {
  stop(
    "The final candidate file is missing columns: ",
    paste(
      setdiff(required_candidate_cols, names(candidates_34)),
      collapse = ", "
    )
  )
}
if (!all(c("peak_id", "vip_max") %in% names(vip_all_34))) {
  stop("PLS_VIP_scores.csv must contain peak_id and vip_max.")
}

if (anyDuplicated(vip_all_34$peak_id) > 0 ||
    anyDuplicated(candidates_34$peak_id) > 0) {
  stop("Duplicated peak IDs found in Figure 3.4 input files.")
}
if (!all(candidates_34$peak_id %in% vip_all_34$peak_id)) {
  stop("Some final candidate IDs are missing from PLS VIP scores.")
}

# Preserve existing Tier 1 / Tier 2 colours.
tier_levels_34 <- c("Tier 1", "Tier 2", "Tier 3")
tier_fill_34 <- c(
  "Tier 1" = "#8ADCF0",
  "Tier 2" = "#CABBE3",
  "Tier 3" = "#B6C9B7",
  "Other" = SCI_GREY_LIGHT
)
tier_edge_34 <- c(
  "Tier 1" = "#43B8E8",
  "Tier 2" = "#9A83C5",
  "Tier 3" = "#80A888",
  "Other" = SCI_GREY_MID
)

# Explicit colours for NPC PATHWAY categories (not ClassyFire classes).
# The palette is intentionally distinct from the two Tier bar colours.
pathway_levels_34 <- c(
  "Fatty acids",
  "Terpenoids",
  "Shikimates and Phenylpropanoids",
  "Amino acids and Peptides",
  "Alkaloids",
  "Unknown"
)

# The SCI baseline has five Prism hue families plus neutral greys.
# Pathway categories require five non-grey groups in addition to
# the blue / purple Tier fills; use a restrained extended palette.
# Fatty acids = Prism red family; Shikimates = warm orange family;
# Amino acids = Prism teal family; green and brown are extensions.
pathway_cols_34 <- c(
  "Fatty acids" = "#E36B83",
  "Terpenoids" = "#58A66A",
  "Shikimates and Phenylpropanoids" = "#D7A13F",
  "Amino acids and Peptides" = "#4B9EA6",
  "Alkaloids" = "#966A4B",
  "Unknown" = "#BFBFBF"
)

# Treat empty/missing NPC annotations as Unknown rather than
# silently replacing them with other chemical classification fields.
candidates_34 <- candidates_34 %>%
  mutate(
    peak_id = tolower(trimws(peak_id)),
    final_tier_short = trimws(final_tier_short),
    npc_pathway = trimws(as.character(npc_pathway)),
    npc_pathway = if_else(
      is.na(npc_pathway) | npc_pathway == "",
      "Unknown", npc_pathway
    )
  )

unexpected_pathways_34 <- setdiff(
  unique(candidates_34$npc_pathway),
  pathway_levels_34
)
if (length(unexpected_pathways_34) > 0) {
  stop(
    "Unexpected NPC pathway labels, please inspect rather than guess: ",
    paste(unexpected_pathways_34, collapse = "; ")
  )
}
if (!all(candidates_34$final_tier_short %in% tier_levels_34)) {
  stop("Unexpected candidate Tier labels in final candidate table.")
}

vip_all_34 <- vip_all_34 %>%
  mutate(peak_id = tolower(trimws(peak_id)))

# Use VIP values from the same output of the final PLS model.
vip_consistency_34 <- candidates_34 %>%
  left_join(
    vip_all_34 %>%
      select(peak_id, vip_max_reference = vip_max),
    by = "peak_id"
  )

if (any(
  abs(vip_consistency_34$vip_max -
      vip_consistency_34$vip_max_reference) > 0.011,
  na.rm = TRUE
)) {
  stop("Candidate VIP scores disagree with PLS_VIP_scores.csv.")
}

# Panel A: rank all 275 features and highlight only selected candidates.
vip_panel_A_34 <- vip_all_34 %>%
  left_join(
    candidates_34 %>%
      select(peak_id, final_tier_short),
    by = "peak_id"
  ) %>%
  mutate(
    tier = if_else(
      is.na(final_tier_short), "Other", final_tier_short
    ),
    tier = factor(
      tier, levels = c("Tier 1", "Tier 2", "Tier 3", "Other")
    )
  ) %>%
  arrange(desc(vip_max), peak_id) %>%
  mutate(vip_rank = row_number())

# Panel B: same final Tier logic as in the thesis, sorted by VIP score.
candidates_panel_B_34 <- candidates_34 %>%
  arrange(
    factor(final_tier_short, levels = tier_levels_34),
    desc(vip_max),
    peak_id
  ) %>%
  mutate(
    peak_id_plot = factor(
      peak_id,
      levels = rev(peak_id)
    ),
    final_tier_short = factor(
      final_tier_short, levels = tier_levels_34
    ),
    npc_pathway = factor(
      npc_pathway, levels = pathway_levels_34
    )
  )

# Count summary for the figure legend / manuscript Results.
pathway_summary_34 <- candidates_panel_B_34 %>%
  count(npc_pathway, name = "n_features", .drop = FALSE) %>%
  filter(n_features > 0) %>%
  arrange(desc(n_features))

write_csv(
  pathway_summary_34,
  file.path(thesis_table_dir, "Figure_3_4_NPC_pathway_counts.csv")
)
write_csv(
  candidates_panel_B_34 %>%
    select(peak_id, final_tier_short, vip_max, npc_pathway),
  file.path(thesis_table_dir, "Figure_3_4_candidate_VIP_pathways.csv")
)

message(
  "Figure 3.4 inputs: ", nrow(vip_panel_A_34),
  " VIP features; ", nrow(candidates_panel_B_34),
  " inhibition candidates."
)
print(pathway_summary_34)

# ============================================================
# 14. Plot Figure 3.4 A and B — revised compact layout
# Keep the original 2-panel composition; avoid a separate
# pathway panel, which previously caused legend clipping.
# ============================================================

p_vip_34A <- ggplot(
  vip_panel_A_34,
  aes(x = vip_rank, y = vip_max, fill = tier, colour = tier)
) +
  geom_col(width = 0.98, linewidth = 0.12) +
  geom_hline(
    yintercept = c(1.5, 2),
    linetype = c("dotted", "dashed"),
    colour = SCI_COLOUR_GUIDE_LINE,
    linewidth = SCI_AXIS_LINEWIDTH
  ) +
  scale_fill_manual(values = tier_fill_34, guide = "none") +
  scale_colour_manual(values = tier_edge_34, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 50, 100, 150, 200, 250, 275),
    expand = expansion(mult = c(0, 0.004))
  ) +
  labs(x = "LC-MS features ranked by VIP score", y = "VIP score") +
  theme_sci_prism(legend_position = "none")

# Panel B now contains both the VIP bars and small coloured
# pathway squares in ONE plotting coordinate system.
# This guarantees exact row alignment and no second y-axis.
# Pathway uses colour (not fill), keeping Tier bar fill unchanged.
PATHWAY_X_34 <- 3.72
PATHWAY_AXIS_MAX_34 <- 4.05

pathway_legend_labels_34 <- c(
  "Fatty acids" = "Fatty acids (12)",
  "Terpenoids" = "Terpenoids (6)",
  "Shikimates and Phenylpropanoids" = "Shikimates / Phenylpropanoids (5)",
  "Amino acids and Peptides" = "Amino acids / Peptides (3)",
  "Alkaloids" = "Alkaloids (2)",
  "Unknown" = "Unknown (1)"
)

p_vip_34B <- ggplot(
  candidates_panel_B_34,
  aes(x = vip_max, y = peak_id_plot)
) +
  geom_vline(
    xintercept = c(1.5, 2),
    linetype = "dotted",
    colour = SCI_COLOUR_GUIDE_LINE,
    linewidth = 0.3
  ) +
  geom_col(
    aes(fill = final_tier_short),
    width = 0.80,
    colour = "#A6A6A6",
    linewidth = 0.17,
    show.legend = c(fill = TRUE, colour = FALSE)
  ) +
  geom_text(
    aes(label = sprintf("%.2f", vip_max)),
    hjust = -0.28,
    size = 2.4,
    family = SCI_FONT
  ) +
  # Shape 15 is a filled square using the COLOUR aesthetic.
  # A single plot avoids patchwork alignment and legend issues.
  geom_point(
    aes(colour = npc_pathway),
    x = PATHWAY_X_34,
    shape = 15,
    size = 3.0,
    show.legend = c(colour = TRUE, fill = FALSE)
  ) +
  annotate(
    "text",
    x = PATHWAY_X_34,
    y = nrow(candidates_panel_B_34) + 1.2,
    label = "Putative pathway",
    size = 2.5,
    family = SCI_FONT,
    hjust = 0.5
  ) +
  scale_fill_manual(
    values = tier_fill_34,
    breaks = c("Tier 1", "Tier 2"),
    labels = c("Tier 1 (n = 16)", "Tier 2 (n = 13)"),
    name = NULL,
    guide = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(colour = NA, linewidth = 0)
    )
  ) +
  scale_colour_manual(
    values = pathway_cols_34,
    breaks = pathway_levels_34,
    labels = pathway_legend_labels_34[pathway_levels_34],
    name = NULL,
    guide = guide_legend(
      order = 2, nrow = 2, byrow = TRUE,
      override.aes = list(shape = 15, size = 3.3)
    )
  ) +
  scale_x_continuous(
    limits = c(0, PATHWAY_AXIS_MAX_34),
    breaks = c(0, 1, 2, 3),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.55, 1.65))
  ) +
  labs(x = "VIP score", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_sci_prism(legend_position = "bottom") +
  theme(
    axis.text.y = element_text(size = 7.0),
    axis.line.y = element_line(linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.justification = "center",
    legend.key.size = grid::unit(3.6, "mm"),
    legend.spacing.y = grid::unit(1, "mm"),
    legend.text = element_text(size = 7.4),
    legend.margin = margin(2, 0, 0, 0),
    panel.grid.major.x = element_line(
      colour = "grey89", linewidth = 0.2, linetype = "dotted"
    ),
    plot.margin = margin(3, 8, 3, 5)
  )

# Use one shared-width plot per row and a compact lower-panel ratio.
# 'A' and 'B' labels placed by patchwork outside the plot bodies.
p_figure_3_4 <- (p_vip_34A / p_vip_34B) +
  patchwork::plot_layout(heights = c(1.0, 2.45)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

# ============================================================
# 15. Save Figure 3.4
# ============================================================

save_sci_plot(
  plot_object = p_vip_34A,
  output_dir = thesis_plot_dir,
  file_stem = "Figure_3_4A_VIP_distribution",
  width_mm = 180,
  height_mm = 65
)

# New filename: preserves the previous combined Figure 3.4 output.
save_sci_plot(
  plot_object = p_figure_3_4,
  output_dir = thesis_plot_dir,
  file_stem = "Figure_3_4_VIP_pathway_compact_v3",
  width_mm = 180,
  height_mm = 215
)

message("Figure 3.4 v3: Tier legend separated; pathway palette updated.")
