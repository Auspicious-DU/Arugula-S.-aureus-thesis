# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# ============================================================
# THESIS_phenotype_figures.R
# Figure 3.2: phenotype variation in Eruca sativa samples
#
# A: AUC ratio distribution
# B: Delta-T0.10 distribution (12-h capped)
# C: AUC ratio versus Delta-T0.10 (individual level)
#
# IMPORTANT
# - Use the harmonised phenotype table directly (589 individuals).
# - Delta-T capping was already performed upstream; do not recalculate.
# - Only the neutral references AUC ratio = 1 and Delta-T = 0
#   are displayed. No biological-group classification is applied.
# - Do not modify the underlying source data.
# ============================================================

# ------------------------------------------------------------
# 0. Project paths and SCI plotting theme
# ------------------------------------------------------------

DT_project_dir <- repo_dir

DT_phenotype_file <- file.path(repo_dir, "output", "delta_t", "00_preprocessed_data", "DT_phenotype_12h_clean.csv")

# Confirmed location in the uploaded 00_DT_setup.R
SCI_theme_file <- file.path(repo_dir, "scripts", "figures", "plot_theme_prism_like_v2.R")

THESIS_output_dir <- file.path(repo_dir, "output", "thesis_figures")

THESIS_plot_dir <- file.path(THESIS_output_dir, "plots")
THESIS_table_dir <- file.path(THESIS_output_dir, "tables")

if (!file.exists(SCI_theme_file)) {
  stop("SCI theme not found: ", SCI_theme_file)
}
source(SCI_theme_file)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(patchwork)
})

dir.create(THESIS_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(THESIS_table_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read and validate the single authoritative phenotype file
# ------------------------------------------------------------

if (!file.exists(DT_phenotype_file)) {
  stop("Phenotype input not found: ", DT_phenotype_file)
}

phenotype_raw <- readr::read_csv(
  DT_phenotype_file,
  show_col_types = FALSE
)

required_columns <- c(
  "sample_id_clean",
  "auc_ratio_mean",
  "dt0_10_12h",
  "dt0_10_12h_status"
)

missing_columns <- setdiff(required_columns, names(phenotype_raw))
if (length(missing_columns) > 0) {
  stop("Required columns missing: ", paste(missing_columns, collapse = ", "))
}

if (anyDuplicated(phenotype_raw$sample_id_clean) > 0) {
  stop("Duplicated sample_id_clean values: check the phenotype source.")
}

if (!is.numeric(phenotype_raw$auc_ratio_mean) ||
    !is.numeric(phenotype_raw$dt0_10_12h)) {
  stop("AUC ratio or Delta-T0.10 is not numeric.")
}

phenotype_plot <- phenotype_raw %>%
  dplyr::select(
    sample_id_clean,
    dplyr::any_of("accession"),
    auc_ratio_mean,
    dt0_10_12h,
    dt0_10_12h_status
  ) %>%
  dplyr::filter(
    !is.na(auc_ratio_mean),
    !is.na(dt0_10_12h),
    is.finite(auc_ratio_mean),
    is.finite(dt0_10_12h)
  )

if (nrow(phenotype_plot) < 3) {
  stop("Fewer than three complete phenotype observations.")
}

message("Complete individual-level phenotype observations: ", nrow(phenotype_plot))
if (nrow(phenotype_plot) != 589L) {
  warning("Expected 589 rows in the uploaded dataset; found ", nrow(phenotype_plot))
}

readr::write_csv(
  phenotype_plot,
  file.path(THESIS_table_dir, "Figure_3_2_plotting_data.csv")
)

# ------------------------------------------------------------
# 2. Summaries for writing Results (not placed on graph)
# ------------------------------------------------------------

summarise_trait <- function(x, trait_name) {
  tibble::tibble(
    trait = trait_name,
    n = length(x),
    mean = mean(x),
    sd = stats::sd(x),
    min = min(x),
    q1 = unname(stats::quantile(x, 0.25)),
    median = stats::median(x),
    q3 = unname(stats::quantile(x, 0.75)),
    max = max(x),
    n_below_reference = sum(x < if (trait_name == "AUC ratio") 1 else 0),
    n_equal_reference = sum(x == if (trait_name == "AUC ratio") 1 else 0),
    n_above_reference = sum(x > if (trait_name == "AUC ratio") 1 else 0)
  )
}

phenotype_summary <- dplyr::bind_rows(
  summarise_trait(phenotype_plot$auc_ratio_mean, "AUC ratio"),
  summarise_trait(phenotype_plot$dt0_10_12h, "Delta-T0.10 (h), 12-h capped")
)

status_summary <- phenotype_plot %>%
  dplyr::count(dt0_10_12h_status, name = "n") %>%
  dplyr::mutate(percent = 100 * n / sum(n))

readr::write_csv(
  phenotype_summary,
  file.path(THESIS_table_dir, "Figure_3_2_phenotype_summary.csv")
)
readr::write_csv(
  status_summary,
  file.path(THESIS_table_dir, "Figure_3_2_DT_capping_status.csv")
)

# ------------------------------------------------------------
# 3. Correlation: Pearson and Spearman, full individual-level data
# ------------------------------------------------------------

pearson_test <- stats::cor.test(
  phenotype_plot$auc_ratio_mean,
  phenotype_plot$dt0_10_12h,
  method = "pearson"
)

# exact = FALSE avoids the exact-p issue caused by tied Delta-T values.
spearman_test <- suppressWarnings(stats::cor.test(
  phenotype_plot$auc_ratio_mean,
  phenotype_plot$dt0_10_12h,
  method = "spearman",
  exact = FALSE
))

correlation_summary <- tibble::tibble(
  method = c("Pearson", "Spearman"),
  n = nrow(phenotype_plot),
  estimate = c(unname(pearson_test$estimate), unname(spearman_test$estimate)),
  p_value = c(pearson_test$p.value, spearman_test$p.value)
)

readr::write_csv(
  correlation_summary,
  file.path(THESIS_table_dir, "Figure_3_2_correlation_summary.csv")
)

# Display Pearson on the scatter for now; report both in CSV.
# This is a descriptive association between two phenotype metrics.
pearson_r_label <- sprintf("Pearson r = %.3f", unname(pearson_test$estimate))
pearson_p_label <- paste0("p ", if (pearson_test$p.value < 0.001) {
  "< 0.001"
} else {
  paste0("= ", formatC(pearson_test$p.value, digits = 3, format = "f"))
})
cor_annotation <- paste(pearson_r_label, pearson_p_label, sep = "\n")

# ------------------------------------------------------------
# 4. Colour identities and plot-specific settings
# ------------------------------------------------------------

AUC_FILL <- prism_dataset_colour(1, role = "main")   # Purple
AUC_EDGE <- prism_dataset_colour(1, role = "edge")
DT_FILL  <- prism_dataset_colour(2, role = "main")   # Blue
DT_EDGE  <- prism_dataset_colour(2, role = "edge")

# Match the previous phenotype histogram bin widths, but remove
# biological grouping and thresholds. No data are clipped.
AUC_BINWIDTH <- 0.05
DT_BINWIDTH <- 0.25

n_label <- paste0("n = ", nrow(phenotype_plot))

# ------------------------------------------------------------
# 5. Panel A: AUC-ratio distribution (count histogram)
# ------------------------------------------------------------

p_auc <- ggplot2::ggplot(
  phenotype_plot,
  ggplot2::aes(x = auc_ratio_mean)
) +
  ggplot2::geom_histogram(
    binwidth = AUC_BINWIDTH,
    boundary = 0,
    fill = AUC_FILL,
    colour = AUC_EDGE,
    alpha = SCI_ALPHA,
    linewidth = SCI_SYMBOL_BORDER
  ) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = SCI_COLOUR_REF_LINE,
    linewidth = SCI_AXIS_LINEWIDTH
  ) +
  ggplot2::annotate(
    "text", x = Inf, y = Inf,
    label = n_label,
    hjust = 1.08, vjust = 1.3,
    family = SCI_FONT,
    size = pt_to_mm(SCI_AXIS_TEXT_SIZE)
  ) +
  ggplot2::scale_x_continuous(breaks = seq(0, 2, by = 0.5)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::labs(
    x = "AUC ratio",
    y = "Number of samples"
  ) +
  theme_sci_prism(legend_position = "none")

# ------------------------------------------------------------
# 6. Panel B: Delta-T0.10 distribution (12-h capped phenotype)
# ------------------------------------------------------------

p_dt <- ggplot2::ggplot(
  phenotype_plot,
  ggplot2::aes(x = dt0_10_12h)
) +
  ggplot2::geom_histogram(
    binwidth = DT_BINWIDTH,
    boundary = 0,
    fill = DT_FILL,
    colour = DT_EDGE,
    alpha = SCI_ALPHA,
    linewidth = SCI_SYMBOL_BORDER
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = SCI_COLOUR_ZERO_LINE,
    linewidth = SCI_AXIS_LINEWIDTH
  ) +
  ggplot2::annotate(
    "text", x = Inf, y = Inf,
    label = n_label,
    hjust = 1.08, vjust = 1.3,
    family = SCI_FONT,
    size = pt_to_mm(SCI_AXIS_TEXT_SIZE)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::labs(
    x = expression(Delta*T[0.10]~"(h)"),
    y = "Number of samples"
  ) +
  theme_sci_prism(legend_position = "none")

# ------------------------------------------------------------
# 7. Panel C: AUC ratio versus Delta-T0.10
# ------------------------------------------------------------

p_cor <- ggplot2::ggplot(
  phenotype_plot,
  ggplot2::aes(x = auc_ratio_mean, y = dt0_10_12h)
) +
  # A simple linear trend for visual orientation, not a causal model.
  geom_point(
    shape = 21,
    size = SCI_POINT_SIZE,
    stroke = SCI_SYMBOL_BORDER,
    fill = "grey65",
    colour = "grey55",
    alpha = 0.65
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = SCI_CURVE_LINEWIDTH,
    colour = "black"
  ) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = SCI_COLOUR_REF_LINE,
    linewidth = SCI_AXIS_LINEWIDTH
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = SCI_COLOUR_ZERO_LINE,
    linewidth = SCI_AXIS_LINEWIDTH
  ) +
  ggplot2::annotate(
    "text", x = Inf, y = Inf,
    label = cor_annotation,
    hjust = 1.08, vjust = 1.2,
    family = SCI_FONT,
    size = pt_to_mm(SCI_AXIS_TEXT_SIZE),
    lineheight = 1.1
  ) +
  ggplot2::labs(
    x = "AUC ratio",
    y = expression(Delta*T[0.10]~"(h)")
  ) +
  theme_sci_prism(legend_position = "none")

# ------------------------------------------------------------
# 8. Assemble Figure 3.2
# ------------------------------------------------------------

# A and B on the first row, C spans the width of the second row.
# This arrangement remains readable at ~180-mm overall width.

p_combined <- (p_auc | p_dt) / p_cor +
  patchwork::plot_layout(heights = c(1, 1.1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = SCI_FONT,
      face = "bold",
      size = 9
    ),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_auc)
print(p_dt)
print(p_cor)
print(p_combined)

# ------------------------------------------------------------
# 9. Export PDF (vector) + PNG (600 dpi), using existing helper
# ------------------------------------------------------------

save_sci_plot(
  plot_object = p_auc,
  output_dir = THESIS_plot_dir,
  file_stem = "Figure_3_2A_AUC_ratio_distribution",
  width_mm = 89,
  height_mm = 68
)

save_sci_plot(
  plot_object = p_dt,
  output_dir = THESIS_plot_dir,
  file_stem = "Figure_3_2B_DT010_12h_distribution",
  width_mm = 89,
  height_mm = 68
)

save_sci_plot(
  plot_object = p_cor,
  output_dir = THESIS_plot_dir,
  file_stem = "Figure_3_2C_AUC_DT010_relationship",
  width_mm = 180,
  height_mm = 75
)

save_sci_plot(
  plot_object = p_combined,
  output_dir = THESIS_plot_dir,
  file_stem = "Figure_3_2_Phenotypic_variation_combined",
  width_mm = 180,
  height_mm = 150
)

message("THESIS Figure 3.2 outputs saved to: ", THESIS_plot_dir)
message("Summary tables saved to: ", THESIS_table_dir)
message("Dataset: 12-h harmonised phenotype; no re-capping or grouping.")
