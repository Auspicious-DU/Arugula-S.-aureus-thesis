# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# ============================================================
# plot_theme_prism_like.R
# General SCI figure style matched to Prism baseline
# ============================================================

# Required packages
library(ggplot2)
library(scales)
library(grid)

# ============================================================
# 1. FONT
# ============================================================

# Windows 注册 Arial
# 你的电脑目前需要这一句，避免 Arial font warning
if (.Platform$OS.type == "windows") {
  windowsFonts(
    Arial = windowsFont("Arial")
  )
}

SCI_FONT <- "Arial"

SCI_AXIS_TITLE_SIZE  <- 8
SCI_AXIS_TEXT_SIZE   <- 7
SCI_LEGEND_TEXT_SIZE <- 7


# ============================================================
# 2. LINE / SYMBOL SETTINGS
# ============================================================

# Prism 使用 pt；
# ggplot2 的 linewidth 更接近 mm，
# 所以统一先转换。

pt_to_mm <- function(x) {
  x * 25.4 / 72
}

# Prism baseline
SCI_AXIS_LINEWIDTH  <- pt_to_mm(0.5)
SCI_CURVE_LINEWIDTH <- pt_to_mm(1.0)
SCI_ERROR_LINEWIDTH <- pt_to_mm(0.5)
SCI_SYMBOL_BORDER   <- pt_to_mm(0.5)

# Prism Symbol size = 2
# 与 ggplot2 size 不是完全相同单位。
# 当前经过视觉比较，先固定为 1.0。
SCI_POINT_SIZE <- 1.0

# Prism transparency = 49%
# R alpha = opacity = 51%
SCI_ALPHA <- 0.51

# Prism standard tick length 的近似值
SCI_TICK_LENGTH <- grid::unit(
  3,
  "pt"
)


# ============================================================
# 3. PRISM COLOUR PALETTE
# ============================================================

# Light colours
# Purple RGB 164,138,211
# Blue   RGB 28,197,254
# Teal   RGB 111,199,207
# Orange RGB 251,162,125
# Red    RGB 251,125,128

PRISM_LIGHT <- c(
  Purple = "#A48AD3",
  Blue   = "#1CC5FE",
  Teal   = "#6FC7CF",
  Orange = "#FBA27D",
  Red    = "#FB7D80"
)


# Dark colours
# Purple RGB 44,20,83
# Blue   RGB 17,76,232
# Teal   RGB 14,111,124
# Orange RGB 251,79,6
# Red    RGB 251,0,5

PRISM_DARK <- c(
  Purple = "#2C1453",
  Blue   = "#114CE8",
  Teal   = "#0E6F7C",
  Orange = "#FB4F06",
  Red    = "#FB0005"
)


# Prism Dataset 1–10 的主颜色：
# Dataset 1–5 = Light
# Dataset 6–10 = Dark

PRISM_MAIN_10 <- c(
  unname(PRISM_LIGHT),
  unname(PRISM_DARK)
)

names(PRISM_MAIN_10) <-
  paste0(
    "Dataset",
    1:10
  )


# Prism Dataset 1–10 的边框 / error bar 颜色：
# Dataset 1–5 = Dark
# Dataset 6–10 = Light

PRISM_EDGE_10 <- c(
  unname(PRISM_DARK),
  unname(PRISM_LIGHT)
)

names(PRISM_EDGE_10) <-
  paste0(
    "Dataset",
    1:10
  )

# ============================================================
# 3B. NEUTRAL GREY COLOURS
#
# Reusable neutral colours for non-dataset elements:
# - non-highlighted / background categories
# - reference lines
# - threshold lines
# - zero / identity lines
#
# These colours do NOT belong to the Prism Dataset 1-10
# palette and should not replace dataset-specific colours.
# ============================================================

SCI_GREY_DARK  <- "grey35"
SCI_GREY_MID   <- "grey55"
SCI_GREY_LIGHT <- "grey75"
SCI_GREY_PALE  <- "grey88"


# ------------------------------------------------------------
# Semantic aliases
# ------------------------------------------------------------

# Neutral data categories
SCI_COLOUR_OTHER_FILL <- SCI_GREY_LIGHT
SCI_COLOUR_OTHER_EDGE <- SCI_GREY_MID

# Reference / guide lines
SCI_COLOUR_REF_LINE   <- SCI_GREY_MID
SCI_COLOUR_GUIDE_LINE <- SCI_GREY_LIGHT

# Common analytical reference lines
SCI_COLOUR_ZERO_LINE  <- SCI_GREY_MID
SCI_COLOUR_THRESHOLD  <- SCI_GREY_MID
SCI_COLOUR_IDENTITY   <- SCI_GREY_MID
# ============================================================
# 4. HELPER: GET PRISM DATASET COLOURS
# ============================================================

prism_dataset_colour <- function(
    dataset_number,
    role = c(
      "main",
      "edge"
    )
) {
  
  role <- match.arg(role)
  
  if (
    any(
      !dataset_number %in% 1:10
    )
  ) {
    stop(
      "dataset_number must be between 1 and 10."
    )
  }
  
  palette <- if (
    role == "main"
  ) {
    PRISM_MAIN_10
  } else {
    PRISM_EDGE_10
  }
  
  unname(
    palette[
      dataset_number
    ]
  )
}


# ============================================================
# 5. PRISM-LIKE SCI THEME
# ============================================================

theme_sci_prism <- function(
    legend_position = "right",
    plot_margin_pt = 4
) {
  
  ggplot2::theme_classic(
    base_family = SCI_FONT,
    base_size = SCI_AXIS_TEXT_SIZE
  ) +
    
    ggplot2::theme(
      
      # ---------- Global text ----------
      text = ggplot2::element_text(
        family = SCI_FONT,
        face = "plain",
        colour = "black"
      ),
      
      
      # ---------- Axis titles ----------
      axis.title = ggplot2::element_text(
        family = SCI_FONT,
        face = "plain",
        size = SCI_AXIS_TITLE_SIZE,
        colour = "black"
      ),
      
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(
          t = 4,
          unit = "pt"
        )
      ),
      
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(
          r = 4,
          unit = "pt"
        )
      ),
      
      
      # ---------- Tick labels ----------
      axis.text = ggplot2::element_text(
        family = SCI_FONT,
        face = "plain",
        size = SCI_AXIS_TEXT_SIZE,
        colour = "black"
      ),
      
      
      # ---------- Axes ----------
      axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = SCI_AXIS_LINEWIDTH
      ),
      
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = SCI_AXIS_LINEWIDTH
      ),
      
      axis.ticks.length =
        SCI_TICK_LENGTH,
      
      
      # ---------- No grid ----------
      panel.grid =
        ggplot2::element_blank(),
      
      
      # ---------- Background ----------
      panel.background =
        ggplot2::element_blank(),
      
      plot.background =
        ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
      
      # ---------- Facet strips ----------
      strip.background =
        ggplot2::element_blank(),
      
      strip.text =
        ggplot2::element_text(
          family = SCI_FONT,
          face = "plain",
          size = SCI_AXIS_TEXT_SIZE,
          colour = "black"
        ),
      
      # ---------- Legend ----------
      legend.position =
        legend_position,
      
      legend.title =
        ggplot2::element_blank(),
      
      legend.text =
        ggplot2::element_text(
          family = SCI_FONT,
          face = "plain",
          size = SCI_LEGEND_TEXT_SIZE,
          colour = "black"
        ),
      
      legend.background =
        ggplot2::element_blank(),
      
      legend.key =
        ggplot2::element_blank(),
      
      legend.key.width =
        grid::unit(
          14,
          "pt"
        ),
      
      legend.key.height =
        grid::unit(
          8,
          "pt"
        ),
      
      
      # ---------- No internal figure title ----------
      plot.title =
        ggplot2::element_blank(),
      
      plot.subtitle =
        ggplot2::element_blank(),
      
      
      # ---------- Outer margin ----------
      plot.margin =
        ggplot2::margin(
          t = plot_margin_pt,
          r = plot_margin_pt,
          b = plot_margin_pt,
          l = plot_margin_pt,
          unit = "pt"
        )
    )
}


# ============================================================
# 6. COMMON SHAPES
# ============================================================

# 21 circle
# 22 square
# 24 triangle up
# 25 triangle down
# 23 diamond
#
# Dataset 6–10 repeats the same shape series.

PRISM_SHAPES_10 <- c(
  21, 22, 24, 25, 23,
  21, 22, 24, 25, 23
)


# ============================================================
# 7. EXPORT FUNCTION
# ============================================================
#
# 一次调用即可输出：
# 1. Vector PDF —— AI / manuscript master
# 2. 600-dpi PNG —— PPT / preview / raster backup
#
# 注意：
# 根据你的电脑设置，这里显式使用 "\\" 拼接路径，
# 不使用 file.path()。
# ============================================================

save_sci_plot <- function(
    plot_object,
    output_dir,
    file_stem,
    width_mm,
    height_mm,
    save_pdf = TRUE,
    save_png = TRUE,
    png_dpi = 600
) {
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (save_pdf) {
    
    pdf_path <- file.path(output_dir, paste0(file_stem, ".pdf"))
    
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot_object,
      width = width_mm,
      height = height_mm,
      units = "mm",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
    
    message(
      "PDF saved to: ",
      pdf_path
    )
  }
  
  
  if (save_png) {
    
    png_path <- file.path(output_dir, paste0(file_stem, ".png"))
    
    ggplot2::ggsave(
      filename = png_path,
      plot = plot_object,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = png_dpi,
      bg = "white"
    )
    
    message(
      "PNG saved to: ",
      png_path
    )
  }
  
  invisible(NULL)
}