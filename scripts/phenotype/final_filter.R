# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

library(dplyr)
library(stringr)
library(readr)
library(readxl)

# =========================
# 0. 路径
# =========================
validation_dir <- file.path(repo_dir, "input", "phenotype", "validation list")
output_dir <- file.path(repo_dir, "output", "phenotype")

rank1_file   <- file.path(validation_dir, "Arugula LCMS+Seqence samples list_22122025_for Yuhan.xlsx")
tracker_file <- file.path(validation_dir, "Sample Status list.xlsx")

result_file  <- file.path(output_dir, "Rank1_summary_with_thresholds_final_wide_dedup_earliest.csv")

# =========================
# 1. 小工具
# =========================
clean_text <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "Na", "N/A", "NULL")] <- NA
  x
}

# =========================
# 2. 读数据
# =========================
rank1_raw   <- read_excel(rank1_file, sheet = "Sheet1")
tracker_raw <- read_excel(tracker_file, sheet = "Samples")
result_wide <- read_csv(result_file, show_col_types = FALSE)

# =========================
# 3. 整理 Rank1 和 tracker
# =========================
rank1_ids <- rank1_raw %>%
  transmute(
    SampleID = clean_text(ID),
    Accession_rank1 = clean_text(Accession)
  ) %>%
  filter(!is.na(SampleID)) %>%
  distinct()

tracker_ids <- tracker_raw %>%
  transmute(
    SampleID = clean_text(`Sample ID`),
    Tracker_Status = clean_text(Status),
    Tracker_TestDate = as.character(`Test Date`),
    Tracker_AssayBatch = clean_text(`Assay Batch`)
  ) %>%
  filter(!is.na(SampleID)) %>%
  distinct()

# =========================
# 4. 给 result 加上 Rank1 / tracker 信息
# =========================
result_annotated <- result_wide %>%
  mutate(SampleID = clean_text(SampleID)) %>%
  left_join(rank1_ids, by = "SampleID") %>%
  left_join(tracker_ids, by = "SampleID") %>%
  mutate(
    In_Rank1 = !is.na(Accession_rank1),
    Tracker_Tested = Tracker_Status == "Tested"
  )

# =========================
# 5. 最终 clean 表
#    只保留：Rank1 + Tested
# =========================
final_clean_summary <- result_annotated %>%
  filter(In_Rank1, Tracker_Tested)

# 被过滤掉的记录也保留一份
filtered_out_records <- result_annotated %>%
  filter(!(In_Rank1 & Tracker_Tested)) %>%
  mutate(
    FilterReason = case_when(
      !In_Rank1 ~ "Not in Rank1 list",
      In_Rank1 & is.na(Tracker_Status) ~ "No tracker status",
      In_Rank1 & Tracker_Status == "Out of stock" ~ "Tracker status = Out of stock",
      In_Rank1 & Tracker_Status == "Exclude" ~ "Tracker status = Exclude",
      TRUE ~ "Other"
    )
  )

# =========================
# 6. 简单汇总
# =========================
filter_summary <- tibble(
  Category = c(
    "Input rows (dedup result)",
    "Final clean rows (Rank1 + Tested)",
    "Filtered out rows"
  ),
  N = c(
    nrow(result_wide),
    nrow(final_clean_summary),
    nrow(filtered_out_records)
  )
)

filtered_out_by_reason <- filtered_out_records %>%
  count(FilterReason, sort = TRUE)

# =========================
# 7. 输出
# =========================
write_csv(
  final_clean_summary,
  file.path(output_dir, "Rank1_summary_with_thresholds_final_wide_clean.csv")
)

write_csv(
  filtered_out_records,
  file.path(output_dir, "Rank1_summary_filtered_out_records.csv")
)

write_csv(
  filter_summary,
  file.path(output_dir, "Rank1_summary_filter_summary.csv")
)

write_csv(
  filtered_out_by_reason,
  file.path(output_dir, "Rank1_summary_filtered_out_by_reason.csv")
)

# =========================
# 8. 控制台查看
cat("\n================ QC SUMMARY: Final Clean Summary ================\n")
# =========================
print(filter_summary)
print(filtered_out_by_reason)
cat("================================================================\n\n")