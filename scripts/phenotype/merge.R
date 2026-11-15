# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# =========================
# 1. 路径
# =========================
project_dir <- file.path(repo_dir, "input", "phenotype")
output_dir <- file.path(repo_dir, "output", "phenotype")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_clean_dir <- file.path(project_dir, "data_clean")

# 如果你已经把输出整理到子文件夹里，请改成对应路径
output_2rep_dir <- file.path(output_dir, "2_replicate")
output_3rep_dir <- file.path(output_dir, "3_replicate")


# =========================
# 2. 读取 2rep 宽表 / 长表
# =========================
summary_2rep <- read.csv(
  file.path(output_2rep_dir, "summary_with_thresholds_all.csv"),
  stringsAsFactors = FALSE
)

long_2rep <- read.csv(
  file.path(output_2rep_dir, "threshold_summary_long_all.csv"),
  stringsAsFactors = FALSE
)

# =========================
# 3. 读取 3rep 宽表 / 长表
# =========================
summary_3rep <- read.csv(
  file.path(output_3rep_dir, "summary_with_thresholds_all_3rep.csv"),
  stringsAsFactors = FALSE
)

long_3rep <- read.csv(
  file.path(output_3rep_dir, "threshold_summary_long_all_3rep.csv"),
  stringsAsFactors = FALSE
)

# =========================
# 4. 加来源标签
# =========================
summary_2rep <- summary_2rep %>%
  mutate(ReplicateType = "2rep")

summary_3rep <- summary_3rep %>%
  mutate(ReplicateType = "3rep")

long_2rep <- long_2rep %>%
  mutate(ReplicateType = "2rep")

long_3rep <- long_3rep %>%
  mutate(ReplicateType = "3rep")

# =========================
# 5. 合并宽表
# =========================
rank1_wide <- bind_rows(summary_2rep, summary_3rep)

# 可选：按日期/样本排序
rank1_wide <- rank1_wide %>%
  arrange(TestDate, SampleID)

# =========================
# 6. 合并长表
#    先把两批 threshold long 拼起来
# =========================
rank1_long_raw <- bind_rows(long_2rep, long_3rep)

# 只保留 sample，不保留 MeOH（如果你画 sample phenotype 图更方便）
rank1_long_raw <- rank1_long_raw %>%
  filter(SampleID != "MeOH")

# =========================
# 7. 给长表补上 Accession / AUC / ratio 等信息
#    从宽表中抽取基础信息再 join
# =========================
meta_cols <- rank1_wide %>%
  select(
    FileName, TestDate, AssayBatch, BatchID, SampleID, ReplicateType,
    Accession, AUC_mean, AUC_sd, AUC_ratio_mean, AUC_ratio_sd
  ) %>%
  distinct()

rank1_long <- rank1_long_raw %>%
  left_join(
    meta_cols,
    by = c("FileName", "TestDate", "AssayBatch", "BatchID", "SampleID", "ReplicateType")
  ) %>%
  arrange(TestDate, SampleID, threshold_name)

# =========================
# 8. 可选：再做一个“超长表”
#    把 wide 里的 threshold 宽列也展开成长表（备用）
# =========================
rank1_wide_to_long <- rank1_wide %>%
  pivot_longer(
    cols = matches("^T0\\.(05|1|15)_(T_mean|T_sd|dT)$"),
    names_to = c("threshold_name", ".value"),
    names_pattern = "(T0\\.(?:05|1|15))_(T_mean|T_sd|dT)"
  ) %>%
  arrange(TestDate, SampleID, threshold_name)

# =========================
# 9. 保存最终表
# =========================
write.csv(
  rank1_wide,
  file.path(output_dir, "Rank1_summary_with_thresholds_final_wide.csv"),
  row.names = FALSE
)

write.csv(
  rank1_long,
  file.path(output_dir, "Rank1_thresholds_final_long.csv"),
  row.names = FALSE
)

write.csv(
  rank1_wide_to_long,
  file.path(output_dir, "Rank1_summary_with_thresholds_final_long_from_wide.csv"),
  row.names = FALSE
)

# =========================
# 10. 快速检查
# =========================
dim(rank1_wide)
dim(rank1_long)

rank1_wide %>% count(ReplicateType)
rank1_long %>% count(ReplicateType, threshold_name)

head(rank1_wide)
head(rank1_long)
