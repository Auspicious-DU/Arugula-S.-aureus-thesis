# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

library(dplyr)
library(readr)

result_wide_dedup_earliest <- result_wide %>%
  mutate(
    TestDate = as.Date(TestDate),
    SampleID = trimws(as.character(SampleID))
  ) %>%
  arrange(SampleID, TestDate) %>%
  group_by(SampleID) %>%
  slice(1) %>%
  ungroup()

removed_duplicates <- result_wide %>%
  mutate(
    TestDate = as.Date(TestDate),
    SampleID = trimws(as.character(SampleID))
  ) %>%
  arrange(SampleID, TestDate) %>%
  group_by(SampleID) %>%
  slice(-1) %>%
  ungroup()

write_csv(
  result_wide_dedup_earliest,
  file.path(repo_dir, "output", "phenotype", "Rank1_summary_with_thresholds_final_wide_dedup_earliest.csv")
)

write_csv(
  removed_duplicates,
  file.path(repo_dir, "output", "phenotype", "Rank1_summary_removed_duplicates.csv")
)

cat("\n================ QC SUMMARY: Duplicate Handling ================\n")
cat("Input rows before deduplication :", nrow(result_wide), "\n")
cat("Rows after deduplication        :", nrow(result_wide_dedup_earliest), "\n")
cat("Removed duplicate rows          :", nrow(removed_duplicates), "\n")
cat("Rule applied                    : keep the earliest TestDate for each SampleID\n")
cat("Rationale                       : later V1 batch showed abnormal values\n")
cat("================================================================\n\n")