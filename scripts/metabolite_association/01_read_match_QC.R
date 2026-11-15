# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}

# =========================
# 01_read_match_qc.R
# =========================
# =========================
# 4. Read files
# =========================
metab_raw <- read_excel(metab_file, sheet = 1) %>%
  clean_names()

anno_raw <- read_excel(anno_file, sheet = 1) %>%
  clean_names()

pheno_raw <- read_excel(pheno_file, sheet = 1) %>%
  clean_names()

# =========================
# 5. Check column names
# =========================
names(metab_raw)
names(anno_raw)
names(pheno_raw)

head(metab_raw)
head(anno_raw)
head(pheno_raw)

dim(metab_raw)
dim(anno_raw)
dim(pheno_raw)

# =========================
# 6. Inspect sample ID format
# =========================
metab_raw %>%
  dplyr::select(sample_id) %>%
  slice(1:20) %>%
  print(n = 20)

pheno_raw %>%
  dplyr::select(sample_id) %>%
  slice(1:20) %>%
  print(n = 20)

# =========================
# 7. Clean sample IDs，统一成小写的a和.
# =========================

clean_sample_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("：", ":") %>%          # handle Chinese colon if any
    str_replace_all("[:/\\\\]", ".") %>%   # replace :, /, \ with .
    str_replace_all("\\s+", "") %>%        # remove spaces
    str_replace_all("\\.+", ".") %>%       # collapse multiple dots
    str_replace_all("^\\.|\\.$", "")       # remove leading/trailing dots
}

metab_id_check <- metab_raw %>%
  mutate(sample_id_clean = clean_sample_id(sample_id)) %>%
  dplyr::select(sample_id, sample_id_clean)

pheno_id_check <- pheno_raw %>%
  mutate(sample_id_clean = clean_sample_id(sample_id)) %>%
  dplyr::select(sample_id, sample_id_clean)

metab_id_check %>%
  slice(1:20) %>%
  print(n = 20)

pheno_id_check %>%
  slice(1:20) %>%
  print(n = 20)

# =========================
# 8. Check duplicated sample IDs
# =========================

metab_dup <- metab_id_check %>%
  count(sample_id_clean) %>%
  filter(n > 1)

pheno_dup <- pheno_id_check %>%
  count(sample_id_clean) %>%
  filter(n > 1)

metab_dup
pheno_dup

# =========================
# 9. Sample matching summary
# =========================

metab_ids <- metab_id_check %>%
  distinct(sample_id_clean)

pheno_ids <- pheno_id_check %>%
  distinct(sample_id_clean)

matched_ids <- inner_join(metab_ids, pheno_ids, by = "sample_id_clean")

metab_only <- anti_join(metab_ids, pheno_ids, by = "sample_id_clean")

pheno_only <- anti_join(pheno_ids, metab_ids, by = "sample_id_clean")

matching_summary <- tibble(
  category = c(
    "Metabolite samples",
    "Phenotype samples",
    "Matched samples",
    "Metabolite only",
    "Phenotype only"
  ),
  n = c(
    nrow(metab_ids),
    nrow(pheno_ids),
    nrow(matched_ids),
    nrow(metab_only),
    nrow(pheno_only)
  )
)

matching_summary

# =========================
# 10. Save ID check outputs
# =========================

dir.create(file.path(output_dir, "01_sample_matching"), showWarnings = FALSE)

write_csv(matching_summary,
          file.path(output_dir, "01_sample_matching", "sample_matching_summary.csv"))

write_csv(matched_ids,
          file.path(output_dir, "01_sample_matching", "matched_sample_ids.csv"))

write_csv(metab_only,
          file.path(output_dir, "01_sample_matching", "metabolite_only_sample_ids.csv"))

write_csv(pheno_only,
          file.path(output_dir, "01_sample_matching", "phenotype_only_sample_ids.csv"))

# =========================
# 11. Inspect unmatched phenotype sample
# =========================

pheno_id_check %>%
  filter(sample_id_clean %in% pheno_only$sample_id_clean)

pheno_raw %>%
  mutate(sample_id_clean = clean_sample_id(sample_id)) %>%
  filter(sample_id_clean %in% pheno_only$sample_id_clean)

# =========================
# 12. Check peak matching between matrix and annotation
# =========================

peak_cols <- names(metab_raw) %>%
  keep(~ str_detect(.x, "^peak"))

anno_peak_ids <- anno_raw %>%
  mutate(peak_id_clean = str_to_lower(peak_id)) %>%
  pull(peak_id_clean)

peak_matching_summary <- tibble(
  category = c(
    "Peak columns in metabolite matrix",
    "Peak IDs in annotation",
    "Matched peaks",
    "Matrix only peaks",
    "Annotation only peaks"
  ),
  n = c(
    length(peak_cols),
    length(anno_peak_ids),
    length(intersect(peak_cols, anno_peak_ids)),
    length(setdiff(peak_cols, anno_peak_ids)),
    length(setdiff(anno_peak_ids, peak_cols))
  )
)

peak_matching_summary

# =========================
# 13. Create cleaned data tables
# =========================

metab_clean <- metab_raw %>%
  mutate(sample_id_clean = clean_sample_id(sample_id)) %>%
  rename(
    metab_sample_id_original = sample_id,
    metab_group = group
  ) %>%
  dplyr::select(
    sample_id_clean,
    meta_id,
    metab_group,
    metab_sample_id_original,
    all_of(peak_cols)
  )

pheno_clean <- pheno_raw %>%
  mutate(sample_id_clean = clean_sample_id(sample_id)) %>%
  rename(
    pheno_sample_id_original = sample_id
  )

anno_clean <- anno_raw %>%
  mutate(peak_id_clean = str_to_lower(peak_id)) %>%
  dplyr::select(peak_id_clean, everything())

# =========================
# 14. Inner join phenotype and metabolite data
# =========================

matched_data <- pheno_clean %>%
  inner_join(metab_clean, by = "sample_id_clean")

dim(matched_data)

# =========================
# 15. Save matched dataset
# =========================

dir.create(file.path(output_dir, "02_matched_data"), showWarnings = FALSE)

write_csv(
  matched_data,
  file.path(output_dir, "02_matched_data", "matched_pheno_metab_data.csv")
)

write_csv(
  anno_clean,
  file.path(output_dir, "02_matched_data", "peak_annotation_clean.csv")
)

write_csv(
  peak_matching_summary,
  file.path(output_dir, "02_matched_data", "peak_matching_summary.csv")
)