# Run from the GitHub repository root in RStudio (Session > Set Working Directory > Choose Directory).
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_dir, "README.md"))) {
  stop("Set the R working directory to the repository root (folder containing README.md).")
}



# =========================
# 0. Packages
# =========================
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# =========================
# 1. Project paths
# =========================
project_dir <- file.path(repo_dir, "input", "phenotype")
raw_auto_dir <- file.path(project_dir, "raw_auto", "2_replicate")
data_clean_dir <- file.path(project_dir, "data_clean")
output_dir <- file.path(repo_dir, "output", "phenotype", "2_replicate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 2. Files
# =========================
files <- list.files(
  path = raw_auto_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

# =========================
# 3. Parse file info
# =========================
parse_file_info <- function(file_path) {
  file_name <- basename(file_path)
  
  date_raw <- str_extract(file_name, "\\d{2}\\.\\d{1,2}\\.\\d{1,2}")
  
  test_date <- NA
  if (!is.na(date_raw)) {
    parts <- str_split(date_raw, "\\.", simplify = TRUE)
    yy <- parts[1]
    mm <- parts[2]
    dd <- parts[3]
    test_date <- as.Date(sprintf("20%02d-%02d-%02d",
                                 as.integer(yy), as.integer(mm), as.integer(dd)))
  }
  
  assay_batch <- str_match(file_name, "Sample\\s+(.+)\\.xlsx$")[, 2]
  
  batch_id <- if (!is.na(test_date) && !is.na(assay_batch)) {
    paste0(assay_batch, "-", format(test_date, "%Y%m%d"))
  } else {
    NA_character_
  }
  
  tibble(
    FileName = file_name,
    TestDate = test_date,
    AssayBatch = assay_batch,
    BatchID = batch_id
  )
}

# =========================
# 4. Threshold function
# =========================
calc_threshold_time <- function(time_vec, od_vec, threshold = 0.1) {
  
  df <- data.frame(Time_h = time_vec, OD = od_vec) %>%
    filter(!is.na(Time_h), !is.na(OD)) %>%
    arrange(Time_h)
  
  if (nrow(df) == 0) return(NA_real_)
  if (all(df$OD < threshold)) return(NA_real_)
  
  idx <- which(df$OD >= threshold)[1]
  
  if (idx == 1) {
    return(df$Time_h[idx])
  }
  
  t1 <- df$Time_h[idx - 1]
  t2 <- df$Time_h[idx]
  y1 <- df$OD[idx - 1]
  y2 <- df$OD[idx]
  
  if (y2 == y1) return(t2)
  
  t_cross <- t1 + (threshold - y1) / (y2 - y1) * (t2 - t1)
  return(t_cross)
}

# =========================
# 5. Read one file: sheet1
# =========================
read_curve_sheet_one <- function(file_path) {
  
  raw_num <- read_excel(
    file_path,
    sheet = 1,
    range = "A61:BG100",
    col_names = FALSE
  )
  
  # 去掉全空列
  raw_num <- raw_num %>%
    select(where(~ !all(is.na(.))))
  
  # 临时列名
  colnames(raw_num) <- c("Time_h", paste0("V", 2:ncol(raw_num)))
  
  # 转数值
  raw_num <- raw_num %>%
    mutate(across(everything(), as.numeric))
  
  # 去掉空行
  raw_num <- raw_num %>%
    filter(!is.na(Time_h)) %>%
    filter(rowSums(!is.na(.)) > 1)
  
  raw_num
}

# =========================
# 6. Read one file: sheet2
# =========================
read_summary_sheet_one <- function(file_path) {
  
  auc_sheet <- read_excel(
    file_path,
    sheet = 2,
    col_names = FALSE
  )
  
  # 左侧 A:D
  auc_left <- auc_sheet %>%
    slice(3:n()) %>%
    select(1:4)
  
  colnames(auc_left) <- c("SampleID", "Accession", "AUC_mean", "AUC_sd")
  
  auc_left <- auc_left %>%
    filter(!is.na(SampleID)) %>%
    mutate(
      SampleID = as.character(SampleID),
      Accession = as.character(Accession),
      AUC_mean = as.numeric(AUC_mean),
      AUC_sd = as.numeric(AUC_sd)
    )
  
  # 右侧 F:I
  auc_right <- auc_sheet %>%
    slice(3:n()) %>%
    select(6:9)
  
  colnames(auc_right) <- c("SampleID", "Accession", "AUC_ratio_mean", "AUC_ratio_sd")
  
  auc_right <- auc_right %>%
    filter(!is.na(SampleID)) %>%
    mutate(
      SampleID = as.character(SampleID),
      Accession = as.character(Accession),
      AUC_ratio_mean = as.numeric(AUC_ratio_mean),
      AUC_ratio_sd = as.numeric(AUC_ratio_sd)
    )
  
  summary_tbl <- auc_left %>%
    left_join(auc_right, by = c("SampleID", "Accession"))
  
  summary_tbl
}

# =========================
# 7. Process one file
# =========================
process_one_file <- function(file_path, thresholds = c(0.05, 0.1, 0.15)) {
  
  info <- parse_file_info(file_path)
  file_name <- info$FileName[1]
  
  log_tbl <- tibble(
    FileName = file_name,
    TestDate = info$TestDate[1],
    AssayBatch = info$AssayBatch[1],
    BatchID = info$BatchID[1],
    status = "OK",
    message = NA_character_
  )
  
  # ---- read sheet1 ----
  raw_num <- tryCatch(
    read_curve_sheet_one(file_path),
    error = function(e) e
  )
  
  if (inherits(raw_num, "error")) {
    log_tbl$status <- "FAIL"
    log_tbl$message <- paste("sheet1 read error:", raw_num$message)
    return(list(
      curve_long = NULL,
      summary_tbl = NULL,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  # ---- read sheet2 ----
  summary_tbl <- tryCatch(
    read_summary_sheet_one(file_path),
    error = function(e) e
  )
  
  if (inherits(summary_tbl, "error")) {
    log_tbl$status <- "FAIL"
    log_tbl$message <- paste("sheet2 read error:", summary_tbl$message)
    return(list(
      curve_long = NULL,
      summary_tbl = NULL,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  # ---- 基本结构检查 ----
  if (ncol(raw_num) < 3) {
    log_tbl$status <- "SKIP"
    log_tbl$message <- "curve block has fewer than 3 columns"
    return(list(
      curve_long = NULL,
      summary_tbl = summary_tbl,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  # 后面样本列数必须是偶数（每个样本2列）
  sample_col_n <- ncol(raw_num) - 3
  
  if (sample_col_n < 0 || sample_col_n %% 2 != 0) {
    log_tbl$status <- "SKIP"
    log_tbl$message <- paste0("curve sample columns invalid: ncol(raw_num)=", ncol(raw_num))
    return(list(
      curve_long = NULL,
      summary_tbl = summary_tbl,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  sample_n_curve <- sample_col_n / 2
  sample_n_summary <- nrow(summary_tbl)
  
  if (sample_n_curve != sample_n_summary) {
    log_tbl$status <- "SKIP"
    log_tbl$message <- paste0(
      "sample count mismatch: curve=", sample_n_curve,
      ", summary=", sample_n_summary
    )
    return(list(
      curve_long = NULL,
      summary_tbl = summary_tbl,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  # ---- 命名 curve 列 ----
  sample_ids <- summary_tbl$SampleID
  
  new_names <- c(
    "Time_h",
    "MeOH__rep1",
    "MeOH__rep2",
    as.vector(rbind(
      paste0(sample_ids, "__rep1"),
      paste0(sample_ids, "__rep2")
    ))
  )
  
  if (length(new_names) != ncol(raw_num)) {
    log_tbl$status <- "SKIP"
    log_tbl$message <- "column naming length mismatch"
    return(list(
      curve_long = NULL,
      summary_tbl = summary_tbl,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  colnames(raw_num) <- new_names
  
  # ---- 转成长表 ----
  curve_long <- raw_num %>%
    pivot_longer(
      cols = -Time_h,
      names_to = "SampleCol",
      values_to = "OD"
    ) %>%
    mutate(
      SampleID = str_remove(SampleCol, "__rep\\d+$"),
      Replicate = as.integer(str_extract(SampleCol, "(?<=__rep)\\d+$")),
      Group = if_else(SampleID == "MeOH", "Control", "Sample")
    ) %>%
    select(Time_h, SampleID, Replicate, Group, OD) %>%
    mutate(
      FileName = info$FileName[1],
      TestDate = info$TestDate[1],
      AssayBatch = info$AssayBatch[1],
      BatchID = info$BatchID[1]
    ) %>%
    select(FileName, TestDate, AssayBatch, BatchID, everything())
  
  # ---- 时间点检查：放宽，不要求一定25个 ----
  time_check <- curve_long %>%
    group_by(SampleID, Replicate) %>%
    summarise(
      n_time = n(),
      min_time = min(Time_h, na.rm = TRUE),
      max_time = max(Time_h, na.rm = TRUE),
      is_increasing = all(diff(Time_h) > 0),
      .groups = "drop"
    )
  
  # 至少要有足够时间点，且时间递增
  if (any(time_check$n_time < 20) || any(!time_check$is_increasing)) {
    log_tbl$status <- "SKIP"
    log_tbl$message <- "insufficient time points or non-increasing time order"
    return(list(
      curve_long = curve_long,
      summary_tbl = summary_tbl,
      threshold_summary_long = NULL,
      threshold_summary_wide = NULL,
      summary_with_thresholds = NULL,
      log = log_tbl
    ))
  }
  
  # ---- 给 summary 加文件信息 ----
  summary_tbl <- summary_tbl %>%
    mutate(
      FileName = info$FileName[1],
      TestDate = info$TestDate[1],
      AssayBatch = info$AssayBatch[1],
      BatchID = info$BatchID[1]
    ) %>%
    select(FileName, TestDate, AssayBatch, BatchID, everything())
  
  # ---- 批量算 thresholds ----
  threshold_tbl <- map_dfr(thresholds, function(thr) {
    curve_long %>%
      group_by(FileName, TestDate, AssayBatch, BatchID, SampleID, Replicate, Group) %>%
      summarise(
        T_value = calc_threshold_time(Time_h, OD, threshold = thr),
        .groups = "drop"
      ) %>%
      mutate(
        threshold = thr,
        threshold_name = paste0("T", thr)
      )
  })
  
  threshold_summary_long <- threshold_tbl %>%
    group_by(FileName, TestDate, AssayBatch, BatchID, threshold, threshold_name, SampleID, Group) %>%
    summarise(
      T_mean = mean(T_value, na.rm = TRUE),
      T_sd = sd(T_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- 同板 MeOH reference ----
  meoh_ref <- threshold_summary_long %>%
    filter(SampleID == "MeOH") %>%
    select(FileName, threshold, threshold_name, MeOH_T_mean = T_mean)
  
  threshold_summary_long <- threshold_summary_long %>%
    left_join(meoh_ref, by = c("FileName", "threshold", "threshold_name")) %>%
    mutate(
      dT = T_mean - MeOH_T_mean
    )
  
  # 宽表（不保留 MeOH）
  threshold_summary_wide <- threshold_summary_long %>%
    filter(SampleID != "MeOH") %>%
    select(FileName, TestDate, AssayBatch, BatchID, SampleID, threshold_name, T_mean, T_sd, dT) %>%
    pivot_wider(
      names_from = threshold_name,
      values_from = c(T_mean, T_sd, dT),
      names_glue = "{threshold_name}_{.value}"
    )
  
  summary_with_thresholds <- summary_tbl %>%
    left_join(
      threshold_summary_wide,
      by = c("FileName", "TestDate", "AssayBatch", "BatchID", "SampleID")
    )
  
  list(
    curve_long = curve_long,
    summary_tbl = summary_tbl,
    threshold_summary_long = threshold_summary_long,
    threshold_summary_wide = threshold_summary_wide,
    summary_with_thresholds = summary_with_thresholds,
    log = log_tbl
  )
}

# =========================
# 8. Batch process
# =========================
results <- map(files, ~ process_one_file(.x, thresholds = c(0.05, 0.1, 0.15)))

# =========================
# 9. Combine outputs
# =========================
curve_long_all <- bind_rows(map(results, "curve_long"))
summary_tbl_all <- bind_rows(map(results, "summary_tbl"))
threshold_summary_long_all <- bind_rows(map(results, "threshold_summary_long"))
threshold_summary_wide_all <- bind_rows(map(results, "threshold_summary_wide"))
summary_with_thresholds_all <- bind_rows(map(results, "summary_with_thresholds"))
processing_log <- bind_rows(map(results, "log"))

# =========================
# 10. Save
# =========================
write.csv(
  curve_long_all,
  file.path(output_dir, "curve_long_all.csv"),
  row.names = FALSE
)

write.csv(
  summary_tbl_all,
  file.path(output_dir, "summary_tbl_all.csv"),
  row.names = FALSE
)

write.csv(
  threshold_summary_long_all,
  file.path(output_dir, "threshold_summary_long_all.csv"),
  row.names = FALSE
)

write.csv(
  threshold_summary_wide_all,
  file.path(output_dir, "threshold_summary_wide_all.csv"),
  row.names = FALSE
)

write.csv(
  summary_with_thresholds_all,
  file.path(output_dir, "summary_with_thresholds_all.csv"),
  row.names = FALSE
)

write.csv(
  processing_log,
  file.path(output_dir, "processing_log.csv"),
  row.names = FALSE
)

# 看一下日志
processing_log


threshold_summary_long_all %>%
  filter(SampleID == "MeOH") %>%
  group_by(threshold_name) %>%
  summarise(
    mean_T = mean(T_mean, na.rm = TRUE),
    sd_between_batches = sd(T_mean, na.rm = TRUE),
    .groups = "drop"
  )

