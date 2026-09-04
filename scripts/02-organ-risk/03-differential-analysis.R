# Linear Regression Analysis on Organ Risk Index across Subgroups

library(dplyr)
library(tibble)

# Set file paths and output directory
input_organ_index <- "organ_risk_index.csv"
input_sample_info  <- "sample_metadata.csv"
output_dir         <- "lm-results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
organ_index <- read.csv(input_organ_index, header = TRUE, row.names = 1, check.names = FALSE)
info        <- read.csv(input_sample_info,  header = TRUE, row.names = 1, check.names = FALSE)
if (!all(c("Group", "E4_status", "Range") %in% names(info))) {
  stop("sample_metadata.csv must contain Group, E4_status, and Range columns.")
}

# 2. Preprocess Sample Metadata
info <- info %>%
  rownames_to_column("Sample_ID") %>%
  mutate(
    Group_E4 = case_when(
      Group == "AD" & E4_status == "E4_nocarrier" ~ "AD_E4_nocarrier",
      Group == "AD" & E4_status == "E4_carrier"   ~ "AD_E4_carrier",
      Group == "NC" & E4_status == "E4_nocarrier" ~ "NC_E4_nocarrier",
      Group == "NC" & E4_status == "E4_carrier"   ~ "NC_E4_carrier",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group_E4), !is.na(Range)) %>%
  column_to_rownames("Sample_ID")

# Align samples
common_samples <- intersect(rownames(organ_index), rownames(info))
if (length(common_samples) < 3) stop("At least three shared samples are required.")
organ_index <- organ_index[common_samples, , drop = FALSE]
info        <- info[common_samples, , drop = FALSE]
organ_index[] <- lapply(organ_index, function(x) as.numeric(as.character(x)))

# 3. Helper Function for Fitting Linear Regression Models
run_lm <- function(data, group_vec, group_levels) {
  
  res <- data.frame(
    Estimate  = rep(NA, ncol(data)),
    SE        = rep(NA, ncol(data)),
    T_value   = rep(NA, ncol(data)),
    P_value   = rep(NA, ncol(data)),
    row.names = colnames(data)
  )
  
  group_vec <- factor(group_vec, levels = group_levels)
  
  for (k in seq_len(ncol(data))) {
    df <- data.frame(
      y = data[, k],
      g = group_vec
    )
    
    fit <- try(lm(y ~ g, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    sm <- summary(fit)$coefficients
    if (nrow(sm) < 2) next
    
    res[k, ] <- sm[2, c(1, 2, 3, 4)]
  }
  
  res$FDR <- p.adjust(res$P_value, method = "BH")
  res <- res[order(res$FDR), , drop = FALSE]
  
  return(res)
}

# 4. Perform Comparisons & Save Results

# Comparison 1: AD vs. NC
res_AD_vs_NC <- run_lm(
  data         = organ_index, 
  group_vec    = info$Group, 
  group_levels = c("NC", "AD")
)
write.csv(res_AD_vs_NC, file.path(output_dir, "organ_index_lm_AD_vs_NC.csv"))

# Comparison 2: AD E4 carrier vs. nocarrier
idx_AD <- info$Group %in% "AD"
if (sum(idx_AD) > 0) {
  res_AD_E4 <- run_lm(
    data         = organ_index[idx_AD, , drop = FALSE], 
    group_vec    = info$Group_E4[idx_AD],
    group_levels = c("AD_E4_nocarrier", "AD_E4_carrier")
  )
  write.csv(res_AD_E4, file.path(output_dir, "organ_index_lm_AD_E4_carrier_vs_nocarrier.csv"))
}

# Comparison 3: NC E4 carrier vs. nocarrier
idx_NC <- info$Group %in% "NC"
if (sum(idx_NC) > 0) {
  res_NC_E4 <- run_lm(
    data         = organ_index[idx_NC, , drop = FALSE],
    group_vec    = info$Group_E4[idx_NC],
    group_levels = c("NC_E4_nocarrier", "NC_E4_carrier")
  )
  write.csv(res_NC_E4, file.path(output_dir, "organ_index_lm_NC_E4_carrier_vs_nocarrier.csv"))
}
