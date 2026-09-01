# Linear Model Differential Analysis of Organ-Specific ssGSEA Pathway Scores

library(dplyr)
library(tibble)
library(tidyr)

# File paths and configurations
ssgsea_file <- "organ_pathway_scores.csv"
info_file   <- "sample_metadata.csv"
output_dir  <- "ssgsea-lm-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Load Data & Align Samples
ssgsea_mat <- as.matrix(read.csv(ssgsea_file, row.names = 1, check.names = FALSE))
storage.mode(ssgsea_mat) <- "numeric"
info       <- read.csv(info_file, header = TRUE, row.names = 1, check.names = FALSE)
info       <- tibble::rownames_to_column(info, "Sample")

info$Group <- factor(info$Group, levels = c("NC", "AD"))
info$Range <- as.character(info$Range)

common_samples <- intersect(colnames(ssgsea_mat), info$Sample)
ssgsea_mat     <- ssgsea_mat[, common_samples, drop = FALSE]
info           <- info %>% filter(Sample %in% common_samples)
info           <- info[match(colnames(ssgsea_mat), info$Sample), ]

# Helper function for LM modeling
run_lm_analysis <- function(mat, metadata) {
  pathway_names <- rownames(mat)
  
  res_df <- data.frame(
    Pathway   = pathway_names,
    Estimate  = NA_real_,
    SE        = NA_real_,
    T_value   = NA_real_,
    P_value   = NA_real_,
    AD_median = NA_real_,
    NC_median = NA_real_,
    stringsAsFactors = FALSE
  )
  rownames(res_df) <- pathway_names
  
  for (i in seq_along(pathway_names)) {
    pw <- pathway_names[i]
    y  <- as.numeric(mat[i, ])
    
    tryCatch({
      fit   <- lm(y ~ Group, data = metadata)
      coefs <- summary(fit)$coefficients
      
      if ("GroupAD" %in% rownames(coefs)) {
        res_df[pw, c("Estimate", "SE", "T_value", "P_value")] <- coefs["GroupAD", 1:4]
      }
      
      res_df[pw, "AD_median"] <- median(y[metadata$Group == "AD"], na.rm = TRUE)
      res_df[pw, "NC_median"] <- median(y[metadata$Group == "NC"], na.rm = TRUE)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  res_df <- res_df %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    arrange(FDR)
  
  return(res_df)
}


# Part 1: Stratified Analysis by Range

all_ranges  <- sort(unique(info$Range))
range_list  <- list()

for (rg in all_ranges) {
  info_sub   <- info %>% filter(Range == rg)
  mat_sub    <- ssgsea_mat[, info_sub$Sample, drop = FALSE]
  info_sub   <- info_sub[match(colnames(mat_sub), info_sub$Sample), ]
  
  res_sub       <- run_lm_analysis(mat_sub, info_sub)
  res_sub$Range <- rg
  
  range_list[[rg]] <- res_sub
  
  write.csv(res_sub, file.path(output_dir, paste0("Organ_Pathway_lm.AD_vs_NC.Range_", rg, ".csv")), row.names = FALSE)
}

# Merge all stratified results and separate Organ/Pathway
all_range_res <- bind_rows(range_list) %>%
  separate(Pathway, into = c("Organ", "Pathway"), sep = "__", remove = FALSE, extra = "merge", fill = "right")

write.csv(all_range_res, file.path(output_dir, "Organ_Pathway_lm.AD_vs_NC.Stratified.csv"), row.names = FALSE)


# Part 2: All Samples Combined Analysis

all_res <- run_lm_analysis(ssgsea_mat, info)

all_res_split <- all_res %>%
  separate(Pathway, into = c("Organ", "Pathway"), sep = "__", remove = FALSE, extra = "merge", fill = "right")

write.csv(all_res_split, file.path(output_dir, "Organ_Pathway_lm.AD_vs_NC.All.csv"), row.names = FALSE)
