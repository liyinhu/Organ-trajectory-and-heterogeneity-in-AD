# Shrinkage Partial Correlation Analysis for Subtype-Specific Significant Proteins

library(corpcor)
library(dplyr)

# File paths and configurations
protein_expr_file <- "protein_adjusted.csv"
cluster_info_file <- "cluster_membership.csv"
diff_dir          <- "diff-results"                         # Directory containing cluster_X_vs_others.csv
output_dir        <- "partial-correlation-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cluster_names <- c("1", "2", "3")

# 1. Load Data
protein_expr <- read.csv(protein_expr_file, header = TRUE, row.names = 1, check.names = FALSE)
cluster_info <- read.csv(cluster_info_file, header = TRUE, row.names = 1, check.names = FALSE)

# Align samples
common_samples <- intersect(rownames(protein_expr), rownames(cluster_info))
protein_expr   <- protein_expr[common_samples, , drop = FALSE]
cluster_info   <- cluster_info[common_samples, , drop = FALSE]
group          <- as.character(cluster_info$Cluster)

# 2. Function for Shrinkage Partial Correlation
get_partial_cor_shrinkage <- function(cls_name) {
  diff_path <- file.path(diff_dir, paste0("cluster-", cls_name, "-vs-others.csv"))
  if (!file.exists(diff_path)) {
    return(NULL)
  }
  
  diff_tab <- read.csv(diff_path, row.names = 1, check.names = FALSE)
  if (!"BH_adj_P" %in% colnames(diff_tab)) {
    return(NULL)
  }
  
  sig_prots <- rownames(diff_tab[!is.na(diff_tab$BH_adj_P) & diff_tab$BH_adj_P < 0.05, , drop = FALSE])
  sig_prots <- intersect(sig_prots, colnames(protein_expr))
  
  if (length(sig_prots) < 2) {
    return(NULL)
  }
  
  sub_data  <- protein_expr[group == cls_name, sig_prots, drop = FALSE]
  n_sample  <- nrow(sub_data)
  
  if (n_sample <= 3) {
    return(NULL)
  }
  
  # Compute Shrinkage Correlation
  sh_cor   <- cor.shrink(as.matrix(sub_data), verbose = FALSE)
  pcor_mat <- cor2pcor(sh_cor)
  
  rownames(pcor_mat) <- colnames(pcor_mat) <- sig_prots
  diag(pcor_mat)     <- 0
  
  df_used <- n_sample - 3
  
  full_res <- as.data.frame(as.table(pcor_mat)) %>%
    rename(Protein1 = Var1, Protein2 = Var2, Partial_r = Freq) %>%
    filter(as.character(Protein1) < as.character(Protein2)) %>%
    mutate(
      t_stat   = Partial_r * sqrt(df_used / (1 - Partial_r^2 + 1e-16)),
      P_value  = 2 * pt(abs(t_stat), df = df_used, lower.tail = FALSE),
      Cluster  = cls_name
    ) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    arrange(P_value)
  
  return(full_res)
}

# 3. Main Loop
all_subtype_results <- list()

for (cls in cluster_names) {
  res <- tryCatch({
    get_partial_cor_shrinkage(cls)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(res)) {
    all_subtype_results[[cls]] <- res
  }
}

# 4. Save Combined Outputs
if (length(all_subtype_results) > 0) {
  final_df <- bind_rows(all_subtype_results)
  
  write.csv(final_df, file.path(output_dir, "Subtype_Partial_Correlation.csv"), row.names = FALSE)
  saveRDS(all_subtype_results, file.path(output_dir, "Subtype_Partial_Correlation.rds"))
}
