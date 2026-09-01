# Differential Protein Expression Analysis (One-vs-Rest per Cluster)

library(dplyr)

# Inputs and output directories
input_protein_data <- "protein_adjusted.csv"
input_cluster_info  <- "cluster_membership.csv"
output_diff_dir     <- "diff-results"

if (!dir.exists(output_diff_dir)) {
  dir.create(output_diff_dir, recursive = TRUE)
}

# 1. Load Data
protein_data <- read.csv(input_protein_data, header = TRUE, row.names = 1, check.names = FALSE)
cluster_info <- read.csv(input_cluster_info,  header = TRUE, row.names = 1, check.names = FALSE)

# Align samples
common_samples <- intersect(rownames(protein_data), rownames(cluster_info))
protein_data   <- protein_data[common_samples, , drop = FALSE]
cluster_info   <- cluster_info[common_samples, , drop = FALSE]

cluster_labels  <- as.character(cluster_info$Cluster)
unique_clusters <- sort(unique(cluster_labels))

# 2. Linear Regression per Cluster
for (current_cluster in unique_clusters) {
  
  temp_data  <- protein_data
  temp_group <- ifelse(cluster_labels == current_cluster, "Target", "Others")
  temp_data$Comparison_Group <- factor(temp_group, levels = c("Others", "Target"))
  
  num_proteins <- ncol(protein_data)
  test_summary <- data.frame(
    Estimate  = rep(NA, num_proteins),
    SE        = rep(NA, num_proteins),
    t_value   = rep(NA, num_proteins),
    P_value   = rep(NA, num_proteins),
    row.names = colnames(protein_data)
  )
  
  for (i in seq_len(num_proteins)) {
    tryCatch({
      fit <- lm(temp_data[, i] ~ Comparison_Group, data = temp_data)
      test_summary[i, 1:4] <- summary(fit)$coefficients[2, 1:4]
    }, error = function(e) {})
  }
  
  test_summary$BH_adj_P <- p.adjust(test_summary$P_value, method = "BH")
  test_summary          <- test_summary[order(test_summary$P_value), , drop = FALSE]
  
  diff_file_name <- file.path(output_diff_dir, paste0("cluster-", current_cluster, "-vs-others.csv"))
  write.csv(test_summary, diff_file_name, row.names = TRUE)
}
