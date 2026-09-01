# Calculate Subtype-Specific Driver Scores and Ranks

library(dplyr)
library(readr)
library(tidyr)

# File paths and configurations
net_file          <- "protein_network.csv"
protein_organ_file<- "protein_to_organ.tsv"
diff_dir          <- "diff-results"                            # Directory with cluster_X_vs_others.csv
output_dir        <- "driver-score-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cluster_names <- c("1", "2", "3")

# Filtering Thresholds
filter_network_by_p <- TRUE
network_p_cutoff     <- 0.05
filter_network_by_r <- TRUE
network_r_cutoff     <- 0.15
eps                 <- 1e-300

# Load Mapping & Network Files
protein2organ_df <- read.table(protein_organ_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("Organ", "Protein") %in% names(protein2organ_df))) {
  stop("protein_to_organ.tsv must contain Organ and Protein columns.")
}
protein2organ_df <- protein2organ_df %>%
  transmute(Protein_ID = Protein, Organ = Organ)

all_net_df <- read_csv(net_file, show_col_types = FALSE)

# Safe Z-Score Normalization
zscore_safe <- function(x) {
  if (length(x) <= 1) return(rep(0, length(x)))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  return((x - m) / s)
}

# Main Subtype Analysis Loop
for (cls in cluster_names) {
  
  # A. Read Differential Expression File
  diff_path <- file.path(diff_dir, paste0("cluster-", cls, "-vs-others.csv"))
  if (!file.exists(diff_path)) next
  
  diff_df <- read.csv(diff_path, header = TRUE, check.names = FALSE)
  colnames(diff_df)[1] <- "Protein_ID"
  
  diff_df <- diff_df %>%
    dplyr::select(Protein_ID, Estimate, BH_adj_P) %>%
    mutate(
      Estimate = as.numeric(Estimate),
      BH_adj_P = as.numeric(BH_adj_P),
      BH_adj_P = ifelse(is.na(BH_adj_P), 1, pmax(BH_adj_P, eps))
    ) %>%
    filter(!is.na(Protein_ID), !is.na(Estimate)) %>%
    mutate(Differential_score = abs(Estimate) * (-log10(BH_adj_P)))
  
  # B. Extract Subtype-Specific Partial Correlation Network
  net_sub <- all_net_df %>%
    filter(as.character(Cluster) == as.character(cls)) %>%
    mutate(Partial_r = as.numeric(Partial_r), P_value = as.numeric(P_value)) %>%
    filter(!is.na(Protein1), !is.na(Protein2), Protein1 != Protein2)
  
  if (filter_network_by_p) net_sub <- net_sub %>% filter(P_value < network_p_cutoff)
  if (filter_network_by_r) net_sub <- net_sub %>% filter(abs(Partial_r) >= network_r_cutoff)
  
  # C. Compute Network Score (Sum of Absolute Partial Correlations)
  net_long <- bind_rows(
    net_sub %>% transmute(Protein_ID = Protein1, abs_r = abs(Partial_r)),
    net_sub %>% transmute(Protein_ID = Protein2, abs_r = abs(Partial_r))
  )
  
  net_score_df <- net_long %>%
    group_by(Protein_ID) %>%
    summarise(
      Network_score = sum(abs_r, na.rm = TRUE),
      Degree        = n(),
      .groups       = "drop"
    )
  
  # D. Merge Differential & Network Scores to Derive Driver Score & Ranks
  score_df <- diff_df %>%
    dplyr::left_join(net_score_df, by = "Protein_ID") %>%
    dplyr::mutate(
      Network_score = tidyr::replace_na(Network_score, 0),
      Degree        = tidyr::replace_na(Degree, 0)
    ) %>%
    dplyr::mutate(
      Z_diff       = zscore_safe(Differential_score),
      Z_net        = zscore_safe(Network_score),
      Driver_score = Z_diff + Z_net
    ) %>%
    dplyr::mutate(
      Driver_rank = rank(-Driver_score, ties.method = "min"),
      Diff_rank   = rank(-Differential_score, ties.method = "min"),
      Net_rank    = rank(-Network_score, ties.method = "min")
    ) %>%
    dplyr::arrange(Driver_rank)
  
  # Export Results
  out_prefix <- file.path(output_dir, paste0("Subtype", cls))
  
  write_csv(score_df, paste0(out_prefix, "_DriverScores.csv"))
  
  score_df_with_organ <- score_df %>%
    dplyr::left_join(protein2organ_df %>% dplyr::select(Protein_ID, Organ), by = "Protein_ID")
  
  write_csv(score_df_with_organ, paste0(out_prefix, "_DriverScores.Organ.csv"))
}
