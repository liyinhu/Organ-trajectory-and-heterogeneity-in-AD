# Organ Risk Network Analysis using Glasso and Bootstrap Resampling

library(bootnet)
library(dplyr)

# Input and output file path configurations
input_organ_index <- "organ_risk_index.csv"
input_cluster_info<- "cluster_membership.csv"
output_dir        <- "network-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Parameter settings
cluster_names <- c("1", "2", "3")
B_SAMPLES     <- 1000
TUNING        <- 0.5
n_cores       <- max(1, parallel::detectCores(logical = TRUE) - 1)

# 1. Load Data
organ <- read.csv(input_organ_index, header = TRUE, row.names = 1, check.names = FALSE)
info  <- read.csv(input_cluster_info, header = TRUE, row.names = 1, check.names = FALSE)

# Align samples
common_samples <- intersect(rownames(organ), rownames(info))
if (length(common_samples) < 6) stop("At least six shared samples are required for network analysis.")
organ <- organ[common_samples, , drop = FALSE]
info  <- info[common_samples, , drop = FALSE]

group <- info$Cluster

# 2. Data Normalization (Z-score Scaling)
organ_z <- scale(organ)
organ_z <- as.data.frame(organ_z)

# Remove zero-variance / low-variance features
organ_z    <- organ_z[, apply(organ_z, 2, var) > 1e-6, drop = FALSE]
node_names <- colnames(organ_z)

# 3. Subsetting Data by Clusters
data_C1 <- organ_z[group == cluster_names[1], , drop = FALSE]
data_C2 <- organ_z[group == cluster_names[2], , drop = FALSE]
data_C3 <- organ_z[group == cluster_names[3], , drop = FALSE]

# 4. Bootstrap EBICglasso Network Inference
set.seed(123)

boot_C1 <- bootnet(data_C1, default = "EBICglasso", nBoots = B_SAMPLES, 
                   type = "nonparametric", tuning = TUNING, nCores = n_cores, threshold = TRUE)
boot_C2 <- bootnet(data_C2, default = "EBICglasso", nBoots = B_SAMPLES, 
                   type = "nonparametric", tuning = TUNING, nCores = n_cores, threshold = TRUE)
boot_C3 <- bootnet(data_C3, default = "EBICglasso", nBoots = B_SAMPLES, 
                   type = "nonparametric", tuning = TUNING, nCores = n_cores, threshold = TRUE)

# 5. Extract Edge Properties Function
get_final_edges_robust <- function(boot_obj, name) {
  adj_mat <- as.matrix(boot_obj$sample$graph)
  
  orig_df <- as.data.frame(as.table(adj_mat)) %>%
    rename(Node1 = Var1, Node2 = Var2, Weight = Freq) %>%
    filter(Weight != 0) %>%
    mutate(
      Node1 = as.character(Node1),
      Node2 = as.character(Node2),
      tmp1  = pmin(Node1, Node2),
      tmp2  = pmax(Node1, Node2)
    ) %>%
    filter(Node1 < Node2) %>%
    select(-Node1, -Node2) %>%
    rename(Node1 = tmp1, Node2 = tmp2)
  
  boot_stats <- boot_obj$bootTable %>%
    filter(type == "edge") %>%
    mutate(
      node1 = as.character(node1),
      node2 = as.character(node2),
      tmp1  = pmin(node1, node2),
      tmp2  = pmax(node1, node2)
    ) %>%
    group_by(tmp1, tmp2) %>%
    summarise(
      Boot_Mean = mean(value, na.rm = TRUE),
      CI_lower  = quantile(value, 0.025, na.rm = TRUE),
      CI_upper  = quantile(value, 0.975, na.rm = TRUE),
      p_value   = mean(value == 0 | sign(value) != sign(mean(value, na.rm = TRUE)), na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    rename(Node1 = tmp1, Node2 = tmp2)
  
  final_df <- orig_df %>%
    left_join(boot_stats, by = c("Node1", "Node2")) %>%
    mutate(
      Cluster      = name,
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  return(final_df)
}

# Combine and Export Edges
final_edges_all <- bind_rows(
  get_final_edges_robust(boot_C1, paste0("Cluster", cluster_names[1])),
  get_final_edges_robust(boot_C2, paste0("Cluster", cluster_names[2])),
  get_final_edges_robust(boot_C3, paste0("Cluster", cluster_names[3]))
)

write.csv(final_edges_all, file.path(output_dir, "Organ_Network_Edges.csv"), row.names = FALSE)

# 6. Extract Node Centrality Metrics Function
get_final_nodes_safe <- function(boot_obj, name, nodes) {
  adj_mat  <- as.matrix(boot_obj$sample$graph)
  strength <- rowSums(abs(adj_mat))
  
  node_df <- data.frame(
    Cluster  = name,
    Node     = nodes,
    Strength = strength,
    row.names = NULL
  )
  
  return(node_df)
}

# Combine and Export Nodes
final_nodes_all <- bind_rows(
  get_final_nodes_safe(boot_C1, cluster_names[1], node_names),
  get_final_nodes_safe(boot_C2, cluster_names[2], node_names),
  get_final_nodes_safe(boot_C3, cluster_names[3], node_names)
)

write.csv(final_nodes_all, file.path(output_dir, "Organ_Network_Node.csv"), row.names = FALSE)
