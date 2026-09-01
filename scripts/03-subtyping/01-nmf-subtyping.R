# NMF Subtyping and Pattern Analysis on Organ Risk Index

library(NMF)
library(tibble)
library(Cairo)

# Set file paths and output directory
input_organ_index <- "organ_risk_index.csv"
input_sample_info  <- "sample_metadata.csv"
  output_dir         <- "nmf-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Load Data & Preprocess
organ_index <- read.csv(input_organ_index, header = TRUE, row.names = 1, check.names = FALSE)
info        <- read.csv(input_sample_info,  header = TRUE, row.names = 1, check.names = FALSE)

# Filter for AD group samples
ad_info <- info[info$Group == "AD", , drop = FALSE]
common_ad_samples <- intersect(rownames(organ_index), rownames(ad_info))
if (length(common_ad_samples) < 3) {
  stop("At least three AD samples shared by the organ index and metadata are required.")
}

organ_index_aligned <- organ_index[common_ad_samples, , drop = FALSE]
normalized_data     <- t(organ_index_aligned) # Rows: Organs, Columns: Samples

# 2. Rank Evaluation
max_rank <- min(dim(normalized_data))
if (max_rank < 2) stop("The input matrix is too small for NMF rank evaluation.")
ranks <- 2:min(15, max_rank)
eval_result <- nmf(normalized_data, ranks, nrun = 5)

# Save Evaluation Metrics
sink(file.path(output_dir, "summary_rank_evaluation.txt"))
print(summary(eval_result))
sink()

saveRDS(eval_result, file = file.path(output_dir, "nmf_rank_evaluation.rds"))

# 3. Perform NMF decomposition for a selected rank
r <- min(3, max(ranks))
seed <- 123
res_rank <- nmf(normalized_data, rank = r, nrun = 50, seed = seed, method = "brunet")

saveRDS(res_rank, file = file.path(output_dir, paste0("nmf-rank-", r, ".rds")))

# Extract W (Basis) and H (Coefficient) Matrices
W <- basis(res_rank) # Organs x Ranks
H <- coef(res_rank)  # Ranks x Samples
sample_clusters <- predict(res_rank)

# Color Palette Configuration
cluster_colors   <- c("1" = "#F7C97E", "2" = "#A4BD9C", "3" = "#80B1D3", "4" = "#ECA8A9", "5" = "#CFAFD4")
consensus_colors <- c("1" = "#b58980", "2" = "#d7caa0", "3" = "#C2ABCB", "4" = "#a8b9a7", "5" = "#93B3B3")

# 4. Save Clustering Maps
CairoPDF(file.path(output_dir, paste0("consensus-map-rank-", r, ".pdf")), width = 10, height = 8)
consensusmap(res_rank, color = "-RdBu:50",
             annColors = list(basis = cluster_colors, consensus = consensus_colors))
dev.off()

CairoPDF(file.path(output_dir, paste0("basis-map-rank-", r, ".pdf")), width = 10, height = 8)
basismap(res_rank, color = "Reds:50", scale = "r1",
         annColors = list(basis = cluster_colors))
dev.off()

# 5. Export Processed Data & Matrices
all_ad_data <- t(rbind(normalized_data, Cluster = sample_clusters))
all_ad_merged <- merge(all_ad_data, ad_info["Range"], by = "row.names")
all_ad_merged <- column_to_rownames(all_ad_merged, var = "Row.names")

write.csv(W, file.path(output_dir, paste0("w-matrix-rank-", r, ".csv")), row.names = TRUE)
write.csv(H, file.path(output_dir, paste0("h-matrix-rank-", r, ".csv")), row.names = TRUE)
write.csv(all_ad_merged, file.path(output_dir, paste0("all-ad-clusters-rank-", r, ".csv")), row.names = TRUE)
