# Soft Clustering Analysis of Time-Series/Longitudinal Data using Mfuzz

library(Mfuzz)
library(Biobase)

# Create output directory
output_dir <- "mfuzz-results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 1. Load Expression Data
expr <- read.csv("protein_time_series.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert data to ExpressionSet and standardize
expr_mat <- as.matrix(expr)
eset <- new("ExpressionSet", exprs = expr_mat)
eset_std <- standardise(eset)

# 2. Estimate Fuzzification Parameter (m)
m_est <- mestimate(eset_std)
cat("Estimated fuzzification parameter m:", m_est, "\n")

# Evaluate minimum centroid distance across potential cluster numbers
dmin <- Dmin(eset_std, m = m_est, crange = 2:5, repeats = 3, visu = FALSE)

pdf(file.path(output_dir, "cluster_number_selection.pdf"), width = 6, height = 5)
plot(
  2:5, dmin, type = "b", pch = 19,
  xlab = "Number of clusters",
  ylab = "Minimum centroid distance",
  main = "Selection of Cluster Number"
)
dev.off()

# 3. Perform Soft Clustering
num_clusters <- 4
set.seed(1234)
cl <- mfuzz(eset_std, c = num_clusters, m = m_est)

# Define time-point labels for plotting
time_labels <- colnames(expr)

# 4. Plot Clusters
pdf(file.path(output_dir, paste0("mfuzz_clusters_c", num_clusters, ".pdf")), width = 7, height = 6)
mfuzz.plot2(
  eset_std, cl = cl, mfrow = c(2, 2), time.labels = time_labels, x11 = FALSE,
  colo = "fancy", centre = TRUE, centre.col = "#8B0000",
  xlab = "Time Range"
)
dev.off()

# 5. Export Feature Membership Scores
cluster_membership <- acore(eset_std, cl = cl, min.acore = 0)
for (i in seq_len(num_clusters)) {
  write.csv(
    cluster_membership[[i]], 
    file = file.path(output_dir, paste0("cluster_", i, "_members.csv")),
    row.names = FALSE
  )
}
