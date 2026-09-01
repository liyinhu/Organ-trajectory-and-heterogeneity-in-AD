# Differential Expression Analysis between Two Groups

# Create output directory
output_dir <- "diff"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
protein_data <- read.csv("protein_adjusted.csv", header = TRUE, row.names = 1,
                         check.names = FALSE)
info <- read.csv("sample_metadata.csv", header = TRUE, row.names = 1,
                 check.names = FALSE)
if (!"Group" %in% names(info)) stop("sample_metadata.csv must contain a Group column.")

# Align samples and match group info
inte <- intersect(rownames(protein_data), rownames(info))
if (length(inte) < 3) stop("At least three shared samples are required.")
protein_data <- protein_data[inte, , drop = FALSE]
info <- info[inte, , drop = FALSE]

protein_data[] <- lapply(protein_data, function(x) as.numeric(as.character(x)))
group <- factor(info$Group, levels = c("NC", "AD"))
if (nlevels(droplevels(group)) < 2) stop("Both NC and AD groups are required.")

# 2. Linear Regression for Differential Analysis
n_features <- ncol(protein_data)
feature_names <- colnames(protein_data)

test_summary <- data.frame(
  Estimate = numeric(n_features),
  SE = numeric(n_features),
  T_value = numeric(n_features),
  P_value = numeric(n_features),
  row.names = feature_names
)

for (i in seq_len(n_features)) {
  tryCatch({
    fit <- lm(protein_data[, i] ~ group)
    coefficients <- summary(fit)$coefficients
    if (nrow(coefficients) >= 2) {
      test_summary[i, 1:4] <- coefficients[2, 1:4]
    }
  }, error = function(e) {
    cat("Error fitting model for feature", feature_names[i], ":", conditionMessage(e), "\n")
  })
}

# 3. Multiple Testing Correction (Benjamini-Hochberg)
test_summary$FDR_BH <- p.adjust(test_summary$P_value, method = "BH")
test_summary <- test_summary[order(test_summary$FDR_BH), ]

# 4. Save Differential Analysis Results
write.csv(test_summary, file.path(output_dir, "differential_results.csv"), row.names = TRUE)
