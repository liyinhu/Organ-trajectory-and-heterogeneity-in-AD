# Organ Feature Selection via Lasso Regression

library(glmnet)
library(caret)
library(dplyr)

# Create output directory
output_dir <- "lasso-results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
protein_df <- read.csv("protein_adjusted.csv", header = TRUE, row.names = 1)
protein2organ_df <- read.csv("protein_to_organ.tsv", header = TRUE, sep = "\t")
organ_index_df <- read.csv("organ_risk_index.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Parameters
top_n <- 10
min_features <- 5

# Align samples across datasets
common_samples <- Reduce(intersect, list(
  rownames(organ_index_df),
  rownames(protein_df)
))

organ_index_df <- organ_index_df[common_samples, , drop = FALSE]
protein_df <- protein_df[common_samples, , drop = FALSE]
protein2organ_df <- protein2organ_df[, c("Protein", "Organ")]
protein_df[] <- lapply(protein_df, function(x) as.numeric(as.character(x)))
organ_index_df[] <- lapply(organ_index_df, function(x) as.numeric(as.character(x)))

all_feature_results <- list()  # Top N features per organ
all_coef_results <- list()     # All non-zero coefficient features

# 2. Iterate Over Each Organ
for (organ in colnames(organ_index_df)) {
  message("Processing organ: ", organ)
  
  # Retrieve mapped proteins for current organ
  protein_features <- protein2organ_df %>%
    filter(Organ == organ) %>%
    pull(Protein) %>%
    unique() %>%
    intersect(colnames(protein_df))
  
  # Filter out organs with insufficient feature count
  if (length(protein_features) < min_features) {
    message("Skipping ", organ, ": feature count less than ", min_features)
    next
  }
  
  # Build feature matrix and standardize
  X <- scale(protein_df[, protein_features, drop = FALSE])
  y <- organ_index_df[[organ]]

  complete <- complete.cases(X, y)
  X <- X[complete, , drop = FALSE]
  y <- y[complete]
  if (length(y) < 6) {
    message("Skipping ", organ, ": at least six complete samples are required.")
    next
  }
  
  # Partition training and testing sets
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, , drop = FALSE]
  y_test <- y[-train_idx]
  
  # Fit Lasso model with 5-fold cross-validation
  nfolds <- min(5, nrow(X_train))
  if (nfolds < 2) next
  cvfit <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = nfolds)
  best_lambda <- cvfit$lambda.min
  final_model <- glmnet(X_train, y_train, alpha = 1, lambda = best_lambda)
  
  # Extract non-zero coefficients
  coef_df <- coef(final_model)
  coef_table <- data.frame(
    feature = rownames(coef_df),
    coefficient = as.vector(coef_df)
  ) %>%
    filter(feature != "(Intercept)", coefficient != 0)
  
  if (nrow(coef_table) == 0) {
    message("No non-zero features selected for organ: ", organ)
    next
  }
  
  # Store all non-zero features
  all_coef_results[[organ]] <- coef_table %>%
    mutate(organ = organ) %>%
    select(organ, feature, coefficient)
  
  # Extract Top N features by absolute coefficient magnitude
  top_features <- coef_table %>%
    mutate(abs_coef = abs(coefficient)) %>%
    arrange(desc(abs_coef)) %>%
    slice_head(n = top_n) %>%
    mutate(organ = organ) %>%
    select(organ, feature, coefficient)
  
  all_feature_results[[organ]] <- top_features
}

# 3. Export Results
if (length(all_feature_results) > 0) {
  final_top_df <- do.call(rbind, all_feature_results)
  final_all_df <- do.call(rbind, all_coef_results)
  
  write.csv(final_top_df, file.path(output_dir, "lasso_features_top.csv"), row.names = FALSE)
  write.csv(final_all_df, file.path(output_dir, "lasso_features_all.csv"), row.names = FALSE)
  message("Lasso feature selection completed successfully.")
} else {
  warning("No features selected for any organ.")
}
