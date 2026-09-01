# PERMANOVA Analysis and Covariate Residual Adjustment

library(vegan)

# Create output directory
output_dir <- "adjusted"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
data <- read.csv("protein_expression.csv", header = TRUE, row.names = 1)
info <- read.csv("sample_metadata.csv", header = TRUE, row.names = 1)
if (!all(c("Group", "Sex", "Age", "TDI", "HQ", "E4_status") %in% names(info))) {
  stop("sample_metadata.csv is missing one or more required covariates.")
}

# Align sample IDs
inte <- intersect(rownames(data), rownames(info))
if (length(inte) < 3) stop("At least three shared samples are required.")
data <- data[inte, , drop = FALSE]
info <- info[inte, , drop = FALSE]
data[] <- lapply(data, function(x) as.numeric(as.character(x)))

# 2. Variable Encoding
qual_map <- c(
  "None_of_any_above" = 0, 
  "junior_high_school" = 1, 
  "senior_high_school" = 2, 
  "bachelor_degree_or_higher" = 3,
  "professional_qualifications" = 4, 
  "Prefer_not_to_answer" = 5
)
info$Highest.Qualification <- factor(qual_map[as.character(info$HQ)])

e4_map <- c("E4_carrier" = 1, "E4_nocarrier" = 0)
info$E4 <- factor(e4_map[as.character(info$E4_status)])

# 3. PERMANOVA on Unadjusted Raw Data
adonis_original <- adonis2(
  data ~ Group + Sex + Age + TDI + Highest.Qualification + E4,
  data = info, 
  permutations = 999, 
  method = "euclidean"
)
write.csv(adonis_original, file.path(output_dir, "permanova-original.csv"))

# Residual Correction Function
residual_correct <- function(y, covariates) {
  model <- lm(y ~ ., data = covariates, na.action = na.exclude)
  residuals(model) + mean(y, na.rm = TRUE)
}

# 4. Adjustment 1: Adjusting All Covariates
covariates_full <- info[, c("Sex", "Age", "TDI", "Highest.Qualification", "E4"), drop = FALSE]
covariates_full$Sex <- as.numeric(factor(covariates_full$Sex))

adjusted_data_full <- sapply(colnames(data), function(x) {
  residual_correct(data[, x], covariates_full)
})
rownames(adjusted_data_full) <- rownames(data)

adonis_adjusted_full <- adonis2(
  adjusted_data_full ~ Group + Sex + Age + TDI + Highest.Qualification + E4,
  data = info, 
  permutations = 999, 
  method = "euclidean"
)

write.csv(adonis_adjusted_full, file.path(output_dir, "permanova-adjusted-full.csv"))
write.csv(adjusted_data_full, file.path(output_dir, "protein-adjusted-full.csv"))

# 5. Adjustment 2: Adjusting Covariates Excluding E4
covariates_no_e4 <- info[, c("Sex", "Age", "TDI", "Highest.Qualification"), drop = FALSE]
covariates_no_e4$Sex <- as.numeric(factor(covariates_no_e4$Sex))

adjusted_data_no_e4 <- sapply(colnames(data), function(x) {
  residual_correct(data[, x], covariates_no_e4)
})
rownames(adjusted_data_no_e4) <- rownames(data)

adonis_adjusted_no_e4 <- adonis2(
  adjusted_data_no_e4 ~ Group + Sex + Age + TDI + Highest.Qualification + E4,
  data = info, 
  permutations = 999, 
  method = "euclidean"
)

write.csv(adonis_adjusted_no_e4, file.path(output_dir, "permanova-adjusted-no-e4.csv"))
write.csv(adjusted_data_no_e4, file.path(output_dir, "protein-adjusted-no-e4.csv"))
