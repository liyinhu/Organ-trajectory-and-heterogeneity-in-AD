# Construct Organ-Specific Pathways & Calculate ssGSEA Scores

library(msigdbr)
library(GSVA)
library(data.table)
library(dplyr)

# File paths and configurations
protein_organ_file <- "protein_to_organ.tsv"
protein_expr_file  <- "protein_adjusted.csv"
output_dir         <- "organ-ssgsea-results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Intermediate and Final Output Files
gmt_output_file <- file.path(output_dir, "c5.organ_specific_pathways.Hs.gmt")
rds_output_file <- file.path(output_dir, "c5.organ_specific_pathways.Hs.rds")
ssgsea_rds_file <- file.path(output_dir, "organ_pathway_scores.rds")
ssgsea_csv_file <- file.path(output_dir, "organ_pathway_scores.csv")

# Thresholds
min_gene_cutoff <- 5

# 1.1 Fetch C5 gene sets (GO/HP) for Homo sapiens
c5_df <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5"),
  error = function(e) msigdbr(species = "Homo sapiens", collection = "C5")
)
pathway_list <- split(c5_df$gene_symbol, c5_df$gs_name)
pathway_list <- lapply(pathway_list, unique)

# 1.2 Load Protein-Organ Mapping
organ_df <- fread(protein_organ_file, sep = "\t") %>%
  mutate(Protein = as.character(Protein), Organ = as.character(Organ))

organ_protein_list <- split(organ_df$Protein, organ_df$Organ)
organ_protein_list <- lapply(organ_protein_list, unique)

# 1.3 Construct Organ-Specific Gene Sets
organ_pathway_sets <- list()

for (org in names(organ_protein_list)) {
  org_proteins <- organ_protein_list[[org]]
  
  for (pw in names(pathway_list)) {
    genes <- intersect(org_proteins, pathway_list[[pw]])
    
    if (length(genes) >= min_gene_cutoff) {
      set_name <- paste(org, pw, sep = "__")
      organ_pathway_sets[[set_name]] <- genes
    }
  }
}

# Helper function to save GMT format
write_gmt <- function(gene_sets, file_path) {
  con <- file(file_path, "w")
  for (nm in names(gene_sets)) {
    line <- paste(c(nm, "NA", gene_sets[[nm]]), collapse = "\t")
    writeLines(line, con)
  }
  close(con)
}

# Save intermediate GMT and RDS
saveRDS(organ_pathway_sets, file = rds_output_file)
write_gmt(organ_pathway_sets, gmt_output_file)

# 2.1 Load Expression Data
expr_df <- fread(protein_expr_file, data.table = FALSE)

sample_col <- colnames(expr_df)[1]
expr_mat   <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- expr_df[[sample_col]]

expr_mat <- t(expr_mat)
mode(expr_mat) <- "numeric"

# Remove proteins with all NA values
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) < ncol(expr_mat), , drop = FALSE]

# 2.2 Intersect Pathways with Expression Universe
protein_universe   <- rownames(expr_mat)
organ_pathway_sets <- lapply(
  organ_pathway_sets,
  function(v) intersect(v, protein_universe)
)

# Filter out sets with fewer than min_gene_cutoff proteins in the matrix
organ_pathway_sets <- organ_pathway_sets[sapply(organ_pathway_sets, length) >= min_gene_cutoff]
if (length(organ_pathway_sets) == 0) {
  stop("No organ-specific gene sets overlap the expression matrix.")
}

# 2.3 Calculate ssGSEA Pathway Scores
param_organ <- ssgseaParam(
  exprData = expr_mat,
  geneSets = organ_pathway_sets,
  alpha    = 0.25,
  normalize= TRUE,
  minSize  = min_gene_cutoff
)

ssgsea_pathway_organ <- gsva(param_organ, verbose = FALSE)

# 2.4 Save Final ssGSEA Score Matrix
saveRDS(ssgsea_pathway_organ, file = ssgsea_rds_file)
write.csv(ssgsea_pathway_organ, file = ssgsea_csv_file, row.names = TRUE)
