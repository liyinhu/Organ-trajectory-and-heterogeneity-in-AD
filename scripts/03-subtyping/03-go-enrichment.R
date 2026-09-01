# GO Enrichment Analysis for Differential Proteins per Cluster

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

cluster_names   <- c("1", "2", "3")
diff_dir        <- "diff-results"
output_go_dir   <- "go-results"

if (!dir.exists(output_go_dir)) {
  dir.create(output_go_dir, recursive = TRUE)
}

# Helper function for GO enrichment
do_enrich <- function(genes, cls_name, type) {
  if (length(genes) < 5) {
    return(NULL)
  }
  
  gene_convert <- tryCatch({
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(gene_convert) || nrow(gene_convert) == 0) {
    return(NULL)
  }
  
  go_res <- tryCatch({
    enrichGO(
      gene          = gene_convert$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(go_res) && nrow(go_res) > 0) {
    out_file <- file.path(output_go_dir, paste0("go-cluster-", cls_name, "-", tolower(type), ".csv"))
    write.csv(as.data.frame(go_res), out_file, row.names = FALSE)
    return(go_res)
  }
  
  return(NULL)
}

# Iterate through clusters and run enrichment
for (cls_name in cluster_names) {
  diff_path <- file.path(diff_dir, paste0("cluster-", cls_name, "-vs-others.csv"))
  if (!file.exists(diff_path)) next
  
  diff_tab <- read.csv(diff_path, row.names = 1, check.names = FALSE)
  
  all_prots  <- rownames(diff_tab[!is.na(diff_tab$BH_adj_P) & diff_tab$BH_adj_P < 0.05, , drop = FALSE])
  up_prots   <- rownames(diff_tab[!is.na(diff_tab$BH_adj_P) & diff_tab$BH_adj_P < 0.05 & diff_tab$Estimate > 0, , drop = FALSE])
  down_prots <- rownames(diff_tab[!is.na(diff_tab$BH_adj_P) & diff_tab$BH_adj_P < 0.05 & diff_tab$Estimate < 0, , drop = FALSE])
  
  do_enrich(all_prots,  cls_name, "All")
  do_enrich(up_prots,   cls_name, "Up")
  do_enrich(down_prots, cls_name, "Down")
}
