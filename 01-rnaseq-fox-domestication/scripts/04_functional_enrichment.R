# 04_functional_enrichment.R ----------------------------------------------
# GO over-representation analysis on significant DEGs, using the
# Canis familiaris annotation (the silver fox genome is annotated against
# the dog reference -- ENSCAFG IDs).
#
# Inputs : data_out/tables/de_results_significant.csv
#          data_out/tables/de_results_full.csv
# Outputs: data_out/tables/go_enrichment_{BP,MF,CC}.csv
#          data_out/figures/07_go_dotplot.png
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Cf.eg.db)
  library(ggplot2)
})

root   <- rprojroot::find_root(rprojroot::has_file("README.md"))
tables <- file.path(root, "data_out", "tables")
figs   <- file.path(root, "data_out", "figures")

sig  <- read.csv(file.path(tables, "de_results_significant.csv"),
                 stringsAsFactors = FALSE)
full <- read.csv(file.path(tables, "de_results_full.csv"),
                 stringsAsFactors = FALSE)

# --- Gene universe (all tested genes with a mappable Entrez ID) ---------
universe <- unique(na.omit(full$entrez_id))
gene_sig <- unique(na.omit(sig$entrez_id))

message("Significant genes with Entrez ID: ", length(gene_sig),
        "  |  Universe: ", length(universe))

run_go <- function(ont) {
  enrichGO(gene          = gene_sig,
           universe      = universe,
           OrgDb         = org.Cf.eg.db,
           keyType       = "ENTREZID",
           ont           = ont,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.2,
           readable      = TRUE)
}

go_list <- list(BP = run_go("BP"),
                MF = run_go("MF"),
                CC = run_go("CC"))

for (ont in names(go_list)) {
  go <- go_list[[ont]]
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    write.csv(as.data.frame(go),
              file.path(tables, paste0("go_enrichment_", ont, ".csv")),
              row.names = FALSE)
    message("GO ", ont, ": ", nrow(as.data.frame(go)), " terms")
  } else {
    message("GO ", ont, ": no enriched terms at these thresholds.")
  }
}

# --- Dotplot of top BP terms --------------------------------------------
go_bp <- go_list$BP
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p <- dotplot(go_bp, showCategory = 15) +
    ggtitle("Top GO Biological Process terms",
            subtitle = "Over-representation among tame-vs-aggressive DEGs") +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9),
          axis.text.y = element_text(size = 8))
  ggsave(file.path(figs, "07_go_dotplot.png"), p,
         width = 9, height = 7, dpi = 150)
}

# --- Session info -------------------------------------------------------
sink(file.path(root, "data_out", "session_info.txt"))
sessionInfo()
sink()

message("Enrichment done. See data_out/tables/go_enrichment_*.csv")
