# 04_marker_enrichment.R -------------------------------------------------
# Unbiased marker discovery per cluster via FindAllMarkers, heatmap of the
# top markers, and GO Biological Process enrichment per cluster to validate
# the biological identity of each cluster.
#
# Inputs : data_out/rds/pbmc_annotated.rds
# Outputs: data_out/tables/cluster_markers.csv
#          data_out/tables/cluster_markers_top10.csv
#          data_out/tables/go_enrichment_by_cluster.csv
#          data_out/figures/08_top_markers_heatmap.png
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")
rds_dir <- file.path(root, "data_out", "rds")

pbmc <- readRDS(file.path(rds_dir, "pbmc_annotated.rds"))

# Run marker discovery on the original Louvain cluster labels (numeric),
# which are more granular than the collapsed cell-type labels.
Idents(pbmc) <- "seurat_clusters"

# Seurat's FindMarkers uses `future` parallelism; the default global-export
# cap (500 MB) is too small for SCT assays. Raise it to 2 GB.
options(future.globals.maxSize = 2 * 1024^3)

# --- 1. FindAllMarkers --------------------------------------------------
# Wilcoxon rank-sum, only positive markers, require expression in >=25% of
# cluster cells and >=0.25 log2FC vs. rest. DefaultAssay set to SCT so
# we're testing SCT-normalised data.
DefaultAssay(pbmc) <- "SCT"
pbmc <- PrepSCTFindMarkers(pbmc, verbose = FALSE)

markers <- FindAllMarkers(pbmc,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          verbose = FALSE)
write.csv(markers, file.path(tables, "cluster_markers.csv"),
          row.names = FALSE)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()
write.csv(top10, file.path(tables, "cluster_markers_top10.csv"),
          row.names = FALSE)

message("Found ", nrow(markers), " cluster markers across ",
        length(unique(markers$cluster)), " clusters.")

# --- 2. Top-markers heatmap --------------------------------------------
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

# Downsample to avoid an impossibly wide heatmap
set.seed(42)
pbmc_sub <- subset(pbmc,
                   downsample = min(50, min(table(pbmc$seurat_clusters))))

png(file.path(figs, "08_top_markers_heatmap.png"),
    width = 2400, height = 2000, res = 220)
print(
  DoHeatmap(pbmc_sub, features = top5$gene, size = 3) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 6))
)
dev.off()

# --- 3. GO enrichment per cluster --------------------------------------
# Use top 50 markers per cluster for enrichment. Map symbols -> Entrez.
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup()

all_universe <- unique(markers$gene)
universe_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                         keys = all_universe,
                                         keytype = "SYMBOL",
                                         column = "ENTREZID",
                                         multiVals = "first")
universe_entrez <- unique(na.omit(universe_entrez))

go_all <- list()
for (cl in sort(unique(top50$cluster))) {
  gsyms <- top50$gene[top50$cluster == cl]
  gids  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = gsyms,
                                 keytype = "SYMBOL",
                                 column = "ENTREZID",
                                 multiVals = "first")
  gids <- unique(na.omit(gids))
  if (length(gids) < 5) next
  ego <- enrichGO(gene = gids, universe = universe_entrez,
                  OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    df <- as.data.frame(ego)
    df$cluster <- cl
    go_all[[as.character(cl)]] <- df
  }
}

if (length(go_all) > 0) {
  go_df <- do.call(rbind, go_all)
  go_df <- go_df[, c("cluster", "ID", "Description", "GeneRatio",
                     "pvalue", "p.adjust", "Count", "geneID")]
  write.csv(go_df, file.path(tables, "go_enrichment_by_cluster.csv"),
            row.names = FALSE)
  message("GO enrichment: ", nrow(go_df), " terms across ",
          length(go_all), " clusters.")
} else {
  message("No GO terms enriched at these thresholds.")
}

# Session info
sink(file.path(root, "data_out", "session_info.txt"))
sessionInfo()
sink()

message("Marker + enrichment step done.")
