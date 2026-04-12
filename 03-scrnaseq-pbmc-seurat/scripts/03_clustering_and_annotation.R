# 03_clustering_and_annotation.R -----------------------------------------
# PCA -> UMAP -> Louvain clustering -> cell-type annotation via canonical
# markers and module scoring.
#
# Inputs : data_out/rds/pbmc_sct.rds
# Outputs: data_out/rds/pbmc_annotated.rds
#          data_out/figures/03_elbow_plot.png
#          data_out/figures/04_umap_clusters.png
#          data_out/figures/05_umap_annotated.png
#          data_out/figures/06_canonical_markers_violin.png
#          data_out/figures/07_canonical_markers_umap.png
#          data_out/figures/09_cluster_module_scores.png
#          data_out/tables/cluster_annotations.csv
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")
rds_dir <- file.path(root, "data_out", "rds")

pbmc <- readRDS(file.path(rds_dir, "pbmc_sct.rds"))
set.seed(42)

# --- 1. PCA + elbow -----------------------------------------------------
pbmc <- RunPCA(pbmc, verbose = FALSE)

p_el <- ElbowPlot(pbmc, ndims = 50) +
  ggtitle("PCA elbow — standard deviation explained per PC") +
  theme_minimal(base_size = 11)
ggsave(file.path(figs, "03_elbow_plot.png"), p_el,
       width = 7, height = 4, dpi = 150)

# 10 PCs is the standard call for PBMC 3k; confirmed by the elbow.
n_pcs <- 10

# --- 2. UMAP + clustering ----------------------------------------------
pbmc <- RunUMAP(pbmc,       dims = 1:n_pcs, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs, verbose = FALSE)
pbmc <- FindClusters(pbmc,  resolution = 0.5, verbose = FALSE)

message("Clusters found: ", length(unique(Idents(pbmc))))

p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE,
                  label.size = 4, repel = TRUE) +
  ggtitle("UMAP — unsupervised Louvain clusters (resolution 0.5)")
ggsave(file.path(figs, "04_umap_clusters.png"), p_umap,
       width = 7, height = 6, dpi = 150)

# --- 3. Canonical markers -----------------------------------------------
canonical <- c("IL7R","CCR7","S100A4",     # CD4 T
               "CD14","LYZ",                # CD14 monocyte
               "MS4A1","CD79A",              # B cell
               "CD8A",                       # CD8 T
               "FCGR3A","MS4A7",             # FCGR3A monocyte
               "GNLY","NKG7",                # NK
               "FCER1A","CST3",              # DC
               "PPBP")                       # Platelet
canonical <- canonical[canonical %in% rownames(pbmc)]

p_vln <- VlnPlot(pbmc, features = canonical, pt.size = 0, stack = TRUE,
                 flip = TRUE) +
  NoLegend() +
  ggtitle("Canonical PBMC markers across clusters")
ggsave(file.path(figs, "06_canonical_markers_violin.png"), p_vln,
       width = 9, height = 8, dpi = 150)

p_feat <- FeaturePlot(pbmc,
                      features = c("IL7R","CD14","MS4A1","CD8A",
                                   "FCGR3A","GNLY","FCER1A","PPBP"),
                      ncol = 4, pt.size = 0.3)
ggsave(file.path(figs, "07_canonical_markers_umap.png"), p_feat,
       width = 14, height = 7, dpi = 150)

# --- 4. Module scoring for cell types -----------------------------------
# AddModuleScore computes the average expression of a gene set vs. a
# random control set. Useful as an unbiased cross-check of cluster identity.
modules <- list(
  CD4_T        = c("IL7R", "CCR7", "LEF1", "TCF7"),
  CD8_T        = c("CD8A", "CD8B", "GZMK"),
  B_cell       = c("MS4A1", "CD79A", "CD79B", "CD19"),
  NK           = c("GNLY", "NKG7", "GZMB", "PRF1"),
  CD14_mono    = c("CD14", "LYZ", "S100A8", "S100A9"),
  FCGR3A_mono  = c("FCGR3A", "MS4A7", "CDKN1C"),
  DC           = c("FCER1A", "CST3", "HLA-DQA1"),
  Platelet     = c("PPBP", "PF4", "GP1BA")
)
for (nm in names(modules)) {
  feats <- intersect(modules[[nm]], rownames(pbmc))
  pbmc  <- AddModuleScore(pbmc, features = list(feats),
                          name = paste0(nm, "_score"))
}
score_cols <- grep("_score1$", colnames(pbmc@meta.data), value = TRUE)

# Per-cluster mean module score -> used for annotation call
score_mat <- pbmc@meta.data %>%
  as_tibble() %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean)) %>%
  as.data.frame()
write.csv(score_mat, file.path(tables, "cluster_module_scores.csv"),
          row.names = FALSE)

# Heatmap of module scores by cluster
score_long <- tidyr::pivot_longer(score_mat,
                                   cols = -seurat_clusters,
                                   names_to = "module",
                                   values_to = "score")
score_long$module <- sub("_score1$", "", score_long$module)
p_score <- ggplot(score_long,
                  aes(x = module, y = seurat_clusters, fill = score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2980b9", mid = "white", high = "#c0392b",
                       midpoint = 0) +
  labs(x = "Cell-type module", y = "Cluster",
       title = "Mean module score per cluster") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figs, "09_cluster_module_scores.png"), p_score,
       width = 8, height = 5, dpi = 150)

# --- 5. Automatic annotation from module scores -------------------------
# Assign each cluster to the module with the highest mean score.
annot_df <- data.frame(
  cluster = score_mat$seurat_clusters,
  best_module = NA_character_,
  best_score = NA_real_
)
for (i in seq_len(nrow(score_mat))) {
  row <- score_mat[i, score_cols]
  mx  <- which.max(as.numeric(row))
  annot_df$best_module[i] <- sub("_score1$", "", score_cols[mx])
  annot_df$best_score[i]  <- as.numeric(row)[mx]
}
print(annot_df)
write.csv(annot_df, file.path(tables, "cluster_annotations.csv"),
          row.names = FALSE)

# Create a named vector to rename clusters on the Seurat object
new_ids <- annot_df$best_module
names(new_ids) <- annot_df$cluster
# Seurat's $<- needs an unnamed vector aligned to cells (names on RHS
# get interpreted as cell barcodes and confuse the merge)
pbmc$cell_type <- unname(new_ids[as.character(pbmc$seurat_clusters)])

Idents(pbmc) <- "cell_type"

p_umap2 <- DimPlot(pbmc, reduction = "umap", label = TRUE,
                   label.size = 4, repel = TRUE) +
  ggtitle("UMAP — annotated cell types (module-score assignment)")
ggsave(file.path(figs, "05_umap_annotated.png"), p_umap2,
       width = 8, height = 6, dpi = 150)

saveRDS(pbmc, file.path(rds_dir, "pbmc_annotated.rds"))
message("Annotation done. Saved pbmc_annotated.rds")
