# 02_qc_exploration.R ------------------------------------------------------
# QC and exploratory analysis of the silver fox prefrontal cortex RNA-seq
# data. We focus on the prefrontal cortex (primary tissue of interest in
# Wang et al. 2018) and set up the DESeq2 dataset used downstream.
#
# Inputs : data_in/GSE76517_fox_cortex_read_counts_by_gene.txt.gz
#          data_in/sample_metadata.csv
# Outputs: data_out/figures/01_library_sizes.png
#          data_out/figures/02_pca.png
#          data_out/figures/03_sample_distance_heatmap.png
#          data_out/tables/normalised_counts.csv
#          data_out/tables/vst_counts.csv
#          data_out/rds/dds.rds   (DESeq2 dataset, for script 03)
#          data_out/rds/vsd.rds   (VST object, for script 03)
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

root      <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in   <- file.path(root, "data_in")
figs      <- file.path(root, "data_out", "figures")
tables    <- file.path(root, "data_out", "tables")
rds_dir   <- file.path(root, "data_out", "rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load counts and metadata -----------------------------------------
counts_file <- file.path(data_in, "GSE76517_fox_cortex_read_counts_by_gene.txt.gz")
counts_raw  <- read.table(gzfile(counts_file), header = TRUE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$dog_gene_id
counts_mat <- as.matrix(counts_raw[, -1, drop = FALSE])

# RSEM outputs non-integer counts; DESeq2 needs integers
counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

meta_all <- read.csv(file.path(data_in, "sample_metadata.csv"),
                     stringsAsFactors = FALSE)

# Keep only prefrontal cortex samples and match to count columns.
# Count cols look like "fox181"; titles contain "ID181".
meta_pfc <- meta_all[meta_all$tissue == "prefrontal_cortex", ]
fox_num  <- sub(".*ID", "", meta_pfc$title)
meta_pfc$sample_id <- paste0("fox", fox_num)

# Reorder meta to match count matrix column order
common <- intersect(colnames(counts_mat), meta_pfc$sample_id)
counts_mat <- counts_mat[, common, drop = FALSE]
meta_pfc   <- meta_pfc[match(common, meta_pfc$sample_id), ]
stopifnot(all(colnames(counts_mat) == meta_pfc$sample_id))

meta_pfc$group <- factor(meta_pfc$group, levels = c("aggressive", "tame"))
# aggressive is the reference -> log2FC > 0 means "up in tame"
message("Samples per group:")
print(table(meta_pfc$group))

# --- 2. Build DESeq2 dataset ---------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = meta_pfc,
                              design    = ~ group)

# Pre-filter: remove genes with < 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message("Genes retained after filtering: ", nrow(dds))

# Run DESeq2 (normalisation + dispersion + Wald)
dds <- DESeq(dds)

# Variance stabilising transform for QC / visualisation
vsd <- vst(dds, blind = TRUE)

# --- 3. Library sizes ----------------------------------------------------
lib_df <- data.frame(sample     = colnames(dds),
                     group      = colData(dds)$group,
                     total_reads = colSums(counts(dds)))

p_lib <- ggplot(lib_df, aes(x = reorder(sample, total_reads),
                            y = total_reads / 1e6, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(aggressive = "#c0392b", tame = "#2980b9")) +
  labs(x = NULL, y = "Total reads (millions)",
       title = "Library size per sample",
       subtitle = "Prefrontal cortex, GSE76517") +
  theme_minimal(base_size = 11)

ggsave(file.path(figs, "01_library_sizes.png"), p_lib,
       width = 7, height = 6, dpi = 150)

# --- 4. PCA --------------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = group, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 30,
                           show.legend = FALSE) +
  scale_colour_manual(values = c(aggressive = "#c0392b", tame = "#2980b9")) +
  labs(x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       title = "PCA of VST-transformed counts",
       subtitle = "Prefrontal cortex, tame vs aggressive foxes") +
  theme_minimal(base_size = 11)

# ggrepel is optional; fall back to geom_text if missing
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = group, label = name)) +
    geom_point(size = 3) +
    geom_text(size = 3, vjust = -0.8, show.legend = FALSE) +
    scale_colour_manual(values = c(aggressive = "#c0392b", tame = "#2980b9")) +
    labs(x = paste0("PC1: ", pct_var[1], "% variance"),
         y = paste0("PC2: ", pct_var[2], "% variance"),
         title = "PCA of VST-transformed counts",
         subtitle = "Prefrontal cortex, tame vs aggressive foxes") +
    theme_minimal(base_size = 11)
}

ggsave(file.path(figs, "02_pca.png"), p_pca,
       width = 7, height = 5.5, dpi = 150)

# --- 5. Sample-to-sample distance heatmap --------------------------------
sample_dists     <- dist(t(assay(vsd)))
sample_dist_mat  <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- paste0(colData(vsd)$sample_id, " (",
                                    colData(vsd)$group, ")")
colnames(sample_dist_mat) <- NULL

png(file.path(figs, "03_sample_distance_heatmap.png"),
    width = 2000, height = 1800, res = 200)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample Euclidean distance (VST)")
dev.off()

# --- 6. Save normalised / VST counts and the DESeq2 object ---------------
write.csv(counts(dds, normalized = TRUE),
          file.path(tables, "normalised_counts.csv"))
write.csv(assay(vsd),
          file.path(tables, "vst_counts.csv"))

saveRDS(dds, file.path(rds_dir, "dds.rds"))
saveRDS(vsd, file.path(rds_dir, "vsd.rds"))

message("QC done. Figures in data_out/figures/, tables in data_out/tables/.")
