# 02_qc_and_normalization.R ----------------------------------------------
# Load 10x matrix into a Seurat object, compute QC metrics, filter low-quality
# cells, run SCTransform normalisation.
#
# Inputs : data_in/filtered_gene_bc_matrices/hg19/
# Outputs: data_out/rds/pbmc_sct.rds
#          data_out/figures/01_qc_violin_prefilter.png
#          data_out/figures/02_qc_scatter.png
#          data_out/tables/qc_summary.csv
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in <- file.path(root, "data_in", "filtered_gene_bc_matrices", "hg19")
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")
rds_dir <- file.path(root, "data_out", "rds")

# --- 1. Load + create Seurat object -------------------------------------
counts <- Read10X(data.dir = data_in)
pbmc   <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                              min.cells = 3, min.features = 200)
message("Initial: ", ncol(pbmc), " cells x ", nrow(pbmc), " genes")

# --- 2. Mitochondrial QC ------------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Violin of QC metrics BEFORE filtering
p_qc <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1) &
  theme(plot.title = element_text(size = 10))
ggsave(file.path(figs, "01_qc_violin_prefilter.png"), p_qc,
       width = 10, height = 4, dpi = 150)

# Scatter: nCount vs nFeature, coloured by %mt
p_sc1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p_sc2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(figs, "02_qc_scatter.png"), p_sc1 + p_sc2,
       width = 10, height = 4, dpi = 150)

# --- 3. Filtering -------------------------------------------------------
# Thresholds chosen from the distributions above (standard PBMC 3k cutoffs).
pre_n <- ncol(pbmc)
pbmc  <- subset(pbmc,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 2500 &
                         percent.mt < 5)
post_n <- ncol(pbmc)
message("After QC: ", post_n, " cells (removed ", pre_n - post_n, ")")

qc_df <- data.frame(
  stage          = c("pre_filter", "post_filter"),
  n_cells        = c(pre_n, post_n),
  median_nFeature = c(NA, median(pbmc$nFeature_RNA)),
  median_nCount   = c(NA, median(pbmc$nCount_RNA)),
  median_pct_mt   = c(NA, round(median(pbmc$percent.mt), 2))
)
write.csv(qc_df, file.path(tables, "qc_summary.csv"), row.names = FALSE)

# --- 4. SCTransform normalisation ---------------------------------------
# SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData.
# glmGamPoi accelerates the fit.
pbmc <- SCTransform(pbmc,
                    method = "glmGamPoi",
                    vars.to.regress = "percent.mt",
                    verbose = FALSE)

saveRDS(pbmc, file.path(rds_dir, "pbmc_sct.rds"))
message("SCTransform done. Saved SCT-normalised Seurat object.")
