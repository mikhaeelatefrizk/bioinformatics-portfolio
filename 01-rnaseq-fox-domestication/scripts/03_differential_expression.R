# 03_differential_expression.R --------------------------------------------
# Wald test for tame vs aggressive in the prefrontal cortex, with apeglm
# log-fold-change shrinkage. Annotates Ensembl dog gene IDs with SYMBOL and
# ENTREZID via org.Cf.eg.db, and writes a publication-ready volcano plot.
#
# Inputs : data_out/rds/dds.rds, data_out/rds/vsd.rds
# Outputs: data_out/tables/de_results_full.csv
#          data_out/tables/de_results_significant.csv
#          data_out/figures/04_volcano.png
#          data_out/figures/05_ma_plot.png
#          data_out/figures/06_top_genes_heatmap.png
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(biomaRt)
  library(ggplot2)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
rds_dir <- file.path(root, "data_out", "rds")
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")

dds <- readRDS(file.path(rds_dir, "dds.rds"))
vsd <- readRDS(file.path(rds_dir, "vsd.rds"))

# --- 1. Results + apeglm LFC shrinkage -----------------------------------
# resultsNames(dds)[2] is the tame_vs_aggressive contrast
coef_name <- grep("tame", resultsNames(dds), value = TRUE)
message("Using contrast: ", coef_name)

res_raw <- results(dds, name = coef_name, alpha = 0.05)
res     <- lfcShrink(dds, coef = coef_name, type = "apeglm")

summary(res)

# --- 2. Annotate via Ensembl 100 archive (last CanFam3.1 release) --------
# The counts use old CanFam3.1 Ensembl IDs (ENSCAFG00000...). Current
# org.Cf.eg.db has switched to the UU_Cfam_GSD assembly (ENSCAFG00845...),
# so we query Ensembl's Apr 2020 archive which still hosts CanFam3.1.
ens <- rownames(res)
message("Mapping ", length(ens), " gene IDs via Ensembl 100 (CanFam3.1) ...")

mart <- useEnsembl(biomart = "genes",
                   dataset = "clfamiliaris_gene_ensembl",
                   version = 100)
map <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "entrezgene_id",
                            "description"),
             filters = "ensembl_gene_id",
             values  = ens,
             mart    = mart)

# De-duplicate: one row per ensembl ID (keep first non-NA entrez)
map <- map[order(map$ensembl_gene_id, is.na(map$entrezgene_id)), ]
map <- map[!duplicated(map$ensembl_gene_id), ]

# Save the mapping table so script 04 doesn't need to hit biomaRt again
write.csv(map, file.path(tables, "gene_annotation_canfam31.csv"),
          row.names = FALSE)

sym    <- map$external_gene_name[match(ens, map$ensembl_gene_id)]
entrez <- map$entrezgene_id     [match(ens, map$ensembl_gene_id)]
descr  <- map$description       [match(ens, map$ensembl_gene_id)]

res_df <- data.frame(
  ensembl_id  = ens,
  symbol      = unname(sym),
  entrez_id   = unname(entrez),
  description = unname(descr),
  baseMean    = res$baseMean,
  log2FC      = res$log2FoldChange,
  lfcSE       = res$lfcSE,
  pvalue      = res$pvalue,
  padj        = res$padj,
  stringsAsFactors = FALSE
)
# Sort by padj (NAs last), then by abs(log2FC)
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FC), na.last = TRUE), ]

write.csv(res_df, file.path(tables, "de_results_full.csv"), row.names = FALSE)

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
write.csv(sig, file.path(tables, "de_results_significant.csv"), row.names = FALSE)

message("DE results: ", nrow(res_df), " genes total, ",
        nrow(sig), " significant at padj < 0.05.")

# --- 3. Volcano plot -----------------------------------------------------
# Use gene symbols where available, else Ensembl ID
labels <- ifelse(is.na(res_df$symbol) | res_df$symbol == "",
                 res_df$ensembl_id, res_df$symbol)
res_for_plot <- res_df
rownames(res_for_plot) <- make.unique(labels)

png(file.path(figs, "04_volcano.png"), width = 2200, height = 1800, res = 220)
print(
  EnhancedVolcano(res_for_plot,
    lab            = rownames(res_for_plot),
    x              = "log2FC",
    y              = "padj",
    pCutoff        = 0.05,
    FCcutoff       = 1,
    pointSize      = 2,
    labSize        = 3.2,
    title          = "Tame vs Aggressive foxes (prefrontal cortex)",
    subtitle       = "DESeq2 + apeglm LFC shrinkage",
    caption        = paste0("padj < 0.05 & |log2FC| > 1  |  ",
                            "Data: Wang et al. 2018, GSE76517"),
    legendPosition = "right",
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    colAlpha       = 0.7)
)
dev.off()

# --- 4. MA plot ----------------------------------------------------------
png(file.path(figs, "05_ma_plot.png"), width = 1800, height = 1400, res = 200)
plotMA(res, ylim = c(-4, 4),
       main = "MA plot: tame vs aggressive (apeglm-shrunken)")
dev.off()

# --- 5. Heatmap of top 30 DEGs -------------------------------------------
top_n <- 30
top_genes <- head(sig$ensembl_id, top_n)
if (length(top_genes) >= 4) {
  vst_mat <- assay(vsd)[top_genes, , drop = FALSE]
  # z-score across samples
  vst_z   <- t(scale(t(vst_mat)))
  lbl     <- sig$symbol[match(top_genes, sig$ensembl_id)]
  lbl[is.na(lbl) | lbl == ""] <- top_genes[is.na(lbl) | lbl == ""]
  rownames(vst_z) <- make.unique(lbl)

  anno_col <- data.frame(group = colData(vsd)$group,
                         row.names = colnames(vst_z))
  ann_colours <- list(group = c(aggressive = "#c0392b", tame = "#2980b9"))

  png(file.path(figs, "06_top_genes_heatmap.png"),
      width = 2200, height = 2400, res = 220)
  pheatmap(vst_z,
           annotation_col    = anno_col,
           annotation_colors = ann_colours,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           main = paste0("Top ", nrow(vst_z),
                         " DEGs (padj < 0.05), VST z-score"))
  dev.off()
} else {
  message("Not enough significant genes for a top-genes heatmap.")
}

message("DE done. Wrote volcano, MA plot, top-genes heatmap, and DE tables.")
