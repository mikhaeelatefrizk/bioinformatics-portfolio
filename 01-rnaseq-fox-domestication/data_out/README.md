# data_out/

Everything produced by the analysis pipeline — cleaned, versioned, interview-ready.

## tables/

| File | Produced by | Contents |
|------|-------------|----------|
| `normalised_counts.csv` | `02_qc_exploration.R` | DESeq2 size-factor-normalised counts |
| `vst_counts.csv` | `02_qc_exploration.R` | Variance-stabilising-transformed counts for PCA/heatmaps |
| `de_results_full.csv` | `03_differential_expression.R` | All genes, log2FC, p-value, padj, base mean |
| `de_results_significant.csv` | `03_differential_expression.R` | padj < 0.05 only |
| `go_enrichment_BP.csv` | `04_functional_enrichment.R` | GO Biological Process over-representation |
| `go_enrichment_MF.csv` | `04_functional_enrichment.R` | GO Molecular Function |
| `go_enrichment_CC.csv` | `04_functional_enrichment.R` | GO Cellular Component |

## figures/

| File | Produced by | What it shows |
|------|-------------|---------------|
| `01_library_sizes.png` | QC | Per-sample total counts (sanity check for coverage) |
| `02_pca.png` | QC | PC1 vs PC2 of VST-transformed counts, coloured by group |
| `03_sample_distance_heatmap.png` | QC | Euclidean distance between samples — do replicates cluster? |
| `04_volcano.png` | DE | log2FC vs. −log10(padj), with top genes labelled |
| `05_ma_plot.png` | DE | log2FC vs. mean expression (after apeglm shrinkage) |
| `06_top_genes_heatmap.png` | DE | Z-scored VST counts of top 30 DEGs across samples |
| `07_go_dotplot.png` | Enrichment | Top enriched GO BP terms |

## session_info.txt

Output of `sessionInfo()` from the final run — pins R version and every package
version used, so results are reproducible.
