# data_out/

## tables/
| File | Contents |
|------|----------|
| `qc_summary.csv` | Pre- and post-filter cell counts, median nFeature/nCount/percent.mt |
| `cluster_markers.csv` | All cluster-defining genes (Wilcoxon, FindAllMarkers) |
| `cluster_markers_top10.csv` | Top 10 markers per cluster, ranked by avg_log2FC |
| `cluster_annotations.csv` | Cluster ID → cell type label with supporting markers |
| `go_enrichment_by_cluster.csv` | GO BP enrichment on each cluster's top markers |

## figures/
| File | Contents |
|------|----------|
| `01_qc_violin_prefilter.png` | nFeature / nCount / percent.mt before QC |
| `02_qc_scatter.png` | nCount vs nFeature, coloured by percent.mt |
| `03_elbow_plot.png` | PCA elbow plot for choosing # PCs |
| `04_umap_clusters.png` | UMAP coloured by unsupervised cluster |
| `05_umap_annotated.png` | UMAP with biological cell-type labels |
| `06_canonical_markers_violin.png` | Violin plots of canonical markers across clusters |
| `07_canonical_markers_umap.png` | Feature plots of canonical markers on UMAP |
| `08_top_markers_heatmap.png` | Heatmap of top 5 markers per cluster |
| `09_cluster_module_scores.png` | Cell-type module scores per cluster |

## rds/
- `pbmc_annotated.rds` — fully processed, clustered, annotated Seurat object.
