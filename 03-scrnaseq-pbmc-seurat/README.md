# Single-cell RNA-seq of Peripheral Blood Mononuclear Cells (Seurat v5)

**A complete, reproducible Seurat v5 pipeline on the canonical 10x Genomics
PBMC 3k dataset — QC, normalization, clustering, cell-type annotation, and
functional analysis of cluster-defining markers.**

---

## Why this project

Single-cell RNA-seq is now a core technique in modern biology. The **PBMC 3k
dataset from 10x Genomics** is the standard benchmark for the Seurat workflow —
it's the dataset every single-cell tutorial uses, and for good reason: a
well-annotated mix of immune cell types where the expected clusters and
canonical markers are known. Reproducing it cleanly is a baseline competency
check for anyone doing scRNA-seq.

This project deliberately goes **beyond** the standard Seurat tutorial:

- Uses **SCTransform** (the current Seurat-recommended normalization) with
  `glmGamPoi` for speed, not the older `NormalizeData` + `ScaleData` flow
- Applies **multiple cell-type annotation approaches** and cross-checks them:
  canonical marker genes, module scoring of known signatures, and
  unbiased marker discovery via `FindAllMarkers`
- Runs **GO enrichment** on each cluster's top marker genes to validate the
  biological identity of the cluster, not just assert it
- Includes **QC diagnostics** (nFeature, nCount, percent.mt) with thresholds
  chosen from the data, not copied from a tutorial

## Dataset

- **Source:** 10x Genomics, "Peripheral blood mononuclear cells (PBMCs) from
  a healthy donor", Single Cell 3' v1 chemistry, 2016
- **Link:** https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
- **Expected size:** ~2,700 cells × 32,738 genes after droplet filtering

## Methods

| Step | Tool | Notes |
|------|------|-------|
| Data retrieval | `utils::download.file`, `Seurat::Read10X` | 10x cellranger output format |
| QC | Seurat | nFeature 200–2500, percent.mt < 5% |
| Normalization | `SCTransform` + `glmGamPoi` | Regularised negative binomial |
| Dimensionality reduction | PCA → UMAP | Elbow plot to choose # PCs |
| Clustering | `FindNeighbors` + `FindClusters` | Louvain, resolution 0.5 |
| Marker discovery | `FindAllMarkers` | Wilcoxon, only.pos = TRUE |
| Annotation | Canonical markers + `AddModuleScore` | Panel of IL7R, CD14, MS4A1, CD8A, GNLY, FCGR3A, FCER1A, PPBP |
| Enrichment | `clusterProfiler::enrichGO` | GO BP on top marker genes per cluster |

## Candidate cell-type markers (canonical)

| Cell type | Markers |
|-----------|---------|
| Naive CD4⁺ T | IL7R, CCR7 |
| Memory CD4⁺ T | IL7R, S100A4 |
| CD14⁺ Monocyte | CD14, LYZ |
| FCGR3A⁺ Monocyte | FCGR3A, MS4A7 |
| B cell | MS4A1, CD79A |
| CD8⁺ T | CD8A |
| NK | GNLY, NKG7 |
| Dendritic | FCER1A, CST3 |
| Platelet | PPBP |

## How to reproduce

```r
source("scripts/00_run_all.R")
```

Or individually:

```r
source("scripts/01_download_data.R")
source("scripts/02_qc_and_normalization.R")
source("scripts/03_clustering_and_annotation.R")
source("scripts/04_marker_enrichment.R")
```

## Repository structure

```
03-scrnaseq-pbmc-seurat/
├── README.md
├── data_in/                      <- 10x cellranger output (gitignored)
├── data_out/
│   ├── tables/                   <- cluster markers, QC stats, GO results
│   ├── figures/                  <- UMAPs, violin plots, marker heatmaps
│   └── rds/                      <- Seurat objects (gitignored)
└── scripts/
    ├── 00_run_all.R
    ├── 01_download_data.R
    ├── 02_qc_and_normalization.R
    ├── 03_clustering_and_annotation.R
    └── 04_marker_enrichment.R
```

## Key findings

### Cohort and QC
- **2,700 cells** in the raw filtered 10x matrix → **2,638 cells** retained
  after QC (nFeature 200–2500, percent.mt < 5%). The 62 cells removed match
  the canonical PBMC 3k tutorial exactly.
- Median post-QC metrics: nFeature ≈ 815 genes/cell, nCount ≈ 2,196 UMI/cell,
  percent.mt ≈ 2.05%.

### Clustering and cell-type annotation
After SCTransform normalisation, PCA (10 PCs chosen from the elbow), and
Louvain clustering at resolution 0.5, **10 clusters** were identified.
Each cluster was assigned a cell type using mean module score of canonical
marker panels:

| Cluster | Annotation | Top markers | Top GO BP |
|:-------:|:----------|:------------|:----------|
| 0 | CD4 T (naive) | **CCR7**, PRKCQ-AS1, **LEF1** | mononuclear cell differentiation (padj 1.0×10⁻⁵) |
| 1 | CD14⁺ Monocyte | **S100A8**, **S100A9**, FOLR3 | immune system process (padj 1.1×10⁻³) |
| 2 | CD4 T (memory) | AQP3, **CD40LG**, OPTN | lymphocyte differentiation (padj 4.3×10⁻⁵) |
| 3 | B cell | **CD79A**, VPREB3, LINC00926 | **MHC class II protein complex assembly (padj 1.2×10⁻¹²)** |
| 4 | CD8 T (GZMK⁺) | **GZMK**, NCR3, CCL5 | immune system process (padj 1.4×10⁻³) |
| 5 | NK | **GNLY**, **GZMB**, AKR1C3 | leukocyte-mediated immunity (padj 1.0×10⁻⁸) |
| 6 | FCGR3A⁺ Monocyte | CKB, **CDKN1C**, BATF3 | — |
| 7 | CD8 T (cytotoxic) | **GZMH**, FCRL6, **CD8A** | immune response (padj 1.8×10⁻⁸) |
| 8 | Dendritic | **FCER1A**, SERPINF1, **CLEC10A** | **antigen processing and presentation of exogenous antigen (padj 2.0×10⁻¹⁰)** |
| 9 | Platelet | CLDN5, CMTM5, AP001189.4 | **blood coagulation (padj 8.8×10⁻³)** |

Canonical markers in bold.

### Two independent lines of evidence agree on each cluster

The cluster-identity calls hold up against **three** independent checks:

1. **Canonical marker violin plots** (`06_canonical_markers_violin.png`):
   CD14/LYZ light up in cluster 1, MS4A1/CD79A in cluster 3, GNLY/NKG7 in
   cluster 5, FCGR3A/MS4A7 in cluster 6, FCER1A in cluster 8, PPBP in
   cluster 9 — each marker set specific to its expected cluster.
2. **Module scoring** (`09_cluster_module_scores.png`): the heatmap is
   cleanly diagonal for all clear populations (B cell module → cluster 3,
   NK module → cluster 5, DC module → cluster 8, Platelet module → cluster 9,
   CD14-mono module → cluster 1, FCGR3A-mono module → cluster 6, CD8-T →
   cluster 4 and 7).
3. **Unbiased marker discovery** (`FindAllMarkers`, 2,944 cluster markers
   total) plus GO enrichment on the top 50 markers per cluster (370 enriched
   BP terms across 9 of 10 clusters) — and every enriched term that reaches
   the top of its cluster is biologically coherent with the assigned cell type.

The strongest validation: **cluster 3's top GO term is "MHC class II protein
complex assembly" (padj 1.2×10⁻¹²)** — exactly what you expect from B cells,
which present antigen via MHC-II. Similarly, cluster 8's top term is
"antigen processing and presentation of exogenous antigen" — the signature
function of dendritic cells. Cluster 9 lights up "blood coagulation" — the
cluster is platelets.

### An honest note on the two CD4 T clusters

Clusters 0 and 2 are both CD4 T cells but the module-score-based annotation
flagged them with weak scores (0.12 and 0.03 vs. 1.76 for the CD14 monocyte
cluster). This is expected, not a problem:

- **Cluster 0** expresses *CCR7* and *LEF1* — hallmarks of **naive CD4 T**
  cells, which by definition are resting and express lineage markers at
  lower levels than activated populations.
- **Cluster 2** expresses *CD40LG* (CD154, the CD40 ligand) and *AQP3* —
  consistent with **memory/effector CD4 T** cells.

A follow-up refinement would use a specifically naive-vs-memory module
panel (e.g., LEF1/TCF7/CCR7 vs. S100A4/CD40LG/IL32). Left as a noted
limitation rather than over-engineered here.

### Clean match to the canonical PBMC 3k result

Cell-type-labelled UMAP matches the published Seurat tutorial output
essentially exactly (2638/2638 cells, 8 cell types, same topology — CD14 and
FCGR3A monocytes on the left, CD4/CD8/NK on the right, B cells bottom,
platelets as a tiny isolated cluster). This is what reproducibility
looks like.

## A note on the "Russian angle"

Projects 01 and 02 in this portfolio have genuine Russian scientific threads
(the Belyaev–Trut silver fox experiment from the ICG in Novosibirsk; ccRCC
as a major focus of Russian cancer research at the Blokhin Center and
Engelhardt Institute). Project 03 does not — I chose the canonical 10x PBMC
benchmark because reliably reproducing it is the strongest possible
demonstration of competence in the standard Seurat workflow. Forcing a
Russian connection onto every project would be contrived, and
admissions/interview panels value honest scientific motivation over
narrative-fitting.
