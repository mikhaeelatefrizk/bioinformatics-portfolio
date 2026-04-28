# Bioinformatics Portfolio

Three end-to-end reproducible bioinformatics analyses on public datasets —
bulk RNA-seq differential expression, clinical survival analysis, and
single-cell RNA-seq — each built from raw data to publication-quality figures.

**Author:** Mikhaeel Wahba · Rotterdam · [GitHub](https://github.com/mikhaeelatefrizk) · [mikhaeelatefrizk@proton.me](mailto:mikhaeelatefrizk@proton.me)

---

## The three projects

### 01 · RNA-seq differential expression
**Silver fox domestication — Belyaev / Trut experiment, Novosibirsk**

A reanalysis of the 1959-present Russian silver fox domestication experiment
([Wang et al. 2018, *PNAS*](https://doi.org/10.1073/pnas.1800889115)).
Twelve tame vs. twelve aggressive foxes, prefrontal-cortex RNA-seq,
DESeq2 + apeglm differential expression, clusterProfiler GO enrichment.

- **679 DE genes** at padj < 0.05 (584 up-regulated in the tame line, 95 down)
- **PCDHGA1** (padj 2.4×10⁻⁶) and **DKKL1** (padj 7.6×10⁻⁶) — the two genes
  Wang et al. specifically called out as their top hits — recovered as
  ranks 3 and 5 here with DESeq2 instead of edgeR (independent methodology
  cross-check)
- **16 enriched GO Biological Process terms** including gliogenesis and
  vasculature development
- Dataset: [GEO GSE76517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76517)

📂 [`01-rnaseq-fox-domestication/`](01-rnaseq-fox-domestication/)

### 02 · Clinical survival analysis
**TCGA kidney renal clear cell carcinoma**

Survival modelling of 529 ccRCC patients from The Cancer Genome Atlas
(TCGA-KIRC), testing a hypothesis-driven panel of eight canonical ccRCC
markers for prognostic signal. Kaplan-Meier stratification, univariate and
multivariate Cox proportional-hazards regression, Schoenfeld-residual PH
diagnostics.

- **EPAS1 (HIF-2α) is the only gene that remains a significant independent
  predictor** of overall survival after adjustment for age, sex, and AJCC stage
  (HR 0.73 [0.61–0.89], p = 0.0013, padj = 0.011)
- Clinically meaningful: EPAS1 is the direct target of **belzutifan**, the
  FDA-approved HIF-2α inhibitor for ccRCC
- Model concordance jumps from 0.64 → 0.77 when clinical covariates are added
- PH diagnostics flagged honestly: stage violates PH in every model
  (a known KIRC issue), EPAS1 itself satisfies PH

📂 [`02-tcga-survival-kidney-cancer/`](02-tcga-survival-kidney-cancer/)

### 03 · Single-cell RNA-seq
**Peripheral blood mononuclear cells — Seurat v5 workflow**

A clean, reproducible scRNA-seq pipeline on the canonical 10x Genomics PBMC 3k
benchmark: QC → SCTransform → PCA → UMAP → Louvain clustering → cell-type
annotation by module scoring → unbiased marker discovery → per-cluster GO
enrichment.

- **2,638 cells** retained from 2,700 after QC (matches the Seurat tutorial exactly)
- **10 clusters** cleanly resolved into 8 immune populations:
  CD4 T (naive + memory), CD8 T (GZMK⁺ + cytotoxic), NK, B cell,
  CD14⁺ and FCGR3A⁺ monocytes, dendritic cells, platelets
- Cluster identities validated by three independent lines of evidence:
  canonical markers, module scoring, and unbiased GO enrichment
- **2,944 cluster markers** and **370 GO terms** across 9 clusters
- Strongest validation: B-cell cluster's top GO term is
  "MHC class II protein complex assembly" (padj 1.2×10⁻¹²)

📂 [`03-scrnaseq-pbmc-seurat/`](03-scrnaseq-pbmc-seurat/)

## Project structure

Every project follows the same layout:

```
<project-name>/
├── README.md           # biological context, methods, findings
├── data_in/            # raw inputs (gitignored if large)
├── data_out/
│   ├── tables/         # every CSV referenced in the README
│   ├── figures/        # every figure referenced in the README
│   └── rds/            # large serialised R objects (gitignored)
└── scripts/
    ├── 00_run_all.R
    ├── 01_...
    ├── 02_...
    └── ...
```

Scripts are numbered and meant to be run in order. `00_run_all.R` orchestrates
the whole pipeline from raw download to final figures.

## How to reproduce

Clone the repo, then from an R session at the repo root:

```r
# One-time: install every dependency for every project (~10-30 min)
source("install_dependencies.R")

# Run any project end-to-end
setwd("01-rnaseq-fox-domestication")  ;  source("scripts/00_run_all.R")
setwd("../02-tcga-survival-kidney-cancer")  ;  source("scripts/00_run_all.R")
setwd("../03-scrnaseq-pbmc-seurat")  ;  source("scripts/00_run_all.R")
```

Each project writes a `session_info.txt` after a full run so package versions
are pinned and reproducible.

## Environment

- **OS:** Windows 11
- **R:** 4.5.3
- **Bioconductor:** 3.22
- **Key packages:** DESeq2, Seurat 5, TCGAbiolinks, clusterProfiler,
  survival, survminer, EnhancedVolcano

## A note on methodology

Every analytical choice in this portfolio is documented in its project README:
why DESeq2 rather than edgeR for project 1, why Ensembl 100 for the fox
annotation (genome-assembly mismatch), why module scoring is cross-checked
against unbiased markers in project 3, why stage violates PH in the KIRC
Cox models. Caveats and limitations are flagged inline rather than hidden —
the PCA outliers in project 1, the weak naive-CD4 module scores in project 3,
the expression-vs-mutation-status distinction for *PBRM1* / *SETD2* in
project 2.

This is the approach an admissions committee or technical interviewer will
actually stress-test you on — not whether the volcano plot is pretty, but
whether you understand what you did and why.

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgements

- **Project 1:** Wang, Pipes, Trut, Herbeck, Vladimirova, Gulevich,
  Kharlamova, Johnson, Acland, Kukekova & Clark (2018, *PNAS*). And everyone
  at the Institute of Cytology and Genetics in Novosibirsk who has kept
  the silver fox experiment running for 66 years.
- **Project 2:** The Cancer Genome Atlas Research Network; NCI Genomic Data
  Commons.
- **Project 3:** 10x Genomics for the public PBMC 3k benchmark; the Satija
  lab for the Seurat tutorial that made this dataset the community reference.
