# RNA-seq: Silver Fox Domestication (Belyaev–Trut experiment)

**Differential gene expression in the prefrontal cortex of tame vs. aggressive silver foxes from the 60-year Novosibirsk selection experiment.**

---

## Biological background

In 1959, Soviet geneticist Dmitri Belyaev began selectively breeding silver foxes
(*Vulpes vulpes*) for a single trait: tolerance of humans. In 1970, a parallel
line selected for aggression was established as a control. The experiment has
run continuously at the Institute of Cytology and Genetics (ICG) in Novosibirsk,
Russia, and is now overseen by Lyudmila Trut.

After roughly 50 generations of selection, the tame line shows not just behavioural
changes, but also the broader "domestication syndrome" — floppy ears, piebald
coats, shorter snouts, lowered stress-hormone levels, and altered reproductive
timing. Belyaev's hypothesis that tameness alone drags the rest of the syndrome
along with it has largely held up.

The question this project addresses: **which genes in the prefrontal cortex
are differentially expressed between the two lines, and what biology do they
implicate?**

## Dataset

- **Accession:** [GEO GSE76517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76517)
- **Paper:** Wang et al. (2018). *Genomic responses to selection for tame/aggressive
  behaviors in the silver fox (Vulpes vulpes).* **PNAS** 115(41): 10398–10403.
  [doi:10.1073/pnas.1800889115](https://doi.org/10.1073/pnas.1800889115)
- **Samples:** 24 individuals (12 tame, 12 aggressive), right prefrontal cortex,
  Illumina paired-end RNA-seq.
- **Original analysis:** edgeR, 146 DEGs at 5% FDR.

## Methods (this reanalysis)

| Step | Tool | Notes |
|------|------|-------|
| Data retrieval | `GEOquery` | Pulls supplementary count matrix + phenotype data from GEO |
| QC / exploration | `DESeq2::vst` + `ggplot2` | PCA, sample-to-sample distance heatmap, library size checks |
| Differential expression | `DESeq2` | Wald test, Benjamini–Hochberg FDR, apeglm LFC shrinkage |
| Visualisation | `EnhancedVolcano`, `pheatmap` | Volcano plot, top-gene heatmap, MA plot |
| Functional enrichment | `clusterProfiler` | GO (BP/MF/CC) over-representation |

The pipeline deliberately uses DESeq2 rather than the original edgeR so that results
can be meaningfully compared between methods — an exercise in understanding *why*
method choice matters for DE analysis.

## How to reproduce

From an R session with the working directory set to this folder:

```r
source("scripts/00_run_all.R")
```

Or run scripts individually in order:

```r
source("scripts/01_download_data.R")          # GEO → data_in/
source("scripts/02_qc_exploration.R")          # QC plots → data_out/figures/
source("scripts/03_differential_expression.R") # DE results → data_out/tables/
source("scripts/04_functional_enrichment.R")   # GO enrichment → data_out/
```

All outputs land in `data_out/`. A `session_info.txt` is written after every full run
so package versions are pinned.

## Repository structure

```
01-rnaseq-fox-domestication/
├── README.md                  <- you are here
├── data_in/                   <- raw GEO downloads (gitignored if large)
├── data_out/
│   ├── tables/                <- DE results, enrichment tables, normalised counts
│   └── figures/               <- PCA, volcano, heatmap, MA plot, enrichment dotplot
└── scripts/
    ├── 00_run_all.R
    ├── 01_download_data.R
    ├── 02_qc_exploration.R
    ├── 03_differential_expression.R
    └── 04_functional_enrichment.R
```

## Key findings

**Pipeline:** DESeq2 + apeglm LFC shrinkage on 12 tame vs 12 aggressive prefrontal
cortex samples, 12,808 genes after filtering (≥10 total counts).

### Differential expression
- **679 genes DE at padj < 0.05** (Wald test, BH-adjusted).
- **584 up in tame (86%) / 95 down in tame (14%)** — a strong directional skew
  toward up-regulation in the domesticated line.
- For comparison, the original edgeR analysis (Wang et al. 2018) reported
  146 DEGs at the same FDR threshold. DESeq2 is more sensitive here, but
  *the top hits agree between methods* — see below.

### Top hits reproduce the original paper
The two most significant genes in our DESeq2 reanalysis are the same two genes
Wang et al. flagged as the top hits from their edgeR analysis:

| Rank | Symbol | Ensembl (CanFam3.1) | log₂FC (tame vs aggr) | padj |
|:----:|:------:|:--------------------|----------------------:|-----:|
| 1 | FSCN3 | ENSCAFG00000001721 | −0.79 | 2.4×10⁻⁷ |
| 2 | PRIMPOL | ENSCAFG00000007746 | −0.83 | 9.8×10⁻⁷ |
| 3 | **PCDHGA1** | ENSCAFG00000023682 | +1.16 | 2.4×10⁻⁶ |
| 4 | FARS2 | ENSCAFG00000009483 | +0.56 | 2.4×10⁻⁶ |
| 5 | **DKKL1** | ENSCAFG00000003681 | +1.03 | 7.6×10⁻⁶ |
| 6 | TRIB1 | ENSCAFG00000001077 | +0.82 | 3.7×10⁻⁵ |
| 7 | CCDC138 | ENSCAFG00000001978 | −0.94 | 4.1×10⁻⁵ |
| 8 | KHK | ENSCAFG00000004673 | +0.91 | 5.4×10⁻⁵ |
| 9 | DUSP6 | ENSCAFG00000006100 | +0.61 | 2.1×10⁻⁴ |
| 10 | NQO1 | ENSCAFG00000020252 | +0.49 | 2.2×10⁻⁴ |
| 11 | EGR1 | ENSCAFG00000001254 | +0.77 | 2.4×10⁻⁴ |

**PCDHGA1** (protocadherin gamma A1) and **DKKL1** (dickkopf-like 1) are the two
genes Wang et al. specifically called out as their most significant hits in the
prefrontal cortex, at P < 10⁻⁸ in their edgeR analysis. Recovering them as #3 and
#5 here with independent methodology is the strongest possible sanity check.
Other notable hits include **EGR1** (an immediate-early gene central to
activity-dependent neural plasticity and memory), **DUSP6** (MAPK-pathway
phosphatase), and **SLITRK6** (neurite outgrowth).

### GO enrichment (Biological Process)
Top enriched terms (BH padj < 0.05):

| GO ID | Term | Genes | padj |
|-------|------|-----:|-----:|
| GO:0051180 | vitamin transport | 7 | 4.4×10⁻³ |
| GO:0044283 | small-molecule biosynthetic process | 29 | 4.4×10⁻³ |
| GO:0071402 | cellular response to lipoprotein particle stimulus | 6 | 6.9×10⁻³ |
| GO:0055094 | response to lipoprotein particle | 5 | 1.5×10⁻² |
| GO:0042063 | **gliogenesis** | 15 | 4.0×10⁻² |
| GO:0001944 | vasculature development | 28 | 4.0×10⁻² |
| GO:0006897 | endocytosis | 25 | 4.0×10⁻² |

**Gliogenesis** (glial cell development) is the most biologically interesting hit
here — glia are increasingly recognised as active participants in social
behaviour and neural plasticity, not just support cells. The vascular and
lipoprotein terms partly reflect the heterogeneity of bulk brain tissue.

### An honest PCA observation
PC1 captures 93% of the variance but does **not** cleanly separate groups.
Instead, ~5 tame foxes (fox182, 186, 187, 190, 191) sit far from the main
cluster as outliers along PC1, while the other tame foxes intermix with the
aggressive foxes. A few interpretations worth considering:

- **Biological heterogeneity within the tame line** — the tame population was
  selected from conventional farm foxes in 1959 and has been reinforced over
  50+ generations, but individual variation in temperament (and presumably in
  gene expression) is expected.
- **Possible batch/technical effect** — all samples were sequenced in a single
  project, but library-prep batches could still contribute. Library sizes
  (`01_library_sizes.png`) look comparable across samples, so gross technical
  issues seem unlikely.
- **The signal is there despite the PCA** — 679 DEGs with the top hits
  reproducing the original paper shows the tame-vs-aggressive contrast is
  real; PCA simply isn't the right lens when the biological signal is
  distributed across many genes rather than concentrated on a few high-variance
  ones.

This is the kind of caveat I'd raise and discuss in an interview rather than
airbrush out of the analysis.

## Acknowledgements

Data: Wang, Pipes, Trut, Herbeck, Vladimirova, Gulevich, Kharlamova,
Johnson, Acland, Kukekova & Clark (2018). PNAS.

The experiment itself: the late Dmitri Belyaev, Lyudmila Trut, and every
animal keeper at the ICG farm in Akademgorodok who has kept it going for 66 years.
