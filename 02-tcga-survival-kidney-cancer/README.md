# TCGA Survival Analysis: Kidney Renal Clear Cell Carcinoma (KIRC)

**Survival analysis of clear-cell renal cell carcinoma using TCGA-KIRC, testing
whether a panel of known and candidate prognostic genes stratify overall survival.**

---

## Background

Clear-cell renal cell carcinoma (ccRCC) is the most common kidney cancer,
accounting for roughly 75% of renal cancers. It is driven largely by biallelic
loss of the *VHL* tumour-suppressor gene, which stabilises HIF-1α/HIF-2α and
rewires cell metabolism and angiogenesis. ccRCC is an important model cancer for
metabolic dysregulation, and has long been a focus of Russian oncology research
— notably at the N.N. Blokhin National Medical Research Center of Oncology in
Moscow and the Engelhardt Institute of Molecular Biology.

The question: **can expression levels of individual genes (measured in the primary
tumour at resection) predict how long a patient survives?**

## Dataset

- **Cohort:** TCGA-KIRC (The Cancer Genome Atlas, Kidney Renal Clear Cell Carcinoma)
- **Access:** Programmatic via `TCGAbiolinks` from the NCI Genomic Data Commons
- **Samples:** ~530 primary tumours with matched clinical follow-up
- **Data types used:** STAR-aligned gene counts (RNA-seq), clinical (overall survival, stage, age, sex)

## Candidate genes

Hypothesis-driven panel of canonical ccRCC markers:

| Gene | Rationale |
|------|-----------|
| **CA9** | HIF target, classical ccRCC marker used clinically |
| **VEGFA** | Driver of ccRCC angiogenesis, target of anti-angiogenic therapy |
| **MKI67** | Proliferation marker (Ki-67) — clinically used prognostic marker |
| **BAP1** | Common tumour-suppressor loss in ccRCC, associated with poor prognosis |
| **SETD2** | Chromatin modifier, loss associated with poor prognosis |
| **PBRM1** | PBAF complex member, loss common in ccRCC |
| **HIF1A** | Master regulator of hypoxic response, central to ccRCC biology |
| **EPAS1** | HIF-2α, the dominant HIF in ccRCC |

## Methods

| Step | Tool | Notes |
|------|------|-------|
| Data retrieval | `TCGAbiolinks` | Queries GDC, downloads + assembles SummarizedExperiment |
| Preprocessing | `DESeq2::vst` | VST-normalised expression for downstream modelling |
| Kaplan-Meier | `survival`, `survminer` | Median-split on each candidate gene, log-rank test |
| Multivariate modelling | `survival::coxph` | Cox PH adjusting for age, sex, pathologic stage |
| PH diagnostics | `survival::cox.zph` | Schoenfeld residual test of proportional-hazards assumption |

## How to reproduce

```r
source("scripts/00_run_all.R")
```

Or individually:

```r
source("scripts/01_download_tcga.R")       # GDC → data_in/
source("scripts/02_preprocess.R")           # merged survival+expression → data_out/rds/
source("scripts/03_km_analysis.R")          # KM curves + log-rank per gene
source("scripts/04_cox_regression.R")       # univariate + multivariate Cox
```

The first run downloads ~1–2 GB of TCGA data into `data_in/GDCdata/`; this is
gitignored. Subsequent runs use the local cache.

## Repository structure

```
02-tcga-survival-kidney-cancer/
├── README.md
├── data_in/                     <- GDC download cache (gitignored)
├── data_out/
│   ├── tables/                  <- KM log-rank, Cox coefficients, Schoenfeld tests
│   ├── figures/                 <- KM curves per gene + forest plot
│   └── rds/                     <- merged analysis-ready data
└── scripts/
    ├── 00_run_all.R
    ├── 01_download_tcga.R
    ├── 02_preprocess.R
    ├── 03_km_analysis.R
    └── 04_cox_regression.R
```

## Key findings

### Cohort
529 ccRCC primary tumours with usable survival data, 173 deaths (32.7% event
rate), median follow-up 1,191 days (~3.3 years). Stage distribution: I (230),
II (53), III (102), IV (56), missing (88). Male/female: 344/185. Median age
61 years.

### Univariate Cox regression (gene alone, HR per unit VST)

| Gene | HR | 95% CI | p (Wald) | padj (BH) | Concordance |
|------|----:|:---:|-------:|------:|:---:|
| **EPAS1** | **0.66** | 0.58–0.75 | **6.3×10⁻¹¹** | **5.1×10⁻¹⁰** | 0.64 |
| **MKI67** | 1.40 | 1.21–1.63 | 6.8×10⁻⁶ | 2.7×10⁻⁵ | 0.60 |
| PBRM1 | 0.71 | 0.55–0.92 | 0.011 | 0.028 | 0.55 |
| SETD2 | 0.74 | 0.55–1.00 | 0.053 | 0.11 | 0.55 |
| BAP1 | 0.86 | 0.61–1.22 | 0.41 | 0.65 | 0.50 |
| HIF1A | 1.07 | 0.85–1.34 | 0.57 | 0.69 | 0.52 |
| CA9 | 0.98 | 0.91–1.05 | 0.60 | 0.69 | 0.52 |
| VEGFA | 1.01 | 0.88–1.15 | 0.92 | 0.92 | 0.51 |

### Multivariate Cox (adjusted for age, sex, AJCC stage)

| Gene | HR | 95% CI | p | padj (BH) | Concordance |
|------|----:|:---:|-------:|------:|:---:|
| **EPAS1** | **0.73** | 0.61–0.89 | **0.0013** | **0.011** | 0.77 |
| HIF1A | 1.22 | 0.96–1.55 | 0.10 | 0.41 | 0.77 |
| VEGFA | 1.13 | 0.94–1.35 | 0.19 | 0.45 | 0.77 |
| MKI67 | 1.12 | 0.94–1.33 | 0.22 | 0.45 | 0.77 |

**EPAS1 (HIF-2α) is the only gene that remains a significant independent
predictor of overall survival after adjustment for age, sex, and AJCC stage.**

### Biological interpretation

- **EPAS1 (HIF-2α)** is the dominant hypoxia-inducible factor in ccRCC — *VHL*
  loss stabilises HIF-2α, which then drives tumour angiogenesis and growth.
  It is the direct target of **belzutifan**, an FDA-approved HIF-2α inhibitor
  for ccRCC (approved Aug 2021, expanded 2023). The protective association
  seen here (high EPAS1 → better survival, HR 0.73 adjusted) may seem
  counterintuitive for a cancer driver, but the clinical and translational
  literature shows that ccRCC tumours split into HIF-1α-high vs HIF-2α-high
  subgroups with distinct prognoses — HIF-2α-high ccRCC tends to be the
  classical, well-differentiated, less aggressive form.

- **MKI67 (Ki-67)** is a proliferation marker and behaves as expected — high
  expression predicts worse survival univariately, but the signal is largely
  absorbed by stage once we adjust, because advanced-stage tumours proliferate
  more aggressively by definition.

- **PBRM1** and **SETD2** are *tumour suppressors* commonly lost in ccRCC.
  Higher expression is associated with better survival univariately, consistent
  with their roles, but the multivariate effect is not significant — here
  expression may be a noisier readout than mutation status (not tested here).

- **Adding clinical covariates dramatically improves model concordance**
  (0.60–0.64 → 0.77 across models) — stage remains the single most powerful
  predictor of ccRCC survival, as expected clinically.

### Proportional-hazards diagnostics

Schoenfeld residual tests (`data_out/tables/schoenfeld_tests.csv`) show that
in the EPAS1 multivariate model the gene itself comfortably satisfies PH
(p = 0.82), as do age (p = 0.28) and sex (p = 0.34). The global test is
mildly violated (p = 0.022), driven entirely by **stage** (p = 0.005) — a
well-known issue in ccRCC, since early-stage patients have nearly-flat hazards
while late-stage hazards rise sharply early. A fully rigorous follow-up would
stratify on stage rather than include it as a covariate. This is left
as a note for an interview discussion rather than airbrushed out of the
analysis.

### Comparison with the original TCGA-KIRC publication
The Cancer Genome Atlas Research Network (2013, *Nature* 499:43–49) reported
that loss-of-function mutations in *PBRM1*, *BAP1*, and *SETD2* are recurrent
in ccRCC, with *BAP1*-mutant tumours having the worst prognosis. Here we test
*expression* rather than mutation status, which is a strictly weaker proxy —
a gene can be inactivated by mutation while still producing transcripts. That
our expression-based survival associations for *PBRM1* and *SETD2* are weaker
than reported for mutation status is expected and consistent with this
distinction.
