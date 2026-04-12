# data_in/

TCGA data from the NCI Genomic Data Commons (GDC), downloaded via
`TCGAbiolinks::GDCdownload()`.

## Structure

```
data_in/
└── GDCdata/
    └── TCGA-KIRC/
        ├── Transcriptome_Profiling/
        │   └── Gene_Expression_Quantification/
        │       └── ...STAR count files (gitignored)
        └── (clinical indexed supplement)
```

The raw GDC downloads (~1–2 GB) are gitignored. Re-run
`scripts/01_download_tcga.R` to reproduce.

After download, the assembled SummarizedExperiment is cached to
`data_in/tcga_kirc_se.rds` for fast reloading.
