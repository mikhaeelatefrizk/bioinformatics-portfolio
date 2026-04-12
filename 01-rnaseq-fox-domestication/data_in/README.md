# data_in/

Raw inputs, as downloaded — nothing in this folder is edited by hand.

## Contents (after running `scripts/01_download_data.R`)

| File | Source | Description |
|------|--------|-------------|
| `GSE76517_series_matrix.txt.gz` | GEO | Phenotype / sample metadata |
| `GSE76517_*_counts.txt.gz` | GEO supplementary | Raw gene-level count matrix |
| `sample_metadata.csv` | derived | Cleaned sample sheet (24 samples, tame vs. aggressive) |

Large files (`*.gz`, `*.tar`) are gitignored — re-run the download script to reproduce.
