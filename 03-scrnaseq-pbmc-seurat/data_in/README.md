# data_in/

The raw 10x Genomics PBMC 3k cellranger output, downloaded via
`scripts/01_download_data.R`.

## Structure after download

```
data_in/
├── pbmc3k_filtered_gene_bc_matrices.tar.gz   (downloaded archive)
└── filtered_gene_bc_matrices/
    └── hg19/
        ├── barcodes.tsv       (cell barcodes, one per line)
        ├── genes.tsv          (Ensembl ID, gene symbol, one gene per line)
        └── matrix.mtx         (sparse count matrix)
```

The archive and unpacked files are gitignored. Re-run the download script
to reproduce.
