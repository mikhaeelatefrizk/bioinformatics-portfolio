# install_dependencies.R -------------------------------------------------
# One-time setup script. Installs every R package required to reproduce
# all three projects in this portfolio. Run once after cloning the repo.
#
# Usage:   source("install_dependencies.R")
# Time:    ~10-30 min on a fresh machine (Bioconductor installs from source)
# ------------------------------------------------------------------------

message("Setting up CRAN + Bioconductor ...")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# CRAN packages
cran_pkgs <- c(
  "rprojroot", "ggplot2", "ggrepel", "patchwork", "dplyr", "tidyr",
  "pheatmap", "RColorBrewer", "Seurat", "survival", "survminer"
)
to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(to_install)) {
  message("Installing CRAN: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_pkgs <- c(
  "DESeq2", "GEOquery", "apeglm", "EnhancedVolcano",       # Project 1
  "clusterProfiler", "AnnotationDbi", "biomaRt",
  "org.Cf.eg.db", "org.Hs.eg.db",                          # annotation
  "TCGAbiolinks", "SummarizedExperiment",                  # Project 2
  "glmGamPoi", "SingleCellExperiment"                      # Project 3
)
to_install <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(to_install)) {
  message("Installing Bioconductor: ", paste(to_install, collapse = ", "))
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

message("\nAll dependencies installed.")
message("Run each project with: source('<project-folder>/scripts/00_run_all.R')")
