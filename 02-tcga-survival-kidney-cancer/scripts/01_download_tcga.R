# 01_download_tcga.R -----------------------------------------------------
# Download TCGA-KIRC (Kidney Renal Clear Cell Carcinoma) RNA-seq counts and
# clinical data from the NCI Genomic Data Commons via TCGAbiolinks.
#
# Outputs: data_in/GDCdata/         (GDC cache, gitignored)
#          data_in/tcga_kirc_se.rds (assembled SummarizedExperiment)
#          data_in/tcga_kirc_clinical.csv
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in <- file.path(root, "data_in")
dir.create(data_in, showWarnings = FALSE, recursive = TRUE)

# Keep GDC cache inside data_in so it's self-contained
setwd(data_in)

se_path <- file.path(data_in, "tcga_kirc_se.rds")

if (!file.exists(se_path)) {
  message("Querying GDC for TCGA-KIRC STAR gene counts ...")
  query_rna <- GDCquery(
    project          = "TCGA-KIRC",
    data.category    = "Transcriptome Profiling",
    data.type        = "Gene Expression Quantification",
    workflow.type    = "STAR - Counts",
    sample.type      = c("Primary Tumor")
  )

  message("Downloading (this takes a while the first time) ...")
  GDCdownload(query_rna, method = "api", files.per.chunk = 50)

  message("Preparing SummarizedExperiment ...")
  se <- GDCprepare(query_rna, summarizedExperiment = TRUE)
  saveRDS(se, se_path)
  message("Saved SummarizedExperiment: ", se_path)
} else {
  message("SummarizedExperiment cache found, skipping download.")
  se <- readRDS(se_path)
}

message("Samples: ", ncol(se), "  |  Genes: ", nrow(se))

# Clinical: pull the indexed clinical supplement (has follow-up data)
clin_path <- file.path(data_in, "tcga_kirc_clinical.csv")
if (!file.exists(clin_path)) {
  message("Querying GDC for clinical data ...")
  clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")
  write.csv(clinical, clin_path, row.names = FALSE)
  message("Saved clinical table: ", clin_path, " (", nrow(clinical), " patients)")
} else {
  message("Clinical cache found.")
}

message("\nDownload step complete.")
