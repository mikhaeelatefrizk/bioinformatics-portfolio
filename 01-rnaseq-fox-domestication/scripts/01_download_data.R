# 01_download_data.R ------------------------------------------------------
# Download the Wang et al. 2018 silver fox prefrontal-cortex RNA-seq data
# from GEO (GSE76517) and assemble a clean sample metadata table.
#
# Inputs : none (pulls from GEO over the network)
# Outputs: data_in/GSE76517_*  (series matrix + supplementary count files)
#          data_in/sample_metadata.csv
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
})

# --- paths ---------------------------------------------------------------
project_root <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in <- file.path(project_root, "data_in")
dir.create(data_in, showWarnings = FALSE, recursive = TRUE)

gse_id <- "GSE76517"

# --- 1. Series matrix (phenotype/sample metadata) ------------------------
message("Downloading series matrix for ", gse_id, " ...")
gse <- getGEO(gse_id, destdir = data_in, GSEMatrix = TRUE, getGPL = FALSE)
# getGEO returns a list (one per platform); silver fox study is single-platform
eset <- gse[[1]]

pheno <- Biobase::pData(eset)

# --- 2. Supplementary files (count matrices) -----------------------------
message("Downloading supplementary files for ", gse_id, " ...")
supp_paths <- getGEOSuppFiles(
  gse_id,
  baseDir = data_in,
  makeDirectory = FALSE,
  fetch_files = TRUE
)
message("Supplementary files saved:\n", paste(rownames(supp_paths), collapse = "\n"))

# --- 3. Build a clean sample metadata table ------------------------------
# GEO phenoData is messy; keep only the columns we need and derive `group`
# from the sample title / characteristics columns.
keep_cols <- c("title", "geo_accession", "source_name_ch1",
               grep("characteristics", colnames(pheno), value = TRUE))
meta <- pheno[, intersect(keep_cols, colnames(pheno)), drop = FALSE]

# Derive group (tame vs aggressive) from title or characteristics
title_lc <- tolower(as.character(meta$title))
meta$group <- ifelse(grepl("^tame", title_lc), "tame",
              ifelse(grepl("^aggr", title_lc), "aggressive", NA_character_))

# Derive tissue (prefrontal cortex vs basal forebrain) from title if present
meta$tissue <- ifelse(grepl("prefrontal|pfc", title_lc), "prefrontal_cortex",
               ifelse(grepl("basal|forebrain|bf", title_lc), "basal_forebrain",
                      NA_character_))

write.csv(meta, file.path(data_in, "sample_metadata.csv"), row.names = FALSE)
message("Wrote sample_metadata.csv with ", nrow(meta), " rows.")

message("Done. Inspect data_in/ to see what GEO supplied.")
