# 01_download_data.R -----------------------------------------------------
# Download the 10x Genomics PBMC 3k dataset (cellranger filtered output).
#
# Outputs: data_in/pbmc3k_filtered_gene_bc_matrices.tar.gz
#          data_in/filtered_gene_bc_matrices/hg19/{barcodes,genes,matrix}
# ------------------------------------------------------------------------

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in <- file.path(root, "data_in")
dir.create(data_in, showWarnings = FALSE, recursive = TRUE)

url      <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
tar_path <- file.path(data_in, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
unpacked <- file.path(data_in, "filtered_gene_bc_matrices", "hg19")

if (!file.exists(file.path(unpacked, "matrix.mtx"))) {
  if (!file.exists(tar_path)) {
    message("Downloading PBMC 3k from 10x Genomics ...")
    options(timeout = 600)
    download.file(url, tar_path, mode = "wb")
  }
  message("Unpacking ...")
  untar(tar_path, exdir = data_in)
} else {
  message("Data already present, skipping.")
}

stopifnot(file.exists(file.path(unpacked, "matrix.mtx")))
message("Done. Matrix at: ", unpacked)
