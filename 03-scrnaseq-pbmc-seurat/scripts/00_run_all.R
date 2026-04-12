# 00_run_all.R -----------------------------------------------------------
root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
scripts <- file.path(root, "scripts")

message("\n==== 01: download data ====")
source(file.path(scripts, "01_download_data.R"))
message("\n==== 02: QC + SCTransform ====")
source(file.path(scripts, "02_qc_and_normalization.R"))
message("\n==== 03: clustering + annotation ====")
source(file.path(scripts, "03_clustering_and_annotation.R"))
message("\n==== 04: markers + enrichment ====")
source(file.path(scripts, "04_marker_enrichment.R"))

message("\nAll done.")
