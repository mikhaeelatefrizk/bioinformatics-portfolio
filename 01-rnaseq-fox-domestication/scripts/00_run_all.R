# 00_run_all.R ------------------------------------------------------------
# Run the whole pipeline in order. Assumes working directory is the project
# root (or any subfolder -- rprojroot will find it).
# --------------------------------------------------------------------------
root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
scripts <- file.path(root, "scripts")

message("\n==== 01: download data ====")
source(file.path(scripts, "01_download_data.R"))
message("\n==== 02: QC / exploration ====")
source(file.path(scripts, "02_qc_exploration.R"))
message("\n==== 03: differential expression ====")
source(file.path(scripts, "03_differential_expression.R"))
message("\n==== 04: functional enrichment ====")
source(file.path(scripts, "04_functional_enrichment.R"))

message("\nAll done.")
