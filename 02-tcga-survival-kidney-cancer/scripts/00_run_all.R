# 00_run_all.R -----------------------------------------------------------
root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
scripts <- file.path(root, "scripts")

message("\n==== 01: download TCGA-KIRC ====")
source(file.path(scripts, "01_download_tcga.R"))
message("\n==== 02: preprocess ====")
source(file.path(scripts, "02_preprocess.R"))
message("\n==== 03: Kaplan-Meier ====")
source(file.path(scripts, "03_km_analysis.R"))
message("\n==== 04: Cox regression ====")
source(file.path(scripts, "04_cox_regression.R"))

message("\nAll done.")
