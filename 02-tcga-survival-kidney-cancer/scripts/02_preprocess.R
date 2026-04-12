# 02_preprocess.R --------------------------------------------------------
# Build an analysis-ready table: VST-normalised expression of candidate
# genes + survival outcomes + clinical covariates, one row per patient.
#
# Inputs : data_in/tcga_kirc_se.rds
# Outputs: data_out/rds/analysis_data.rds
#          data_out/tables/clinical_summary.csv
#          data_out/figures/04_cohort_overview.png
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

root     <- rprojroot::find_root(rprojroot::has_file("README.md"))
data_in  <- file.path(root, "data_in")
tables   <- file.path(root, "data_out", "tables")
figs     <- file.path(root, "data_out", "figures")
rds_dir  <- file.path(root, "data_out", "rds")

se <- readRDS(file.path(data_in, "tcga_kirc_se.rds"))
message("Loaded SE: ", nrow(se), " genes x ", ncol(se), " samples")

# Candidate genes (hypothesis-driven, canonical ccRCC markers)
gene_panel <- c("CA9", "VEGFA", "MKI67", "BAP1", "SETD2",
                "PBRM1", "HIF1A", "EPAS1")

# --- 1. De-duplicate samples to one per patient -------------------------
# TCGA barcodes: TCGA-XX-XXXX-01A-... -> first 12 chars = patient
colData(se)$patient_id <- substr(colnames(se), 1, 12)
# Keep first sample per patient
first_idx <- !duplicated(colData(se)$patient_id)
se <- se[, first_idx]
message("After de-dup: ", ncol(se), " unique patients")

# --- 2. VST transform for modelling -------------------------------------
counts <- assay(se, "unstranded")
mode(counts) <- "integer"
# DESeq2 needs integers and non-negative
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = colData(se),
                              design    = ~ 1)
# Filter very low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message("Genes after filtering: ", nrow(dds))

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

# --- 3. Subset to candidate genes ---------------------------------------
rd <- rowData(se)
gene_to_ensembl <- rd$gene_id[match(gene_panel, rd$gene_name)]
names(gene_to_ensembl) <- gene_panel
# Keep only genes that survived filtering
available <- gene_to_ensembl[gene_to_ensembl %in% rownames(vst_mat)]
message("Candidate genes available: ",
        paste(names(available), collapse = ", "))

expr_panel <- t(vst_mat[available, , drop = FALSE])
colnames(expr_panel) <- names(available)

# --- 4. Build survival + clinical table ---------------------------------
cd <- as.data.frame(colData(se))

# Overall survival: time to death (Dead) or last follow-up (Alive)
cd$os_time   <- ifelse(cd$vital_status == "Dead",
                       cd$days_to_death,
                       cd$days_to_last_follow_up)
cd$os_event  <- ifelse(cd$vital_status == "Dead", 1L, 0L)

# Collapse stage to ordinal (I/II/III/IV)
cd$stage_simple <- gsub("Stage ", "", cd$ajcc_pathologic_stage)
cd$stage_simple <- gsub("A$|B$|C$", "", cd$stage_simple)
cd$stage_simple <- factor(cd$stage_simple, levels = c("I","II","III","IV"))

keep_cols <- c("patient_id", "os_time", "os_event",
               "age_at_index", "gender", "stage_simple", "vital_status")
cd_small  <- cd[, keep_cols]

# Merge expression + clinical
stopifnot(all(rownames(expr_panel) == rownames(cd_small)))
analysis <- cbind(cd_small, as.data.frame(expr_panel))

# Filter to patients with usable survival data
analysis <- analysis[!is.na(analysis$os_time) & analysis$os_time > 0, ]
message("Patients with valid survival: ", nrow(analysis))
message("Events (deaths): ", sum(analysis$os_event),
        "  (", round(100*mean(analysis$os_event), 1), "%)")

saveRDS(analysis, file.path(rds_dir, "analysis_data.rds"))

# --- 5. Cohort summary ---------------------------------------------------
summary_df <- data.frame(
  variable = c("Total patients",
               "Events (deaths)",
               "Censored",
               "Median follow-up (days)",
               "Median age (years)",
               "Male / Female",
               "Stage I",  "Stage II", "Stage III", "Stage IV",
               "Stage missing"),
  value = c(
    nrow(analysis),
    sum(analysis$os_event),
    sum(analysis$os_event == 0),
    round(median(analysis$os_time, na.rm = TRUE)),
    round(median(analysis$age_at_index, na.rm = TRUE), 1),
    paste(sum(analysis$gender == "male", na.rm = TRUE), "/",
          sum(analysis$gender == "female", na.rm = TRUE)),
    sum(analysis$stage_simple == "I",  na.rm = TRUE),
    sum(analysis$stage_simple == "II", na.rm = TRUE),
    sum(analysis$stage_simple == "III", na.rm = TRUE),
    sum(analysis$stage_simple == "IV", na.rm = TRUE),
    sum(is.na(analysis$stage_simple))
  )
)
write.csv(summary_df, file.path(tables, "clinical_summary.csv"),
          row.names = FALSE)
print(summary_df)

# --- 6. Cohort overview figure ------------------------------------------
p_stage <- analysis %>%
  mutate(stage_simple = ifelse(is.na(stage_simple), "Unknown",
                               as.character(stage_simple))) %>%
  mutate(stage_simple = factor(stage_simple,
                               levels = c("I","II","III","IV","Unknown"))) %>%
  count(stage_simple, vital_status) %>%
  ggplot(aes(x = stage_simple, y = n, fill = vital_status)) +
  geom_col() +
  scale_fill_manual(values = c(Alive = "#2980b9", Dead = "#c0392b")) +
  labs(x = "AJCC pathologic stage", y = "Patients",
       title = "TCGA-KIRC cohort overview",
       subtitle = paste0(nrow(analysis), " patients, ",
                         sum(analysis$os_event), " deaths")) +
  theme_minimal(base_size = 11)

ggsave(file.path(figs, "04_cohort_overview.png"), p_stage,
       width = 7, height = 5, dpi = 150)

message("Preprocessing done.")
