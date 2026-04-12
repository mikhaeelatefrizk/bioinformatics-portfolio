# 03_km_analysis.R -------------------------------------------------------
# Kaplan-Meier analysis for each candidate gene. Samples are split at the
# median expression into "high" and "low" groups, and a log-rank test is
# run on overall survival. One KM plot per gene is saved.
#
# Inputs : data_out/rds/analysis_data.rds
# Outputs: data_out/figures/01_km_<gene>.png  (one per gene)
#          data_out/tables/km_logrank_results.csv
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
rds_dir <- file.path(root, "data_out", "rds")
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")

analysis <- readRDS(file.path(rds_dir, "analysis_data.rds"))

# Candidate genes = columns beyond the clinical block
clinical_cols <- c("patient_id", "os_time", "os_event",
                   "age_at_index", "gender", "stage_simple", "vital_status")
gene_cols <- setdiff(colnames(analysis), clinical_cols)
message("Testing ", length(gene_cols), " genes: ",
        paste(gene_cols, collapse = ", "))

# --- Per-gene KM + log-rank ---------------------------------------------
km_results <- data.frame()

for (g in gene_cols) {
  d <- analysis[, c("os_time", "os_event", g)]
  colnames(d)[3] <- "expr"
  med <- median(d$expr, na.rm = TRUE)
  d$group <- factor(ifelse(d$expr >= med, "High", "Low"),
                    levels = c("Low", "High"))

  fit <- survfit(Surv(os_time, os_event) ~ group, data = d)
  lr  <- survdiff(Surv(os_time, os_event) ~ group, data = d)
  # chisq-based p-value
  p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  # HR for High vs Low via univariate Cox (matching this dichotomy)
  cox <- coxph(Surv(os_time, os_event) ~ group, data = d)
  sm  <- summary(cox)
  hr  <- sm$coefficients[1, "exp(coef)"]
  hr_low  <- sm$conf.int[1, "lower .95"]
  hr_high <- sm$conf.int[1, "upper .95"]

  km_results <- rbind(km_results, data.frame(
    gene     = g,
    median_vst = round(med, 3),
    p_logrank  = p_lr,
    HR_high_vs_low = hr,
    HR_lower_95    = hr_low,
    HR_upper_95    = hr_high,
    n_high   = sum(d$group == "High"),
    n_low    = sum(d$group == "Low"),
    events_high = sum(d$os_event[d$group == "High"]),
    events_low  = sum(d$os_event[d$group == "Low"])
  ))

  # KM plot
  plot_title <- sprintf("%s expression and overall survival (TCGA-KIRC)", g)
  plot_sub   <- sprintf("Median split  |  log-rank p = %.2e  |  HR (High vs Low) = %.2f [%.2f-%.2f]",
                        p_lr, hr, hr_low, hr_high)
  km_plot <- ggsurvplot(
    fit,
    data       = d,
    pval       = TRUE,
    risk.table = TRUE,
    conf.int   = TRUE,
    palette    = c(Low = "#2980b9", High = "#c0392b"),
    legend.title = g,
    legend.labs  = c("Low", "High"),
    xlab       = "Time (days)",
    ylab       = "Overall survival",
    title      = plot_title,
    subtitle   = plot_sub,
    ggtheme    = theme_minimal(base_size = 11)
  )

  png(file.path(figs, paste0("01_km_", g, ".png")),
      width = 2000, height = 2000, res = 220)
  print(km_plot)
  dev.off()
}

# Adjust log-rank p-values across the panel (BH)
km_results$p_logrank_BH <- p.adjust(km_results$p_logrank, method = "BH")
km_results <- km_results[order(km_results$p_logrank), ]

write.csv(km_results, file.path(tables, "km_logrank_results.csv"),
          row.names = FALSE)
print(km_results)

message("KM analysis done. ", sum(km_results$p_logrank < 0.05),
        " genes with log-rank p<0.05.")
