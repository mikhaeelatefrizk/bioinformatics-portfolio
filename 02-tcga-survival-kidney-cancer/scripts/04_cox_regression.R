# 04_cox_regression.R ----------------------------------------------------
# Cox proportional-hazards modelling:
#  (a) Univariate Cox per gene (gene alone, expression as continuous)
#  (b) Multivariate Cox per gene (gene + age + sex + stage)
#  (c) Schoenfeld residuals to test the PH assumption
#
# Inputs : data_out/rds/analysis_data.rds
# Outputs: data_out/tables/cox_univariate.csv
#          data_out/tables/cox_multivariate.csv
#          data_out/tables/schoenfeld_tests.csv
#          data_out/figures/02_forest_univariate.png
#          data_out/figures/03_forest_multivariate.png
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

root    <- rprojroot::find_root(rprojroot::has_file("README.md"))
rds_dir <- file.path(root, "data_out", "rds")
figs    <- file.path(root, "data_out", "figures")
tables  <- file.path(root, "data_out", "tables")

analysis <- readRDS(file.path(rds_dir, "analysis_data.rds"))

clinical_cols <- c("patient_id", "os_time", "os_event",
                   "age_at_index", "gender", "stage_simple", "vital_status")
gene_cols <- setdiff(colnames(analysis), clinical_cols)

# --- 1. Univariate Cox (continuous expression) --------------------------
uni_rows <- lapply(gene_cols, function(g) {
  f <- as.formula(sprintf("Surv(os_time, os_event) ~ `%s`", g))
  fit <- coxph(f, data = analysis)
  s   <- summary(fit)
  data.frame(
    gene = g,
    HR_per_unit_vst  = s$coefficients[1, "exp(coef)"],
    HR_lower_95      = s$conf.int[1, "lower .95"],
    HR_upper_95      = s$conf.int[1, "upper .95"],
    p_wald           = s$coefficients[1, "Pr(>|z|)"],
    concordance      = s$concordance[1]
  )
})
uni <- do.call(rbind, uni_rows)
uni$p_wald_BH <- p.adjust(uni$p_wald, method = "BH")
uni <- uni[order(uni$p_wald), ]
write.csv(uni, file.path(tables, "cox_univariate.csv"), row.names = FALSE)

# --- 2. Multivariate Cox (gene + age + sex + stage) ---------------------
multi_rows <- lapply(gene_cols, function(g) {
  f <- as.formula(sprintf(
    "Surv(os_time, os_event) ~ `%s` + age_at_index + gender + stage_simple", g))
  fit <- try(coxph(f, data = analysis), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  s <- summary(fit)
  coef_tab <- s$coefficients
  ci_tab   <- s$conf.int
  # Gene row is the first coefficient
  data.frame(
    gene         = g,
    HR_gene      = coef_tab[1, "exp(coef)"],
    HR_gene_low  = ci_tab[1, "lower .95"],
    HR_gene_high = ci_tab[1, "upper .95"],
    p_gene       = coef_tab[1, "Pr(>|z|)"],
    HR_age       = coef_tab["age_at_index", "exp(coef)"],
    p_age        = coef_tab["age_at_index", "Pr(>|z|)"],
    HR_male      = coef_tab["gendermale", "exp(coef)"],
    p_male       = coef_tab["gendermale", "Pr(>|z|)"],
    n            = s$n,
    events       = s$nevent,
    concordance  = s$concordance[1]
  )
})
multi <- do.call(rbind, Filter(Negate(is.null), multi_rows))
multi$p_gene_BH <- p.adjust(multi$p_gene, method = "BH")
multi <- multi[order(multi$p_gene), ]
write.csv(multi, file.path(tables, "cox_multivariate.csv"), row.names = FALSE)

# --- 3. Schoenfeld test of PH assumption (multivariate models) ----------
sch_rows <- lapply(gene_cols, function(g) {
  f <- as.formula(sprintf(
    "Surv(os_time, os_event) ~ `%s` + age_at_index + gender + stage_simple", g))
  fit <- try(coxph(f, data = analysis), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  zp  <- cox.zph(fit)
  tbl <- as.data.frame(zp$table)
  tbl$term <- rownames(tbl)
  tbl$gene_model <- g
  tbl
})
sch <- do.call(rbind, Filter(Negate(is.null), sch_rows))
write.csv(sch, file.path(tables, "schoenfeld_tests.csv"), row.names = FALSE)

# --- 4. Forest plots -----------------------------------------------------
forest_plot <- function(d, title, out_png) {
  d$gene <- factor(d$gene, levels = rev(d$gene[order(d$HR_gene)]))
  p <- ggplot(d, aes(x = HR_gene, y = gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = HR_gene_low, xmax = HR_gene_high),
                   height = 0.2) +
    geom_point(aes(colour = p_gene < 0.05), size = 3) +
    scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey40"),
                        name = "p < 0.05") +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL,
         title = title,
         subtitle = "Multivariate Cox adjusting for age, sex, stage") +
    theme_minimal(base_size = 11)
  ggsave(out_png, p, width = 8, height = 5, dpi = 150)
}

# Univariate forest plot (renaming for reuse)
uni_fp <- data.frame(
  gene = uni$gene,
  HR_gene      = uni$HR_per_unit_vst,
  HR_gene_low  = uni$HR_lower_95,
  HR_gene_high = uni$HR_upper_95,
  p_gene       = uni$p_wald
)
forest_plot(uni_fp, "Univariate Cox: HR per unit VST expression",
            file.path(figs, "02_forest_univariate.png"))
forest_plot(multi, "Multivariate Cox: HR per unit VST expression",
            file.path(figs, "03_forest_multivariate.png"))

# --- Session info -------------------------------------------------------
sink(file.path(root, "data_out", "session_info.txt"))
sessionInfo()
sink()

message("Cox regression done. ",
        sum(multi$p_gene < 0.05, na.rm = TRUE),
        " genes significant in multivariate model (p<0.05).")
