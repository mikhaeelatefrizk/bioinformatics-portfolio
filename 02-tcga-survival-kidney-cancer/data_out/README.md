# data_out/

Analysis outputs.

## tables/
| File | Contents |
|------|----------|
| `clinical_summary.csv` | Cohort demographics (n, age, sex, stage, events) |
| `km_logrank_results.csv` | Per-gene median-split log-rank test (p-value, HR, CI) |
| `cox_univariate.csv` | Per-gene univariate Cox PH (HR, CI, p, Wald) |
| `cox_multivariate.csv` | Multivariate Cox PH (each gene + age + sex + stage) |
| `schoenfeld_tests.csv` | Proportional-hazards assumption diagnostics |

## figures/
| File | Contents |
|------|----------|
| `01_km_<gene>.png` | Kaplan-Meier curves, high- vs low-expression groups |
| `02_forest_univariate.png` | Forest plot of univariate HRs across candidate genes |
| `03_forest_multivariate.png` | Forest plot of fully-adjusted HRs |
| `04_cohort_overview.png` | Cohort demographics bar chart |

## rds/
- `analysis_data.rds` — merged clinical + expression table ready for modelling.
