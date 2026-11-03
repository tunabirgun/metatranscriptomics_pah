# PAH lung transcriptomic meta-analysis

Work-in-progress analysis code for a registered multi-cohort study of pulmonary arterial hypertension lung tissue. The pipeline covers bulk-expression preprocessing, study-level differential expression, meta-analysis, single-cell localization, enrichment and diagnostic-panel evaluation.

[Registered protocol](https://doi.org/10.5281/zenodo.21399521)

## Environment and inputs

Restore the package environment from the repository root:

```r
install.packages("renv")
renv::restore()
```

`renv.lock` records the package set. `00_config.R` defines paths and the seed. Data directories must be supplied beside the scripts; the legacy `submission/scripts/` layout is also detected. Raw data are not included.

## Data

All input datasets are public in the NCBI Gene Expression Omnibus (GEO). Discovery: GSE113439, GSE117261, GSE53408, GSE15197, GSE48149, GSE254617, GSE272776. External validation: GSE208592. Single-cell reference: GSE293580. No raw data is stored in this repository.

## Pipeline

Run scripts separately in numerical order and inspect outputs between stages. Use `python` for `.py` files and `Rscript` for `.R` files.

| Stage | Scripts |
|---|---|
| Phenotypes, platform-specific preprocessing and QC | `01a`–`01e` |
| Per-study differential expression | `02_differential_expression.R` |
| Meta-analysis and sensitivity | `03a_meta_analysis.R`, `03b_sensitivity_analyses.R` |
| Enrichment, single-cell localization and deconvolution | `04`–`06`, including `04b` and `04c` |
| Classifier evaluation | `07_classifier_validation.R` |
| Figures | `08b`–`08d`, with enrichment outputs from `04c` |
| Reviewer-response analyses | `09a`–`09c` |

Example invocation from the configured project root:

```bash
Rscript 02_differential_expression.R
```

**Protocol distinction:** `07b_insilico_qpcr.R` is a post-hoc analysis, not part of the registered plan. It evaluates the panel relative to TBP/SDHA reference genes in selected cohorts. The `09` analyses are reviewer-response extensions. Preserve these distinctions when reporting results; implementation alone does not demonstrate scientific validation.

See each script and the registered protocol for exact inputs, thresholds, models and execution scope.
