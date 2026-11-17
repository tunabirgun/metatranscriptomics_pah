# PAH lung transcriptomic meta-analysis

Analysis code for a registered multi-cohort study of pulmonary arterial hypertension (PAH) lung tissue: bulk-expression preprocessing, per-study differential expression, meta-analysis, single-cell localization, enrichment and diagnostic-panel evaluation.

> **Status:** work in progress; the code may change.

[Registered protocol](https://doi.org/10.5281/zenodo.21399521)

## Data

All inputs are public in the NCBI Gene Expression Omnibus (GEO). Raw data are not included in this repository.

| Role | Accessions |
|---|---|
| Discovery | GSE113439, GSE117261, GSE53408, GSE15197, GSE48149, GSE254617, GSE272776 |
| External validation | GSE208592 |
| Single-cell reference | GSE293580 |

## Environment

`renv.lock` pins the R and Bioconductor package set. Restore it from the repository root:

```r
install.packages("renv")
renv::restore()
```

The `Dockerfile` builds the same environment. It copies the scripts into the image and expects the GEO data to be mounted at run time.

```bash
docker build -t pah-pipeline .
```

## Usage

Run every script from the project root, the directory that holds `02_acquisition/raw/` with the downloaded GEO series. `00_config.R` defines the relative paths, the output directories and the random seed. The scripts also run when kept in `submission/scripts/` under the project root. Scripts `04b` and `07b` additionally read `manuscript/tables/Table3_diagnostic_panel.csv`, which is not part of this repository.

Run the scripts one at a time in numerical order and inspect each output before continuing. Use `python` for `01a_build_phenotypes.py` and `Rscript` for the `.R` files:

```bash
python 01a_build_phenotypes.py
Rscript 02_differential_expression.R
```

| Stage | Scripts |
|---|---|
| Phenotypes, platform-specific preprocessing and QC | `01a`-`01e` |
| Per-study differential expression | `02_differential_expression.R` |
| Meta-analysis and sensitivity analyses | `03a_meta_analysis.R`, `03b_sensitivity_analyses.R` |
| Enrichment, single-cell localization and deconvolution | `04_functional_enrichment.R`, `04b_enrichment_expanded.R`, `04c_enrichment_figures.R`, `05_singlecell_localisation.R`, `06_deconvolution_ssgsea.R` |
| Classifier evaluation | `07_classifier_validation.R` |
| In silico qPCR (post hoc) | `07b_insilico_qpcr.R` |
| Figures | `08b_figures.R`, `08c_singlecell_figures.R`, `08d_composite_figures.R` |
| Reviewer-response extensions | `09a_composition_adjusted.R`, `09b_leave_one_cohort_out.R`, `09c_panel_diagnostics.R` |

## Interpretation

- `07b_insilico_qpcr.R` is a post hoc analysis outside the registered plan. It expresses the panel genes relative to the TBP and SDHA reference genes in two cohorts, GSE117261 (microarray) and GSE254617 (RNA-seq).
- The `09` scripts are reviewer-response extensions, also outside the registered plan.
- Report these analyses as such. Running the code does not establish scientific validity; exact inputs, thresholds and models are in the scripts and the registered protocol.
