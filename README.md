# Analysis scripts (work in progress)

Reproduces the PAH lung transcriptomic meta-analysis. Run from the project root; paths are relative. Seed 1234.

| Script | What it does |
|--------|--------------|
| `00_config.R` | Paths, seed, colour palette (sourced by the others) |
| `01a_build_phenotypes.py` | Case/control labels from GEO metadata |
| `01b_preprocess_affymetrix.R` | RMA normalisation (Affymetrix) |
| `01c_preprocess_rnaseq.R` | RNA-seq count / FPKM preprocessing |
| `01d_preprocess_agilent_illumina.R` | Agilent + Illumina normalisation |
| `01e_qc_pca.R` | QC PCA |
| `02_differential_expression.R` | Per-study differential expression |
| `03a_meta_analysis.R` | Random-effects meta-analysis + rank aggregation |
| `03b_sensitivity_analyses.R` | Leave-one-out and cohort-collapse robustness |
| `04_functional_enrichment.R` | GO / KEGG / Reactome enrichment |
| `05_singlecell_localisation.R` | Single-cell cell-of-origin |
| `06_deconvolution_ssgsea.R` | Immune/stromal deconvolution |
| `07_classifier_validation.R` | Diagnostic panel + external validation |
| `08b_figures.R`, `08c_singlecell_figures.R`, `08d_composite_figures.R` | Figures |
