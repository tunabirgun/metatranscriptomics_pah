# Shared configuration. Run scripts from the PROJECT ROOT directory
# (e.g., Rscript submission/scripts/01b_preprocess_affymetrix.R). All paths are relative.
# ROOT auto-detects whether you are at the project root or inside submission/scripts/.
ROOT <- if (dir.exists("02_acquisition")) "." else if (dir.exists("../../02_acquisition")) "../.." else "."
.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
suppressWarnings(suppressMessages({ set.seed(1234) }))
PATHS <- list(
  raw   = file.path(ROOT, "02_acquisition/raw"),
  pheno = file.path(ROOT, "03_preprocessing/pheno"),
  proc  = file.path(ROOT, "03_preprocessing/processed"),
  de    = file.path(ROOT, "04_diff_expression"),
  meta  = file.path(ROOT, "05_meta_analysis"),
  func  = file.path(ROOT, "06_functional"),
  sc    = file.path(ROOT, "07_singlecell"),
  decon = file.path(ROOT, "08_deconvolution"),
  clf   = file.path(ROOT, "09_classifier"),
  fig   = file.path(ROOT, "results/figures"),
  tab   = file.path(ROOT, "results/tables")
)
# color-blind-safe palette (Okabe-Ito); reused across all figures
PAL <- c(PAH="#D55E00", control="#0072B2", up="#D55E00", down="#0072B2",
         ns="#999999", accent="#009E73", warn="#E69F00", purple="#CC79A7")
read_pheno <- function(gse) {
  p <- read.delim(file.path(PATHS$pheno, paste0(gse, "_pheno.tsv")),
                  stringsAsFactors=FALSE, colClasses="character")
  p[p$group %in% c("PAH","control"), ]
}
