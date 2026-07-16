# Step 1a: RMA normalisation of Affymetrix Gene 1.0 ST cohorts (GPL6244).
# CEL -> oligo::rma -> transcript-cluster -> gene symbol (max-mean collapse).
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(oligo); library(hugene10sttranscriptcluster.db); library(AnnotationDbi) })

affy_sets <- c("GSE113439","GSE117261","GSE53408")

collapse_by_symbol <- function(mat) {
  # keep, per gene symbol, the probe row with the highest mean expression
  sym <- rownames(mat)
  ord <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym)
  out <- mat[keep, , drop = FALSE]; rownames(out) <- sym[keep]
  out
}

for (gse in affy_sets) {
  
  message("== ", gse, " ==")
  tar <- list.files(file.path(PATHS$raw, "microarray", gse), pattern = "_RAW\\.tar$", full.names = TRUE)
  tmp <- file.path(PATHS$proc, "celtmp", gse); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  untar(tar, exdir = tmp)
  celgz <- list.files(tmp, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
  for (f in celgz) R.utils::gunzip(f, overwrite = TRUE, remove = TRUE)
  cels <- list.files(tmp, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)

  raw <- read.celfiles(cels, verbose = FALSE)
  eset <- rma(raw)                       # gene-level (transcript cluster) RMA, log2
  ex <- exprs(eset)
  colnames(ex) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(ex)))  # CEL -> GSM

  # annotate transcript clusters -> SYMBOL
  ids <- rownames(ex)
  map <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = ids,
                               columns = "SYMBOL", keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID), ]
  ex <- ex[map$PROBEID, , drop = FALSE]; rownames(ex) <- map$SYMBOL
  ex <- collapse_by_symbol(ex)

  # subset to eligible samples
  ph <- read_pheno(gse)
  common <- intersect(colnames(ex), ph$gsm)
  ex <- ex[, common, drop = FALSE]; ph <- ph[match(common, ph$gsm), ]
  stopifnot(all(colnames(ex) == ph$gsm))

  saveRDS(list(expr = ex, pheno = ph, platform = "affy_hugene10st", scale = "log2_rma"),
          file = file.path(PATHS$proc, paste0(gse, ".rds")))
  message(sprintf("  genes=%d  samples=%d (PAH=%d control=%d)",
                  nrow(ex), ncol(ex), sum(ph$group=="PAH"), sum(ph$group=="control")))
  # sanity: NPPB/NPPA and BMPR2 direction (known PAH-associated)
  for (mk in c("NPPA","NPPB","BMPR2")) if (mk %in% rownames(ex)) {
    d <- mean(ex[mk, ph$group=="PAH"]) - mean(ex[mk, ph$group=="control"])
    message(sprintf("  %-6s PAH-control log2 diff = %+.2f", mk, d))
  }
  unlink(tmp, recursive = TRUE)

}
message("DONE affy RMA")
