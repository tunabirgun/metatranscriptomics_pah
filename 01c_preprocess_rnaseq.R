# Step 1b: RNA-seq preprocessing. Produces per-cohort raw-count (or log2-FPKM)
# matrices harmonised to HGNC symbol, with verified case/control labels.
# GSE254617 (counts), GSE272776 (FPKM only -> deviation D1), GSE208592 (validation counts).
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(data.table); library(org.Hs.eg.db); library(AnnotationDbi) })

ens2sym <- function(ens) {
  ens <- sub("\\..*$", "", ens)
  m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens,
        column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  unname(m)
}
collapse_sum <- function(counts, sym) {           # counts: matrix, sym: gene symbols
  ok <- !is.na(sym) & sym != ""
  counts <- counts[ok, , drop = FALSE]; sym <- sym[ok]
  agg <- rowsum(counts, group = sym)              # sum counts per symbol
  agg
}
collapse_maxmean <- function(mat, sym) {          # for FPKM: keep highest-mean row per symbol
  ok <- !is.na(sym) & sym != ""
  mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
  ord <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym); out <- mat[keep, , drop = FALSE]; rownames(out) <- sym[keep]; out
}
title2gsm <- function(gse) {                      # series-matrix title -> GSM
  d <- readRDS_sm(gse); setNames(d$gsm, d$title)
}
readRDS_sm <- function(gse) {
  # minimal series-matrix title/gsm/group via the pheno file + raw titles
  ph <- read.delim(file.path(PATHS$pheno, paste0(gse,"_pheno.tsv")), colClasses="character")
  ph
}

## ---- GSE254617 (counts; columns PHBI_xxx_Sxx) ----
message("== GSE254617 ==")
dt <- fread(file.path(PATHS$raw,"rnaseq/GSE254617/GSE254617_txiCounts.csv.gz"))
genes <- dt[[1]]; mat <- as.matrix(dt[,-1]); rownames(mat) <- genes
mat <- round(mat)                                  # tximport counts -> integer for DESeq2
storage.mode(mat) <- "integer"
ph <- read.delim(file.path(PATHS$pheno,"GSE254617_pheno.tsv"), colClasses="character")
ph <- ph[ph$group %in% c("PAH","control"), ]
colkey <- sub("_S[0-9]+$", "", colnames(mat))      # PHBI_004_S50 -> PHBI_004
keep <- colkey %in% ph$title
mat <- mat[, keep, drop=FALSE]; colkey <- colkey[keep]
# sum sequencing-replicate columns per subject (10 subjects run twice) -> one column/subject
mat <- t(rowsum(t(mat), group = colkey))
storage.mode(mat) <- "integer"
ph <- ph[match(colnames(mat), ph$title), ]
mat <- collapse_sum(mat, rownames(mat))            # symbols already; collapse dups
saveRDS(list(expr=mat, pheno=ph, platform="rnaseq_hiseq4000", scale="counts"),
        file=file.path(PATHS$proc,"GSE254617.rds"))
message(sprintf("  genes=%d samples=%d (PAH=%d control=%d)", nrow(mat), ncol(mat),
                sum(ph$group=="PAH"), sum(ph$group=="control")))


## ---- GSE272776 (FPKM only; NC=control, LUN=PAH) — DEVIATION D1 ----
message("== GSE272776 (FPKM) ==")
dt <- fread(file.path(PATHS$raw,"rnaseq/GSE272776/GSE272776_gene_fpkm.txt.gz"))
samp <- grep("^(NC|LUN)_", names(dt), value=TRUE)
mat <- as.matrix(dt[, ..samp]); rownames(mat) <- dt$gene_name  # gene_name column provided
grp <- ifelse(grepl("^NC_", samp), "control", "PAH")
logfpkm <- log2(mat + 1)
logfpkm <- collapse_maxmean(logfpkm, rownames(logfpkm))
ph <- data.frame(gsm=samp, title=samp, group=grp, subtype=grp, reason="", stringsAsFactors=FALSE)
saveRDS(list(expr=logfpkm, pheno=ph, platform="rnaseq_novaseq_FPKM", scale="log2fpkm"),
        file=file.path(PATHS$proc,"GSE272776.rds"))
message(sprintf("  genes=%d samples=%d (PAH=%d control=%d)", nrow(logfpkm), ncol(logfpkm),
                sum(grp=="PAH"), sum(grp=="control")))
for (mk in c("NPPA","NPPB","COL1A1")) if (mk %in% rownames(logfpkm)) {
  d <- mean(logfpkm[mk, grp=="PAH"]) - mean(logfpkm[mk, grp=="control"])
  message(sprintf("  marker %-6s PAH-control log2 = %+.2f", mk, d)) }


## ---- GSE208592 (VALIDATION counts; Control_Lung_/PAH_Lung_) — processed in isolation ----
message("== GSE208592 (VALIDATION) ==")
dt <- fread(file.path(PATHS$raw,"rnaseq/GSE208592/GSE208592_Counts_Lung.csv.gz"))
genes <- dt[[1]]; mat <- as.matrix(dt[,-1]); rownames(mat) <- genes
mat <- round(mat); storage.mode(mat) <- "integer"
grp <- ifelse(grepl("^Control", colnames(mat)), "control", "PAH")
sym <- ens2sym(rownames(mat))
mat <- collapse_sum(mat, sym)
ph <- data.frame(gsm=colnames(mat), title=colnames(mat), group=grp, subtype=grp,
                 reason="", stringsAsFactors=FALSE)
saveRDS(list(expr=mat, pheno=ph, platform="rnaseq_novaseq", scale="counts", validation=TRUE),
        file=file.path(PATHS$proc,"GSE208592.rds"))
message(sprintf("  genes=%d samples=%d (PAH=%d control=%d)", nrow(mat), ncol(mat),
                sum(grp=="PAH"), sum(grp=="control")))

message("DONE rnaseq")
