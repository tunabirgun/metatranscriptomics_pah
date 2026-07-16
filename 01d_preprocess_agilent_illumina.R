# Step 1c: Agilent one-colour (GSE15197, raw->normalize) and Illumina (GSE48149, series matrix).
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(limma); library(GEOquery); library(data.table) })

collapse_maxmean <- function(mat, sym) {
  ok <- !is.na(sym) & sym != "" & !grepl("^GE_|CONTROL|^ERCC|^\\(\\+\\)|^\\(-\\)", sym, ignore.case=TRUE)
  mat <- mat[ok,,drop=FALSE]; sym <- sym[ok]
  ord <- order(rowMeans(mat), decreasing=TRUE); mat <- mat[ord,,drop=FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym); out <- mat[keep,,drop=FALSE]; rownames(out) <- sym[keep]; out
}

## ---- GSE15197: Agilent 4x44K one-colour, raw -> normalize (limma) ----
message("== GSE15197 (Agilent raw) ==")
tmp <- file.path(PATHS$proc,"agtmp"); dir.create(tmp, showWarnings=FALSE, recursive=TRUE)
untar(file.path(PATHS$raw,"microarray/GSE15197/GSE15197_RAW.tar"), exdir=tmp)
gz <- list.files(tmp, pattern="\\.txt\\.gz$", full.names=TRUE)
for (f in gz) R.utils::gunzip(f, overwrite=TRUE, remove=TRUE)
files <- list.files(tmp, pattern="\\.txt$", full.names=TRUE)
gsm <- sub("^(GSM[0-9]+).*","\\1", basename(files))
x <- read.maimages(files, source="agilent.median", green.only=TRUE,
                   other.columns="gIsWellAboveBG", verbose=FALSE)
colnames(x) <- gsm
y <- backgroundCorrect(x, method="normexp", offset=16)      # standard Agilent 1-colour
y <- normalizeBetweenArrays(y, method="quantile")           # -> log2 scale
ex <- y$E; rownames(ex) <- y$genes$ProbeName
sym <- y$genes$GeneName
ctl <- y$genes$ControlType != 0
ex <- ex[!ctl,,drop=FALSE]; sym <- sym[!ctl]
ex <- collapse_maxmean(ex, sym)
ph <- read_pheno("GSE15197"); common <- intersect(colnames(ex), ph$gsm)
ex <- ex[,common,drop=FALSE]; ph <- ph[match(common, ph$gsm),]
saveRDS(list(expr=ex, pheno=ph, platform="agilent_4x44k", scale="log2_quantile"),
        file=file.path(PATHS$proc,"GSE15197.rds"))
message(sprintf("  genes=%d samples=%d (PAH=%d control=%d)", nrow(ex), ncol(ex),
                sum(ph$group=="PAH"), sum(ph$group=="control")))
unlink(tmp, recursive=TRUE)


## ---- GSE48149: Illumina custom, processed series matrix ----
message("== GSE48149 (Illumina series matrix) ==")
es <- getGEO(filename=file.path(PATHS$raw,"microarray/GSE48149/GSE48149_series_matrix.txt.gz"),
             getGPL=TRUE)
ex <- Biobase::exprs(es)
fd <- Biobase::fData(es)
symcol <- grep("^(Symbol|GENE_SYMBOL|ILMN_Gene|Gene.?[Ss]ymbol)$", colnames(fd), value=TRUE)[1]
message("  symbol column: ", symcol)
sym <- fd[[symcol]]
ex <- log2(pmax(ex, 1))                                    # linear intensities -> log2
ex <- normalizeBetweenArrays(ex, method="quantile")
ex <- collapse_maxmean(ex, sym)
ph <- read_pheno("GSE48149"); common <- intersect(colnames(ex), ph$gsm)
ex <- ex[,common,drop=FALSE]; ph <- ph[match(common, ph$gsm),]
saveRDS(list(expr=ex, pheno=ph, platform="illumina_gpl16221", scale="log2_quantile"),
        file=file.path(PATHS$proc,"GSE48149.rds"))
message(sprintf("  genes=%d samples=%d (PAH=%d control=%d)", nrow(ex), ncol(ex),
                sum(ph$group=="PAH"), sum(ph$group=="control")))

message("DONE microarray other")
