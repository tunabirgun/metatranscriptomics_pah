# Step 2 (§8.1): per-study differential expression. Each cohort analysed independently;
# emits a standardised table (gene, log2FC, SE, stat, pvalue, FDR). +log2FC = up in PAH.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(limma); library(edgeR); library(DESeq2); library(data.table) })

micro <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149")
covars <- list(GSE117261=c("sex"), GSE254617=c("sex","batch"))  # available covariates (§8.1)

get_pheno_cov <- function(gse, ids, idcol="gsm") {
  ph <- read.delim(file.path(PATHS$pheno, paste0(gse,"_pheno.tsv")), colClasses="character")
  ph[match(ids, ph[[idcol]]), , drop=FALSE]
}
build_design <- function(group, cov=NULL) {
  group <- factor(group, levels=c("control","PAH"))
  df <- data.frame(group=group)
  form <- "~ group"
  if (!is.null(cov)) for (cn in names(cov)) {
    v <- cov[[cn]]
    if (length(unique(v[!is.na(v) & v!=""])) > 1) { df[[cn]] <- factor(v); form <- paste(form,"+",cn) }
  }
  list(mm=model.matrix(as.formula(form), data=df), df=df, form=form)
}
save_de <- function(gse, tab) {
  fwrite(tab, file.path(PATHS$de, paste0(gse,"_de.tsv")), sep="\t")
  sig <- sum(tab$FDR < 0.05, na.rm=TRUE)
  up <- sum(tab$FDR<0.05 & tab$log2FC>0, na.rm=TRUE); dn <- sum(tab$FDR<0.05 & tab$log2FC<0, na.rm=TRUE)
  message(sprintf("  %s: genes=%d  FDR<0.05=%d (up=%d dn=%d)", gse, nrow(tab), sig, up, dn))
}

## ---- microarray: limma ----
for (gse in micro) {
  o <- readRDS(file.path(PATHS$proc, paste0(gse,".rds")))
  ex <- o$expr; ph <- o$pheno
  cov <- NULL
  if (!is.null(covars[[gse]])) {
    pc <- get_pheno_cov(gse, ph$gsm); cov <- setNames(lapply(covars[[gse]], function(k) pc[[k]]), covars[[gse]])
  }
  des <- build_design(ph$group, cov)
  fit <- eBayes(lmFit(ex, des$mm))
  cf <- "groupPAH"
  tt <- topTable(fit, coef=cf, number=Inf, sort.by="none")
  se <- tt$logFC / tt$t
  tab <- data.table(gene=rownames(ex), log2FC=tt$logFC, SE=abs(se),
                    stat=tt$t, pvalue=tt$P.Value, FDR=tt$adj.P.Val)
  attr(tab,"model") <- des$form
  save_de(gse, tab); message("     model: ", des$form)

}

## ---- GSE254617: DESeq2 (counts) + sex + batch ----
{ gse <- "GSE254617"; o <- readRDS(file.path(PATHS$proc, paste0(gse,".rds")))
  ex <- o$expr; ph <- o$pheno
  pc <- get_pheno_cov(gse, ph$title, "title")
  cd <- data.frame(group=factor(ph$group, levels=c("control","PAH")),
                   sex=factor(pc$sex), batch=factor(make.names(pc$batch)))
  keep <- edgeR::filterByExpr(ex, group=cd$group)             # low-count filter (§7)
  dds <- DESeqDataSetFromMatrix(ex[keep,], colData=cd, design=~ sex + batch + group)
  dds <- DESeq(dds, quiet=TRUE)
  res <- as.data.frame(results(dds, contrast=c("group","PAH","control")))
  tab <- data.table(gene=rownames(res), log2FC=res$log2FoldChange, SE=res$lfcSE,
                    stat=res$stat, pvalue=res$pvalue, FDR=res$padj)
  tab <- tab[!is.na(log2FC) & !is.na(SE)]
  attr(tab,"model") <- "~ sex + batch + group"
  save_de(gse, tab); message("     model: ~ sex + batch + group (filtered ", sum(keep)," genes)")
 }

## ---- GSE272776: limma-trend on log2-FPKM (DEVIATION D1) ----
{ gse <- "GSE272776"; o <- readRDS(file.path(PATHS$proc, paste0(gse,".rds")))
  ex <- o$expr; ph <- o$pheno
  keep <- rowMeans(ex) > 1                                    # drop near-zero FPKM
  des <- build_design(ph$group)
  fit <- eBayes(lmFit(ex[keep,], des$mm), trend=TRUE)
  tt <- topTable(fit, coef="groupPAH", number=Inf, sort.by="none")
  se <- tt$logFC/tt$t
  tab <- data.table(gene=rownames(ex[keep,]), log2FC=tt$logFC, SE=abs(se),
                    stat=tt$t, pvalue=tt$P.Value, FDR=tt$adj.P.Val)
  attr(tab,"model") <- "~ group (limma-trend, log2FPKM)"
  save_de(gse, tab); message("     model: limma-trend on log2FPKM [D1]")
 }

message("DONE per-study DE")
