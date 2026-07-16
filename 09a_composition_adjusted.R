# Composition-adjusted meta-analysis: add ssGSEA lineage scores as covariates in each
# per-study model, re-run DE, re-meta, and report which consensus genes survive.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(limma); library(edgeR); library(DESeq2); library(GSVA); library(metafor); library(data.table) })
set.seed(1234)

sets <- list(
  Endothelial=c("PECAM1","CDH5","VWF","CLDN5"),
  Fibroblast =c("COL1A1","PDGFRA","DCN","LUM"),
  SMC        =c("ACTA2","MYH11","TAGLN","CNN1"),
  Myeloid    =c("CD68","CD163","LYZ","MARCO","CSF1R") )   # main shifting lineages
disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
to_log <- function(o){ if (o$scale=="counts"){ d<-DGEList(o$expr); d<-calcNormFactors(d); cpm(d,log=TRUE,prior.count=1) } else o$expr }
lin_scores <- function(logexpr){
  par <- ssgseaParam(exprData=as.matrix(logexpr), geneSets=lapply(sets,function(s) intersect(s,rownames(logexpr))))
  t(gsva(par, verbose=FALSE))   # samples x lineages
}

FC <- list(); SE <- list()
for (g in disc){
  o <- readRDS(file.path(PATHS$proc, paste0(g,".rds"))); ph <- o$pheno
  le <- to_log(o); sc <- lin_scores(le)
  cov <- as.data.frame(scale(sc)); grp <- factor(ph$group, levels=c("control","PAH"))
  if (o$scale=="counts"){
    keep <- edgeR::filterByExpr(o$expr, group=grp)
    cd <- cbind(data.frame(group=grp), cov)
    dds <- DESeqDataSetFromMatrix(o$expr[keep,], colData=cd,
             design=as.formula(paste("~", paste(c(colnames(cov),"group"),collapse="+"))))
    dds <- DESeq(dds, quiet=TRUE); r <- as.data.frame(results(dds, name="group_PAH_vs_control"))
    FC[[g]] <- setNames(r$log2FoldChange, rownames(r)); SE[[g]] <- setNames(r$lfcSE, rownames(r))
  } else {
    trend <- (o$scale=="log2fpkm")
    des <- model.matrix(as.formula(paste("~ group +", paste(colnames(cov),collapse="+"))), data=cbind(cov, group=grp))
    ex <- if (trend) le[rowMeans(le)>1,,drop=FALSE] else le
    fit <- eBayes(lmFit(ex, des), trend=trend)
    tt <- topTable(fit, coef="groupPAH", number=Inf, sort.by="none")
    FC[[g]] <- setNames(tt$logFC, rownames(ex)); SE[[g]] <- setNames(abs(tt$logFC/tt$t), rownames(ex))
  }
  message("adjusted DE done: ", g)
}

genes <- Reduce(intersect, lapply(FC, names))
FCm <- sapply(disc, function(g) FC[[g]][genes]); SEm <- sapply(disc, function(g) SE[[g]][genes])
rownames(FCm)<-genes; rownames(SEm)<-genes; SEm[SEm<=0|!is.finite(SEm)]<-NA

b<-p<-nc<-rep(NA_real_,length(genes))
for (i in seq_along(genes)){ yi<-FCm[i,]; si<-SEm[i,]; ok<-is.finite(yi)&is.finite(si); if(sum(ok)<3) next
  m<-tryCatch(rma(yi[ok],sei=si[ok],method="REML",control=list(maxiter=200)),error=function(e) tryCatch(rma(yi[ok],sei=si[ok],method="DL"),error=function(e2) NULL))
  if(is.null(m)) next; b[i]<-as.numeric(m$b); p[i]<-m$pval; nc[i]<-sum(sign(yi[ok])==sign(b[i])) }
adj <- data.table(gene=genes, adj_log2FC=round(b,3), adj_FDR=p.adjust(p,"BH"), n_concord=nc)

cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))
out <- merge(cons[,.(gene, orig_log2FC=round(pooled_log2FC,3), orig_FDR=FDR, direction)], adj, by="gene", all.x=TRUE)
out[, survives_adjustment := !is.na(adj_FDR) & adj_FDR<0.05 & sign(adj_log2FC)==sign(orig_log2FC)]
setorder(out, -survives_adjustment, adj_FDR)
fwrite(out, file.path(ROOT,"results/tables/review_composition_adjusted_consensus.tsv"), sep="\t")
cat(sprintf("\nConsensus genes surviving composition adjustment: %d / %d\n", sum(out$survives_adjustment), nrow(out)))
cat("Surviving:", paste(out[survives_adjustment==TRUE, gene], collapse=", "), "\n")
cat("Attenuated:", paste(out[survives_adjustment==FALSE, gene], collapse=", "), "\n")
