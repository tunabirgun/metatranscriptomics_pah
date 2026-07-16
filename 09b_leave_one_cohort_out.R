# Leave-one-cohort-out validation: train LASSO on 6 discovery cohorts, test on the 7th.
# Extra generalisation evidence beyond the single external cohort.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(glmnet); library(pROC); library(edgeR); library(data.table) })
set.seed(1234)
disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))$gene
to_log <- function(o){ if (o$scale=="counts"){ d<-DGEList(o$expr); d<-calcNormFactors(d); cpm(d,log=TRUE,prior.count=1) } else o$expr }
zrows <- function(m){ (m - rowMeans(m)) / (apply(m,1,sd)+1e-9) }

O <- lapply(disc, function(g) readRDS(file.path(PATHS$proc, paste0(g,".rds")))); names(O)<-disc
pool <- Reduce(intersect, lapply(O, function(o) rownames(o$expr)))
feat <- intersect(cons, pool)
Xc <- lapply(disc, function(g){ t(zrows(to_log(O[[g]])[feat,,drop=FALSE])) }); names(Xc)<-disc
yc <- lapply(disc, function(g) ifelse(O[[g]]$pheno$group=="PAH",1,0)); names(yc)<-disc

res <- data.table(held_out=character(), n=integer(), n_PAH=integer(), AUC=numeric(), CI=character())
for (h in disc){
  tr <- setdiff(disc, h)
  Xtr <- do.call(rbind, Xc[tr]); ytr <- unlist(yc[tr])
  set.seed(1234); cv <- cv.glmnet(Xtr, ytr, family="binomial", alpha=1, nfolds=10)
  p <- as.numeric(predict(cv, Xc[[h]], s="lambda.1se", type="response"))
  r <- pROC::roc(yc[[h]], p, quiet=TRUE); ci <- as.numeric(pROC::ci.auc(r))
  res <- rbind(res, data.table(held_out=h, n=length(yc[[h]]), n_PAH=sum(yc[[h]]),
                               AUC=round(as.numeric(r$auc),3), CI=sprintf("%.2f-%.2f",ci[1],ci[3])))
}
fwrite(res, file.path(ROOT,"results/tables/review_loco_validation.tsv"), sep="\t")
print(res)
cat(sprintf("\nMean leave-one-cohort-out AUC: %.3f (range %.3f-%.3f); external GSE208592 AUC 0.90\n",
            mean(res$AUC), min(res$AUC), max(res$AUC)))
