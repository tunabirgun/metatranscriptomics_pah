# Reviewer-response diagnostics: RF nested-CV, panel robustness to KRT4/ITLN2,
# single-cohort artifact check for LILRA2/GIMAP6, and the two-method intersection.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(glmnet); library(randomForest); library(pROC); library(edgeR); library(data.table) })
set.seed(1234)
dir.create(file.path(ROOT,"results/tables"), showWarnings=FALSE, recursive=TRUE)

disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))$gene
to_log <- function(o){ if (o$scale=="counts"){ d<-DGEList(o$expr); d<-calcNormFactors(d); cpm(d,log=TRUE,prior.count=1) } else o$expr }
zrows <- function(m){ (m - rowMeans(m)) / (apply(m,1,sd)+1e-9) }

val <- readRDS(file.path(PATHS$proc,"GSE208592.rds"))
pool <- Reduce(intersect, c(lapply(disc,function(g) rownames(readRDS(file.path(PATHS$proc,paste0(g,".rds")))$expr)), list(rownames(val$expr))))
feat_all <- intersect(cons, pool)

build_XY <- function(feat){
  Xl<-list(); yl<-c()
  for (g in disc){ o<-readRDS(file.path(PATHS$proc,paste0(g,".rds"))); m<-zrows(to_log(o)[feat,,drop=FALSE]); Xl[[g]]<-t(m); yl<-c(yl, ifelse(o$pheno$group=="PAH",1,0)) }
  Xtr<-do.call(rbind,Xl); Xte<-t(zrows(to_log(val)[feat,,drop=FALSE])); yte<-ifelse(val$pheno$group=="PAH",1,0)
  list(Xtr=Xtr,ytr=yl,Xte=Xte,yte=yte)
}
nestedCV <- function(Xtr,ytr,K=10,method="lasso"){
  set.seed(1234); folds<-sample(rep(1:K,length.out=nrow(Xtr))); oof<-rep(NA,nrow(Xtr))
  for(k in 1:K){ tr<-folds!=k; te<-folds==k
    if(method=="lasso"){ cv<-cv.glmnet(Xtr[tr,],ytr[tr],family="binomial",alpha=1,nfolds=10); oof[te]<-predict(cv,Xtr[te,],s="lambda.1se",type="response")[,1] }
    else { rf<-randomForest(Xtr[tr,],factor(ytr[tr]),ntree=1000); oof[te]<-predict(rf,Xtr[te,],type="prob")[,"1"] } }
  as.numeric(pROC::auc(ytr,oof,quiet=TRUE))
}
finalLASSO_val <- function(D){
  set.seed(1234); cvf<-cv.glmnet(D$Xtr,D$ytr,family="binomial",alpha=1,nfolds=10)
  co<-as.matrix(coef(cvf,s="lambda.1se")); panel<-setdiff(rownames(co)[co[,1]!=0],"(Intercept)")
  p<-predict(cvf,D$Xte,s="lambda.1se",type="response")[,1]; r<-pROC::roc(D$yte,p,quiet=TRUE)
  list(panel=panel, auc=as.numeric(r$auc), ci=as.numeric(pROC::ci.auc(r)))
}

out <- list()
# (E) RF nested-CV on full 28 features (for Table 3b symmetry)
D <- build_XY(feat_all)
out$rf_nestedCV_auc <- round(nestedCV(D$Xtr,D$ytr,method="rf"),3)
out$lasso_nestedCV_auc <- round(nestedCV(D$Xtr,D$ytr,method="lasso"),3)

# (3) panel robustness: drop KRT4, drop ITLN2, drop both
for (drop in list(c("KRT4"), c("ITLN2"), c("KRT4","ITLN2"))){
  f <- setdiff(feat_all, drop); r <- finalLASSO_val(build_XY(f))
  out[[paste0("drop_",paste(drop,collapse="_"))]] <- list(val_auc=round(r$auc,3),
      ci=sprintf("%.2f-%.2f",r$ci[1],r$ci[3]), panel_size=length(r$panel), krt4_in_panel="KRT4" %in% r$panel)
}

# (artifact) LILRA2 / GIMAP6: per-cohort log2FC + single-cohort dominance
M <- readRDS(file.path(PATHS$meta,"meta_objects.rds"))
res <- fread(file.path(PATHS$meta,"meta_rma_all_genes.tsv"))
art <- lapply(c("LILRA2","GIMAP6","SOSTDC1"), function(g){
  fc<-M$FC[g,]; se<-M$SE[g,]
  loo_z <- sapply(seq_along(fc), function(i){ w<-1/se[-i]^2; b<-sum(fc[-i]*w)/sum(w); b/sqrt(1/sum(w)) })
  list(gene=g, per_cohort_log2FC=round(fc,2), min_loo_absZ=round(min(abs(loo_z)),1),
       full_z=round(res[gene==g,z],1), I2=round(res[gene==g,I2],1), all_same_sign=all(sign(fc)==sign(fc[1])))
})
out$artifact_check <- art

# (intersection) metafor consensus  vs  RRA FDR<0.05
rra_up<-fread(file.path(PATHS$meta,"rra_up.tsv")); rra_dn<-fread(file.path(PATHS$meta,"rra_dn.tsv"))
rra_sig<-union(rra_up[FDR<0.05,Name], rra_dn[FDR<0.05,Name])
inter <- intersect(cons, rra_sig)
out$two_method_intersection <- list(n=length(inter), genes=sort(inter))

saveRDS(out, file.path(ROOT,"results/tables/review_diagnostics.rds"))
str(out, max.level=2)
