# Step 6 (§9): compact diagnostic panel from the consensus signature.
# Per-cohort gene-wise z-standardisation (D6); nested-CV LASSO (+RF); ONE-SHOT validation on GSE208592.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(glmnet); library(randomForest); library(pROC); library(edgeR)
                   library(data.table); library(ggplot2) })
set.seed(1234)

disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))$gene

to_log <- function(o){ if (o$scale=="counts"){ d<-DGEList(o$expr); d<-calcNormFactors(d); cpm(d,log=TRUE,prior.count=1) } else o$expr }
zrows  <- function(m){ (m - rowMeans(m)) / (apply(m,1,sd)+1e-9) }   # z per gene across samples

# feature set = consensus genes present in ALL discovery + validation
val <- readRDS(file.path(PATHS$proc,"GSE208592.rds"))
gene_pool <- Reduce(intersect, c(lapply(disc, function(g) rownames(readRDS(file.path(PATHS$proc,paste0(g,".rds")))$expr)),
                                 list(rownames(val$expr))))
feat <- intersect(cons, gene_pool)
message("panel candidate features (consensus ∩ all cohorts): ", length(feat))

# assemble discovery matrix (samples x genes), z per cohort
Xl <- list(); yl <- c()
for (g in disc) {
  o <- readRDS(file.path(PATHS$proc, paste0(g,".rds")))
  m <- to_log(o)[feat,,drop=FALSE]; m <- zrows(m)
  Xl[[g]] <- t(m); yl <- c(yl, ifelse(o$pheno$group=="PAH",1,0))
}
Xtr <- do.call(rbind, Xl); ytr <- yl
Xte <- t(zrows(to_log(val)[feat,,drop=FALSE])); yte <- ifelse(val$pheno$group=="PAH",1,0)
message(sprintf("train: %d samples (%d PAH); test GSE208592: %d (%d PAH)",
                nrow(Xtr), sum(ytr), nrow(Xte), sum(yte)))

# ---- nested CV on discovery (outer 10-fold; inner cv.glmnet lambda) ----
set.seed(1234); K <- 10; folds <- sample(rep(1:K, length.out=nrow(Xtr)))
oof <- rep(NA_real_, nrow(Xtr))
for (k in 1:K) {
  tr <- folds!=k; te <- folds==k
  cvm <- cv.glmnet(Xtr[tr,], ytr[tr], family="binomial", alpha=1, nfolds=10)
  oof[te] <- as.numeric(predict(cvm, Xtr[te,], s="lambda.1se", type="response"))
}
auc_cv <- as.numeric(pROC::auc(response=ytr, predictor=oof, quiet=TRUE))
message(sprintf("nested-CV discovery AUC (LASSO, lambda.1se): %.3f", auc_cv))

# ---- final LASSO on all discovery -> compact panel ----
set.seed(1234); cvf <- cv.glmnet(Xtr, ytr, family="binomial", alpha=1, nfolds=10)
co <- as.matrix(coef(cvf, s="lambda.1se")); panel <- rownames(co)[co[,1]!=0]; panel <- setdiff(panel,"(Intercept)")
message(sprintf("compact panel size (lambda.1se nonzero): %d", length(panel)))
panel_dt <- data.table(gene=rownames(co), coef=co[,1])[coef!=0]
fwrite(panel_dt, file.path(PATHS$clf,"lasso_panel_coefficients.tsv"), sep="\t")

# ---- ONE-SHOT external validation ----
p_lasso <- as.numeric(predict(cvf, Xte, s="lambda.1se", type="response"))
roc_l <- pROC::roc(yte, p_lasso, quiet=TRUE); auc_l <- as.numeric(roc_l$auc)
ci_l <- as.numeric(pROC::ci.auc(roc_l))
brier <- mean((p_lasso-yte)^2)

# RF comparison
set.seed(1234); rf <- randomForest(x=Xtr, y=factor(ytr), ntree=1000, importance=TRUE)
p_rf <- predict(rf, Xte, type="prob")[,"1"]; auc_rf <- as.numeric(pROC::auc(yte, p_rf, quiet=TRUE))
imp <- data.table(gene=rownames(rf$importance), MeanDecreaseGini=rf$importance[,"MeanDecreaseGini"])[order(-MeanDecreaseGini)]
fwrite(imp, file.path(PATHS$clf,"rf_importance.tsv"), sep="\t")

val_out <- data.table(model=c("LASSO","RandomForest"),
                      discovery_nestedCV_AUC=c(round(auc_cv,3),NA),
                      validation_AUC=c(round(auc_l,3),round(auc_rf,3)),
                      AUC_CI95=c(sprintf("%.3f-%.3f",ci_l[1],ci_l[3]),NA),
                      Brier=c(round(brier,3),NA), panel_size=c(length(panel),length(feat)))
fwrite(val_out, file.path(PATHS$clf,"validation_performance.tsv"), sep="\t")
print(val_out)

# ---- figures: ROC + calibration ----
rocdf <- data.table(fpr=1-roc_l$specificities, tpr=roc_l$sensitivities)[order(fpr)]
pr <- ggplot(rocdf, aes(fpr,tpr)) + geom_line(color=PAL["PAH"], linewidth=1) +
  geom_abline(linetype=2, color="grey60") +
  annotate("text", x=.6, y=.15, label=sprintf("AUC = %.3f (95%% CI %.3f-%.3f)",auc_l,ci_l[1],ci_l[3]), size=3.5) +
  labs(x="False positive rate", y="True positive rate") + theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank())
ggsave(file.path(PATHS$fig,"validation_ROC.png"), pr, width=5.5, height=5, dpi=300)
ggsave(file.path(PATHS$fig,"validation_ROC.svg"), pr, width=5.5, height=5)

cal <- data.table(p=p_lasso, y=yte)[, bin:=cut(p, seq(0,1,0.2), include.lowest=TRUE)][
  , .(mean_pred=mean(p), obs=mean(y), n=.N), by=bin][order(mean_pred)]
pc <- ggplot(cal, aes(mean_pred, obs)) + geom_abline(linetype=2,color="grey60") +
  geom_point(aes(size=n), color=PAL["control"]) + geom_line(color=PAL["control"]) +
  xlim(0,1)+ylim(0,1) + labs(x="Mean predicted probability", y="Observed PAH fraction", size="n") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(PATHS$fig,"validation_calibration.png"), pc, width=5.5, height=5, dpi=300)
ggsave(file.path(PATHS$fig,"validation_calibration.svg"), pc, width=5.5, height=5)

saveRDS(list(cvf=cvf, rf=rf, feat=feat, panel=panel, auc_cv=auc_cv, auc_l=auc_l, auc_rf=auc_rf,
             ci=ci_l, brier=brier, p_lasso=p_lasso, yte=yte), file.path(PATHS$clf,"classifier_objects.rds"))
message("DONE classifier")
