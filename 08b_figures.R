# Regenerate all figures in Arial with italic gene symbols.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(matrixStats); library(edgeR) })
FIG <- PATHS$fig
save3 <- function(p,name,w,h){ ggsave(file.path(FIG,paste0(name,".png")),p,width=w,height=h,dpi=300)
                               ggsave(file.path(FIG,paste0(name,".svg")),p,width=w,height=h) }

M   <- readRDS(file.path(PATHS$meta,"meta_objects.rds"))
res <- fread(file.path(PATHS$meta,"meta_rma_all_genes.tsv"))
cons<- fread(file.path(PATHS$meta,"consensus_signature.tsv"))
core<- readLines(file.path(PATHS$meta,"robust_core_genes.txt"))

## Fig 2 volcano (italic gene labels)
res[, cat:="ns"]; res[gene %in% cons$gene & pooled_log2FC>0,cat:="up"]; res[gene %in% cons$gene & pooled_log2FC<0,cat:="down"]
lab <- res[gene %in% c(core,"VIPR1","CX3CR1","COL3A1")]
p2 <- ggplot(res, aes(pooled_log2FC,-log10(FDR))) +
  geom_point(aes(color=cat), size=ifelse(res$cat=="ns",0.5,1.6), alpha=ifelse(res$cat=="ns",0.3,0.9)) +
  scale_color_manual(values=c(up=PAL[["up"]],down=PAL[["down"]],ns=PAL[["ns"]]),guide="none") +
  geom_vline(xintercept=c(-0.585,0.585),linetype=2,color="grey55") +
  geom_hline(yintercept=-log10(0.05),linetype=2,color="grey55") +
  geom_text_repel(data=lab, aes(label=gene), fontface="italic", size=3.4,
                  max.overlaps=30, box.padding=0.4, seg.color="grey60") +
  labs(x=expression(Random-effects~pooled~log[2]~fold~change~(PAH~vs~control)), y=expression(-log[10]~meta-FDR)) +
  theme(panel.grid.minor=element_blank())
save3(p2,"meta_volcano",8,6)

## Fig 3 cohort consistency heatmap as ggplot (italic gene rows)
FC <- M$FC[cons$gene,]; ord <- order(cons$direction,-abs(cons$pooled_log2FC))
dt <- as.data.table(as.table(FC[ord,])); setnames(dt,c("Gene","Cohort","log2FC"))
dt[, Gene:=factor(Gene,levels=rev(cons$gene[ord]))]
dt[, core:=ifelse(Gene %in% core,"•","")]
lim <- quantile(abs(dt$log2FC),0.98,na.rm=TRUE)
p3 <- ggplot(dt, aes(Cohort,Gene,fill=log2FC)) + geom_tile(color="grey92") +
  scale_fill_gradient2(low=PAL[["down"]],mid="white",high=PAL[["up"]],limits=c(-lim,lim),oob=scales::squish,name=expression(log[2]~FC)) +
  labs(x=NULL,y=NULL) +
  theme(axis.text.y=element_text(face="italic",size=8),
        axis.text.x=element_text(angle=45,hjust=1,size=9),
        panel.grid=element_blank())
save3(p3,"consensus_cohort_heatmap",7,8)

## Fig 4 forest (italic facet strips)
disc <- colnames(M$FC)
rows <- rbindlist(lapply(core,function(g){ fc<-M$FC[g,];se<-M$SE[g,]
  s<-data.table(gene=g,study=disc,est=fc,lo=fc-1.96*se,hi=fc+1.96*se,type="study")
  pr<-res[gene==g]; pooled<-data.table(gene=g,study="Pooled (RE)",est=pr$pooled_log2FC,
    lo=pr$pooled_log2FC-1.96*pr$pooled_SE,hi=pr$pooled_log2FC+1.96*pr$pooled_SE,type="pooled")
  rbind(s,pooled)}))
rows[,study:=factor(study,levels=rev(c(disc,"Pooled (RE)")))]; rows[,gene:=factor(gene,levels=core)]
p4 <- ggplot(rows, aes(est,study,color=type)) + geom_vline(xintercept=0,linetype=2,color="grey60") +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=0.25,linewidth=0.4) +
  geom_point(aes(size=type,shape=type)) + facet_wrap(~gene,scales="free_x",ncol=3) +
  scale_color_manual(values=c(study=PAL[["control"]],pooled=PAL[["PAH"]]),guide="none") +
  scale_shape_manual(values=c(study=16,pooled=18),guide="none") +
  scale_size_manual(values=c(study=1.8,pooled=3.2),guide="none") +
  labs(x=expression(log[2]~fold~change~(PAH~vs~control)),y=NULL) +
  theme(strip.text=element_text(face="italic"),panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7))
save3(p4,"Figure4_forest_robustcore",11,7)

## Fig 5 functional ORA (Times, GO terms not italic)
ego <- fread(file.path(PATHS$func,"ORA_GO_BP.tsv"))[order(p.adjust)][1:15]
ego[, Description:=factor(Description,levels=rev(Description))]
p5 <- ggplot(ego, aes(-log10(p.adjust),Description,size=Count,color=-log10(p.adjust))) + geom_point() +
  scale_color_viridis_c() + labs(x=expression(-log[10]~adjusted~p),y=NULL,size="Genes",color=NULL) +
  theme(panel.grid.minor=element_blank())
save3(p5,"functional_ORA_GO",9,5.5)

## Fig 7 deconvolution as ggplot tile (Times, lineages not genes)
eff <- fread(file.path(PATHS$decon,"lineage_effect_cohensd.tsv"))
dm <- melt(eff, id.vars="lineage", variable.name="Cohort", value.name="d")
dm[, lineage:=factor(lineage,levels=eff$lineage[order(rowMeans(as.matrix(eff[,-1])))])]
lim7 <- max(abs(dm$d),na.rm=TRUE)
p7 <- ggplot(dm, aes(Cohort,lineage,fill=d)) + geom_tile(color="grey92") +
  geom_text(aes(label=sprintf("%.2f",d)),size=2.6) +
  scale_fill_gradient2(low=PAL[["control"]],mid="white",high=PAL[["PAH"]],limits=c(-lim7,lim7),name="Cohen's d") +
  labs(x=NULL,y=NULL) + theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
save3(p7,"deconvolution_heatmap",8,6)

## Fig 8 validation ROC + calibration (Times)
C <- readRDS(file.path(PATHS$clf,"classifier_objects.rds"))
library(pROC); roc_l <- pROC::roc(C$yte,C$p_lasso,quiet=TRUE)
rocdf <- data.table(fpr=1-roc_l$specificities,tpr=roc_l$sensitivities)[order(fpr)]
p8a <- ggplot(rocdf,aes(fpr,tpr)) + geom_line(color=PAL[["PAH"]],linewidth=1) +
  geom_abline(linetype=2,color="grey60") +
  annotate("text",x=.6,y=.15,label=sprintf("AUC = %.3f (95%% CI %.3f-%.3f)",C$auc_l,C$ci[1],C$ci[3]),size=3.5) +
  labs(x="False positive rate",y="True positive rate") + theme(panel.grid.minor=element_blank())
save3(p8a,"validation_ROC",5.5,5)
cal <- data.table(p=C$p_lasso,y=C$yte)[, bin:=cut(p,seq(0,1,0.2),include.lowest=TRUE)][,.(mp=mean(p),obs=mean(y),n=.N),by=bin][order(mp)]
p8b <- ggplot(cal,aes(mp,obs)) + geom_abline(linetype=2,color="grey60") +
  geom_point(aes(size=n),color=PAL[["control"]]) + geom_line(color=PAL[["control"]]) +
  xlim(0,1)+ylim(0,1)+labs(x="Mean predicted probability",y="Observed PAH fraction",size="n") + theme(panel.grid.minor=element_blank())
save3(p8b,"validation_calibration",5.5,5)

## Fig S1 QC PCA (Times)
sets <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776","GSE208592")
pca <- rbindlist(lapply(sets,function(g){ o<-readRDS(file.path(PATHS$proc,paste0(g,".rds")));ex<-o$expr
  if(o$scale=="counts"){d<-edgeR::DGEList(ex);d<-edgeR::calcNormFactors(d);ex<-edgeR::cpm(d,log=TRUE,prior.count=1)}
  v<-matrixStats::rowVars(ex);top<-head(order(v,decreasing=TRUE),2000);pc<-prcomp(t(ex[top,]),scale.=TRUE)
  ve<-round(100*pc$sdev^2/sum(pc$sdev^2),1)
  data.table(PC1=pc$x[,1],PC2=pc$x[,2],group=o$pheno$group,lab=sprintf("%s (%.0f%%, %.0f%%)",g,ve[1],ve[2]))}))
pS1 <- ggplot(pca,aes(PC1,PC2,color=group)) + geom_point(size=1.6,alpha=0.85) +
  facet_wrap(~lab,scales="free",ncol=4) + scale_color_manual(values=PAL) +
  labs(color=NULL) + theme(strip.text=element_text(face="bold",size=8),panel.grid.minor=element_blank(),legend.position="top")
save3(pS1,"QC_pca_overview",12,6.5)
message("DONE regen ggplot figures")
