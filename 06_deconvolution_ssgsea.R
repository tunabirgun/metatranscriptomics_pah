# Step 7 (§8.5 additional): immune/stromal deconvolution via ssGSEA (GSVA), per discovery cohort.
# Curated canonical lineage marker sets (D7). Compares PAH vs control; effect per cohort.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(GSVA); library(edgeR); library(data.table); library(ggplot2); library(pheatmap) })


sets <- list(
  T_cell      = c("CD3D","CD3E","CD2","CD8A","IL7R"),
  B_cell      = c("CD19","MS4A1","CD79A","CD79B"),
  NK_cell     = c("NKG7","KLRD1","GNLY","NCAM1"),
  Myeloid     = c("CD68","CD163","LYZ","MARCO","CSF1R"),
  Neutrophil  = c("FCGR3B","CSF3R","S100A8","S100A9"),
  Mast        = c("TPSAB1","CPA3","MS4A2"),
  Endothelial = c("PECAM1","CDH5","VWF","CLDN5"),
  Fibroblast  = c("COL1A1","PDGFRA","DCN","LUM"),
  SMC         = c("ACTA2","MYH11","TAGLN","CNN1"),
  Cytotoxic   = c("GZMB","PRF1","GZMK","IFNG") )

disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
to_log <- function(o){ if (o$scale=="counts"){ d<-DGEList(o$expr); d<-calcNormFactors(d); cpm(d,log=TRUE,prior.count=1) } else o$expr }
cohend <- function(x,y){ (mean(x)-mean(y))/sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2)) }

eff <- matrix(NA, length(sets), length(disc), dimnames=list(names(sets), disc))
pval <- eff
for (g in disc) {
  o <- readRDS(file.path(PATHS$proc, paste0(g,".rds"))); m <- to_log(o)
  par <- ssgseaParam(exprData=as.matrix(m), geneSets=lapply(sets, function(s) intersect(s, rownames(m))))
  sc <- gsva(par, verbose=FALSE)                         # lineages x samples
  grp <- o$pheno$group
  for (ln in rownames(sc)) {
    x <- sc[ln, grp=="PAH"]; y <- sc[ln, grp=="control"]
    eff[ln,g] <- cohend(x,y); pval[ln,g] <- wilcox.test(x,y)$p.value
  }
}
out <- as.data.table(eff, keep.rownames="lineage")
fwrite(out, file.path(PATHS$decon,"lineage_effect_cohensd.tsv"), sep="\t")
fwrite(as.data.table(pval, keep.rownames="lineage"), file.path(PATHS$decon,"lineage_wilcox_p.tsv"), sep="\t")

# summary: mean effect across cohorts + n cohorts significant & concordant
summ <- data.table(lineage=names(sets),
  mean_d = round(rowMeans(eff),3),
  n_up_PAH = rowSums(eff>0),
  n_sig = rowSums(pval<0.05, na.rm=TRUE),
  n_sig_concordant = sapply(seq_len(nrow(eff)), function(i) sum(pval[i,]<0.05 & sign(eff[i,])==sign(rowMeans(eff)[i]), na.rm=TRUE)))
setorder(summ, -mean_d)
fwrite(summ, file.path(PATHS$decon,"deconv_summary.tsv"), sep="\t")
print(summ)

# heatmap lineage x cohort (Cohen's d, PAH vs control)
png(file.path(PATHS$fig,"deconvolution_heatmap.png"), width=8, height=6, units="in", res=300)
pheatmap(eff, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.2f",
         color=colorRampPalette(c(PAL[["control"]],"white",PAL[["PAH"]]))(51),
         breaks=seq(-max(abs(eff),na.rm=TRUE),max(abs(eff),na.rm=TRUE),length.out=52),
         main="", fontsize=9)
dev.off()
svg(file.path(PATHS$fig,"deconvolution_heatmap.svg"), width=8, height=6)
pheatmap(eff, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.2f",
         color=colorRampPalette(c(PAL[["control"]],"white",PAL[["PAH"]]))(51),
         breaks=seq(-max(abs(eff),na.rm=TRUE),max(abs(eff),na.rm=TRUE),length.out=52), fontsize=9)
dev.off()

message("DONE deconvolution")
