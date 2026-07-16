# Step 1d: QC — per-cohort PCA (top-variable genes) coloured by group; outlier scan.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(ggplot2); library(matrixStats); library(edgeR) })

sets <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149",
          "GSE254617","GSE272776","GSE208592")
pca_df <- list(); qc_rows <- list()
for (gse in sets) {
  o <- readRDS(file.path(PATHS$proc, paste0(gse, ".rds")))
  ex <- o$expr; ph <- o$pheno
  if (o$scale == "counts") {                       # log-CPM for visualisation only
    d <- edgeR::DGEList(ex); d <- edgeR::calcNormFactors(d)
    ex <- edgeR::cpm(d, log=TRUE, prior.count=1)
  }
  v <- matrixStats::rowVars(ex); top <- head(order(v, decreasing=TRUE), 2000)
  pc <- prcomp(t(ex[top,,drop=FALSE]), scale.=TRUE)
  ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
  df <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], group=ph$group, gse=gse,
                   lab=sprintf("%s (PC1 %.0f%%, PC2 %.0f%%)", gse, ve[1], ve[2]))
  pca_df[[gse]] <- df
  # outlier scan: |robust z| > 5 on PC1 or PC2 within cohort
  z <- function(x){ (x-median(x))/ (mad(x)+1e-9) }
  out <- ph$gsm[abs(z(df$PC1))>5 | abs(z(df$PC2))>5]
  qc_rows[[gse]] <- data.frame(gse=gse, n=ncol(ex), genes=nrow(o$expr),
    scale=o$scale, n_outlier=length(out),
    outliers=ifelse(length(out), paste(out, collapse=";"), ""))
}
allpca <- do.call(rbind, pca_df)
qc <- do.call(rbind, qc_rows)
write.table(qc, file.path(PATHS$tab,"QC_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(qc)

p <- ggplot(allpca, aes(PC1, PC2, color=group)) +
  geom_point(size=1.8, alpha=0.85) +
  facet_wrap(~lab, scales="free", ncol=4) +
  scale_color_manual(values=PAL) +
  labs(x="PC1", y="PC2", color=NULL) +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold", size=8),
        panel.grid.minor=element_blank(), legend.position="top")
ggsave(file.path(PATHS$fig,"QC_pca_overview.png"), p, width=12, height=6.5, dpi=300)
ggsave(file.path(PATHS$fig,"QC_pca_overview.svg"), p, width=12, height=6.5)
message("QC done -> results/figures/QC_pca_overview.{png,svg}")
