# Step 4 (§8.5 primary): single-cell cell-of-origin localisation of the consensus signature.
# Reference rebuilt from GSE293580 raw 10X (D5). Cluster -> canonical lung-lineage annotation.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(Seurat); library(data.table); library(ggplot2) })
set.seed(1234)


# --- extract raw 10X ---
tar <- file.path(PATHS$raw,"scrna/GSE293580/GSE293580_RAW.tar")
tmp <- file.path(PATHS$sc,"tenx"); dir.create(tmp, showWarnings=FALSE, recursive=TRUE)
untar(tar, exdir=tmp)
mtx <- list.files(tmp, pattern="matrix\\.mtx\\.gz$", full.names=TRUE)
samp <- sub("_matrix\\.mtx\\.gz$","",basename(mtx))

# condition map from series matrix titles ("BA-044, Failed Donor, scRNA-seq")
sm <- readLines(gzfile(file.path(PATHS$raw,"scrna/GSE293580/GSE293580_series_matrix.txt.gz")))
tl <- grep("^!Sample_title", sm, value=TRUE)
titles <- strsplit(sub("^!Sample_title\\t","",tl), "\t")[[1]]; titles <- gsub('"','',titles)
cond_map <- sapply(titles, function(x){ p<-strsplit(x,",")[[1]]; setNames(trimws(p[2]), trimws(p[1])) })
id_of <- function(s) sub(".*_(([A-Z]{2})-[0-9]+)$","\\1", s)

objs <- list()
for (i in seq_along(mtx)) {
  d <- dirname(mtx[i]); pre <- samp[i]
  m <- ReadMtx(mtx=mtx[i],
               cells=file.path(d, paste0(pre,"_barcodes.tsv.gz")),
               features=file.path(d, paste0(pre,"_features.tsv.gz")), feature.column=2)
  o <- CreateSeuratObject(m, project=pre, min.cells=3, min.features=200)
  sid <- id_of(pre); o$sample <- sid
  o$condition <- ifelse(sid %in% names(cond_map), cond_map[sid], NA)
  objs[[pre]] <- o
}
sc <- merge(objs[[1]], objs[-1], add.cell.ids=samp)
sc <- JoinLayers(sc)
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern="^MT-")
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
message("cells after QC: ", ncol(sc), " | samples: ", length(unique(sc$sample)))

sc <- NormalizeData(sc, verbose=FALSE)
sc <- FindVariableFeatures(sc, nfeatures=2000, verbose=FALSE)
sc <- ScaleData(sc, verbose=FALSE)
sc <- RunPCA(sc, npcs=30, verbose=FALSE)
sc <- RunUMAP(sc, dims=1:30, verbose=FALSE)
sc <- FindNeighbors(sc, dims=1:30, verbose=FALSE)
sc <- FindClusters(sc, resolution=0.5, verbose=FALSE)

# canonical lung-lineage markers
markers <- list(
  Endothelial=c("PECAM1","CLDN5","VWF","CDH5"),
  SMC_Pericyte=c("ACTA2","MYH11","TAGLN","PDGFRB","NOTCH3"),
  Fibroblast=c("COL1A1","COL1A2","PDGFRA","DCN","LUM"),
  Epithelial=c("EPCAM","SFTPC","SFTPB","AGER","SCGB1A1"),
  Myeloid=c("CD68","LYZ","MARCO","ITGAM","FCGR3A"),
  T_NK=c("CD3E","CD3D","NKG7","CD8A","IL7R"),
  B_Plasma=c("MS4A1","CD79A","MZB1"),
  Mast=c("TPSAB1","CPA3","MS4A2") )
for (ct in names(markers)) {
  g <- intersect(markers[[ct]], rownames(sc))
  sc <- AddModuleScore(sc, features=list(g), name=paste0("sc_",ct), seed=1234)
}
scorecols <- paste0("sc_",names(markers),"1")
cl_scores <- sapply(scorecols, function(sccol) tapply(sc[[sccol]][,1], sc$seurat_clusters, mean))
cl_assign <- names(markers)[apply(cl_scores, 1, which.max)]
names(cl_assign) <- rownames(cl_scores)
sc$celltype <- factor(unname(cl_assign[as.character(sc$seurat_clusters)]))
message("cell-type counts:"); print(table(sc$celltype))

# --- localise consensus signature ---
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))
genes <- intersect(cons$gene, rownames(sc))
avg <- AverageExpression(sc, features=genes, group.by="celltype", assays="RNA")$RNA
# fraction of cells (within cell type) expressing each gene
pct <- sapply(levels(sc$celltype), function(ct){
  cells <- colnames(sc)[sc$celltype==ct]
  Matrix::rowMeans(GetAssayData(sc, layer="counts")[genes, cells, drop=FALSE] > 0) })
origin <- data.table(gene=genes,
  cell_of_origin=colnames(avg)[apply(avg,1,which.max)],
  max_avg_expr=round(apply(avg,1,max),3),
  pct_in_origin=round(sapply(seq_along(genes), function(i) pct[i, which.max(avg[i,])]),3))
origin <- merge(origin, cons[,.(gene,pooled_log2FC,direction)], by="gene")
setorder(origin, cell_of_origin, -max_avg_expr)
fwrite(origin, file.path(PATHS$sc,"signature_cell_of_origin.tsv"), sep="\t")
fwrite(as.data.table(avg, keep.rownames="gene"), file.path(PATHS$sc,"signature_avgexpr_by_celltype.tsv"), sep="\t")
message("\ncell-of-origin distribution:"); print(table(origin$cell_of_origin))

# --- figures ---
p1 <- DimPlot(sc, group.by="celltype", label=TRUE, repel=TRUE) +
      ggtitle(NULL) + theme_minimal(base_size=10)
ggsave(file.path(PATHS$fig,"sc_umap_celltype.png"), p1, width=8, height=6, dpi=300)
ggsave(file.path(PATHS$fig,"sc_umap_celltype.svg"), p1, width=8, height=6)
p2 <- DotPlot(sc, features=genes, group.by="celltype") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7)) +
      labs(x=NULL,y=NULL)
ggsave(file.path(PATHS$fig,"sc_dotplot_signature.png"), p2, width=12, height=5, dpi=300)
ggsave(file.path(PATHS$fig,"sc_dotplot_signature.svg"), p2, width=12, height=5)

saveRDS(sc, file.path(PATHS$sc,"GSE293580_reference.rds"))

message("DONE single-cell")
