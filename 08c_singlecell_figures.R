source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(Seurat); library(data.table); library(ggplot2) })
FIG <- PATHS$fig
sc <- readRDS(file.path(PATHS$sc,"GSE293580_reference.rds"))
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))
genes <- intersect(cons$gene, rownames(sc))

# UMAP: cell types (Times, not italic)
p1 <- DimPlot(sc, group.by="celltype", label=TRUE, repel=TRUE, label.size=3) +
  ggtitle(NULL) + theme_minimal(base_size=10) +
  theme(text=element_text(family="Arial"), legend.text=element_text(family="Arial"))
ggsave(file.path(FIG,"sc_umap_celltype.png"), p1, width=8, height=6, dpi=300)
ggsave(file.path(FIG,"sc_umap_celltype.svg"), p1, width=8, height=6)

# Dotplot: genes on x-axis italic Times
p2 <- DotPlot(sc, features=genes, group.by="celltype") +
  theme_minimal(base_size=10) +
  theme(text=element_text(family="Arial"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7, face="italic"),
        axis.text.y=element_text(family="Arial")) + labs(x=NULL,y=NULL)
ggsave(file.path(FIG,"sc_dotplot_signature.png"), p2, width=12, height=5, dpi=300)
ggsave(file.path(FIG,"sc_dotplot_signature.svg"), p2, width=12, height=5)
message("DONE sc figures")
