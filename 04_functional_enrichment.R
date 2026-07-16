# Step 5 (§8.4): functional enrichment. ORA (GO/KEGG/Reactome) on consensus genes;
# GSEA on the full ranked meta list. Universe = 8,386 tested genes.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db)
                   library(data.table); library(ggplot2) })

res  <- fread(file.path(PATHS$meta,"meta_rma_all_genes.tsv"))
cons <- fread(file.path(PATHS$meta,"consensus_signature.tsv"))
universe <- res[!is.na(pval), gene]
sym2eg <- function(s) bitr(s, "SYMBOL","ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
uni_eg <- sym2eg(universe)

# ---- ORA on the 28 consensus genes ----
cons_eg <- sym2eg(cons$gene)
ego <- enrichGO(cons_eg, OrgDb=org.Hs.eg.db, ont="BP", universe=uni_eg,
                pAdjustMethod="BH", qvalueCutoff=0.1, readable=TRUE)
ekegg <- tryCatch(enrichKEGG(cons_eg, universe=uni_eg, pAdjustMethod="BH"), error=function(e) NULL)
ereact <- tryCatch(enrichPathway(cons_eg, universe=uni_eg, pAdjustMethod="BH",
                                 qvalueCutoff=0.1, readable=TRUE), error=function(e) NULL)
if (!is.null(ego))    fwrite(as.data.table(ego@result),   file.path(PATHS$func,"ORA_GO_BP.tsv"), sep="\t")
if (!is.null(ekegg))  fwrite(as.data.table(ekegg@result), file.path(PATHS$func,"ORA_KEGG.tsv"), sep="\t")
if (!is.null(ereact)) fwrite(as.data.table(ereact@result),file.path(PATHS$func,"ORA_Reactome.tsv"), sep="\t")
message("ORA GO BP terms (q<0.1): ", if(!is.null(ego)) sum(ego@result$qvalue<0.1, na.rm=TRUE) else 0)

# ---- GSEA on full ranked list (signed z) ----
rk <- res[!is.na(z)]
map <- bitr(rk$gene, "SYMBOL","ENTREZID", OrgDb=org.Hs.eg.db)
rk <- merge(rk, as.data.table(map), by.x="gene", by.y="SYMBOL")
rk <- rk[order(-z)]
gl <- setNames(rk$z, rk$ENTREZID); gl <- gl[!duplicated(names(gl))]
gse_go <- tryCatch(gseGO(sort(gl,decreasing=TRUE), OrgDb=org.Hs.eg.db, ont="BP",
                         pAdjustMethod="BH", seed=TRUE, verbose=FALSE), error=function(e) NULL)
gse_re <- tryCatch(gsePathway(sort(gl,decreasing=TRUE), pAdjustMethod="BH",
                              seed=TRUE, verbose=FALSE), error=function(e) NULL)
if (!is.null(gse_go)) fwrite(as.data.table(gse_go@result), file.path(PATHS$func,"GSEA_GO_BP.tsv"), sep="\t")
if (!is.null(gse_re)) fwrite(as.data.table(gse_re@result), file.path(PATHS$func,"GSEA_Reactome.tsv"), sep="\t")
message("GSEA GO BP sig (p.adj<0.05): ", if(!is.null(gse_go)) sum(gse_go@result$p.adjust<0.05) else 0)

# ---- figure: top ORA GO terms ----
if (!is.null(ego) && nrow(ego@result) > 0) {
  d <- as.data.table(ego@result)[order(p.adjust)][1:min(15,.N)]
  d[, Description := factor(Description, levels=rev(Description))]
  p <- ggplot(d, aes(-log10(p.adjust), Description, size=Count, color=-log10(p.adjust))) +
    geom_point() + scale_color_viridis_c() +
    labs(x=expression(-log[10]~adjusted~p), y=NULL, size="Genes", color=NULL) +
    theme_minimal(base_size=10) + theme(panel.grid.minor=element_blank())
  ggsave(file.path(PATHS$fig,"functional_ORA_GO.png"), p, width=9, height=5.5, dpi=300)
  ggsave(file.path(PATHS$fig,"functional_ORA_GO.svg"), p, width=9, height=5.5)
}
message("DONE functional")
