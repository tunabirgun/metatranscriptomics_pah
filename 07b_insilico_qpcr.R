# In silico qPCR for the 13-gene diagnostic panel (adapting Kayalar et al. 2024, Viruses).
# Per-sample expression relative to the geometric mean of two reference genes (TBP, SDHA),
# on the log2 scale. Microarray intensities are used on the same footing as RNA-seq counts.
# One microarray cohort (GSE117261) and one RNA-seq cohort (GSE254617). GAPDH is avoided
# because it is differentially expressed between PAH and control here; TBP and SDHA are not.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t0 <- Sys.time()

panel <- fread(file.path(ROOT, "manuscript/tables/Table3_diagnostic_panel.csv"))
panel <- panel[Gene != ""]; genes <- panel$Gene
REF <- c("TBP", "SDHA")
roles <- c(GSE117261="microarray", GSE254617="RNA-seq")

log2expr <- function(o) if (o$scale == "counts") log2(o$expr + 0.5) else as.matrix(o$expr)  # arrays already log2
map_group <- function(o) {                        # colnames map to gsm or title depending on cohort
  cn <- colnames(o$expr)
  for (k in c("title","gsm")) if (k %in% names(o$pheno) && all(cn %in% o$pheno[[k]]))
    return(o$pheno$group[match(cn, o$pheno[[k]])])
  o$pheno$group
}

long <- list(); stats <- list()
for (g in names(roles)) {
  o <- readRDS(file.path(PATHS$proc, paste0(g, ".rds")))
  L <- log2expr(o); grp <- map_group(o)
  keep <- grp %in% c("PAH","control"); L <- L[, keep, drop=FALSE]; grp <- grp[keep]
  refmean <- colMeans(L[REF, , drop=FALSE])       # log2 geometric mean of the two references
  for (gene in genes) {
    if (!gene %in% rownames(L)) next
    rel <- L[gene, ] - refmean
    long[[paste(g,gene)]] <- data.table(cohort=g, gene=gene, group=grp, rel=rel)
    p <- tryCatch(wilcox.test(rel[grp=="PAH"], rel[grp=="control"])$p.value, error=function(e) NA)
    stats[[paste(g,gene)]] <- data.table(cohort=g, gene=gene,
      delta=mean(rel[grp=="PAH"]) - mean(rel[grp=="control"]), p=p,
      n_PAH=sum(grp=="PAH"), n_ctrl=sum(grp=="control"))
  }
}
long <- rbindlist(long); st <- rbindlist(stats)
st[, padj := p.adjust(p, "BH"), by=cohort]
st[, stars := fifelse(padj<1e-4,"****", fifelse(padj<1e-3,"***",
             fifelse(padj<1e-2,"**", fifelse(padj<0.05,"*","ns"))))]
st[, direction := panel$Direction[match(gene, panel$Gene)]]
fwrite(st, file.path(PATHS$clf, "insilico_qpcr_stats.tsv"), sep="\t")
fwrite(long, file.path(PATHS$clf, "insilico_qpcr_long.tsv"), sep="\t")

# ---- figure: per-gene violin panels (article style), one block per cohort ----
glev <- panel[order(match(Direction,c("down","up")), -abs(LASSO_coef)), Gene]
long[, gene := factor(gene, levels=glev)]; st[, gene := factor(gene, levels=glev)]
long[, group := factor(group, levels=c("control","PAH"))]

grid_plot <- function(gse) {
  d <- long[cohort==gse]; s <- st[cohort==gse]
  yr <- d[, .(ymax=max(rel), ymin=min(rel)), by=gene]
  s <- merge(s, yr, by="gene"); s[, ypos := ymax + 0.16*(ymax - ymin + 1e-6)]
  ggplot(d, aes(group, rel)) +
    geom_violin(aes(fill=group), scale="width", trim=FALSE, linewidth=0.3, alpha=0.85) +
    stat_summary(fun=median, geom="crossbar", width=0.55, linewidth=0.22, colour="grey15") +
    geom_jitter(aes(shape=group), width=0.12, height=0, size=0.45, alpha=0.5, stroke=0) +
    geom_text(data=s, aes(x=1.5, y=ypos, label=stars), inherit.aes=FALSE, size=2.9) +
    facet_wrap(~gene, scales="free_y", ncol=7) +
    scale_fill_manual(values=c(control=unname(PAL["control"]), PAH=unname(PAL["PAH"])),
                      labels=c("Control","PAH"), name=NULL) +
    scale_shape_manual(values=c(control=16, PAH=17), guide="none") +
    scale_x_discrete(labels=c(control="Ctrl", PAH="PAH")) +
    labs(x=NULL, y="log2 expression relative to TBP, SDHA",
         subtitle=sprintf("%s  (%s; PAH n=%d, control n=%d)", gse, roles[[gse]], s$n_PAH[1], s$n_ctrl[1])) +
    theme_minimal(base_size=9) +
    theme(strip.text=element_text(face="italic", size=8.5),
          panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          legend.position="top", plot.subtitle=element_text(face="bold", size=9.5))
}
p <- grid_plot("GSE117261") / grid_plot("GSE254617") +
     plot_annotation(tag_levels="A") + plot_layout(guides="collect") &
     theme(legend.position="top")
ggsave(file.path(PATHS$fig, "insilico_qpcr_panel.png"), p, width=10, height=9, dpi=300)
ggsave(file.path(PATHS$fig, "insilico_qpcr_panel.svg"), p, width=10, height=9)

FIGD <- file.path(ROOT, "manuscript/figures")
for (ext in c("png","svg"))
  file.copy(file.path(PATHS$fig, paste0("insilico_qpcr_panel.",ext)),
            file.path(FIGD, paste0("Figure9_insilico_qpcr.",ext)), overwrite=TRUE)

cat("\n=== significant genes (BH-FDR<0.05) and direction concordance per cohort ===\n")
print(st[, .(sig=sum(padj<0.05), concordant=sum(sign(delta)==ifelse(direction=="up",1,-1))), by=cohort])
if (exists("log_step")) log_step("09c_insilico_qpcr", t0)
message("DONE in silico qPCR")
