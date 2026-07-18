# Figures for the expanded ORA: faceted GO(BP/CC/MF)+KEGG dotplots for the LASSO panel
# and consensus signature, and an up-vs-down comparison for the 613 meta-significant genes.
# Okabe-Ito/viridis, minimal theme, no titles; PNG (300 dpi) + SVG.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(clusterProfiler); library(data.table)
                   library(ggplot2); library(patchwork) })
t0 <- Sys.time()
objs <- readRDS(file.path(PATHS$func, "expanded", "enrich_objects.rds"))
FDR <- 0.05; TOPN <- 8
ONTS <- c(BP="GO: Biological process", CC="GO: Cellular component",
          MF="GO: Molecular function", KEGG="KEGG pathway")

getobj <- function(set, ont) {  # simplified GO; plain KEGG
  x <- objs[[set]][[ont]]
  if (ont == "KEGG") x else x$simp
}
tidy <- function(set, ont, topn=TOPN) {
  e <- getobj(set, ont)
  if (is.null(e)) return(NULL)
  d <- as.data.table(as.data.frame(e))
  if (!nrow(d)) return(NULL)
  d <- d[p.adjust < FDR]
  if (!nrow(d)) return(NULL)
  d[, ratio := sapply(strsplit(GeneRatio, "/"), function(z) as.numeric(z[1])/as.numeric(z[2]))]
  d <- d[order(p.adjust)][1:min(topn, .N)]
  d[, .(ontology = factor(ONTS[[ont]], levels=unname(ONTS)),
        Description, Count, ratio, p.adjust, neglog = -log10(p.adjust))]
}
gather <- function(set) rbindlist(lapply(names(ONTS), function(o) tidy(set, o)), fill=TRUE)

# ---- single-set faceted dotplot ----
dotfig <- function(set, file, h=8.5, w=8.2) {
  d <- gather(set)
  if (is.null(d) || !nrow(d)) { message("no significant terms for ", set); return(invisible()) }
  d[, Description := factor(Description, levels=rev(unique(Description)))]
  p <- ggplot(d, aes(ratio, Description, size=Count, colour=neglog)) +
    geom_point() +
    scale_colour_viridis_c(option="D", name=expression(-log[10]~FDR)) +
    scale_size_continuous(range=c(2,7), name="Genes") +
    facet_grid(ontology ~ ., scales="free_y", space="free_y", switch="y",
               labeller=label_wrap_gen(width=18)) +
    labs(x="Gene ratio", y=NULL) +
    theme_minimal(base_size=11) +
    theme(panel.grid.minor=element_blank(),
          strip.text.y.left=element_text(angle=0, face="bold", hjust=0),
          strip.placement="outside",
          axis.text.y=element_text(size=9),
          legend.position="right")
  ggsave(paste0(file, ".png"), p, width=w, height=h, dpi=300)
  ggsave(paste0(file, ".svg"), p, width=w, height=h)
  message("wrote ", file, " (", nrow(d), " terms)")
}

# ---- up vs down comparison for the 613 set ----
updown_fig <- function(file, h=7, w=8.5) {
  du <- gather("AllSig_up");   if (!is.null(du) && nrow(du)) du[, Direction := "Up (166)"]
  dd <- gather("AllSig_down"); if (!is.null(dd) && nrow(dd)) dd[, Direction := "Down (447)"]
  d <- rbindlist(list(du, dd), fill=TRUE)
  if (!nrow(d)) { message("no up/down terms"); return(invisible()) }
  d[, Direction := factor(Direction, levels=c("Up (166)","Down (447)"))]
  d[, Description := factor(Description, levels=rev(unique(Description[order(ontology, p.adjust)])))]
  p <- ggplot(d, aes(Direction, Description, size=Count, colour=neglog)) +
    geom_point() +
    scale_colour_viridis_c(option="D", name=expression(-log[10]~FDR)) +
    scale_size_continuous(range=c(2,7), name="Genes") +
    facet_grid(ontology ~ ., scales="free_y", space="free_y", switch="y",
               labeller=label_wrap_gen(width=18)) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=11) +
    theme(panel.grid.minor=element_blank(),
          strip.text.y.left=element_text(angle=0, face="bold", hjust=0),
          strip.placement="outside", axis.text.y=element_text(size=9))
  ggsave(paste0(file, ".png"), p, width=w, height=h, dpi=300)
  ggsave(paste0(file, ".svg"), p, width=w, height=h)
  message("wrote ", file, " (", nrow(d), " terms)")
}

FIG <- PATHS$fig
dotfig("Consensus",   file.path(FIG, "enrichment_consensus"))
dotfig("LASSO_panel", file.path(FIG, "enrichment_lasso_panel"))
updown_fig(           file.path(FIG, "enrichment_allsig_updown"))

# copy to manuscript/figures under embed names (Figure 5 upgraded; S2/S3 supplementary)
FIGD <- file.path(ROOT, "manuscript/figures")
embed_map <- list(c("enrichment_consensus",     "Figure5_functional_enrichment"),
                  c("enrichment_lasso_panel",   "FigureS2_lasso_enrichment"),
                  c("enrichment_allsig_updown", "FigureS3_allsig_updown"))
for (cc in embed_map) for (ext in c("png","svg")) {
  src <- file.path(FIG, paste0(cc[1], ".", ext))
  if (file.exists(src)) file.copy(src, file.path(FIGD, paste0(cc[2], ".", ext)), overwrite=TRUE)
}
if (exists("log_step")) log_step("06c_enrichment_figures", t0)
message("DONE figures")
