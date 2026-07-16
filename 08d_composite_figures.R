# Composite multi-panel figures with bold A/B panel labels (Arial). Figure 6 (single-cell)
# and Figure 8 (validation). Individual panels are archived under figures/raw_figures/.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(Seurat); library(data.table); library(ggplot2); library(patchwork); library(pROC) })
FIG <- PATHS$fig
tagtheme <- theme(plot.tag = element_text(face = "bold", size = 16))

## ---- Figure 6: single-cell (A UMAP, B dotplot) ----
sc   <- readRDS(file.path(PATHS$sc, "GSE293580_reference.rds"))
cons <- fread(file.path(PATHS$meta, "consensus_signature.tsv"))
genes <- intersect(cons$gene, rownames(sc))
p_umap <- DimPlot(sc, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle(NULL) + theme_minimal(base_size = 10) +
  theme(text = element_text(family = "Arial"), legend.text = element_text(family = "Arial"))
p_dot <- DotPlot(sc, features = genes, group.by = "celltype") +
  theme_minimal(base_size = 10) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "italic"),
        axis.text.y = element_text(family = "Arial")) + labs(x = NULL, y = NULL)
fig6 <- (p_umap / p_dot) + plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") & tagtheme
ggsave(file.path(FIG, "Figure6_singlecell.png"), fig6, width = 10, height = 11, dpi = 300)
ggsave(file.path(FIG, "Figure6_singlecell.svg"), fig6, width = 10, height = 11)

## ---- Figure 8: validation (A ROC, B calibration) ----
C <- readRDS(file.path(PATHS$clf, "classifier_objects.rds"))
roc_l <- pROC::roc(C$yte, C$p_lasso, quiet = TRUE)
rocdf <- data.table(fpr = 1 - roc_l$specificities, tpr = roc_l$sensitivities)[order(fpr)]
p_roc <- ggplot(rocdf, aes(fpr, tpr)) + geom_line(color = PAL[["PAH"]], linewidth = 1) +
  geom_abline(linetype = 2, color = "grey60") +
  annotate("text", x = .6, y = .15, size = 3.5,
           label = sprintf("AUC = %.3f (95%% CI %.3f-%.3f)", C$auc_l, C$ci[1], C$ci[3])) +
  labs(x = "False positive rate", y = "True positive rate") + theme(panel.grid.minor = element_blank())
cal <- data.table(p = C$p_lasso, y = C$yte)[, bin := cut(p, seq(0,1,0.2), include.lowest = TRUE)][
  , .(mp = mean(p), obs = mean(y), n = .N), by = bin][order(mp)]
p_cal <- ggplot(cal, aes(mp, obs)) + geom_abline(linetype = 2, color = "grey60") +
  geom_point(aes(size = n), color = PAL[["control"]]) + geom_line(color = PAL[["control"]]) +
  xlim(0,1) + ylim(0,1) + labs(x = "Mean predicted probability", y = "Observed PAH fraction", size = "n") +
  theme(panel.grid.minor = element_blank())
fig8 <- (p_roc + p_cal) + plot_annotation(tag_levels = "A") & tagtheme
ggsave(file.path(FIG, "Figure8_validation.png"), fig8, width = 11, height = 5.2, dpi = 300)
ggsave(file.path(FIG, "Figure8_validation.svg"), fig8, width = 11, height = 5.2)
message("DONE composites")
