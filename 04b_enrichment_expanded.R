# Expanded ORA: GO (BP, CC, MF) + KEGG for three gene sets against the 8,386-gene
# tested universe. Sets: 13-gene LASSO panel, 28-gene consensus, 613 meta-significant
# (FDR<0.05; also split into up/down). BH adjustment; GO tables simplified for redundancy.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db)
                   library(data.table) })
t0 <- Sys.time()
OUT <- file.path(PATHS$func, "expanded"); dir.create(OUT, showWarnings=FALSE)

# ---- inputs ----
res  <- fread(file.path(PATHS$meta, "meta_rma_all_genes.tsv"))
cons <- fread(file.path(PATHS$meta, "consensus_signature.tsv"))$gene
lasso <- fread(file.path(ROOT, "manuscript/tables/Table3_diagnostic_panel.csv"))$Gene
sig  <- res[FDR < 0.05]
universe <- res[!is.na(pval), gene]

sym2eg <- function(s) {
  m <- bitr(unique(s), "SYMBOL", "ENTREZID", OrgDb=org.Hs.eg.db)
  message(sprintf("  mapped %d/%d symbols", nrow(m), length(unique(s))))
  m$ENTREZID
}
uni_eg <- sym2eg(universe)

sets <- list(
  LASSO_panel   = lasso,
  Consensus     = cons,
  AllSig_pooled = sig$gene,
  AllSig_up     = sig[direction == "up",   gene],
  AllSig_down   = sig[direction == "down", gene]
)
cat("Gene-set sizes:\n"); print(sapply(sets, length))

# ---- ORA driver: GO BP/CC/MF (+ simplified) and KEGG ----
run_go <- function(eg, ont) {
  e <- tryCatch(enrichGO(eg, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                         universe=uni_eg, pAdjustMethod="BH",
                         pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE),
                error=function(err) NULL)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(list(raw=e, simp=e))
  s <- tryCatch(simplify(e, cutoff=0.7, by="p.adjust", select_fun=min),
                error=function(err) e)
  list(raw=e, simp=s)
}
run_kegg <- function(eg) tryCatch(
  enrichKEGG(eg, organism="hsa", universe=uni_eg, pAdjustMethod="BH",
             pvalueCutoff=0.05, qvalueCutoff=0.1),
  error=function(err) NULL)

objs <- list(); summ <- list()
for (nm in names(sets)) {
  eg <- sym2eg(sets[[nm]])
  message("== ", nm, " (", length(eg), " mapped genes) ==")
  go <- lapply(c(BP="BP", CC="CC", MF="MF"), function(o) run_go(eg, o))
  kg <- run_kegg(eg)
  # readable KEGG symbols
  if (!is.null(kg) && nrow(as.data.frame(kg)) > 0)
    kg <- setReadable(kg, org.Hs.eg.db, keyType="ENTREZID")
  objs[[nm]] <- list(BP=go$BP, CC=go$CC, MF=go$MF, KEGG=kg, genes_eg=eg)
  # write full result tables (simplified GO + KEGG)
  wr <- function(obj, tag) if (!is.null(obj) && nrow(as.data.frame(obj)) > 0)
    fwrite(as.data.table(as.data.frame(obj)),
           file.path(OUT, sprintf("ORA_%s_%s.tsv", nm, tag)), sep="\t")
  wr(go$BP$simp, "GO_BP"); wr(go$CC$simp, "GO_CC"); wr(go$MF$simp, "GO_MF"); wr(kg, "KEGG")
  # also keep raw (unsimplified) GO for the record
  wrraw <- function(obj, tag) if (!is.null(obj) && nrow(as.data.frame(obj)) > 0)
    fwrite(as.data.table(as.data.frame(obj)),
           file.path(OUT, sprintf("ORA_%s_%s_raw.tsv", nm, tag)), sep="\t")
  wrraw(go$BP$raw, "GO_BP"); wrraw(go$CC$raw, "GO_CC"); wrraw(go$MF$raw, "GO_MF")
  sig_n <- function(obj) if (is.null(obj)) 0L else sum(as.data.frame(obj)$p.adjust < 0.05, na.rm=TRUE)
  summ[[nm]] <- c(BP=sig_n(go$BP$simp), CC=sig_n(go$CC$simp),
                  MF=sig_n(go$MF$simp), KEGG=sig_n(kg))
}
saveRDS(objs, file.path(OUT, "enrich_objects.rds"))

# ---- provenance: KEGG release ----
kegg_rel <- tryCatch(paste(readLines("https://rest.kegg.jp/info/kegg", warn=FALSE)[1:2], collapse=" | "),
                     error=function(e) "unavailable")
writeLines(c(paste("run_date:", format(Sys.time(), "%Y-%m-%d")),
             paste("KEGG:", kegg_rel),
             paste("clusterProfiler:", as.character(packageVersion("clusterProfiler"))),
             paste("org.Hs.eg.db:", as.character(packageVersion("org.Hs.eg.db"))),
             paste("universe_genes:", length(universe), "mapped:", length(uni_eg))),
           file.path(OUT, "provenance.txt"))

cat("\n=== significant terms (p.adjust<0.05), simplified GO ===\n")
print(do.call(rbind, summ))
if (exists("log_step")) log_step("06b_enrichment_expanded", t0)
message("DONE expanded enrichment -> ", OUT)
