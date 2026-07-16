# Step 3 (§8.2): meta-analytic synthesis over the 7 discovery cohorts.
# (i) metafor random-effects REML on per-gene log2FC + SE (pooled est, I2, tau2)
# (ii) RobustRankAggreg on signed-significance ranks (up / down)
# Consensus (locked): meta-FDR<0.05 AND |pooled log2FC|>=0.585 AND >=6/7 sign-concordant (D3).
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(metafor); library(RobustRankAggreg); library(data.table) })

disc <- c("GSE113439","GSE117261","GSE53408","GSE15197","GSE48149","GSE254617","GSE272776")
de <- lapply(disc, function(g) fread(file.path(PATHS$de, paste0(g,"_de.tsv"))))
names(de) <- disc

# strict intersection of gene symbols across all 7 (§7 / D4)
genes <- Reduce(intersect, lapply(de, function(d) d$gene))
genes <- sort(unique(genes))
message("intersection genes across 7 cohorts: ", length(genes))

# assemble per-gene matrices of log2FC and SE
FC <- sapply(disc, function(g){ d<-de[[g]]; d[match(genes,d$gene), log2FC] })
SE <- sapply(disc, function(g){ d<-de[[g]]; d[match(genes,d$gene), SE] })
rownames(FC) <- genes; rownames(SE) <- genes
SE[SE<=0 | !is.finite(SE)] <- NA

# ---- (i) random-effects REML per gene ----

res <- data.table(gene=genes, k=NA_integer_, pooled_log2FC=NA_real_, pooled_SE=NA_real_,
                  z=NA_real_, pval=NA_real_, I2=NA_real_, tau2=NA_real_, n_concord=NA_integer_)
for (i in seq_along(genes)) {
  yi <- FC[i,]; sei <- SE[i,]; ok <- is.finite(yi) & is.finite(sei)
  if (sum(ok) < 3) next
  m <- tryCatch(rma(yi=yi[ok], sei=sei[ok], method="REML", control=list(maxiter=200)),
                error=function(e) tryCatch(rma(yi=yi[ok], sei=sei[ok], method="DL"), error=function(e2) NULL))
  if (is.null(m)) next
  b <- as.numeric(m$b)
  res[i, `:=`(k=sum(ok), pooled_log2FC=b, pooled_SE=as.numeric(m$se), z=as.numeric(m$zval),
              pval=as.numeric(m$pval), I2=m$I2, tau2=m$tau2,
              n_concord=sum(sign(yi[ok])==sign(b)))]
}
res[, FDR := p.adjust(pval, "BH")]


# consensus flags
res[, consistent6 := n_concord >= 6]
res[, consistent7 := n_concord == k]
res[, consensus := FDR < 0.05 & abs(pooled_log2FC) >= 0.585 & consistent6 & k == 7]
res[, direction := fifelse(pooled_log2FC > 0, "up", "down")]
res <- res[order(FDR, -abs(pooled_log2FC))]
fwrite(res, file.path(PATHS$meta, "meta_rma_all_genes.tsv"), sep="\t")
cons <- res[consensus == TRUE]
fwrite(cons, file.path(PATHS$meta, "consensus_signature.tsv"), sep="\t")

message(sprintf("meta genes tested: %d", sum(!is.na(res$pval))))
message(sprintf("CONSENSUS (FDR<0.05 & |log2FC|>=0.585 & >=6/7 concordant & k=7): %d  (up=%d down=%d)",
                nrow(cons), sum(cons$direction=="up"), sum(cons$direction=="down")))
message(sprintf("  strict all-7-concordant subset: %d", nrow(res[consensus==TRUE & consistent7])))
message(sprintf("  median I2 among consensus: %.1f%%", median(cons$I2, na.rm=TRUE)))

# ---- (ii) RobustRankAggreg (signed) ----

mk_list <- function(sign_dir) {
  lapply(disc, function(g){
    d <- de[[g]]; d <- d[gene %in% genes]
    d <- if (sign_dir>0) d[log2FC>0] else d[log2FC<0]
    setorder(d, pvalue); d$gene
  })
}
rra_up <- aggregateRanks(glist=mk_list(+1), N=length(genes))
rra_dn <- aggregateRanks(glist=mk_list(-1), N=length(genes))
setDT(rra_up); setDT(rra_dn)
rra_up[, FDR := p.adjust(Score,"BH")]; rra_dn[, FDR := p.adjust(Score,"BH")]
fwrite(rra_up, file.path(PATHS$meta,"rra_up.tsv"), sep="\t")
fwrite(rra_dn, file.path(PATHS$meta,"rra_dn.tsv"), sep="\t")


# agreement metafor vs RRA (top consensus genes appear in RRA significant lists)
rra_sig <- union(rra_up[FDR<0.05, Name], rra_dn[FDR<0.05, Name])
ov <- intersect(cons$gene, rra_sig)
message(sprintf("RRA FDR<0.05 genes: up=%d down=%d; overlap with metafor consensus: %d/%d (%.0f%%)",
                sum(rra_up$FDR<0.05), sum(rra_dn$FDR<0.05), length(ov), nrow(cons),
                100*length(ov)/max(1,nrow(cons))))
saveRDS(list(FC=FC, SE=SE, res=res, cons=cons, rra_up=rra_up, rra_dn=rra_dn),
        file=file.path(PATHS$meta,"meta_objects.rds"))
message("DONE meta-analysis")
