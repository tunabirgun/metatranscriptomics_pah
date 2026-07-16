# Step 3b (§8.3): mandatory sensitivity analyses on the consensus signature.
# Leave-one-out, Toronto-collapse (113439+53408), US-pair-collapse (117261+254617), PVOD-exclusion.
source(if (file.exists("00_config.R")) "00_config.R" else "submission/scripts/00_config.R")
suppressMessages({ library(metafor); library(DESeq2); library(edgeR); library(data.table) })

M <- readRDS(file.path(PATHS$meta,"meta_objects.rds"))
FC <- M$FC; SE <- M$SE; disc <- colnames(FC)
cons0 <- M$cons$gene                                   # original 28
genes <- rownames(FC)

meta_consensus <- function(FCm, SEm) {                 # random-effects; return consensus gene set
  k <- ncol(FCm); b <- p <- nc <- rep(NA_real_, nrow(FCm))
  for (i in seq_len(nrow(FCm))) {
    yi <- FCm[i,]; sei <- SEm[i,]; ok <- is.finite(yi) & is.finite(sei)
    if (sum(ok) < 3) next
    m <- tryCatch(rma(yi=yi[ok], sei=sei[ok], method="REML", control=list(maxiter=200)),
                  error=function(e) tryCatch(rma(yi[ok], sei=sei[ok], method="DL"), error=function(e2) NULL))
    if (is.null(m)) next
    b[i] <- as.numeric(m$b); p[i] <- m$pval; nc[i] <- sum(sign(yi[ok])==sign(b[i]))
  }
  fdr <- p.adjust(p,"BH")
  present <- rowSums(is.finite(FCm) & is.finite(SEm))
  genes[which(fdr<0.05 & abs(b)>=0.585 & nc>=ceiling(0.85*k) & present==k)]
}
fe_pool <- function(cols) {                            # fixed-effect collapse of a study pair -> one column
  w <- 1/SE[,cols]^2; yc <- rowSums(FC[,cols]*w)/rowSums(w); sc <- sqrt(1/rowSums(w))
  list(fc=yc, se=sc)
}
jac <- function(a,b) length(intersect(a,b))/length(union(a,b))
report <- data.table(scenario=character(), n_consensus=integer(), retained_of_28=integer(), jaccard=numeric())
scen_sets <- list()   # collect consensus gene sets per scenario for per-gene robustness

# leave-one-out
for (d in disc) {
  cols <- setdiff(disc, d)
  cs <- meta_consensus(FC[,cols], SE[,cols])
  scen_sets[[paste0("LOO_drop_",d)]] <- cs
  report <- rbind(report, data.table(scenario=paste0("LOO_drop_",d), n_consensus=length(cs),
                  retained_of_28=length(intersect(cs,cons0)), jaccard=round(jac(cs,cons0),3)))
}
# Toronto-collapse
p1 <- fe_pool(c("GSE113439","GSE53408"))
FCt <- cbind(FC[,setdiff(disc,c("GSE113439","GSE53408"))], Toronto=p1$fc)
SEt <- cbind(SE[,setdiff(disc,c("GSE113439","GSE53408"))], Toronto=p1$se)
cs <- meta_consensus(FCt, SEt); scen_sets[["Toronto_collapsed"]] <- cs
report <- rbind(report, data.table(scenario="Toronto_collapsed", n_consensus=length(cs),
                retained_of_28=length(intersect(cs,cons0)), jaccard=round(jac(cs,cons0),3)))
# US-pair-collapse
p2 <- fe_pool(c("GSE117261","GSE254617"))
FCu <- cbind(FC[,setdiff(disc,c("GSE117261","GSE254617"))], USpair=p2$fc)
SEu <- cbind(SE[,setdiff(disc,c("GSE117261","GSE254617"))], USpair=p2$se)
cs <- meta_consensus(FCu, SEu); scen_sets[["USpair_collapsed"]] <- cs
report <- rbind(report, data.table(scenario="USpair_collapsed", n_consensus=length(cs),
                retained_of_28=length(intersect(cs,cons0)), jaccard=round(jac(cs,cons0),3)))

# PVOD-exclusion: re-DE GSE254617 without the 6 PVOD cases, swap its column
o <- readRDS(file.path(PATHS$proc,"GSE254617.rds")); ex<-o$expr; ph<-o$pheno
pc <- read.delim(file.path(PATHS$pheno,"GSE254617_pheno.tsv"), colClasses="character")
pc <- pc[match(ph$title, pc$title),]
drop <- pc$subtype=="PVOD"; keepS <- !drop
cd <- data.frame(group=factor(ph$group[keepS],levels=c("control","PAH")),
                 sex=factor(pc$sex[keepS]), batch=factor(make.names(pc$batch[keepS])))
kg <- edgeR::filterByExpr(ex[,keepS], group=cd$group)
dds <- DESeqDataSetFromMatrix(ex[kg,keepS], colData=cd, design=~sex+batch+group)
dds <- DESeq(dds, quiet=TRUE); rr <- as.data.frame(results(dds, contrast=c("group","PAH","control")))
fc_new <- setNames(rr$log2FoldChange, rownames(rr)); se_new <- setNames(rr$lfcSE, rownames(rr))
FCp <- FC; SEp <- SE
FCp[,"GSE254617"] <- fc_new[genes]; SEp[,"GSE254617"] <- se_new[genes]
cs <- meta_consensus(FCp, SEp); scen_sets[["PVOD_excluded"]] <- cs
report <- rbind(report, data.table(scenario="PVOD_excluded", n_consensus=length(cs),
                retained_of_28=length(intersect(cs,cons0)), jaccard=round(jac(cs,cons0),3)))

fwrite(report, file.path(PATHS$meta,"sensitivity_summary.tsv"), sep="\t")
print(report)

# per-gene robustness: fraction of the 9 scenarios in which each original consensus gene survives
nsc <- length(scen_sets)
surv <- data.table(gene=cons0,
  n_scenarios_survived = sapply(cons0, function(g) sum(sapply(scen_sets, function(s) g %in% s))),
  n_scenarios = nsc)
surv[, frac := round(n_scenarios_survived/n_scenarios, 3)]
surv <- surv[order(-frac)]
fwrite(surv, file.path(PATHS$meta,"consensus_robustness.tsv"), sep="\t")
core <- surv[n_scenarios_survived == nsc, gene]   # survives ALL scenarios = robust core
writeLines(core, file.path(PATHS$meta,"robust_core_genes.txt"))
message(sprintf("robust core (consensus in original + all %d sensitivity scenarios): %d genes", nsc, length(core)))
message(paste(core, collapse=", "))
message("DONE sensitivity")
