while (sink.number() > 0) sink()
suppressPackageStartupMessages({library(GEOquery); library(limma)})
stopifnot(basename(getwd()) == "alzheimers-pipeline-clone")
dir.create("results", showWarnings = FALSE); dir.create("old_results_v1", showWarnings = FALSE)
options(timeout = 1200)
sink("docs/02_console.txt", split = TRUE)

## 0. Keep old results (copy only, never overwrite)
w <- "/mnt/c/Users/DELL/Documents"
for (f in c("GSE44770_all_DE_results.csv", "GSE44770_DEGs_annotated.csv"))
  file.copy(file.path(w, f), file.path("old_results_v1", f), overwrite = FALSE, copy.date = TRUE)

## 1. Data and covariates
eset <- getGEO("GSE44770", GSEMatrix = TRUE, destdir = "data")[[1]]
x <- exprs(eset); pd <- pData(eset)
cvd <- data.frame(
  dx = factor(pd[["disease:ch2"]], levels = c("N", "A")),
  age = as.numeric(pd[["age:ch2"]]), sex = factor(pd[["gender:ch2"]]),
  pmi = as.numeric(pd[["pmi:ch2"]]), ph = as.numeric(pd[["ph:ch2"]]),
  rin = as.numeric(pd[["rin:ch2"]]), batch = factor(pd[["batch:ch2"]]),
  pres = factor(pd[["pres:ch2"]]), row.names = colnames(x))
stopifnot(!anyNA(cvd))
print(table(cvd$dx, cvd$pres))

## 2. Which PCs track which covariates
keep <- rowMeans(is.na(x)) < 0.1; xp <- x[keep, ]
idx <- which(is.na(xp), arr.ind = TRUE); rmed <- apply(xp, 1, median, na.rm = TRUE); xp[idx] <- rmed[idx[, 1]]
pc <- prcomp(t(xp[order(apply(xp, 1, var), decreasing = TRUE)[1:5000], ]), scale. = TRUE)
print(round(summary(pc)$importance[, 1:5], 3))
res <- sapply(c("dx", "age", "sex", "pmi", "ph", "rin", "batch"), function(v) sapply(1:5, function(k) {
  f <- cvd[[v]]; if (is.numeric(f)) cor.test(pc$x[, k], f)$p.value else kruskal.test(pc$x[, k] ~ f)$p.value }))
rownames(res) <- paste0("PC", 1:5); print(signif(res, 2))

## 3. Models (logFC is log10 ratio; 0.3 = 2-fold, 0.176 = 1.5-fold)
fit_tt <- function(form, dat, y = x) {
  d <- model.matrix(form, dat); stopifnot(qr(d)$rank == ncol(d))
  topTable(eBayes(lmFit(y, d)), coef = "dxA", number = Inf, adjust.method = "BH") }
f1 <- ~ dx + age + sex + pmi + ph + rin + batch
tt0 <- fit_tt(~ dx, cvd)
tt1 <- fit_tt(f1, cvd)
s <- cvd$age >= 55
tt2 <- fit_tt(f1, droplevels(cvd[s, ]), x[, s]); cat("age>=55 subset:\n"); print(table(cvd$dx[s]))
old <- read.csv("old_results_v1/GSE44770_all_DE_results.csv", row.names = 1, colClasses = c(probe_id = "character"))
cat("max |logFC| diff, rerun vs old:", max(abs(old[rownames(tt0), "logFC"] - tt0$logFC)), "\n")
cnt <- function(tt, l) sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > l)
tab <- sapply(list(old = old, rerun_unadj = tt0, adjusted = tt1, adj_age55plus = tt2), function(tt)
  c(probes = nrow(tt), FDR05 = sum(tt$adj.P.Val < 0.05), lfc0.3 = cnt(tt, 0.3), lfc0.176 = cnt(tt, 0.176)))
print(tab)
S <- function(tt, l = 0.3) rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) > l]
cat("old DEGs", length(S(old)), "| adjusted", length(S(tt1)), "| overlap", length(intersect(S(old), S(tt1))), "\n")
top50 <- rownames(old)[order(old$adj.P.Val)][1:50]
cat("old top-50 surviving adjustment:", sum(top50 %in% S(tt1)), "\n")
cat("logFC correlation old vs adjusted:", round(cor(old[rownames(tt1), "logFC"], tt1$logFC), 3), "\n")

## 4. Sex-chromosome genes among old vs adjusted results
fd <- fData(eset); print(names(fd))
sym <- grep("symbol|gene.?name", names(fd), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(sym)) { g <- toupper(as.character(fd[[sym]]))
  ix <- which(g %in% c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")); pr <- rownames(fd)[ix]
  print(data.frame(gene = g[ix], old_lfc = round(old[pr, "logFC"], 2), old_FDR = signif(old[pr, "adj.P.Val"], 2),
                   adj_lfc = round(tt1[pr, "logFC"], 2), adj_FDR = signif(tt1[pr, "adj.P.Val"], 2))) }

## 5. Save (tsv, because .gitignore excludes *.csv)
write.table(tt1, "results/de_adjusted.tsv", sep = "\t", quote = FALSE)
write.table(tt0, "results/de_rerun_unadjusted.tsv", sep = "\t", quote = FALSE)
writeLines(capture.output(print(tab), sessionInfo()), "docs/02_de_summary_and_session.txt")
sink()
