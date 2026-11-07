suppressPackageStartupMessages(library(GEOquery))
stopifnot(basename(getwd()) == "alzheimers-pipeline-clone")
eset <- getGEO(filename = "data/GSE44770_series_matrix.txt.gz")
x <- exprs(eset); pd <- pData(eset)
meta <- data.frame(sample = colnames(x),
  label = as.integer(pd[["disease:ch2"]] == "A"),
  age = as.numeric(pd[["age:ch2"]]), sex_male = as.integer(pd[["gender:ch2"]] == "M"),
  pmi = as.numeric(pd[["pmi:ch2"]]), ph = as.numeric(pd[["ph:ch2"]]), rin = as.numeric(pd[["rin:ch2"]]),
  batch2 = as.integer(pd[["batch:ch2"]] == "2"), batch3 = as.integer(pd[["batch:ch2"]] == "3"))
stopifnot(!anyNA(meta))
write.table(meta, "data/ml_meta.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(t(x), "data/ml_expr.tsv", sep = "\t", quote = FALSE, col.names = NA)
print(table(meta$label))
