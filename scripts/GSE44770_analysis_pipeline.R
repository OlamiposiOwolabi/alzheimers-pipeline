# =============================================================================
# GSE44770 Alzheimer's Disease Transcriptomic Analysis Pipeline
# Dataset: NCBI GEO - GSE44770 (Prefrontal Cortex, 230 samples)
# Platform: GPL4372 (Rosetta/Merck Human 44k 1.1 microarray)
# Author: [Owolabi Olamiposi Peter]
# Date: May 2026
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 1: Install and Load Required Packages
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "Biobase",
                        "clusterProfiler", "org.Hs.eg.db",
                        "AnnotationDbi", "DOSE", "enrichplot"),
                     ask = FALSE)

install.packages(c("pheatmap"))

library(GEOquery)
library(limma)
library(Biobase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)


# -----------------------------------------------------------------------------
# SECTION 2: Data Acquisition
# -----------------------------------------------------------------------------

# Download GSE44770 series matrix from NCBI GEO
gse <- getGEO("GSE44770", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Extract expression matrix and sample metadata
data <- exprs(gse[[1]])       # genes x samples expression matrix
meta <- pData(gse[[1]])       # sample phenotype metadata

# Inspect data dimensions and metadata
dim(data)
head(colnames(meta))
table(meta$characteristics_ch2)   # disease: A (Alzheimer) vs N (Normal)

# Check array platform
annotation(gse[[1]])   # GPL4372


# -----------------------------------------------------------------------------
# SECTION 3: Differential Expression Analysis (limma)
# -----------------------------------------------------------------------------

# Extract and clean group labels from metadata
group <- factor(gsub("disease: ", "", meta$characteristics_ch2))
table(group)   # A = 129, N = 101

# Build design matrix (no intercept)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(data, design)

# Define contrast: Alzheimer's vs Normal
contrast <- makeContrasts(A - N, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract all results, BH-adjusted p-values
results <- topTable(fit2, adjust = "BH", number = Inf)

# Summary of logFC distribution (microarray: compressed range expected)
summary(results$logFC)

# Apply thresholds: |logFC| > 0.3 and adj.P < 0.05
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.3, ]
cat("Total significant DEGs:", nrow(sig_genes), "\n")

# Split into up and downregulated
upregulated   <- sig_genes[sig_genes$logFC >  0.3, ]
downregulated <- sig_genes[sig_genes$logFC < -0.3, ]
cat("Upregulated:", nrow(upregulated), "\n")
cat("Downregulated:", nrow(downregulated), "\n")

# Add significance labels for plotting
results$sig <- "Not Significant"
results$sig[results$adj.P.Val < 0.05 & results$logFC >  0.3] <- "Up"
results$sig[results$adj.P.Val < 0.05 & results$logFC < -0.3] <- "Down"


# -----------------------------------------------------------------------------
# SECTION 4: Volcano Plot
# -----------------------------------------------------------------------------

png("GSE44770_volcano.png", width = 2400, height = 1800, res = 300)
plot(results$logFC, -log10(results$adj.P.Val),
     col = ifelse(results$sig == "Up", "red",
            ifelse(results$sig == "Down", "blue", "grey")),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change",
     ylab = "-log10 adjusted P-value",
     main = "Volcano Plot - GSE44770 (AD vs Control)")
abline(v = c(-0.3, 0.3), lty = 2)
abline(h = -log10(0.05), lty = 2)
legend("topright",
       legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c("red", "blue", "grey"), pch = 20, cex = 0.7)
dev.off()


# -----------------------------------------------------------------------------
# SECTION 5: Gene Annotation (Probe ID → Gene Symbol + Entrez ID)
# -----------------------------------------------------------------------------

# Extract feature annotation embedded in GEO object
fData_df <- fData(gse[[1]])
colnames(fData_df)

# Select relevant annotation columns
probe_to_entrez <- fData_df[, c("ID", "Gene ID", "Gene symbol")]

# Add probe IDs as column for merging
results$probe_id   <- rownames(results)
sig_genes$probe_id <- rownames(sig_genes)

# Merge DEG results with annotation
sig_annotated <- merge(sig_genes, probe_to_entrez,
                       by.x = "probe_id", by.y = "ID")

# Remove probes without Entrez IDs
sig_annotated <- sig_annotated[
  !is.na(sig_annotated$`Gene ID`) & sig_annotated$`Gene ID` != "", ]

cat("Annotated DEGs:", nrow(sig_annotated), "\n")


# -----------------------------------------------------------------------------
# SECTION 6: Heatmap of Top 50 DEGs
# -----------------------------------------------------------------------------

# Get top 50 DEGs by adjusted p-value
top50 <- head(sig_annotated[order(sig_annotated$adj.P.Val), ], 50)

# Extract expression data for top 50 probes
heatmap_data <- data[top50$probe_id, ]
rownames(heatmap_data) <- top50$`Gene symbol`

# Sample annotation for heatmap columns
sample_annotation <- data.frame(
  Disease = gsub("disease: ", "", meta$characteristics_ch2),
  row.names = colnames(heatmap_data)
)

# Save heatmap
png("GSE44770_heatmap.png", width = 2400, height = 3000, res = 300)
pheatmap(heatmap_data,
         annotation_col = sample_annotation,
         scale = "row",
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 DEGs - GSE44770 (AD vs Control)")
dev.off()


# -----------------------------------------------------------------------------
# SECTION 7: KEGG Pathway Enrichment Analysis
# -----------------------------------------------------------------------------

# Background gene universe (all annotated probes)
all_annotated <- merge(results, probe_to_entrez, by.x = "probe_id", by.y = "ID")
background <- as.character(all_annotated$`Gene ID`)
background <- background[!is.na(background) & background != ""]

# DEG Entrez IDs
entrez_ids <- as.character(sig_annotated$`Gene ID`)

# Run KEGG enrichment (all DEGs)
kegg_result <- enrichKEGG(gene         = entrez_ids,
                           organism     = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           universe     = background)

summary_kegg <- as.data.frame(kegg_result)
cat("Enriched KEGG pathways:", nrow(summary_kegg), "\n")
head(summary_kegg[, c("Description", "GeneRatio", "pvalue", "p.adjust")], 10)

# Separate upregulated and downregulated KEGG enrichment
up_entrez   <- as.character(sig_annotated$`Gene ID`[sig_annotated$logFC >  0.3])
down_entrez <- as.character(sig_annotated$`Gene ID`[sig_annotated$logFC < -0.3])

kegg_up   <- enrichKEGG(gene = up_entrez,   organism = "hsa",
                         pvalueCutoff = 0.05, universe = background)
kegg_down <- enrichKEGG(gene = down_entrez, organism = "hsa",
                         pvalueCutoff = 0.05, universe = background)

cat("Upregulated pathways:", nrow(as.data.frame(kegg_up)), "\n")
cat("Downregulated pathways:", nrow(as.data.frame(kegg_down)), "\n")

# KEGG visualizations
pdf("GSE44770_KEGG_dotplot.pdf", width = 10, height = 8)
dotplot(kegg_result, showCategory = 20,
        title = "KEGG Pathway Enrichment - GSE44770 (AD vs Control)")
dev.off()

pdf("GSE44770_KEGG_barplot.pdf", width = 10, height = 8)
barplot(kegg_result, showCategory = 20,
        title = "KEGG Pathway Enrichment - GSE44770 (AD vs Control)")
dev.off()


# -----------------------------------------------------------------------------
# SECTION 8: Export Results
# -----------------------------------------------------------------------------

write.csv(results,        "GSE44770_all_DE_results.csv",       row.names = TRUE)
write.csv(sig_annotated,  "GSE44770_DEGs_annotated.csv",       row.names = FALSE)
write.csv(summary_kegg,   "GSE44770_KEGG_all.csv",             row.names = FALSE)
write.csv(as.data.frame(kegg_up),   "GSE44770_KEGG_up.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_down), "GSE44770_KEGG_down.csv", row.names = FALSE)

cat("All results saved to:", getwd(), "\n")


# -----------------------------------------------------------------------------
# SECTION 9: STRING PPI Network Analysis
# -----------------------------------------------------------------------------

library(STRINGdb)

# Extract gene symbols of significant DEGs
deg_symbols <- sig_annotated$`Gene symbol`
deg_symbols <- deg_symbols[deg_symbols != "" & !is.na(deg_symbols)]
cat("Total genes for STRING:", length(deg_symbols), "\n")

# Save gene list
write.table(deg_symbols,
            "GSE44770_DEG_symbols_for_STRING.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Initialize STRING database (human = 9606, score threshold = 400)
string_db <- STRINGdb$new(version = "11.5",
                           species = 9606,
                           score_threshold = 400,
                           input_directory = "")

# Map gene symbols to STRING IDs
deg_df     <- data.frame(gene = deg_symbols)
deg_mapped <- string_db$map(deg_df, "gene", removeUnmappedRows = TRUE)
cat("Genes mapped to STRING:", nrow(deg_mapped), "\n")

# Get interactions between mapped genes
interactions <- string_db$get_interactions(deg_mapped$STRING_id)
cat("Total interactions:", nrow(interactions), "\n")

# Plot and save PPI network
png("GSE44770_PPI_network.png", width = 3000, height = 3000, res = 300)
string_db$plot_network(deg_mapped$STRING_id)
dev.off()

# -----------------------------------------------------------------------------
# SECTION 10: Hub Gene Identification
# -----------------------------------------------------------------------------

# Calculate degree (number of connections) for each gene
all_ids      <- c(interactions$from, interactions$to)
degree_table <- as.data.frame(table(all_ids))
colnames(degree_table) <- c("STRING_id", "degree")

# Merge with gene symbols and remove duplicates
hub_genes       <- merge(degree_table, deg_mapped, by = "STRING_id")
hub_genes_clean <- hub_genes[!duplicated(hub_genes$gene), ]
hub_genes_clean <- hub_genes_clean[order(-hub_genes_clean$degree), ]

# Top 20 hub genes
top20_clean <- head(hub_genes_clean, 20)
print(top20_clean[, c("gene", "degree")])

# Save hub genes
write.csv(hub_genes_clean, "GSE44770_PPI_hub_genes.csv", row.names = FALSE)

# Bar plot of top 20 hub genes
png("GSE44770_hub_genes.png", width = 2400, height = 1800, res = 300)
barplot(top20_clean$degree,
        names.arg = top20_clean$gene,
        las = 2,
        col = "steelblue",
        main = "Top 20 Hub Genes - GSE44770 PPI Network",
        ylab = "Number of Interactions (Degree)",
        cex.names = 0.8)
dev.off()

cat("All PPI results saved to:", getwd(), "\n")
