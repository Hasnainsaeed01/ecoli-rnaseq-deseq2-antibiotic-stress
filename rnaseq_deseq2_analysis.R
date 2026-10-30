#!/usr/bin/env Rscript
# ==============================================================================
# RNA-seq Differential Expression Analysis with DESeq2
# Dataset: GSE220559, Escherichia coli K-12 MG1655 under antibiotic stress
# Contrast: kanamycin-treated vs water control (3 replicates each)
# Reference: Bie et al. 2023, Microbiology Spectrum, doi:10.1128/spectrum.00317-23
# Author: Muhammad Hasnain Saeed
#
# Run from the repository root with:
#   Rscript rnaseq_deseq2_analysis.R
#
# Inputs (data/):
#   Per-sample read-count tables from GEO (GSE220559 supplementary files)
#   ecoli_k12_gene_annotation.csv : b-number to Entrez ID and gene symbol map
# Outputs:
#   figures/ : publication style QC and DEG plots (PNG)
#   results/ : DEG tables, enrichment tables, session info (CSV/TXT)
# ==============================================================================

# --- 0. Setup -----------------------------------------------------------------
set.seed(42)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.EcK12.eg.db)
  library(enrichplot)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
})

# --- 1. Load count files and build metadata -----------------------------------
cat("== 1. Loading count data ==\n")
# GEO sample -> file -> treatment group (verified from GSE220559 series matrix)
sample_map <- data.frame(
  geo_accession = c("GSM6806878", "GSM6806879", "GSM6806880",
                    "GSM6806890", "GSM6806891", "GSM6806892"),
  file = c("data/GSM6806878_KAN_2.txt.gz", "data/GSM6806879_KAN_3.txt.gz",
           "data/GSM6806880_KAN_4.txt.gz", "data/GSM6806890_Wu_2.txt.gz",
           "data/GSM6806891_Wu_3.txt.gz", "data/GSM6806892_Wu_4.txt.gz"),
  sample = c("kanamycin_1", "kanamycin_2", "kanamycin_3",
             "water_1", "water_2", "water_3"),
  treatment = c("kanamycin", "kanamycin", "kanamycin",
                "water", "water", "water"),
  stringsAsFactors = FALSE
)
write.table(sample_map, "results/sample_group_assignment.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("Sample to group assignment logged in results/sample_group_assignment.txt\n")
print(sample_map[, c("geo_accession", "sample", "treatment")], row.names = FALSE)
stopifnot(all(file.exists(sample_map$file)))

read_counts <- function(f) {
  d <- read.table(gzfile(f), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  setNames(d[[2]], d[[1]])  # readcount column, named by b-number locus tag
}
count_list <- lapply(setNames(sample_map$file, sample_map$sample), read_counts)

# Inspect: sample names, genes per sample, identical gene order
cat("Samples:", paste(names(count_list), collapse = ", "), "\n")
cat("Genes per sample:", paste(sapply(count_list, length), collapse = ", "), "\n")
same_order <- all(sapply(count_list, function(x) identical(names(x), names(count_list[[1]]))))
cat("Identical gene order in all files:", same_order, "\n")
stopifnot(same_order)

counts_mat <- do.call(cbind, count_list)
storage.mode(counts_mat) <- "integer"

# Metadata: experimental groups and replicates
metadata <- data.frame(
  sample = colnames(counts_mat),
  treatment = factor(sample_map$treatment, levels = c("water", "kanamycin")),
  replicate = rep(1:3, 2),
  row.names = colnames(counts_mat),
  stringsAsFactors = FALSE
)
cat("Experimental groups:\n"); print(table(metadata$treatment))
cat("Missing values in metadata:", sum(is.na(metadata)), "\n")
cat("Count-matrix columns match metadata rows:",
    all(colnames(counts_mat) == rownames(metadata)), "\n\n")

# --- 2. Low-count filtering ----------------------------------------------------
cat("== 2. Low-count filtering ==\n")
cat("Genes before filtering:", nrow(counts_mat), "\n")
keep <- rowSums(counts_mat) >= 10
counts_mat <- counts_mat[keep, ]
cat("Genes retained (rowSums >= 10):", nrow(counts_mat), "\n\n")

# --- 3. DESeq2 object, normalization, VST --------------------------------------
cat("== 3. DESeq2 normalization ==\n")
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = metadata,
                              design = ~ treatment)
dds <- DESeq(dds)
cat("Size factors:\n"); print(round(sizeFactors(dds), 3))

norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)
cat("VST matrix dimensions:", nrow(assay(vsd)), "x", ncol(assay(vsd)), "\n\n")

# --- 4. Quality control plots ---------------------------------------------------
cat("== 4. QC plots ==\n")

# 4a. Library size barplot
lib_size <- data.frame(sample = colnames(dds),
                       total_counts = colSums(counts(dds)),
                       treatment = metadata$treatment)
p_lib <- ggplot(lib_size, aes(x = sample, y = total_counts / 1e6, fill = treatment)) +
  geom_bar(stat = "identity") +
  labs(title = "Library size per sample (E. coli K-12, GSE220559)",
       x = "Sample", y = "Total counts (millions)", fill = "Treatment") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/01_library_size.png", p_lib, width = 8, height = 5, dpi = 150)

# 4b. Count distribution: boxplots of log10(counts + 1) per sample,
#     raw counts vs DESeq2-normalized counts
dist_df <- rbind(
  data.frame(sample = rep(colnames(counts_mat), each = nrow(counts_mat)),
             treatment = rep(metadata$treatment, each = nrow(counts_mat)),
             value = as.vector(log10(counts_mat + 1)),
             type = "raw counts"),
  data.frame(sample = rep(colnames(norm_counts), each = nrow(norm_counts)),
             treatment = rep(metadata$treatment, each = nrow(norm_counts)),
             value = as.vector(log10(norm_counts + 1)),
             type = "normalized counts")
)
dist_df$type <- factor(dist_df$type, levels = c("raw counts", "normalized counts"))
p_dist <- ggplot(dist_df, aes(x = sample, y = value, fill = treatment)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ type, ncol = 2, scales = "free_y") +
  labs(title = "Count distribution per sample (log10 scale)",
       x = "Sample", y = "log10(count + 1)", fill = "Treatment") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/02_count_distribution.png", p_dist, width = 10, height = 5, dpi = 150)

# 4c. PCA on VST counts
p_pca <- plotPCA(vsd, intgroup = "treatment") +
  ggtitle("PCA of VST-transformed counts") + theme_bw()
ggsave("figures/03_pca.png", p_pca, width = 7, height = 5, dpi = 150)

# 4d. Hierarchical clustering dendrogram (VST, Euclidean distance, complete linkage)
sample_dists <- dist(t(assay(vsd)))
hc <- hclust(sample_dists, method = "complete")
png("figures/04_sample_dendrogram.png", width = 1050, height = 700, res = 150)
plot(hc, main = "Hierarchical clustering of samples (VST, Euclidean, complete linkage)",
     xlab = "", sub = "", cex = 0.9)
dev.off()

# 4e. Sample-to-sample Euclidean distance heatmap
dmat <- as.matrix(sample_dists)
rownames(dmat) <- colnames(dmat) <- paste(colnames(vsd), metadata$treatment, sep = ":")
png("figures/05_sample_distance_heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap(dmat, clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample Euclidean distances (VST)")
dev.off()

# 4f. Sample correlation heatmap
cor_mat <- cor(assay(vsd))
png("figures/06_sample_correlation_heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap(cor_mat, main = "Sample-to-sample correlation (VST)",
         col = colorRampPalette(brewer.pal(9, "YlOrRd"))(255),
         display_numbers = TRUE, number_format = "%.2f", fontsize_number = 10)
dev.off()
cat("QC plots saved to figures/\n\n")

# --- 5. Differential expression: kanamycin vs water ------------------------------
cat("== 5. Differential expression (kanamycin vs water) ==\n")
res <- results(dds, contrast = c("treatment", "kanamycin", "water"))
res <- res[order(res$padj), ]
n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
cat("Genes tested:", nrow(res), "\n")
cat("Significant DEGs (padj < 0.05):", n_sig, "\n\n")

# 5a. MA plot
png("figures/07_ma_plot.png", width = 1050, height = 750, res = 150)
plotMA(res, ylim = c(-4, 4), main = "MA plot: kanamycin vs water")
dev.off()

# 5b. Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$sig <- "not significant"
res_df$sig[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "significant"
p_volc <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("not significant" = "grey60", "significant" = "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano plot: kanamycin vs water",
       x = "log2 fold change", y = "-log10 adjusted p-value", color = "") +
  theme_bw()
ggsave("figures/08_volcano_plot.png", p_volc, width = 8, height = 6, dpi = 150)

# 5c. Top-50 DEG heatmap (VST, row-scaled)
top50 <- head(rownames(res)[!is.na(res$padj) & res$padj < 0.05], 50)
mat50 <- assay(vsd)[top50, ]
png("figures/09_top50_deg_heatmap.png", width = 1200, height = 1500, res = 150)
pheatmap(mat50, scale = "row", show_rownames = FALSE,
         annotation_col = metadata[, "treatment", drop = FALSE],
         main = "Top 50 DEGs (VST counts, row-scaled)")
dev.off()
cat("DEG plots saved to figures/\n\n")

# --- 6. Annotation and DEG tables -------------------------------------------------
cat("== 6. Gene annotation and DEG tables ==\n")
ann <- read.csv("data/ecoli_k12_gene_annotation.csv", stringsAsFactors = FALSE)
rownames(ann) <- ann$b_number

sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
up <- sig[sig$log2FoldChange > 0, ]
down <- sig[sig$log2FoldChange < 0, ]
cat("Upregulated (padj < 0.05, log2FC > 0):", nrow(up), "\n")
cat("Downregulated (padj < 0.05, log2FC < 0):", nrow(down), "\n")

annotate_df <- function(df) {
  df$b_number <- rownames(df)
  df$entrez_id <- ann[rownames(df), "entrez_id"]
  df$symbol <- ann[rownames(df), "symbol"]
  df[, c("b_number", "symbol", "entrez_id", "baseMean",
         "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}
write.csv(annotate_df(as.data.frame(res)), "results/all_results_annotated.csv", row.names = FALSE)
write.csv(annotate_df(as.data.frame(up)), "results/upregulated_degs.csv", row.names = FALSE)
write.csv(annotate_df(as.data.frame(down)), "results/downregulated_degs.csv", row.names = FALSE)
cat("CSVs written to results/\n\n")

# --- 7. GO and KEGG enrichment ------------------------------------------------------
cat("== 7. Enrichment analysis ==\n")
entrez_ids <- na.omit(unique(ann[rownames(sig), "entrez_id"]))
bnumbers_sig <- rownames(sig)

if (length(entrez_ids) >= 10) {
  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(gene = as.character(entrez_ids), OrgDb = org.EcK12.eg.db,
                    ont = ont, pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, readable = TRUE)
    etab <- as.data.frame(ego)
    tag <- tolower(ont)
    write.csv(etab, paste0("results/go_", tag, "_enrichment.csv"), row.names = FALSE)
    cat("GO", ont, "significant terms:", nrow(etab), "\n")
    if (nrow(etab) > 0) {
      png(paste0("figures/10_go_", tag, "_dotplot.png"), width = 1200, height = 900, res = 150)
      print(dotplot(ego, showCategory = 10) + ggtitle(paste("GO", ont, "dot plot: kanamycin vs water")))
      dev.off()
      png(paste0("figures/11_go_", tag, "_barplot.png"), width = 1200, height = 900, res = 150)
      print(barplot(ego, showCategory = 10) + ggtitle(paste("GO", ont, "bar plot: kanamycin vs water")))
      dev.off()
    }
  }
} else {
  cat("Too few mapped ENTREZ IDs for GO enrichment; skipping.\n")
}

# KEGG enrichment (organism eco, b-number locus tags; needs the KEGG web service)
kegg_ok <- tryCatch({
  ekegg <- enrichKEGG(gene = bnumbers_sig, organism = "eco", pvalueCutoff = 0.05)
  ktab <- as.data.frame(ekegg)
  write.csv(ktab, "results/kegg_enrichment.csv", row.names = FALSE)
  cat("KEGG significant pathways:", nrow(ktab), "\n")
  if (nrow(ktab) > 0) {
    png("figures/12_kegg_dotplot.png", width = 1200, height = 900, res = 150)
    print(dotplot(ekegg, showCategory = 10) + ggtitle("KEGG pathway dot plot: kanamycin vs water"))
    dev.off()
    png("figures/13_kegg_barplot.png", width = 1200, height = 900, res = 150)
    print(barplot(ekegg, showCategory = 10) + ggtitle("KEGG pathway bar plot: kanamycin vs water"))
    dev.off()
  }
  TRUE
}, error = function(e) {
  cat("KEGG enrichment skipped (web service unreachable):", conditionMessage(e), "\n")
  FALSE
})
if (!kegg_ok) {
  writeLines(paste0("KEGG enrichment was skipped because the KEGG web service could not be reached ",
                    "during this run. GO results are unaffected."),
             "results/kegg_enrichment_NOTE.txt")
}

# --- 8. Session info ------------------------------------------------------------------
sink("results/session_info.txt")
sessionInfo()
sink()
cat("\nSession info written. Pipeline finished.\n")
