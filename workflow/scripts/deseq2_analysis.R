#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i + 1 <= length(args)) args[[i + 1]] else ""
    out[[gsub("^--", "", key)]] <- val
    i <- i + 2
  }
  out
}

blank_plot <- function(path, label) {
  png(path, width = 1200, height = 900, res = 150)
  plot.new()
  text(0.5, 0.5, label, cex = 1.2)
  dev.off()
}

args <- parse_args()
counts_path <- args[["counts"]]
metadata_path <- args[["metadata"]]
formula_str <- args[["formula"]]
contrast_str <- args[["contrast"]]
padj_threshold <- as.numeric(ifelse(is.null(args[["padj-threshold"]]), "0.05", args[["padj-threshold"]]))
lfc_threshold <- as.numeric(ifelse(is.null(args[["lfc-threshold"]]), "1", args[["lfc-threshold"]]))
out_results <- args[["out-results"]]
out_norm <- args[["out-normalized"]]
out_volcano <- args[["out-volcano"]]
out_heatmap <- args[["out-heatmap"]]

counts_df <- read.delim(counts_path, comment.char = "#", check.names = FALSE)
if (ncol(counts_df) < 7) {
  stop("featureCounts output is malformed; expected at least 7 columns")
}

count_mat <- counts_df[, 7:ncol(counts_df), drop = FALSE]
rownames(count_mat) <- counts_df$Geneid

clean_names <- basename(colnames(count_mat))
clean_names <- sub("\\.sorted\\.bam$", "", clean_names)
clean_names <- sub("\\.bam$", "", clean_names)
colnames(count_mat) <- clean_names

meta <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!("sample" %in% colnames(meta))) {
  stop("metadata must include a 'sample' column")
}

common <- intersect(colnames(count_mat), meta$sample)
if (length(common) < 2) {
  stop("Need at least 2 overlapping samples between counts and metadata")
}

count_mat <- count_mat[, common, drop = FALSE]
meta <- meta[match(common, meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample

if (!("condition" %in% colnames(meta)) || length(unique(meta$condition)) < 2) {
  write.table(
    data.frame(),
    out_results,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(count_mat, out_norm, sep = "\t", quote = FALSE, col.names = NA)
  blank_plot(out_volcano, "No DE test: metadata requires >=2 condition groups")
  blank_plot(out_heatmap, "No heatmap: metadata requires >=2 condition groups")
  quit(status = 0)
}

design_formula <- as.formula(formula_str)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mat)),
  colData = meta,
  design = design_formula
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

contrast_parts <- strsplit(contrast_str, ",")[[1]]
if (length(contrast_parts) != 3) {
  stop("contrast must be formatted as 'factor,levelA,levelB'")
}

res <- results(dds, contrast = contrast_parts)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[order(res_df$padj), ]
write.table(res_df, out_results, sep = "\t", quote = FALSE, row.names = FALSE)

norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$gene <- rownames(norm_counts)
write.table(norm_counts, out_norm, sep = "\t", quote = FALSE, row.names = FALSE)

res_df$sig <- !is.na(res_df$padj) & res_df$padj < padj_threshold & abs(res_df$log2FoldChange) >= lfc_threshold
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.4, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728")) +
  labs(
    title = "Differential Expression Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significant"
  ) +
  theme_minimal(base_size = 12)

ggsave(out_volcano, plot = volcano, width = 8, height = 6, dpi = 200)

sig_genes <- subset(res_df, !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold)
if (nrow(sig_genes) > 2) {
  top_n <- head(sig_genes$gene, 50)
  mat <- as.matrix(norm_counts[norm_counts$gene %in% top_n, common, drop = FALSE])
  rownames(mat) <- norm_counts$gene[norm_counts$gene %in% top_n]
  log_mat <- log2(mat + 1)

  png(out_heatmap, width = 1400, height = 1200, res = 180)
  pheatmap(
    log_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = meta[, setdiff(colnames(meta), "sample"), drop = FALSE],
    main = "Top 50 Differentially Expressed Genes"
  )
  dev.off()
} else {
  blank_plot(out_heatmap, "No significant DE genes for top-50 heatmap")
}
