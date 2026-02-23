#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
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
de_path <- args[["de"]]
genesets_path <- args[["genesets"]]
padj_threshold <- as.numeric(ifelse(is.null(args[["padj-threshold"]]), "0.05", args[["padj-threshold"]]))
lfc_threshold <- as.numeric(ifelse(is.null(args[["lfc-threshold"]]), "1", args[["lfc-threshold"]]))
qvalue_threshold <- as.numeric(ifelse(is.null(args[["qvalue-threshold"]]), "0.05", args[["qvalue-threshold"]]))
out_table <- args[["out-table"]]
out_plot <- args[["out-plot"]]

de <- read.delim(de_path, check.names = FALSE)
required_cols <- c("gene", "log2FoldChange", "padj")
if (!all(required_cols %in% colnames(de))) {
  write.table(data.frame(), out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  blank_plot(out_plot, "DE results missing required columns for enrichment")
  quit(status = 0)
}

sig <- subset(de, !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold)
if (nrow(sig) < 5) {
  write.table(data.frame(), out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  blank_plot(out_plot, "Not enough significant genes for enrichment")
  quit(status = 0)
}

if (is.null(genesets_path) || genesets_path == "" || !file.exists(genesets_path)) {
  write.table(data.frame(), out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  blank_plot(out_plot, "Gene set file not provided. Set paths.gene_sets in config.")
  quit(status = 0)
}

gmt <- read.gmt(genesets_path)
enr <- enricher(
  gene = sig$gene,
  TERM2GENE = gmt,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = qvalue_threshold
)

if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
  write.table(data.frame(), out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  blank_plot(out_plot, "No enriched terms detected")
  quit(status = 0)
}

res <- as.data.frame(enr)
res <- subset(res, !is.na(p.adjust) & p.adjust < qvalue_threshold)
if (nrow(res) == 0) {
  write.table(data.frame(), out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  blank_plot(out_plot, "No enriched terms passed q-value threshold")
  quit(status = 0)
}
write.table(res, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

p <- dotplot(enr, showCategory = 20) +
  ggtitle("Pathway Enrichment (clusterProfiler enricher)")
ggsave(out_plot, plot = p, width = 10, height = 7, dpi = 200)
