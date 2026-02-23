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

write_reason_table <- function(path, reason) {
  out <- data.frame(reason = reason, stringsAsFactors = FALSE)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

args <- parse_args()
genes_path <- args[["genes"]]
genesets_path <- args[["genesets"]]
qvalue_threshold <- as.numeric(ifelse(is.null(args[["qvalue-threshold"]]), "0.05", args[["qvalue-threshold"]]))
out_table <- args[["out-table"]]
out_plot <- args[["out-plot"]]

if (is.null(genes_path) || genes_path == "" || !file.exists(genes_path)) {
  write_reason_table(out_table, "variant_gene_table_missing")
  blank_plot(out_plot, "Variant gene list is missing")
  quit(status = 0)
}

if (is.null(genesets_path) || genesets_path == "" || !file.exists(genesets_path)) {
  write_reason_table(out_table, "gene_set_file_missing")
  blank_plot(out_plot, "Gene set file not provided. Set paths.gene_sets in config.")
  quit(status = 0)
}

genes_df <- read.delim(genes_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!("gene" %in% colnames(genes_df))) {
  write_reason_table(out_table, "variant_gene_table_missing_gene_column")
  blank_plot(out_plot, "Variant gene table missing 'gene' column")
  quit(status = 0)
}

genes <- unique(na.omit(trimws(genes_df$gene)))
if (length(genes) < 5) {
  write_reason_table(out_table, "insufficient_mapped_genes")
  blank_plot(out_plot, "Not enough mapped genes for enrichment")
  quit(status = 0)
}

gmt <- read.gmt(genesets_path)
enr <- enricher(
  gene = genes,
  TERM2GENE = gmt,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = qvalue_threshold
)

if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
  write_reason_table(out_table, "no_terms_detected")
  blank_plot(out_plot, "No enriched terms detected")
  quit(status = 0)
}

res <- as.data.frame(enr)
res <- subset(res, !is.na(p.adjust) & p.adjust < qvalue_threshold)
if (nrow(res) == 0) {
  write_reason_table(out_table, "no_terms_pass_threshold")
  blank_plot(out_plot, "No enriched terms passed q-value threshold")
  quit(status = 0)
}

write.table(res, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

p <- dotplot(enr, showCategory = min(20, nrow(res))) +
  ggtitle("Variant-Gene Pathway Enrichment")
ggsave(out_plot, plot = p, width = 10, height = 7, dpi = 200)
