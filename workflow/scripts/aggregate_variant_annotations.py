#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict

import pandas as pd


FEATURE_PRIORITY = [
    "CDS",
    "UTR",
    "three_prime_utr",
    "five_prime_utr",
    "exon",
    "gene",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate variant-feature overlaps into functional annotations")
    parser.add_argument("--variants", required=True, help="Base variant table TSV")
    parser.add_argument("--overlaps", required=True, help="bedtools intersect output TSV")
    parser.add_argument("--out-annotated", required=True)
    parser.add_argument("--out-genes", required=True)
    parser.add_argument("--out-summary", required=True)
    return parser.parse_args()


def pick_class(features):
    fset = set(features)
    if "CDS" in fset:
        return "CDS"
    if "UTR" in fset or "three_prime_utr" in fset or "five_prime_utr" in fset:
        return "UTR"
    if "exon" in fset:
        return "exon"
    if "gene" in fset:
        return "intronic_or_gene"
    return "intergenic"


def split_semicolon(value):
    if value is None:
        return []
    text = str(value).strip()
    if text in {"", "."}:
        return []
    return [x.strip() for x in text.split(";") if x.strip() and x.strip() != "."]


def main():
    args = parse_args()

    variants = pd.read_csv(args.variants, sep="\t", dtype=str).fillna("")
    if variants.empty:
        variants["functional_class"] = []
        variants["overlapping_features"] = []
        variants["genes"] = []
        variants["n_overlapping_genes"] = []
        variants.to_csv(args.out_annotated, sep="\t", index=False)
        pd.DataFrame(columns=["gene", "variant_count", "functional_classes"]).to_csv(
            args.out_genes, sep="\t", index=False
        )
        pd.DataFrame(columns=["group_type", "group", "count"]).to_csv(args.out_summary, sep="\t", index=False)
        return

    feature_map = defaultdict(set)
    gene_map = defaultdict(set)

    with open(args.overlaps, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 10:
                continue
            key = row[3]
            feature = row[7]
            gene_id = row[8]
            gene_name = row[9]

            if feature:
                feature_map[key].add(feature)
            gene_label = gene_name if gene_name and gene_name != "." else gene_id
            if gene_label and gene_label != ".":
                gene_map[key].add(gene_label)

    functional_classes = []
    overlap_features = []
    overlap_genes = []
    reported_genes = []
    merged_genes = []
    annotation_sources = []

    for key in variants["key"].tolist():
        feats = sorted(feature_map.get(key, set()))
        gs_overlap = sorted(gene_map.get(key, set()))
        functional_classes.append(pick_class(feats))
        overlap_features.append(";".join(feats))
        overlap_genes.append(";".join(gs_overlap))

    reported_gene_series = variants["reported_genes"] if "reported_genes" in variants.columns else pd.Series([""] * len(variants))
    for overlap_text, reported_text in zip(overlap_genes, reported_gene_series.tolist()):
        gs_overlap = set(split_semicolon(overlap_text))
        gs_reported = set(split_semicolon(reported_text))
        gs_merged = sorted(gs_overlap | gs_reported)

        reported_genes.append(";".join(sorted(gs_reported)))
        merged_genes.append(";".join(gs_merged))

        if gs_overlap and gs_reported:
            annotation_sources.append("overlap+reported_info")
        elif gs_overlap:
            annotation_sources.append("overlap")
        elif gs_reported:
            annotation_sources.append("reported_info")
        else:
            annotation_sources.append("none")

    variants["functional_class"] = functional_classes
    variants["overlapping_features"] = overlap_features
    variants["genes_overlap_tmp"] = overlap_genes
    variants["genes_overlap"] = variants["genes_overlap_tmp"]
    variants["genes_reported"] = reported_genes
    variants["genes"] = merged_genes
    variants["annotation_source"] = annotation_sources
    variants["n_overlapping_genes"] = variants["genes_overlap"].apply(lambda x: 0 if not x else len(x.split(";")))
    variants["n_genes_reported"] = variants["genes_reported"].apply(lambda x: 0 if not x else len(x.split(";")))
    variants["n_genes_total"] = variants["genes"].apply(lambda x: 0 if not x else len(x.split(";")))
    variants = variants.drop(columns=["genes_overlap_tmp"])

    variants.to_csv(args.out_annotated, sep="\t", index=False)

    gene_counts = Counter()
    gene_classes = defaultdict(set)
    for _, row in variants.iterrows():
        row_genes = [g for g in str(row["genes"]).split(";") if g]
        for g in row_genes:
            gene_counts[g] += 1
            gene_classes[g].add(row["functional_class"])

    genes_out = pd.DataFrame(
        [
            {
                "gene": gene,
                "variant_count": count,
                "functional_classes": ";".join(sorted(gene_classes[gene])),
            }
            for gene, count in gene_counts.most_common()
        ]
    )
    if genes_out.empty:
        genes_out = pd.DataFrame(columns=["gene", "variant_count", "functional_classes"])
    genes_out.to_csv(args.out_genes, sep="\t", index=False)

    class_counts = variants["functional_class"].value_counts().rename_axis("group").reset_index(name="count")
    class_counts["group_type"] = "functional_class"

    chrom_counts = variants["chrom"].value_counts().rename_axis("group").reset_index(name="count")
    chrom_counts["group_type"] = "chromosome"

    summary = pd.concat([class_counts, chrom_counts], ignore_index=True)
    summary = summary[["group_type", "group", "count"]]
    summary.to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
