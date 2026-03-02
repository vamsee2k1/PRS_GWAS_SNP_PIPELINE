#!/usr/bin/env python3
"""Heuristic AI-style variant prioritization from association + annotation outputs."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Set

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_PRIORITY_GENES = [
    "APOE",
    "TOMM40",
    "TREM2",
    "SORL1",
    "ABCA7",
    "BIN1",
    "PICALM",
    "CR1",
    "CLU",
    "PTK2B",
]


def parse_gene_list(raw: str) -> Set[str]:
    if not raw.strip():
        return {g.upper() for g in DEFAULT_PRIORITY_GENES}
    genes = {g.strip().upper() for g in raw.split(",") if g.strip()}
    return genes or {g.upper() for g in DEFAULT_PRIORITY_GENES}


def to_float_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def functional_class_score(value: str) -> float:
    text = str(value).strip().lower()
    if not text:
        return 0.2
    if any(tok in text for tok in ["splice", "stop", "frameshift", "missense", "nonsense", "cds", "exon"]):
        return 1.0
    if "utr" in text:
        return 0.7
    if "intron" in text or "intronic" in text:
        return 0.4
    if "intergenic" in text:
        return 0.1
    return 0.3


def clinical_score(value: str) -> float:
    text = str(value).strip().lower()
    if not text:
        return 0.0
    if "pathogenic" in text:
        return 1.0
    if "risk" in text:
        return 0.7
    if "uncertain" in text:
        return 0.2
    if "benign" in text:
        return 0.0
    return 0.3


def split_genes(value: str) -> List[str]:
    text = str(value).strip()
    if not text:
        return []
    return [tok.strip().upper() for tok in re.split(r"[;,/|]", text) if tok.strip()]


def load_annotation_subset(annotated_path: str, keys: Set[str]) -> pd.DataFrame:
    usecols = ["key", "functional_class", "genes", "clinical_significance"]
    chunks = []
    for chunk in pd.read_csv(
        annotated_path,
        sep="\t",
        dtype=str,
        usecols=lambda c: c in usecols,
        chunksize=200000,
    ):
        chunk = chunk.fillna("")
        sub = chunk[chunk["key"].isin(keys)]
        if not sub.empty:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=usecols)
    subset = pd.concat(chunks, ignore_index=True)
    return subset.drop_duplicates(subset=["key"], keep="first")


def render_plot(table: pd.DataFrame, top_n: int, out_plot: str) -> None:
    out_path = Path(out_plot)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    top = table.head(top_n).copy()
    if top.empty:
        fig = plt.figure(figsize=(8, 3))
        plt.text(0.5, 0.5, "No variants available for prioritization.", ha="center", va="center")
        plt.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    top["label"] = top.apply(
        lambda r: f"{r.get('id', '.') or '.'} | {r.get('genes', '') or 'NA'}",
        axis=1,
    )
    top = top.iloc[::-1]
    fig = plt.figure(figsize=(12, max(4, 0.3 * len(top) + 1)))
    plt.barh(top["label"], top["priority_score"], color="#0f766e")
    plt.xlabel("AI Priority Score (0-1)")
    plt.ylabel("Variant")
    plt.title("AI Variant Prioritization (Top Ranked Variants)")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def render_report(
    table: pd.DataFrame,
    out_report: str,
    top_n: int,
    trait_context: str,
    priority_genes: Set[str],
) -> None:
    out_path = Path(out_report)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    top = table.head(top_n).copy()
    lines = [
        "# AI Variant Prioritization Report",
        "",
        "## Context",
        "",
        f"- trait_context: `{trait_context or 'unspecified'}`",
        f"- priority_gene_set_size: `{len(priority_genes)}`",
        f"- evaluated_variants: `{len(table)}`",
        "",
        "## Top Ranked Variants",
        "",
    ]
    if top.empty:
        lines.append("- No variants available.")
    else:
        for _, row in top.iterrows():
            lines.append(
                f"- `{row.get('key', '.')}` | score `{row['priority_score']:.3f}` | "
                f"tier `{row['priority_tier']}` | rationale: {row['priority_rationale']}"
            )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prioritize variants using heuristic AI-style scoring.")
    parser.add_argument("--annotated", required=True, help="annotated_variants.tsv path")
    parser.add_argument("--association", required=True, help="variant_association_hits.tsv path")
    parser.add_argument("--out-table", required=True, help="Output prioritized TSV path")
    parser.add_argument("--out-plot", required=True, help="Output prioritized plot path")
    parser.add_argument("--out-report", required=True, help="Output markdown report path")
    parser.add_argument("--top-n", type=int, default=50, help="Top variants count for plot/report")
    parser.add_argument("--trait-context", default="", help="Trait context string")
    parser.add_argument(
        "--priority-genes",
        default="",
        help="Comma-separated high-priority genes used for contextual scoring.",
    )
    args = parser.parse_args()

    assoc = pd.read_csv(args.association, sep="\t", dtype=str).fillna("")
    if "key" not in assoc.columns:
        raise ValueError("Association table must include a 'key' column.")

    keys = set(assoc["key"].astype(str).tolist())
    ann = load_annotation_subset(args.annotated, keys)

    merged = assoc.merge(ann, on="key", how="left", suffixes=("", "_ann"))
    for col in ["functional_class", "genes", "clinical_significance"]:
        if col not in merged.columns:
            merged[col] = ""

    merged["minus_log10_p"] = to_float_series(merged.get("minus_log10_p", pd.Series(dtype=str)))
    merged["pvalue_num"] = to_float_series(merged.get("pvalue", pd.Series(dtype=str)))
    merged["effect_num"] = to_float_series(merged.get("effect", pd.Series(dtype=str)))
    merged["significant_flag"] = (
        merged.get("significant", pd.Series(dtype=str)).astype(str).str.lower().isin(["true", "1", "yes"])
    )

    missing_minus = merged["minus_log10_p"].isna() & merged["pvalue_num"].notna() & (merged["pvalue_num"] > 0)
    merged.loc[missing_minus, "minus_log10_p"] = -merged.loc[missing_minus, "pvalue_num"].map(math.log10)

    merged["sig_score"] = (merged["minus_log10_p"].fillna(0).clip(lower=0, upper=50) / 50.0)
    merged["effect_score"] = (merged["effect_num"].abs().fillna(0).clip(upper=2.0) / 2.0)
    merged["functional_score"] = merged["functional_class"].map(functional_class_score)
    merged["clinical_score"] = merged["clinical_significance"].map(clinical_score)

    priority_genes = parse_gene_list(args.priority_genes)
    merged["gene_tokens"] = merged["genes"].map(split_genes)
    merged["priority_gene_hit"] = merged["gene_tokens"].map(lambda vals: any(g in priority_genes for g in vals))
    merged["gene_score"] = merged["priority_gene_hit"].map(lambda hit: 1.0 if hit else 0.2)
    merged["assoc_hit_score"] = merged["significant_flag"].map(lambda v: 1.0 if v else 0.4)

    merged["priority_score"] = (
        0.35 * merged["sig_score"]
        + 0.20 * merged["effect_score"]
        + 0.15 * merged["functional_score"]
        + 0.10 * merged["clinical_score"]
        + 0.15 * merged["gene_score"]
        + 0.05 * merged["assoc_hit_score"]
    ).clip(lower=0.0, upper=1.0)

    merged["priority_tier"] = pd.cut(
        merged["priority_score"],
        bins=[-1, 0.5, 0.75, 1.01],
        labels=["low", "medium", "high"],
    ).astype(str)

    def rationale(row: pd.Series) -> str:
        reasons: List[str] = []
        if row["priority_gene_hit"]:
            reasons.append("priority_gene_context")
        if row["sig_score"] >= 0.6:
            reasons.append("strong_association_signal")
        if row["functional_score"] >= 0.7:
            reasons.append("coding_or_regulatory_class")
        if row["clinical_score"] >= 0.7:
            reasons.append("clinical_evidence")
        if not reasons:
            reasons.append("composite_signal")
        return ",".join(reasons)

    merged["priority_rationale"] = merged.apply(rationale, axis=1)

    merged = merged.sort_values(
        by=["priority_score", "minus_log10_p", "effect_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)

    ordered_cols = [
        "key",
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "genes",
        "functional_class",
        "clinical_significance",
        "association_mode",
        "pvalue",
        "minus_log10_p",
        "effect",
        "significant",
        "priority_tier",
        "priority_score",
        "priority_rationale",
    ]
    existing_cols = [c for c in ordered_cols if c in merged.columns]
    merged_out = merged[existing_cols + [c for c in merged.columns if c not in existing_cols]]

    out_table = Path(args.out_table)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    merged_out.to_csv(out_table, sep="\t", index=False)

    render_plot(merged_out, top_n=max(1, args.top_n), out_plot=args.out_plot)
    render_report(
        merged_out,
        out_report=args.out_report,
        top_n=max(1, args.top_n),
        trait_context=args.trait_context,
        priority_genes=priority_genes,
    )


if __name__ == "__main__":
    main()

