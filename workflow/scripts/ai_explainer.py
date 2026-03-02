#!/usr/bin/env python3
"""Generate an AI-style narrative report from pipeline output tables."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


def maybe_read_table(path: str) -> Optional[pd.DataFrame]:
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    return pd.read_csv(p, sep="\t", dtype=str).fillna("")


def safe_float(value: str) -> Optional[float]:
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def metric_map(df: Optional[pd.DataFrame]) -> Dict[str, str]:
    if df is None or "metric" not in df.columns or "value" not in df.columns:
        return {}
    out: Dict[str, str] = {}
    for _, row in df.iterrows():
        key = str(row["metric"]).strip()
        val = str(row["value"]).strip()
        if key:
            out[key] = val
    return out


def top_genes_from_assoc(assoc: Optional[pd.DataFrame], limit: int = 8) -> List[str]:
    if assoc is None or assoc.empty:
        return []
    if "genes" not in assoc.columns:
        return []
    counts: Dict[str, int] = {}
    for value in assoc["genes"].astype(str):
        for token in [t.strip() for t in value.replace("/", ";").split(";") if t.strip()]:
            counts[token] = counts.get(token, 0) + 1
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [g for g, _ in ranked[:limit]]


def format_optional_float(value: Optional[float], digits: int = 3) -> str:
    if value is None:
        return "n/a"
    return f"{value:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate AI narrative summary for pipeline outputs.")
    parser.add_argument("--variant-metrics", required=True, help="variant_qc_metrics.tsv")
    parser.add_argument("--association", required=True, help="variant_association_hits.tsv")
    parser.add_argument("--enrichment", required=True, help="variant_enrichment.tsv")
    parser.add_argument("--preflight", default="", help="preflight_checks.tsv (optional)")
    parser.add_argument("--prs-qc", default="", help="prs_qc.tsv (optional)")
    parser.add_argument("--anomaly-table", default="", help="AI anomaly TSV (optional)")
    parser.add_argument("--prioritization-table", default="", help="AI prioritization TSV (optional)")
    parser.add_argument("--mode", default="", help="Pipeline run.mode value")
    parser.add_argument("--analysis-mode", default="", help="Variant analysis mode")
    parser.add_argument("--trait-context", default="", help="Free-text trait context")
    parser.add_argument("--out-report", required=True, help="Output markdown report path")
    args = parser.parse_args()

    variant_metrics_df = maybe_read_table(args.variant_metrics)
    assoc_df = maybe_read_table(args.association)
    enrich_df = maybe_read_table(args.enrichment)
    preflight_df = maybe_read_table(args.preflight)
    prs_qc_df = maybe_read_table(args.prs_qc)
    anomaly_df = maybe_read_table(args.anomaly_table)
    prioritization_df = maybe_read_table(args.prioritization_table)

    vm = metric_map(variant_metrics_df)
    ts_tv = safe_float(vm.get("ts_tv_ratio", ""))
    total_variants = safe_float(vm.get("total_variants", ""))
    sample_count = safe_float(vm.get("sample_count", ""))

    assoc_count = 0 if assoc_df is None else len(assoc_df)
    significant_count = 0
    top_assoc = None
    if assoc_df is not None and not assoc_df.empty:
        sig = assoc_df.get("significant", pd.Series(dtype=str)).astype(str).str.lower().isin(["true", "1", "yes"])
        significant_count = int(sig.sum())
        top_assoc = assoc_df.iloc[0]

    top_genes = top_genes_from_assoc(assoc_df)

    enrich_top = None
    if enrich_df is not None and not enrich_df.empty:
        if "qvalue" in enrich_df.columns:
            qvals = pd.to_numeric(enrich_df["qvalue"], errors="coerce")
            enrich_top = enrich_df.loc[qvals.fillna(999999).idxmin()] if not qvals.empty else enrich_df.iloc[0]
        else:
            enrich_top = enrich_df.iloc[0]

    warning_count = 0
    error_count = 0
    if preflight_df is not None and "status" in preflight_df.columns:
        status = preflight_df["status"].astype(str).str.upper()
        warning_count = int((status == "WARNING").sum())
        error_count = int((status == "ERROR").sum())

    prs_notes = []
    if prs_qc_df is not None and not prs_qc_df.empty:
        prs_map = metric_map(prs_qc_df)
        matched = prs_map.get("matched_model_variants", "n/a")
        model = prs_map.get("model_variants", "n/a")
        rate = prs_map.get("matched_model_rate_pct", "n/a")
        build_match = prs_map.get("build_match", "n/a")
        prs_notes.append(
            f"PRS locus coverage: matched `{matched}` of `{model}` model variants "
            f"(rate `{rate}%`, build_match `{build_match}`)."
        )

    anomaly_notes = []
    if anomaly_df is not None and not anomaly_df.empty and "severity" in anomaly_df.columns:
        sev = anomaly_df["severity"].astype(str).str.lower()
        anomaly_notes.append(
            "AI QC anomaly detector: "
            f"high={int((sev=='high').sum())}, medium={int((sev=='medium').sum())}, "
            f"low={int((sev=='low').sum())}, info={int((sev=='info').sum())}."
        )

    top_priority_note = ""
    if prioritization_df is not None and not prioritization_df.empty:
        row = prioritization_df.iloc[0]
        top_priority_note = (
            f"Top prioritized variant: `{row.get('key', '.')}` "
            f"(score `{row.get('priority_score', 'n/a')}`, tier `{row.get('priority_tier', 'n/a')}`)."
        )

    executive_bullets = [
        f"Run context: mode=`{args.mode or 'n/a'}`, analysis_mode=`{args.analysis_mode or 'n/a'}`, trait_context=`{args.trait_context or 'n/a'}`.",
        f"Retained variants: `{int(total_variants) if total_variants is not None else 'n/a'}`; Ts/Tv: `{format_optional_float(ts_tv, 4)}`; samples: `{int(sample_count) if sample_count is not None else 'n/a'}`.",
        f"Association output size: `{assoc_count}` rows; significant hits: `{significant_count}`.",
        f"Preflight diagnostics: errors=`{error_count}`, warnings=`{warning_count}`.",
    ]
    if top_priority_note:
        executive_bullets.append(top_priority_note)
    executive_bullets.extend(prs_notes)
    executive_bullets.extend(anomaly_notes)

    lines = [
        "# AI Results Explainer",
        "",
        f"Generated (UTC): `{datetime.now(timezone.utc).isoformat()}`",
        "",
        "## Executive Summary",
        "",
    ]
    for bullet in executive_bullets:
        lines.append(f"- {bullet}")

    lines.extend(
        [
            "",
            "## Biological Signal Interpretation",
            "",
            f"- Top genes among association hits: `{', '.join(top_genes) if top_genes else 'n/a'}`.",
        ]
    )

    if top_assoc is not None:
        lines.append(
            "- Lead association row: "
            f"`{top_assoc.get('key', '.')}` with p-value `{top_assoc.get('pvalue', 'n/a')}` "
            f"and effect `{top_assoc.get('effect', 'n/a')}`."
        )

    if enrich_top is not None:
        lines.append(
            "- Top enrichment term: "
            f"`{enrich_top.get('ID', 'n/a')}` "
            f"(qvalue `{enrich_top.get('qvalue', enrich_top.get('p.adjust', 'n/a'))}`)."
        )
    else:
        lines.append("- Enrichment table was empty or unavailable.")

    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- This report is heuristic and intended for prioritization support, not clinical diagnosis.",
            "- Biological conclusions depend on cohort design, reference-build compatibility, and model locus coverage.",
            "- Review primary result tables before making downstream scientific claims.",
            "",
            "## Source Files",
            "",
            f"- variant metrics: `{args.variant_metrics}`",
            f"- association hits: `{args.association}`",
            f"- enrichment: `{args.enrichment}`",
            f"- preflight: `{args.preflight or 'not provided'}`",
            f"- prs qc: `{args.prs_qc or 'not provided'}`",
            f"- anomaly table: `{args.anomaly_table or 'not provided'}`",
            f"- prioritization table: `{args.prioritization_table or 'not provided'}`",
        ]
    )

    out_path = Path(args.out_report)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

