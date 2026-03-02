#!/usr/bin/env python3
"""Heuristic AI-style QC anomaly detector for pipeline outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}


@dataclass
class Finding:
    check: str
    severity: str
    value: str
    expectation: str
    message: str


def safe_float(value: str) -> Optional[float]:
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def load_metric_map(path: str) -> Dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "metric" not in table.columns or "value" not in table.columns:
        raise ValueError(f"{path} must contain 'metric' and 'value' columns")
    return {
        str(row["metric"]).strip(): str(row["value"]).strip()
        for _, row in table.iterrows()
        if str(row["metric"]).strip()
    }


def maybe_read_table(path: str) -> Optional[pd.DataFrame]:
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    return pd.read_csv(p, sep="\t", dtype=str).fillna("")


def add_findings_for_variant_metrics(metrics: Dict[str, str], findings: List[Finding]) -> None:
    total_variants = safe_float(metrics.get("total_variants", ""))
    if total_variants is None:
        findings.append(
            Finding(
                check="variant_count_parse",
                severity="medium",
                value=metrics.get("total_variants", ""),
                expectation="numeric total_variants",
                message="Could not parse total_variants metric.",
            )
        )
    elif total_variants <= 0:
        findings.append(
            Finding(
                check="variant_count_zero",
                severity="high",
                value=str(int(total_variants)),
                expectation="> 0 variants",
                message="No variants were retained; check filtering and inputs.",
            )
        )
    elif total_variants < 100:
        findings.append(
            Finding(
                check="variant_count_low",
                severity="low",
                value=str(int(total_variants)),
                expectation="context-dependent; often > 100",
                message="Low retained variant count may indicate strict filters or sparse test data.",
            )
        )

    ts_tv = safe_float(metrics.get("ts_tv_ratio", ""))
    if ts_tv is None:
        findings.append(
            Finding(
                check="ts_tv_missing",
                severity="medium",
                value=metrics.get("ts_tv_ratio", ""),
                expectation="numeric Ts/Tv ratio",
                message="Ts/Tv ratio is missing or unparsable.",
            )
        )
    elif ts_tv < 1.2 or ts_tv > 3.5:
        findings.append(
            Finding(
                check="ts_tv_outlier",
                severity="high",
                value=f"{ts_tv:.4f}",
                expectation="approx 1.2 to 3.5",
                message="Ts/Tv is outside broad expected ranges for typical human variant sets.",
            )
        )
    elif ts_tv < 1.5 or ts_tv > 3.0:
        findings.append(
            Finding(
                check="ts_tv_borderline",
                severity="low",
                value=f"{ts_tv:.4f}",
                expectation="approx 1.5 to 3.0",
                message="Ts/Tv is borderline; interpret in context of assay, coverage, and filters.",
            )
        )


def add_findings_for_preflight(preflight: Optional[pd.DataFrame], findings: List[Finding]) -> None:
    if preflight is None or "status" not in preflight.columns:
        return
    status = preflight["status"].astype(str).str.upper()
    warning_count = int((status == "WARNING").sum())
    error_count = int((status == "ERROR").sum())
    if error_count > 0:
        findings.append(
            Finding(
                check="preflight_errors",
                severity="high",
                value=str(error_count),
                expectation="0 preflight errors",
                message="Preflight reported errors; run may be invalid for downstream interpretation.",
            )
        )
    if warning_count > 0:
        findings.append(
            Finding(
                check="preflight_warnings",
                severity="low",
                value=str(warning_count),
                expectation="0 preflight warnings preferred",
                message="Preflight warnings present; review preflight_checks.tsv for context.",
            )
        )


def add_findings_for_qc_before_after(qc_table: Optional[pd.DataFrame], findings: List[Finding]) -> None:
    if qc_table is None or qc_table.empty:
        return
    row = qc_table.iloc[0]
    reads_retention = safe_float(row.get("reads_retention_pct", ""))
    if reads_retention is not None:
        if reads_retention < 85:
            findings.append(
                Finding(
                    check="reads_retention_low",
                    severity="high",
                    value=f"{reads_retention:.3f}%",
                    expectation=">= 85%",
                    message="Read retention after trimming is low.",
                )
            )
        elif reads_retention < 92:
            findings.append(
                Finding(
                    check="reads_retention_moderate",
                    severity="medium",
                    value=f"{reads_retention:.3f}%",
                    expectation=">= 92% preferred",
                    message="Read retention is moderate; review trimming thresholds.",
                )
            )

    q30_before = safe_float(row.get("q30_before", ""))
    q30_after = safe_float(row.get("q30_after", ""))
    if q30_before is not None and q30_after is not None and q30_after < q30_before:
        findings.append(
            Finding(
                check="q30_drop_after_trim",
                severity="medium",
                value=f"{q30_before:.4f} -> {q30_after:.4f}",
                expectation="q30_after >= q30_before",
                message="Q30 fraction dropped after trimming; inspect adapter/quality settings.",
            )
        )


def add_findings_for_depth(depth_table: Optional[pd.DataFrame], findings: List[Finding]) -> None:
    if depth_table is None or depth_table.empty:
        return
    mean_depth = safe_float(depth_table.iloc[0].get("mean_depth", ""))
    if mean_depth is None:
        return
    if mean_depth < 0.02:
        findings.append(
            Finding(
                check="mean_depth_very_low",
                severity="high",
                value=f"{mean_depth:.4f}x",
                expectation="context-dependent; often >= 10x for production WGS",
                message="Coverage is extremely low; sensitivity will be limited.",
            )
        )
    elif mean_depth < 0.2:
        findings.append(
            Finding(
                check="mean_depth_low",
                severity="medium",
                value=f"{mean_depth:.4f}x",
                expectation="context-dependent; often >= 10x for production WGS",
                message="Coverage is low; likely suitable for technical/pathway demo, not deep variant sensitivity.",
            )
        )


def write_outputs(findings: List[Finding], out_table: str, out_report: str) -> None:
    out_table_path = Path(out_table)
    out_report_path = Path(out_report)
    out_table_path.parent.mkdir(parents=True, exist_ok=True)
    out_report_path.parent.mkdir(parents=True, exist_ok=True)

    if not findings:
        findings = [
            Finding(
                check="no_major_qc_anomalies",
                severity="info",
                value="none",
                expectation="n/a",
                message="No major QC anomalies were detected by heuristic checks.",
            )
        ]

    result_df = pd.DataFrame([f.__dict__ for f in findings])
    result_df["severity_rank"] = result_df["severity"].map(SEVERITY_ORDER).fillna(99)
    result_df = result_df.sort_values(["severity_rank", "check"]).drop(columns=["severity_rank"])
    result_df.to_csv(out_table_path, sep="\t", index=False)

    severity_counts = result_df["severity"].value_counts().to_dict()
    top_rows = result_df.head(15)
    lines = [
        "# AI QC Anomaly Report",
        "",
        "## Summary",
        "",
        f"- high: {severity_counts.get('high', 0)}",
        f"- medium: {severity_counts.get('medium', 0)}",
        f"- low: {severity_counts.get('low', 0)}",
        f"- info: {severity_counts.get('info', 0)}",
        "",
        "## Findings",
        "",
    ]
    for _, row in top_rows.iterrows():
        lines.append(
            f"- [{row['severity']}] `{row['check']}`: {row['message']} "
            f"(value: `{row['value']}`; expected: `{row['expectation']}`)"
        )

    out_report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Detect QC anomalies from pipeline output tables.")
    parser.add_argument("--variant-metrics", required=True, help="Path to variant_qc_metrics.tsv")
    parser.add_argument("--preflight", default="", help="Optional path to preflight_checks.tsv")
    parser.add_argument("--qc-before-after", default="", help="Optional path to qc_before_after.tsv")
    parser.add_argument("--depth-summary", default="", help="Optional path to depth_summary.tsv")
    parser.add_argument("--out-table", required=True, help="Output anomaly TSV path")
    parser.add_argument("--out-report", required=True, help="Output markdown report path")
    args = parser.parse_args()

    findings: List[Finding] = []
    metrics = load_metric_map(args.variant_metrics)
    add_findings_for_variant_metrics(metrics, findings)
    add_findings_for_preflight(maybe_read_table(args.preflight), findings)
    add_findings_for_qc_before_after(maybe_read_table(args.qc_before_after), findings)
    add_findings_for_depth(maybe_read_table(args.depth_summary), findings)
    write_outputs(findings, args.out_table, args.out_report)


if __name__ == "__main__":
    main()

