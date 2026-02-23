#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_theme(style="whitegrid")


def parse_args():
    parser = argparse.ArgumentParser(description="Create PRS table and distribution plot from PLINK score output")
    parser.add_argument("--sscore", required=True)
    parser.add_argument("--out-table", required=True)
    parser.add_argument("--out-plot", required=True)
    return parser.parse_args()


def find_score_col(columns):
    for suffix in ["_SUM", "_AVG"]:
        matches = [c for c in columns if c.endswith(suffix)]
        if matches:
            return matches[0]
    fallback = [c for c in columns if "SCORE" in c.upper()]
    if fallback:
        return fallback[0]
    raise ValueError("Could not find a PRS score column in PLINK score output")


def main():
    args = parse_args()
    df = pd.read_csv(args.sscore, sep="\t")
    score_col = find_score_col(df.columns)

    sample_col = "IID" if "IID" in df.columns else df.columns[0]
    keep_cols = [sample_col, score_col]
    for extra in ["N_MATCHED_VARIANTS", "MODEL_VARIANTS", "MATCH_RATE"]:
        if extra in df.columns:
            keep_cols.append(extra)
    out = df[keep_cols].copy()
    rename_map = {sample_col: "sample", score_col: "prs"}
    out = out.rename(columns=rename_map)

    if "N_MATCHED_VARIANTS" in out.columns:
        out["N_MATCHED_VARIANTS"] = pd.to_numeric(out["N_MATCHED_VARIANTS"], errors="coerce")
    if "MODEL_VARIANTS" in out.columns:
        out["MODEL_VARIANTS"] = pd.to_numeric(out["MODEL_VARIANTS"], errors="coerce")

    prs_std = out["prs"].std(ddof=0)
    if pd.notna(prs_std) and prs_std > 0:
        out["prs_percentile"] = out["prs"].rank(method="average", pct=True) * 100.0
        out["prs_zscore"] = (out["prs"] - out["prs"].mean()) / prs_std
    else:
        out["prs_percentile"] = 50.0
        out["prs_zscore"] = 0.0

    if "N_MATCHED_VARIANTS" in out.columns:
        out["alzheimers_prs_note"] = out["N_MATCHED_VARIANTS"].apply(
            lambda x: (
                "No matched Alzheimer PRS loci; PRS shown as 0 (no callable AD genetic risk estimate)."
                if pd.notna(x) and x <= 0
                else "Alzheimer PRS computed from matched loci; interpret relative to cohort distribution."
            )
        )
    else:
        out["alzheimers_prs_note"] = "Alzheimer PRS computed; matched-variant coverage column unavailable."

    out.to_csv(args.out_table, sep="\t", index=False)

    plt.figure(figsize=(8, 5))
    sns.histplot(out["prs"], bins=30, kde=True, color="#4C78A8")
    plt.title("Polygenic Risk Score Distribution")
    plt.xlabel("PRS")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=200)


if __name__ == "__main__":
    main()
