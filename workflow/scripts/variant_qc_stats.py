#!/usr/bin/env python3
import argparse
import gzip
import math
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_theme(style="whitegrid")

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_args():
    p = argparse.ArgumentParser(description="Compute biologically meaningful VCF QC and prioritization summaries.")
    p.add_argument("--vcf", required=True)
    p.add_argument("--out-metrics", required=True)
    p.add_argument("--out-missingness", required=True)
    p.add_argument("--out-type-plot", required=True)
    p.add_argument("--out-qual-plot", required=True)
    p.add_argument("--out-dp-plot", required=True)
    p.add_argument("--out-gq-plot", required=True)
    p.add_argument("--out-ab-plot", required=True)
    p.add_argument("--out-missingness-plot", required=True)
    p.add_argument(
        "--max-genotype-samples",
        type=int,
        default=1000,
        help="Skip genotype-level loops if sample count exceeds this threshold. Set <=0 to disable skipping.",
    )
    return p.parse_args()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def blank_plot(path, label):
    plt.figure(figsize=(8, 5))
    plt.text(0.5, 0.5, label, ha="center", va="center", fontsize=11)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def parse_float(token):
    tok = str(token).strip()
    if tok in {"", ".", "NA", "nan", "NAN"}:
        return None
    try:
        val = float(tok)
    except ValueError:
        return None
    if math.isnan(val) or math.isinf(val):
        return None
    return val


def gt_missing(gt):
    return gt in {"", ".", "./.", ".|."}


def is_het(gt):
    g = gt.replace("|", "/").split("/")
    return len(g) == 2 and g[0] != "." and g[1] != "." and g[0] != g[1]


def parse_ad_alt_fraction(ad):
    vals = [parse_float(x) for x in str(ad).split(",")]
    if len(vals) < 2 or vals[0] is None or vals[1] is None:
        return None
    denom = vals[0] + vals[1]
    if denom <= 0:
        return None
    return vals[1] / denom


def main():
    args = parse_args()

    total = 0
    ti = 0
    tv = 0
    type_counts = Counter()
    qual_vals = []
    dp_vals = []
    gq_vals = []
    ab_vals = []

    samples = []
    sample_total = defaultdict(int)
    sample_missing = defaultdict(int)
    process_genotypes = True

    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                if args.max_genotype_samples > 0 and len(samples) > args.max_genotype_samples:
                    process_genotypes = False
                continue
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            total += 1
            ref = parts[3].upper()
            alt_field = parts[4]
            qual = parse_float(parts[5])
            if qual is not None:
                qual_vals.append(qual)

            alts = [a.strip().upper() for a in alt_field.split(",") if a.strip()]
            if len(alts) != 1:
                type_counts["multiallelic"] += 1
            else:
                alt = alts[0]
                if len(ref) == 1 and len(alt) == 1:
                    type_counts["SNP"] += 1
                    if (ref, alt) in TRANSITIONS:
                        ti += 1
                    else:
                        tv += 1
                elif len(ref) != len(alt):
                    type_counts["INDEL"] += 1
                elif len(ref) == len(alt) and len(ref) > 1:
                    type_counts["MNP"] += 1
                else:
                    type_counts["OTHER"] += 1

            if len(parts) < 10 or not samples or not process_genotypes:
                continue

            fmt_keys = parts[8].split(":")
            for sample, sample_field in zip(samples, parts[9:]):
                sample_total[sample] += 1
                vals = sample_field.split(":")
                fdict = {k: vals[i] if i < len(vals) else "" for i, k in enumerate(fmt_keys)}

                gt = fdict.get("GT", "")
                if gt_missing(gt):
                    sample_missing[sample] += 1

                dp = parse_float(fdict.get("DP", ""))
                if dp is not None:
                    dp_vals.append(dp)

                gq = parse_float(fdict.get("GQ", ""))
                if gq is not None:
                    gq_vals.append(gq)

                if is_het(gt):
                    ab = parse_ad_alt_fraction(fdict.get("AD", ""))
                    if ab is not None:
                        ab_vals.append(ab)

    ts_tv_ratio = float(ti) / float(tv) if tv > 0 else np.nan
    metrics = [
        ("total_variants", total),
        ("snp_count", int(type_counts.get("SNP", 0))),
        ("indel_count", int(type_counts.get("INDEL", 0))),
        ("mnp_count", int(type_counts.get("MNP", 0))),
        ("multiallelic_count", int(type_counts.get("multiallelic", 0))),
        ("other_count", int(type_counts.get("OTHER", 0))),
        ("transition_count", ti),
        ("transversion_count", tv),
        ("ts_tv_ratio", "" if np.isnan(ts_tv_ratio) else f"{ts_tv_ratio:.6g}"),
        ("sample_count", len(samples)),
        (
            "genotype_metrics_mode",
            "full"
            if process_genotypes
            else f"skipped_gt_metrics_n>{args.max_genotype_samples}",
        ),
        ("mean_qual", "" if not qual_vals else f"{np.mean(qual_vals):.6g}"),
        ("median_qual", "" if not qual_vals else f"{np.median(qual_vals):.6g}"),
        ("mean_dp", "" if not dp_vals else f"{np.mean(dp_vals):.6g}"),
        ("median_dp", "" if not dp_vals else f"{np.median(dp_vals):.6g}"),
        ("mean_gq", "" if not gq_vals else f"{np.mean(gq_vals):.6g}"),
        ("median_gq", "" if not gq_vals else f"{np.median(gq_vals):.6g}"),
        ("mean_het_allele_balance", "" if not ab_vals else f"{np.mean(ab_vals):.6g}"),
        ("median_het_allele_balance", "" if not ab_vals else f"{np.median(ab_vals):.6g}"),
    ]
    pd.DataFrame(metrics, columns=["metric", "value"]).to_csv(args.out_metrics, sep="\t", index=False)

    if process_genotypes:
        miss_rows = []
        for s in samples:
            denom = sample_total.get(s, 0)
            miss = sample_missing.get(s, 0)
            rate = (miss / denom * 100.0) if denom else np.nan
            miss_rows.append(
                {
                    "sample": s,
                    "missing_variants": miss,
                    "total_variants": denom,
                    "missing_rate_pct": rate,
                }
            )
        miss_df = pd.DataFrame(miss_rows)
    else:
        miss_df = pd.DataFrame(columns=["sample", "missing_variants", "total_variants", "missing_rate_pct"])
    if miss_df.empty:
        miss_df = pd.DataFrame(columns=["sample", "missing_variants", "total_variants", "missing_rate_pct"])
    miss_df.to_csv(args.out_missingness, sep="\t", index=False)

    skipped_label = (
        f"Skipped genotype-level metrics for large cohort (n={len(samples)} > {args.max_genotype_samples})"
        if (not process_genotypes and len(samples) > 0)
        else None
    )

    if type_counts:
        type_df = pd.DataFrame(
            [{"variant_type": k, "count": v} for k, v in sorted(type_counts.items(), key=lambda x: (-x[1], x[0]))]
        )
        plt.figure(figsize=(8, 5))
        sns.barplot(data=type_df, x="variant_type", y="count", color="#3E7CB1")
        plt.title("Variant Type Counts")
        plt.xlabel("Variant Type")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(args.out_type_plot, dpi=200)
        plt.close()
    else:
        blank_plot(args.out_type_plot, "No variant type data")

    if qual_vals:
        plt.figure(figsize=(8, 5))
        sns.histplot(qual_vals, bins=40, color="#2A9D8F")
        plt.title("QUAL Distribution")
        plt.xlabel("QUAL")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(args.out_qual_plot, dpi=200)
        plt.close()
    else:
        blank_plot(args.out_qual_plot, "QUAL not available")

    if dp_vals:
        plt.figure(figsize=(8, 5))
        sns.histplot(dp_vals, bins=50, color="#457B9D")
        plt.title("Depth (DP) Distribution")
        plt.xlabel("DP")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(args.out_dp_plot, dpi=200)
        plt.close()
    else:
        blank_plot(args.out_dp_plot, skipped_label or "DP not available")

    if gq_vals:
        plt.figure(figsize=(8, 5))
        sns.histplot(gq_vals, bins=50, color="#E76F51")
        plt.title("Genotype Quality (GQ) Distribution")
        plt.xlabel("GQ")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(args.out_gq_plot, dpi=200)
        plt.close()
    else:
        blank_plot(args.out_gq_plot, skipped_label or "GQ not available")

    if ab_vals:
        plt.figure(figsize=(8, 5))
        sns.histplot(ab_vals, bins=40, color="#8AB17D")
        plt.title("Allele Balance Distribution (Heterozygous GT)")
        plt.xlabel("Alt Allele Fraction (AD)")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(args.out_ab_plot, dpi=200)
        plt.close()
    else:
        blank_plot(args.out_ab_plot, skipped_label or "AD-based allele balance not available")

    if not miss_df.empty and miss_df["missing_rate_pct"].notna().any():
        top = miss_df.sort_values("missing_rate_pct", ascending=False).head(60)
        plt.figure(figsize=(max(8, 0.3 * len(top)), 5))
        sns.barplot(data=top, x="sample", y="missing_rate_pct", color="#6D597A")
        plt.title("Per-Sample Missingness Rate")
        plt.xlabel("Sample")
        plt.ylabel("Missing Rate (%)")
        plt.xticks(rotation=65, ha="right")
        plt.tight_layout()
        plt.savefig(args.out_missingness_plot, dpi=200)
        plt.close()
    else:
        blank_plot(
            args.out_missingness_plot,
            skipped_label or "Missingness unavailable (no sample genotypes)",
        )


if __name__ == "__main__":
    main()
