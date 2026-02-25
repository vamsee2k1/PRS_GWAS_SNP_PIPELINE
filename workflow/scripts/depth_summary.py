#!/usr/bin/env python3

import csv
import random
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")


def _cfg_depth_value(key, default):
    return snakemake.config.get("depth", {}).get(key, default)


def _safe_int(value, default):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _safe_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


PLOT_MAX_DEPTH = max(1, _safe_int(_cfg_depth_value("summary_plot_max_depth", 200), 200))
FALLBACK_SAMPLE_SIZE = max(0, _safe_int(_cfg_depth_value("summary_plot_sample_size", 50000), 50000))


def sample_from_depth_path(path: Path) -> str:
    return path.name.replace(".depth.tsv", "")


def parse_mosdepth_summary(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return None

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("chrom") != "total":
                continue
            length = _safe_int(row.get("length"), 0)
            depth_sum = _safe_int(row.get("bases"), 0)
            mean_depth = (depth_sum / length) if length > 0 else _safe_float(row.get("mean"), 0.0)
            return {
                "positions": length,
                "depth_sum": depth_sum,
                "mean_depth": float(mean_depth),
                "max_depth": _safe_int(row.get("max"), 0),
            }
    return None


def parse_mosdepth_global_dist(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return None

    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] != "total":
                continue
            depth = _safe_int(parts[1], None)
            frac = _safe_float(parts[2], None)
            if depth is None or frac is None:
                continue
            rows.append((depth, min(1.0, max(0.0, frac))))

    if not rows:
        return None

    rows.sort(key=lambda x: x[0])
    # Make the cumulative proportion monotonic non-increasing across depth.
    fixed = []
    prev = 1.0
    for depth, frac in rows:
        frac = min(frac, prev)
        fixed.append((depth, frac))
        prev = frac
    return fixed


def frac_ge_at_depth(cumulative_rows, threshold: int) -> float:
    if not cumulative_rows:
        return 0.0
    frac_map = {depth: frac for depth, frac in cumulative_rows}
    return float(frac_map.get(threshold, 0.0))


def median_from_cumulative(cumulative_rows) -> float:
    if not cumulative_rows:
        return 0.0

    median_depth = 0
    for depth, frac_ge in cumulative_rows:
        if frac_ge >= 0.5:
            median_depth = depth
        else:
            break
    return float(median_depth)


def reservoir_update(reservoir, value, seen_count, max_size, rng):
    if max_size <= 0:
        return
    if len(reservoir) < max_size:
        reservoir.append(value)
        return
    j = rng.randrange(seen_count)
    if j < max_size:
        reservoir[j] = value


def median_from_histogram(depth_counts: Counter, total_n: int) -> float:
    if total_n == 0:
        return 0.0

    mid1 = (total_n - 1) // 2
    mid2 = total_n // 2
    running = 0
    m1 = None
    m2 = None

    for depth in sorted(depth_counts):
        running += depth_counts[depth]
        if m1 is None and running > mid1:
            m1 = depth
        if running > mid2:
            m2 = depth
            break

    if m1 is None or m2 is None:
        return 0.0
    return float((m1 + m2) / 2.0)


def summarize_from_depth_tsv(depth_file: Path, sample: str):
    total_n = 0
    depth_sum = 0
    max_depth = 0
    ge1 = 0
    ge10 = 0
    ge30 = 0
    ge100 = 0
    depth_counts = Counter()
    sampled_depths = []
    rng = random.Random(7)

    if not depth_file.exists() or depth_file.stat().st_size == 0:
        return {
            "summary": {
                "sample": sample,
                "positions": 0,
                "positions_basis": "depth_rows",
                "mean_depth": 0.0,
                "median_depth": 0.0,
                "pct_cov_ge_1x": 0.0,
                "pct_cov_ge_10x": 0.0,
                "pct_cov_ge_30x": 0.0,
                "pct_cov_ge_100x": 0.0,
                "max_depth": 0,
                "depth_sum": 0,
                "depth_summary_source": "depth_tsv_fallback",
            },
            "sampled_depths": [],
            "plot_source": "fallback",
        }

    with depth_file.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            depth = _safe_int(parts[2], None)
            if depth is None:
                continue

            total_n += 1
            depth_sum += depth
            max_depth = max(max_depth, depth)
            ge1 += int(depth >= 1)
            ge10 += int(depth >= 10)
            ge30 += int(depth >= 30)
            ge100 += int(depth >= 100)
            depth_counts[depth] += 1
            reservoir_update(sampled_depths, float(depth), total_n, FALLBACK_SAMPLE_SIZE, rng)

    if total_n == 0:
        mean_depth = 0.0
        median_depth = 0.0
        pct1 = pct10 = pct30 = pct100 = 0.0
    else:
        mean_depth = depth_sum / total_n
        median_depth = median_from_histogram(depth_counts, total_n)
        pct1 = (ge1 / total_n) * 100.0
        pct10 = (ge10 / total_n) * 100.0
        pct30 = (ge30 / total_n) * 100.0
        pct100 = (ge100 / total_n) * 100.0

    return {
        "summary": {
            "sample": sample,
            "positions": int(total_n),
            "positions_basis": "depth_rows",
            "mean_depth": float(mean_depth),
            "median_depth": float(median_depth),
            "pct_cov_ge_1x": float(pct1),
            "pct_cov_ge_10x": float(pct10),
            "pct_cov_ge_30x": float(pct30),
            "pct_cov_ge_100x": float(pct100),
            "max_depth": int(max_depth),
            "depth_sum": int(depth_sum),
            "depth_summary_source": "depth_tsv_fallback",
        },
        "sampled_depths": sampled_depths,
        "plot_source": "fallback",
    }


def summarize_from_mosdepth(sample: str, summary_path: Path, dist_path: Path):
    summary = parse_mosdepth_summary(summary_path)
    cumulative = parse_mosdepth_global_dist(dist_path)
    if summary is None or cumulative is None:
        return None

    median_depth = median_from_cumulative(cumulative)
    pct1 = frac_ge_at_depth(cumulative, 1) * 100.0
    pct10 = frac_ge_at_depth(cumulative, 10) * 100.0
    pct30 = frac_ge_at_depth(cumulative, 30) * 100.0
    pct100 = frac_ge_at_depth(cumulative, 100) * 100.0

    plot_rows = [
        {
            "sample": sample,
            "depth": depth,
            "pct_cov_ge": frac * 100.0,
        }
        for depth, frac in cumulative
        if depth <= PLOT_MAX_DEPTH
    ]

    return {
        "summary": {
            "sample": sample,
            "positions": int(summary["positions"]),
            "positions_basis": "reference_length",
            "mean_depth": float(summary["mean_depth"]),
            "median_depth": float(median_depth),
            "pct_cov_ge_1x": float(pct1),
            "pct_cov_ge_10x": float(pct10),
            "pct_cov_ge_30x": float(pct30),
            "pct_cov_ge_100x": float(pct100),
            "max_depth": int(summary["max_depth"]),
            "depth_sum": int(summary["depth_sum"]),
            "depth_summary_source": "mosdepth",
        },
        "cumulative_plot_rows": plot_rows,
        "plot_source": "mosdepth",
    }


summary_records = []
mosdepth_plot_records = []
fallback_distribution_records = []

mosdepth_summary_map = {
    Path(str(p)).name.replace(".mosdepth.summary.txt", ""): Path(str(p))
    for p in snakemake.input["mosdepth_summaries"]
}
mosdepth_dist_map = {
    Path(str(p)).name.replace(".mosdepth.global.dist.txt", ""): Path(str(p))
    for p in snakemake.input["mosdepth_dists"]
}

for depth_file in snakemake.input["depths"]:
    path = Path(str(depth_file))
    sample = sample_from_depth_path(path)

    mos_summary_path = mosdepth_summary_map.get(sample)
    mos_dist_path = mosdepth_dist_map.get(sample)
    mos_res = None
    if mos_summary_path is not None and mos_dist_path is not None:
        mos_res = summarize_from_mosdepth(sample, mos_summary_path, mos_dist_path)

    if mos_res is not None:
        summary_records.append(mos_res["summary"])
        mosdepth_plot_records.extend(mos_res["cumulative_plot_rows"])
        continue

    fallback_res = summarize_from_depth_tsv(path, sample)
    summary_records.append(fallback_res["summary"])
    for depth in fallback_res["sampled_depths"]:
        fallback_distribution_records.append(
            {"sample": sample, "depth": min(float(depth), float(PLOT_MAX_DEPTH))}
        )


summary_df = pd.DataFrame(summary_records)
column_order = [
    "sample",
    "positions",
    "positions_basis",
    "mean_depth",
    "median_depth",
    "pct_cov_ge_1x",
    "pct_cov_ge_10x",
    "pct_cov_ge_30x",
    "pct_cov_ge_100x",
    "max_depth",
    "depth_sum",
    "depth_summary_source",
]
summary_df = summary_df.reindex(columns=column_order)
summary_df.to_csv(snakemake.output["table"], sep="\t", index=False)

plt.figure(figsize=(12, 6))

if mosdepth_plot_records:
    plot_df = pd.DataFrame(mosdepth_plot_records)
    for sample, sample_df in plot_df.groupby("sample", sort=True):
        sample_df = sample_df.sort_values("depth")
        plt.plot(sample_df["depth"], sample_df["pct_cov_ge"], label=sample, linewidth=2)
    plt.xlabel(f"Depth (X), clipped at {PLOT_MAX_DEPTH} for display")
    plt.ylabel("Percent of bases with coverage >= X")
    plt.title("Coverage Distribution by Sample (mosdepth global cumulative)")
    plt.ylim(0, 100)
    if plot_df["sample"].nunique() <= 12:
        plt.legend(title="Sample", loc="best")
elif fallback_distribution_records:
    plot_df = pd.DataFrame(fallback_distribution_records)
    sns.violinplot(data=plot_df, x="sample", y="depth", cut=0, inner="quartile")
    plt.yscale("symlog", linthresh=10)
    plt.ylabel(f"Read Depth (symlog, clipped at {PLOT_MAX_DEPTH} for display)")
    plt.xlabel("Sample")
    plt.title("Depth Distribution by Sample (fallback depth.tsv sampling)")
else:
    plt.text(0.5, 0.5, "No depth records found", ha="center", va="center")
    plt.axis("off")

plt.tight_layout()
plt.savefig(snakemake.output["plot"], dpi=200)
