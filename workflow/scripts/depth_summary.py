#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")

summary_records = []
distribution_records = []

for depth_file in snakemake.input["depths"]:
    path = Path(str(depth_file))
    sample = path.name.replace(".depth.tsv", "")

    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "pos", "depth"])
    if df.empty:
        summary_records.append(
            {
                "sample": sample,
                "positions": 0,
                "mean_depth": 0,
                "median_depth": 0,
                "pct_cov_ge_1x": 0,
                "pct_cov_ge_10x": 0,
                "pct_cov_ge_30x": 0,
            }
        )
        continue

    summary_records.append(
        {
            "sample": sample,
            "positions": int(df.shape[0]),
            "mean_depth": float(df["depth"].mean()),
            "median_depth": float(df["depth"].median()),
            "pct_cov_ge_1x": float((df["depth"] >= 1).mean() * 100),
            "pct_cov_ge_10x": float((df["depth"] >= 10).mean() * 100),
            "pct_cov_ge_30x": float((df["depth"] >= 30).mean() * 100),
        }
    )

    sampled = df["depth"].sample(min(50000, len(df)), random_state=7)
    distribution_records.extend(
        [{"sample": sample, "depth": float(x)} for x in sampled.tolist()]
    )

summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(snakemake.output["table"], sep="\t", index=False)

plot_df = pd.DataFrame(distribution_records)
plt.figure(figsize=(12, 6))
if not plot_df.empty:
    sns.violinplot(data=plot_df, x="sample", y="depth", cut=0, inner="quartile")
    plt.yscale("symlog", linthresh=10)
    plt.ylabel("Read Depth (symlog scale)")
    plt.xlabel("Sample")
    plt.title("Depth Distribution by Sample")
else:
    plt.text(0.5, 0.5, "No depth records found", ha="center", va="center")
    plt.axis("off")

plt.tight_layout()
plt.savefig(snakemake.output["plot"], dpi=200)
