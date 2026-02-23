#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")

records = []
for json_path in snakemake.input["jsons"]:
    sample = Path(str(json_path)).name.replace(".fastp.json", "")
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})

    records.append(
        {
            "sample": sample,
            "reads_before": int(before.get("total_reads", 0)),
            "reads_after": int(after.get("total_reads", 0)),
            "bases_before": int(before.get("total_bases", 0)),
            "bases_after": int(after.get("total_bases", 0)),
            "q20_before": float(before.get("q20_rate", 0.0)),
            "q20_after": float(after.get("q20_rate", 0.0)),
            "q30_before": float(before.get("q30_rate", 0.0)),
            "q30_after": float(after.get("q30_rate", 0.0)),
        }
    )

summary = pd.DataFrame(records)
summary["reads_retention_pct"] = (
    (summary["reads_after"] / summary["reads_before"].replace(0, pd.NA)) * 100
).fillna(0)
summary["bases_retention_pct"] = (
    (summary["bases_after"] / summary["bases_before"].replace(0, pd.NA)) * 100
).fillna(0)

summary.to_csv(snakemake.output["table"], sep="\t", index=False)

plot_df = summary.melt(
    id_vars=["sample"],
    value_vars=["reads_before", "reads_after"],
    var_name="stage",
    value_name="reads",
)
plot_df["stage"] = plot_df["stage"].map(
    {"reads_before": "Before QC", "reads_after": "After QC"}
)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sns.barplot(data=plot_df, x="sample", y="reads", hue="stage", ax=axes[0])
axes[0].set_title("Read Counts Before vs After QC")
axes[0].set_ylabel("Read Count")
axes[0].set_xlabel("Sample")
axes[0].tick_params(axis="x", rotation=30)

sns.barplot(data=summary, x="sample", y="reads_retention_pct", color="#4C78A8", ax=axes[1])
axes[1].set_title("Read Retention (%) After QC")
axes[1].set_ylabel("Retention %")
axes[1].set_xlabel("Sample")
axes[1].set_ylim(0, 105)
axes[1].tick_params(axis="x", rotation=30)

plt.tight_layout()
plt.savefig(snakemake.output["plot"], dpi=200)
