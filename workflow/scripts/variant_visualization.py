#!/usr/bin/env python3
import argparse
import gzip
import math
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_theme(style="whitegrid")
NC_CHROM_RE = re.compile(r"^NC_0*([0-9]+)\\.")


def parse_args():
    parser = argparse.ArgumentParser(description="Create variant annotation and association visualizations")
    parser.add_argument("--annotated", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--analysis-mode", default="vcf_interpretation", choices=["vcf_interpretation", "gwas_summary"])
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--effect-threshold", type=float, default=0.0)
    parser.add_argument("--heatmap-top-n", type=int, default=50)
    parser.add_argument("--out-assoc", required=True)
    parser.add_argument("--out-feature-plot", required=True)
    parser.add_argument("--out-chrom-plot", required=True)
    parser.add_argument("--out-volcano", required=True)
    parser.add_argument("--out-manhattan", required=True)
    parser.add_argument("--out-qq", required=True)
    parser.add_argument("--out-heatmap", required=True)
    return parser.parse_args()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def normalize_chrom(chrom):
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]

    m = NC_CHROM_RE.match(c)
    if m:
        num = int(m.group(1))
        if num == 23:
            return "X"
        if num == 24:
            return "Y"
        if num == 12920:
            return "MT"
        return str(num)

    if c in {"M", "MT", "chrM", "chrMT"}:
        return "MT"

    return c


def blank_plot(path, label):
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, label, ha="center", va="center", fontsize=11)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def chrom_sort_key(chrom):
    c = str(chrom)
    if c.isdigit():
        return (0, int(c))
    if c == "X":
        return (1, 23)
    if c == "Y":
        return (1, 24)
    if c in {"MT", "M"}:
        return (1, 25)
    return (2, c)


def gt_alt_dosage(gt):
    if gt in {".", "./.", ".|."}:
        return np.nan
    alleles = gt.replace("|", "/").split("/")
    dose = 0.0
    for a in alleles:
        if a == ".":
            return np.nan
        if a == "0":
            continue
        if a == "1":
            dose += 1.0
        else:
            return np.nan
    return dose


def parse_sample_dose(sample_field, fmt_keys):
    values = sample_field.split(":")
    fdict = {k: values[i] if i < len(values) else "" for i, k in enumerate(fmt_keys)}

    ds = fdict.get("DS", "")
    if ds not in {"", "."}:
        token = ds.split(",")[0]
        try:
            return float(token)
        except ValueError:
            pass

    gt = fdict.get("GT", "")
    if gt:
        return gt_alt_dosage(gt)

    return np.nan


def clinical_sig_to_proxy(raw_sig):
    sig = str(raw_sig or "").strip().lower()
    if not sig or sig in {".", "nan"}:
        return np.nan, np.nan

    if "conflicting" in sig or "uncertain" in sig:
        return 0.2, 0.0
    if "pathogenic" in sig and "likely" not in sig:
        return 1e-4, 1.0
    if "likely_pathogenic" in sig:
        return 1e-3, 0.7
    if "benign" in sig and "likely" not in sig:
        return 1e-4, -1.0
    if "likely_benign" in sig:
        return 1e-3, -0.7
    if "risk_factor" in sig:
        return 5e-3, 0.5
    if "protective" in sig:
        return 5e-3, -0.5

    return 0.05, 0.2


def plot_feature_distribution(df, out_path):
    if "functional_class" not in df.columns or df.empty:
        blank_plot(out_path, "No functional classes available")
        return

    counts = df["functional_class"].fillna("unknown").value_counts().reset_index()
    counts.columns = ["functional_class", "count"]
    plt.figure(figsize=(9, 5))
    sns.barplot(data=counts, x="functional_class", y="count", color="#3E7CB1")
    plt.title("Variant Functional Class Distribution")
    plt.xlabel("Functional Class")
    plt.ylabel("Variant Count")
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_chrom_distribution(df, out_path):
    if "chrom" not in df.columns or df.empty:
        blank_plot(out_path, "No chromosome distribution available")
        return

    counts = df["chrom"].fillna("unknown").value_counts().reset_index()
    counts.columns = ["chrom", "count"]
    counts = counts.sort_values("chrom", key=lambda s: s.map(chrom_sort_key))
    plt.figure(figsize=(10, 5))
    sns.barplot(data=counts, x="chrom", y="count", color="#2A9D8F")
    plt.title("Variant Counts by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Variant Count")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def build_assoc_table(df, p_threshold, effect_threshold):
    assoc = df.copy()
    assoc["pvalue"] = pd.to_numeric(assoc.get("pvalue", pd.Series(dtype=float)), errors="coerce")
    assoc["effect"] = pd.to_numeric(assoc.get("effect", pd.Series(dtype=float)), errors="coerce")
    assoc["qual"] = pd.to_numeric(assoc.get("qual", pd.Series(dtype=float)), errors="coerce")
    assoc["af"] = pd.to_numeric(assoc.get("af", pd.Series(dtype=float)), errors="coerce")

    reported = assoc[
        np.isfinite(assoc["pvalue"])
        & np.isfinite(assoc["effect"])
        & (assoc["pvalue"] > 0)
        & (assoc["pvalue"] <= 1)
    ].copy()

    if not reported.empty:
        reported["association_mode"] = "reported"
        assoc = reported
    else:
        proxy = assoc.copy()
        qual_proxy = np.where(
            np.isfinite(proxy["qual"]) & (proxy["qual"] > 0),
            np.power(10.0, -proxy["qual"].clip(upper=1000) / 10.0),
            np.nan,
        )
        proxy["pvalue"] = qual_proxy

        proxy_effect_from_af = np.where(
            np.isfinite(proxy["af"]) & (proxy["af"] >= 0) & (proxy["af"] <= 1),
            proxy["af"] - 0.5,
            np.nan,
        )
        proxy["effect"] = np.where(np.isfinite(proxy["effect"]), proxy["effect"], proxy_effect_from_af)
        af_proxy_p = np.where(
            np.isfinite(proxy_effect_from_af),
            np.clip(1.0 - (np.abs(proxy_effect_from_af) * 2.0), 1e-6, 1.0),
            np.nan,
        )
        proxy["pvalue"] = np.where(np.isfinite(proxy["pvalue"]), proxy["pvalue"], af_proxy_p)

        proxy = proxy[
            np.isfinite(proxy["pvalue"])
            & np.isfinite(proxy["effect"])
            & (proxy["pvalue"] > 0)
            & (proxy["pvalue"] <= 1)
        ].copy()
        if proxy.empty and "clinical_significance" in assoc.columns:
            clin = assoc.copy()
            pe = clin["clinical_significance"].apply(clinical_sig_to_proxy)
            clin["pvalue"] = pe.apply(lambda x: x[0])
            clin["effect"] = pe.apply(lambda x: x[1])
            clin["pvalue"] = pd.to_numeric(clin["pvalue"], errors="coerce")
            clin["effect"] = pd.to_numeric(clin["effect"], errors="coerce")
            clin = clin[
                np.isfinite(clin["pvalue"])
                & np.isfinite(clin["effect"])
                & (clin["pvalue"] > 0)
                & (clin["pvalue"] <= 1)
            ].copy()
            if not clin.empty:
                clin["association_mode"] = "proxy_clinical_significance"
                proxy = clin

        if proxy.empty:
            return proxy

        has_qual = np.isfinite(proxy["qual"]) & (proxy["qual"] > 0)
        has_af = np.isfinite(proxy["af"]) & (proxy["af"] >= 0) & (proxy["af"] <= 1)
        fallback_mode = np.where(
            has_qual & has_af,
            "proxy_qual_af",
            np.where(has_qual, "proxy_qual_only", np.where(has_af, "proxy_af_only", "proxy")),
        )
        if "association_mode" not in proxy.columns:
            proxy["association_mode"] = fallback_mode
        else:
            missing_mode = proxy["association_mode"].isna() | (proxy["association_mode"] == "")
            proxy.loc[missing_mode, "association_mode"] = fallback_mode[missing_mode]
        assoc = proxy

    assoc["minus_log10_p"] = -np.log10(assoc["pvalue"].clip(lower=1e-300))
    assoc["significant"] = (assoc["pvalue"] <= p_threshold) & (assoc["effect"].abs() >= effect_threshold)
    assoc = assoc.sort_values("pvalue")
    return assoc


def plot_volcano(assoc, out_path, p_threshold, effect_threshold):
    if assoc.empty:
        blank_plot(out_path, "No p-value/effect fields found for variant volcano")
        return

    plot_df = assoc
    if len(plot_df) > 200000:
        plot_df = plot_df.sample(200000, random_state=7)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=plot_df,
        x="effect",
        y="minus_log10_p",
        hue="significant",
        palette={False: "#B0B0B0", True: "#D1495B"},
        alpha=0.55,
        s=14,
        linewidth=0,
    )
    plt.axhline(-math.log10(max(p_threshold, 1e-300)), color="#1D3557", linestyle="--", linewidth=1)
    if effect_threshold > 0:
        plt.axvline(effect_threshold, color="#1D3557", linestyle=":", linewidth=1)
        plt.axvline(-effect_threshold, color="#1D3557", linestyle=":", linewidth=1)
    mode = assoc["association_mode"].iloc[0] if "association_mode" in assoc.columns and not assoc.empty else "reported"
    if mode == "reported":
        plt.title("Variant Association Volcano Plot")
        plt.xlabel("Effect Size (beta or log(OR))")
        plt.ylabel("-log10(p-value)")
    elif mode == "proxy_clinical_significance":
        plt.title("Variant Clinical Significance Volcano (proxy)")
        plt.xlabel("Clinical Significance Effect Proxy")
        plt.ylabel("-log10(proxy confidence)")
    else:
        plt.title("Variant Prioritization Volcano (QUAL/AF proxy)")
        plt.xlabel("Proxy Effect (AF - 0.5)")
        plt.ylabel("-log10(proxy p from QUAL)")
    plt.legend(title="Significant", loc="upper right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_manhattan(assoc, out_path, p_threshold):
    if assoc.empty or "chrom" not in assoc.columns or "pos" not in assoc.columns:
        blank_plot(out_path, "No p-value/chromosome/position data for Manhattan plot")
        return

    man = assoc.copy()
    man["pos"] = pd.to_numeric(man["pos"], errors="coerce")
    man = man.dropna(subset=["pos", "pvalue"])
    if man.empty:
        blank_plot(out_path, "No p-value/chromosome/position data for Manhattan plot")
        return

    if len(man) > 400000:
        man = man.nsmallest(400000, "pvalue")

    chroms = sorted(man["chrom"].unique(), key=chrom_sort_key)

    offsets = {}
    running = 0.0
    tick_pos = []
    tick_lab = []

    for chrom in chroms:
        sub = man[man["chrom"] == chrom]
        if sub.empty:
            continue
        max_pos = sub["pos"].max()
        offsets[chrom] = running
        tick_pos.append(running + sub["pos"].median())
        tick_lab.append(chrom)
        running += max_pos + 1e6

    man["x"] = man.apply(lambda r: r["pos"] + offsets.get(r["chrom"], 0.0), axis=1)

    plt.figure(figsize=(12, 5))
    colors = ["#4C78A8", "#F58518"]
    for idx, chrom in enumerate(chroms):
        sub = man[man["chrom"] == chrom]
        plt.scatter(
            sub["x"],
            sub["minus_log10_p"],
            c=colors[idx % 2],
            s=8,
            alpha=0.6,
            linewidths=0,
        )

    plt.axhline(-math.log10(max(p_threshold, 1e-300)), color="#D1495B", linestyle="--", linewidth=1)
    mode = assoc["association_mode"].iloc[0] if "association_mode" in assoc.columns and not assoc.empty else "reported"
    if mode == "reported":
        title = "Variant Manhattan Plot"
        ylabel = "-log10(p-value)"
    elif mode == "proxy_clinical_significance":
        title = "Variant Clinical Significance Manhattan (proxy)"
        ylabel = "-log10(proxy confidence)"
    else:
        title = "Variant Prioritization Manhattan (QUAL/AF proxy)"
        ylabel = "-log10(proxy p)"
    plt.title(title)
    plt.xlabel("Genomic Position")
    plt.ylabel(ylabel)
    plt.xticks(tick_pos, tick_lab, rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_qq(assoc, out_path, analysis_mode):
    if assoc.empty:
        blank_plot(out_path, "No p-values available for QQ plot")
        return

    mode = assoc["association_mode"].iloc[0] if "association_mode" in assoc.columns and not assoc.empty else "reported"
    if analysis_mode != "gwas_summary" or mode != "reported":
        blank_plot(out_path, "QQ plot is generated only for GWAS-style reported p-values")
        return

    p = pd.to_numeric(assoc["pvalue"], errors="coerce")
    p = p[(p > 0) & (p <= 1)].dropna().sort_values()
    n = len(p)
    if n < 10:
        blank_plot(out_path, "Not enough p-values for QQ plot")
        return

    obs = -np.log10(p.values)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)

    lim = max(obs.max(), exp.max()) * 1.05
    plt.figure(figsize=(6.5, 6.5))
    plt.scatter(exp, obs, s=10, alpha=0.65, color="#2A9D8F", linewidths=0)
    plt.plot([0, lim], [0, lim], linestyle="--", color="#E76F51", linewidth=1)
    plt.title("GWAS QQ Plot")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def extract_heatmap_matrix(vcf_path, selected_keys):
    if not selected_keys:
        return pd.DataFrame()

    target = set(selected_keys)
    samples = []
    rows = {}

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                continue
            if not samples:
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos_raw, _vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if "," in alt:
                continue
            chrom_norm = normalize_chrom(chrom)
            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            key = f"{chrom_norm}:{pos}:{ref}:{alt}"
            if key not in target:
                continue

            fmt_keys = parts[8].split(":")
            doses = [parse_sample_dose(sfield, fmt_keys) for sfield in parts[9:]]
            rows[key] = doses

            if len(rows) == len(target):
                break

    if not rows or not samples:
        return pd.DataFrame()

    mat = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
    return mat


def build_feature_heatmap_matrix(df, keys):
    if not keys:
        return pd.DataFrame()

    sub = df[df["key"].isin(keys)].drop_duplicates("key").set_index("key")
    if sub.empty:
        return pd.DataFrame()

    out = {}
    qual = pd.to_numeric(sub.get("qual", pd.Series(dtype=float)), errors="coerce")
    if qual.notna().any():
        out["qual"] = qual

    af = pd.to_numeric(sub.get("af", pd.Series(dtype=float)), errors="coerce")
    af = af.where((af >= 0) & (af <= 1))
    if af.notna().any():
        out["af"] = af

    pval = pd.to_numeric(sub.get("pvalue", pd.Series(dtype=float)), errors="coerce")
    pval = pval.where((pval > 0) & (pval <= 1))
    if pval.notna().any():
        out["minus_log10_p"] = -np.log10(pval.clip(lower=1e-300))

    effect = pd.to_numeric(sub.get("effect", pd.Series(dtype=float)), errors="coerce")
    if effect.notna().any():
        out["effect"] = effect

    if not out:
        return pd.DataFrame()

    mat = pd.DataFrame(out).reindex(keys).dropna(axis=0, how="all")
    if mat.empty:
        return mat

    # Standardize columns to make mixed scales comparable.
    for col in mat.columns:
        c = mat[col]
        mean = c.mean(skipna=True)
        std = c.std(skipna=True, ddof=0)
        if pd.isna(std) or std == 0:
            mat[col] = c.fillna(mean if not pd.isna(mean) else 0.0)
        else:
            mat[col] = ((c - mean) / std).fillna(0.0)
    return mat


def plot_heatmap(df, vcf_path, out_path, top_n):
    if df.empty:
        blank_plot(out_path, "No variants available for heatmap")
        return

    keys = []
    tmp = df.copy()
    tmp["pvalue"] = pd.to_numeric(tmp.get("pvalue", pd.Series(dtype=float)), errors="coerce")
    tmp["qual"] = pd.to_numeric(tmp.get("qual", pd.Series(dtype=float)), errors="coerce")
    tmp["af"] = pd.to_numeric(tmp.get("af", pd.Series(dtype=float)), errors="coerce")

    ranked = tmp[(tmp["pvalue"] > 0) & (tmp["pvalue"] <= 1)].sort_values("pvalue")
    keys = ranked["key"].head(top_n).tolist()

    if not keys:
        ranked = tmp.sort_values("qual", ascending=False)
        keys = ranked["key"].head(top_n).tolist()

    if not keys:
        ranked = tmp.assign(af_dev=(tmp["af"] - 0.5).abs()).sort_values("af_dev", ascending=False)
        keys = ranked["key"].head(top_n).tolist()

    keys = [k for k in keys if isinstance(k, str) and k]
    if not keys:
        blank_plot(out_path, "No ranked variants available for heatmap")
        return

    mat = extract_heatmap_matrix(vcf_path, keys)
    if mat.empty or mat.shape[1] < 1:
        feature_mat = build_feature_heatmap_matrix(df, keys)
        if feature_mat.empty:
            blank_plot(out_path, "No genotype columns or numeric variant features available for heatmap")
            return

        plt.figure(figsize=(max(7, feature_mat.shape[1] * 1.1), max(6, feature_mat.shape[0] * 0.25)))
        sns.heatmap(feature_mat, cmap="mako", cbar_kws={"label": "Standardized Feature Value"})
        plt.title(f"Variant Feature Heatmap (Top {min(top_n, feature_mat.shape[0])})")
        plt.xlabel("Variant Features")
        plt.ylabel("Variant")
        plt.tight_layout()
        plt.savefig(out_path, dpi=200)
        plt.close()
        return

    mat = mat.reindex(keys).dropna(axis=0, how="all")
    if mat.empty:
        feature_mat = build_feature_heatmap_matrix(df, keys)
        if feature_mat.empty:
            blank_plot(out_path, "Selected variants had no usable genotype values")
            return
        plt.figure(figsize=(max(7, feature_mat.shape[1] * 1.1), max(6, feature_mat.shape[0] * 0.25)))
        sns.heatmap(feature_mat, cmap="mako", cbar_kws={"label": "Standardized Feature Value"})
        plt.title(f"Variant Feature Heatmap (Top {min(top_n, feature_mat.shape[0])})")
        plt.xlabel("Variant Features")
        plt.ylabel("Variant")
        plt.tight_layout()
        plt.savefig(out_path, dpi=200)
        plt.close()
        return

    mat = mat.apply(pd.to_numeric, errors="coerce")
    mat = mat.fillna(mat.mean(axis=1), axis=0).fillna(0.0)

    plt.figure(figsize=(max(8, mat.shape[1] * 0.45), max(6, mat.shape[0] * 0.25)))
    sns.heatmap(mat, cmap="mako", vmin=0, vmax=2, cbar_kws={"label": "Alt Allele Dosage"})
    plt.title(f"Variant Genotype Heatmap (Top {min(top_n, mat.shape[0])})")
    plt.xlabel("Sample")
    plt.ylabel("Variant")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    args = parse_args()

    use_cols = [
        "key",
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "qual",
        "pvalue",
        "effect",
        "af",
        "clinical_significance",
        "functional_class",
        "genes",
    ]
    df = pd.read_csv(args.annotated, sep="\t", dtype=str, usecols=lambda c: c in use_cols).fillna("")

    plot_feature_distribution(df, args.out_feature_plot)
    plot_chrom_distribution(df, args.out_chrom_plot)

    assoc = build_assoc_table(df, args.p_threshold, args.effect_threshold)
    if assoc.empty:
        pd.DataFrame(
            columns=[
                "key",
                "chrom",
                "pos",
                "id",
                "ref",
                "alt",
                "functional_class",
                "genes",
                "pvalue",
                "effect",
                "qual",
                "af",
                "clinical_significance",
                "association_mode",
                "minus_log10_p",
                "significant",
            ]
        ).to_csv(args.out_assoc, sep="\t", index=False)
    else:
        sig = assoc[assoc["significant"]].copy()
        if sig.empty:
            sig = assoc.nsmallest(min(200, len(assoc)), "pvalue").copy()
        cols = [
            "key",
            "chrom",
            "pos",
            "id",
            "ref",
            "alt",
            "functional_class",
            "genes",
            "pvalue",
            "effect",
            "qual",
            "af",
            "clinical_significance",
            "association_mode",
            "minus_log10_p",
            "significant",
        ]
        sig = sig[[c for c in cols if c in sig.columns]]
        sig.to_csv(args.out_assoc, sep="\t", index=False)

    plot_volcano(assoc, args.out_volcano, args.p_threshold, args.effect_threshold)
    plot_manhattan(assoc, args.out_manhattan, args.p_threshold)
    plot_qq(assoc, args.out_qq, args.analysis_mode)
    plot_heatmap(df, args.vcf, args.out_heatmap, max(1, args.heatmap_top_n))


if __name__ == "__main__":
    main()
