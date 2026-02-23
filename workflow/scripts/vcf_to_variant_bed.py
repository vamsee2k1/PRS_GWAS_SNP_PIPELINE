#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import re


PVAL_KEYS = [
    "P",
    "PVALUE",
    "PVAL",
    "P_VALUE",
    "PVAL_U",
    "P_BOLT_LMM_INF",
    "P_BOLT_LMM",
    "PVAL_FIXED",
    "PVAL_RANDOM",
]
MLOG10_PVAL_KEYS = [
    "LP",
    "MLOG10P",
    "MLOGP",
    "NEG_LOG10_P",
    "MINUS_LOG10_P",
    "LOG10P",
]
EFFECT_KEYS = [
    "BETA",
    "EFFECT",
    "EFFECTSIZE",
    "ES",
    "LOG_ODDS",
    "B",
    "ESTIMATE",
    "BETA_FIXED",
    "BETA_RANDOM",
]
OR_KEYS = ["OR", "ODDSRATIO", "ODDS_RATIO"]
AF_KEYS = ["AF", "EAF", "AAF", "ALT_AF", "MAF", "AF_TGP", "AF_EXAC", "AF_ESP"]
REPORTED_GENE_KEYS = ["GENEINFO", "GENE", "GENES", "GENE_SYMBOL", "SYMBOL"]
CONSEQUENCE_KEYS = ["MC", "CONSEQUENCE", "SO", "ANN", "CSQ"]
CLIN_SIG_KEYS = ["CLNSIG", "CLIN_SIG", "CLINICAL_SIGNIFICANCE", "SIGNIFICANCE"]


NC_CHROM_RE = re.compile(r"^NC_0*([0-9]+)\\.")


def parse_args():
    parser = argparse.ArgumentParser(description="Convert VCF into BED + base variant table")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ")
    parser.add_argument("--out-bed", required=True, help="Output BED for feature overlap")
    parser.add_argument("--out-table", required=True, help="Output base variant TSV")
    return parser.parse_args()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def normalize_chrom(chrom):
    c = str(chrom).strip()
    if not c:
        return c

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


def info_to_dict(info_str):
    out = {}
    if info_str in {"", "."}:
        return out
    for token in info_str.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            out[key.upper()] = value
        else:
            out[token.upper()] = "1"
    return out


def first_float(raw):
    if raw is None:
        return None
    token = str(raw).split(",")[0].strip()
    if token in {"", ".", "NA", "nan", "NAN"}:
        return None
    try:
        val = float(token)
    except ValueError:
        return None
    if math.isnan(val) or math.isinf(val):
        return None
    return val


def split_multi_value(raw):
    if raw is None:
        return []
    token = str(raw).strip()
    if token in {"", ".", "NA", "nan", "NAN"}:
        return []
    return [x.strip() for x in token.split(",")]


def alt_index_float(raw, alt_index):
    vals = split_multi_value(raw)
    if not vals:
        return None
    idx = alt_index if alt_index < len(vals) else 0
    return first_float(vals[idx])


def gt_alt_dosage(gt):
    if gt in {"", ".", "./.", ".|."}:
        return None
    dose = 0.0
    for allele in gt.replace("|", "/").split("/"):
        if allele == ".":
            return None
        if allele == "0":
            continue
        if allele == "1":
            dose += 1.0
        else:
            # Multiallelic genotype dosage cannot be assigned to a single ALT.
            return None
    return dose


def parse_sample_alt_dosage(sample_field, fmt_keys):
    values = sample_field.split(":")
    fdict = {k.upper(): values[i] if i < len(values) else "" for i, k in enumerate(fmt_keys)}

    ds = fdict.get("DS", "")
    if ds not in {"", "."}:
        token = ds.split(",")[0]
        val = first_float(token)
        if val is not None:
            return val, "DS"

    gt = fdict.get("GT", "")
    if gt:
        dose = gt_alt_dosage(gt)
        if dose is not None:
            return dose, "GT"

    return None, ""


def estimate_af_from_samples(format_raw, sample_fields):
    total_samples = len(sample_fields)
    if total_samples == 0 or not format_raw or format_raw == ".":
        return None, "", 0, total_samples

    fmt_keys = [x.strip() for x in str(format_raw).split(":")]
    alt_sum = 0.0
    called = 0
    source = ""

    for sample_field in sample_fields:
        dose, src = parse_sample_alt_dosage(sample_field, fmt_keys)
        if dose is None:
            continue
        alt_sum += dose
        called += 1
        if not source:
            source = src

    if called == 0:
        return None, "", 0, total_samples

    af = alt_sum / (2.0 * called)
    af = max(0.0, min(1.0, af))
    return af, f"sample_{source.lower()}" if source else "sample", called, total_samples


def parse_af(info, alt_index):
    for key in AF_KEYS:
        if key in info:
            val = alt_index_float(info[key], alt_index)
            if val is not None:
                return max(0.0, min(1.0, val)), key

    if "AC" in info and "AN" in info:
        ac = alt_index_float(info["AC"], alt_index)
        an = first_float(info["AN"])
        if ac is not None and an is not None and an > 0:
            af = ac / an
            return max(0.0, min(1.0, af)), "AC/AN"

    return None, ""


def parse_reported_genes(info):
    genes = set()
    for key in REPORTED_GENE_KEYS:
        if key not in info:
            continue
        raw = str(info[key]).strip()
        if raw in {"", "."}:
            continue

        if key == "GENEINFO":
            chunks = raw.replace(",", "|").split("|")
            for chunk in chunks:
                chunk = chunk.strip()
                if not chunk:
                    continue
                gene = chunk.split(":", 1)[0].strip()
                if gene and gene != ".":
                    genes.add(gene)
        else:
            chunks = raw.replace("|", ",").split(",")
            for chunk in chunks:
                gene = chunk.strip()
                if gene and gene != ".":
                    genes.add(gene)

    return ";".join(sorted(genes))


def parse_reported_consequence(info):
    terms = set()
    for key in CONSEQUENCE_KEYS:
        if key not in info:
            continue
        raw = str(info[key]).strip()
        if raw in {"", "."}:
            continue

        if key == "MC":
            entries = raw.split(",")
            for entry in entries:
                toks = [x.strip() for x in entry.split("|") if x.strip()]
                if len(toks) >= 2:
                    terms.add(toks[1])
                elif toks:
                    terms.add(toks[0])
        elif key in {"ANN", "CSQ"}:
            entries = raw.split(",")
            for entry in entries:
                toks = [x.strip() for x in entry.split("|")]
                if len(toks) >= 2 and toks[1]:
                    terms.add(toks[1])
                elif toks and toks[0]:
                    terms.add(toks[0])
        else:
            entries = raw.replace("|", ",").split(",")
            for entry in entries:
                term = entry.strip()
                if term and term != ".":
                    terms.add(term)

    return ";".join(sorted(terms))


def parse_clinical_significance(info):
    for key in CLIN_SIG_KEYS:
        if key not in info:
            continue
        raw = str(info[key]).strip()
        if raw in {"", "."}:
            continue
        token = raw.split("|")[0].split(",")[0].strip()
        if not token:
            continue
        token = token.replace(" ", "_").replace("/", "_").lower()
        return token
    return ""


def parse_pvalue(info):
    for key in PVAL_KEYS:
        if key in info:
            val = first_float(info[key])
            if val is None:
                continue
            if val <= 0:
                return 1e-300, key
            if val <= 1:
                return val, key

    for key in MLOG10_PVAL_KEYS:
        if key in info:
            mlog = first_float(info[key])
            if mlog is None:
                continue
            if mlog < 0:
                continue
            try:
                pval = 10 ** (-mlog)
            except OverflowError:
                pval = 0.0
            if pval <= 0:
                pval = 1e-300
            return pval, f"10^-{key}"
    return None, ""


def parse_effect(info):
    for key in EFFECT_KEYS:
        if key in info:
            val = first_float(info[key])
            if val is not None:
                return val, key
    for key in OR_KEYS:
        if key in info:
            val = first_float(info[key])
            if val is not None and val > 0:
                return math.log(val), f"log({key})"
    return None, ""


def format_to_dict(format_raw, sample_raw):
    if not format_raw or not sample_raw or format_raw == "." or sample_raw == ".":
        return {}
    keys = str(format_raw).split(":")
    values = str(sample_raw).split(":")
    out = {}
    for idx, key in enumerate(keys):
        if not key:
            continue
        out[key.upper()] = values[idx] if idx < len(values) else ""
    return out


def fmt_float(x):
    if x is None:
        return ""
    return f"{x:.10g}"


def main():
    args = parse_args()

    with open_text(args.vcf) as fh, open(args.out_bed, "w", encoding="utf-8") as bed_fh, open(
        args.out_table, "w", newline="", encoding="utf-8"
    ) as table_fh:
        writer = csv.writer(table_fh, delimiter="\t")
        writer.writerow(
            [
                "key",
                "chrom",
                "pos",
                "id",
                "ref",
                "alt",
                "qual",
                "filter",
                "pvalue",
                "pvalue_source",
                "effect",
                "effect_source",
                "af",
                "af_source",
                "called_samples",
                "total_samples",
                "reported_genes",
                "reported_consequence",
                "clinical_significance",
            ]
        )

        sample_count = 0
        for line in fh:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_parts = line.rstrip("\n").split("\t")
                sample_count = max(len(header_parts) - 9, 0)
                continue
            if line.startswith("#"):
                continue

            # Keep split bounded unless we explicitly need all sample fields.
            parts = line.rstrip("\n").split("\t", 9)
            if len(parts) < 8:
                continue

            chrom_raw, pos_raw, var_id, ref, alt_field, qual_raw, filt, info_raw = parts[:8]
            format_raw = parts[8] if len(parts) >= 9 else ""
            sample_blob = parts[9] if len(parts) >= 10 else ""
            first_sample_raw = sample_blob.split("\t", 1)[0] if sample_blob else ""
            chrom = normalize_chrom(chrom_raw)

            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            qual = first_float(qual_raw)
            info = info_to_dict(info_raw)
            pval, p_src = parse_pvalue(info)
            effect, e_src = parse_effect(info)
            alts = [x.strip() for x in alt_field.split(",") if x.strip()]

            sample_af = None
            sample_af_source = ""
            called_samples = 0
            total_samples = sample_count
            prefetched_single_af = None
            prefetched_single_af_src = ""

            if len(alts) == 1:
                prefetched_single_af, prefetched_single_af_src = parse_af(info, 0)

            # AF is often already present in INFO (AF or AC/AN), which avoids expensive per-sample parsing.
            if sample_blob and len(alts) == 1 and prefetched_single_af is None:
                sample_fields = sample_blob.split("\t")
                sample_af, sample_af_source, called_samples, total_samples = estimate_af_from_samples(
                    format_raw, sample_fields
                )

            if format_raw and first_sample_raw:
                fmt_data = format_to_dict(format_raw, first_sample_raw)
                if pval is None:
                    pval, p_src = parse_pvalue(fmt_data)
                if effect is None:
                    effect, e_src = parse_effect(fmt_data)

            reported_genes = parse_reported_genes(info)
            reported_consequence = parse_reported_consequence(info)
            clinical_significance = parse_clinical_significance(info)
            for alt_index, alt in enumerate(alts):
                alt = alt.strip()
                if not alt:
                    continue

                if len(alts) == 1 and prefetched_single_af is not None:
                    af, af_src = prefetched_single_af, prefetched_single_af_src
                else:
                    af, af_src = parse_af(info, alt_index)
                if af is None and sample_af is not None and len(alts) == 1:
                    af = sample_af
                    af_src = sample_af_source

                key = f"{chrom}:{pos}:{ref}:{alt}"
                start0 = max(pos - 1, 0)
                bed_fh.write(f"{chrom}\t{start0}\t{pos}\t{key}\n")
                writer.writerow(
                    [
                        key,
                        chrom,
                        pos,
                        var_id,
                        ref,
                        alt,
                        fmt_float(qual),
                        filt,
                        fmt_float(pval),
                        p_src,
                        fmt_float(effect),
                        e_src,
                        fmt_float(af),
                        af_src,
                        called_samples,
                        total_samples,
                        reported_genes,
                        reported_consequence,
                        clinical_significance,
                    ]
                )


if __name__ == "__main__":
    main()
