#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


REQUIRED_CANONICAL = ["CHROM", "POS", "REF", "ALT"]

ALIASES = {
    "CHROM": ["CHROM", "#CHROM", "chrom", "chr", "chromosome", "#dbSNP_hg38_chr"],
    "POS": ["POS", "position", "bp", "dbSNP_hg38_position"],
    "REF": ["REF", "ref", "reference", "ref_allele", "other_allele", "allele2"],
    "ALT": ["ALT", "alt", "alternate", "alt_allele", "nonref_allele", "effect_allele", "allele1"],
    "ID": ["ID", "id", "rsid", "rsID", "SNP", "Top SNP", "variant_id", "markername"],
    "QUAL": ["QUAL", "qual"],
    "FILTER": ["FILTER", "filter"],
    "INFO": ["INFO", "info"],
    "GT": ["GT", "gt", "genotype"],
    "SAMPLE": ["SAMPLE", "sample", "sample_id"],
    "P": ["P", "PVALUE", "PVAL", "P-value", "pvalue", "p_value", "pval"],
    "BETA": ["BETA", "beta", "EFFECT", "effect", "nonref_effect"],
    "OR": ["OR", "or", "OR_nonref"],
    "EAF": ["EAF", "eaf", "AF", "af", "effect_allele_frequency"],
    "SE": ["SE", "se", "stderr", "standard_error"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert variant CSV/TSV tables to VCF with column alias mapping and diagnostics."
    )
    parser.add_argument("--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("--output", required=True, help="Output VCF file")
    parser.add_argument(
        "--ref-fasta",
        default="",
        help="Optional reference FASTA used to derive REF for SNP-only raw tables missing REF",
    )
    parser.add_argument(
        "--report",
        default="",
        help="Optional TSV report describing detected column mappings and missing fields",
    )
    parser.add_argument(
        "--normalized-table",
        default="",
        help="Optional standardized TSV with canonical columns used for conversion",
    )
    return parser.parse_args()


def norm_name(text):
    return re.sub(r"[^a-z0-9]+", "", str(text).strip().lower())


def clean_text(val):
    if pd.isna(val):
        return ""
    text = str(val).strip()
    if text in {"", ".", "NA", "NaN", "nan", "NAN", "None", "none"}:
        return ""
    return text


def pick_column(df, aliases):
    # Exact alias matches first, then normalized-name matches.
    cols = list(df.columns)
    for a in aliases:
        if a in cols:
            return a
    norm_to_col = {norm_name(c): c for c in cols}
    for a in aliases:
        hit = norm_to_col.get(norm_name(a))
        if hit:
            return hit
    return None


def build_mapping(df):
    mapping = {}
    for canonical, aliases in ALIASES.items():
        col = pick_column(df, aliases)
        if col:
            mapping[canonical] = col
    return mapping


def write_report(path, rows):
    if not path:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["item", "status", "detail"])
        writer.writerows(rows)


def to_info_token(row, key):
    return clean_text(row.get(key, ""))


def parse_pos(text):
    txt = clean_text(text)
    if not txt:
        return None
    try:
        return int(txt)
    except ValueError:
        try:
            fval = float(txt)
        except ValueError:
            return None
        if not fval.is_integer():
            return None
        return int(fval)


def is_single_base_allele(allele):
    a = clean_text(allele).upper()
    return len(a) == 1 and a in {"A", "C", "G", "T", "N"}


def fetch_reference_bases(ref_fasta, regions):
    if not regions:
        return {}
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False) as tfh:
        region_file = tfh.name
        for reg in regions:
            tfh.write(reg + "\n")
    try:
        cmd = ["samtools", "faidx", "-r", region_file, ref_fasta]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "samtools faidx failed while deriving REF alleles: "
                + (proc.stderr.strip() or "unknown error")
            )
        seqs = {}
        current = None
        chunks = []
        for line in proc.stdout.splitlines():
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks).upper()
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            seqs[current] = "".join(chunks).upper()
        return seqs
    finally:
        Path(region_file).unlink(missing_ok=True)


def load_fasta_contigs(ref_fasta):
    fai = Path(f"{ref_fasta}.fai")
    if not fai.exists():
        return set()
    contigs = set()
    with open(fai, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            contigs.add(line.split("\t", 1)[0].strip())
    return contigs


def resolve_contig_name(chrom, fasta_contigs):
    c = clean_text(chrom)
    if not c or not fasta_contigs:
        return c
    if c in fasta_contigs:
        return c
    if c.startswith("chr") and c[3:] in fasta_contigs:
        return c[3:]
    if (not c.startswith("chr")) and f"chr{c}" in fasta_contigs:
        return f"chr{c}"
    if c in {"MT", "M"} and "chrM" in fasta_contigs:
        return "chrM"
    if c == "chrM" and "MT" in fasta_contigs:
        return "MT"
    return c


def derive_missing_ref(df_norm, ref_fasta, report_rows):
    if "REF" not in df_norm.columns:
        df_norm["REF"] = ""
    missing_idx = []
    unique_regions = []
    seen_regions = set()
    bad_rows = []
    fasta_contigs = load_fasta_contigs(ref_fasta) if ref_fasta else set()

    for idx, row in df_norm.iterrows():
        ref = clean_text(row.get("REF", ""))
        if ref:
            continue
        chrom = clean_text(row.get("CHROM", ""))
        pos = parse_pos(row.get("POS", ""))
        alt = clean_text(row.get("ALT", "")).upper()
        if not chrom or pos is None:
            bad_rows.append((idx, "cannot derive REF: missing/invalid CHROM or POS"))
            continue
        if not is_single_base_allele(alt):
            bad_rows.append((idx, f"cannot derive REF for non-SNP/ambiguous ALT '{alt or '.'}'"))
            continue
        chrom_for_ref = resolve_contig_name(chrom, fasta_contigs)
        if chrom_for_ref != chrom:
            df_norm.at[idx, "CHROM"] = chrom_for_ref
        region = f"{chrom_for_ref}:{pos}-{pos}"
        missing_idx.append((idx, region))
        if region not in seen_regions:
            seen_regions.add(region)
            unique_regions.append(region)

    if not missing_idx:
        report_rows.append(("ref_derivation", "OK", "No REF derivation needed"))
        return []

    if not ref_fasta:
        report_rows.append(
            ("ref_derivation", "ERROR", "REF missing but no --ref-fasta provided to derive it")
        )
        return bad_rows + [(idx, "REF missing and no reference FASTA provided") for idx, _ in missing_idx]

    seq_by_region = fetch_reference_bases(ref_fasta, unique_regions)
    unresolved = list(bad_rows)
    derived = 0
    for idx, region in missing_idx:
        seq = clean_text(seq_by_region.get(region, ""))
        if is_single_base_allele(seq):
            df_norm.at[idx, "REF"] = seq
            derived += 1
        else:
            unresolved.append((idx, f"reference lookup failed for {region}"))
    report_rows.append(("ref_derivation", "OK", f"Derived REF for {derived} row(s)"))
    if unresolved:
        report_rows.append(
            ("ref_derivation_unresolved", "WARNING", f"Could not derive REF for {len(unresolved)} row(s)")
        )
    return unresolved


def build_normalized_table(df, mapping, report_rows):
    df_norm = pd.DataFrame(index=df.index)
    for canonical, src in mapping.items():
        df_norm[canonical] = df[src].astype(str)

    missing_required = [k for k in REQUIRED_CANONICAL if k not in df_norm.columns]
    if missing_required:
        status = "WARNING" if missing_required == ["REF"] else "ERROR"
        report_rows.append(
            (
                "required_fields",
                status,
                "Missing required fields after alias mapping: " + ", ".join(missing_required),
            )
        )
    else:
        report_rows.append(("required_fields", "OK", "All required fields mapped (REF may still be empty in rows)"))

    return df_norm, missing_required


def normalize_rows_for_vcf(df_norm, report_rows):
    keep_rows = []
    row_errors = []
    for idx, row in df_norm.iterrows():
        chrom = clean_text(row.get("CHROM", ""))
        pos = parse_pos(row.get("POS", ""))
        ref = clean_text(row.get("REF", "")).upper()
        alt = clean_text(row.get("ALT", "")).upper()
        if not chrom or pos is None or not ref or not alt:
            row_errors.append((idx, "missing CHROM/POS/REF/ALT after normalization"))
            continue
        if "," in alt:
            alt = ",".join([x.strip().upper() for x in alt.split(",") if x.strip()])
        if not alt:
            row_errors.append((idx, "ALT empty after normalization"))
            continue
        if ref == alt:
            row_errors.append((idx, "REF equals ALT"))
            continue
        new_row = row.copy()
        new_row["CHROM"] = chrom
        new_row["POS"] = str(pos)
        new_row["REF"] = ref
        new_row["ALT"] = alt
        keep_rows.append(new_row)

    if row_errors and keep_rows:
        report_rows.append(("row_validation", "WARNING", f"Dropping invalid rows: {len(row_errors)}"))
    elif row_errors:
        report_rows.append(("row_validation", "ERROR", f"Invalid rows: {len(row_errors)}"))
    else:
        report_rows.append(("row_validation", "OK", "All rows valid for VCF conversion"))
    return pd.DataFrame(keep_rows), row_errors


def print_diagnostics(report_rows, unresolved_rows, row_errors):
    for item, status, detail in report_rows:
        print(f"[{status}] {item}: {detail}")
    for idx, msg in unresolved_rows[:10]:
        print(f"[WARNING] row {idx + 2}: {msg}")  # +2 for 1-based with header
    if len(unresolved_rows) > 10:
        print(f"[WARNING] ... {len(unresolved_rows) - 10} more REF derivation row issues")
    for idx, msg in row_errors[:10]:
        print(f"[WARNING] row {idx + 2}: {msg}")
    if len(row_errors) > 10:
        print(f"[WARNING] ... {len(row_errors) - 10} more row validation issues")


def main():
    args = parse_args()
    in_path = str(args.input)
    sep = "\t" if in_path.lower().endswith(".tsv") else ","
    df = pd.read_csv(in_path, sep=sep, dtype=str).fillna("")

    report_rows = [
        ("input_file", "OK", in_path),
        ("input_rows", "OK", str(len(df))),
        ("input_columns", "OK", ", ".join([str(c) for c in df.columns])),
    ]

    mapping = build_mapping(df)
    for key in sorted(mapping):
        report_rows.append(("column_mapping", "OK", f"{key} <- {mapping[key]}"))

    unmapped_required = [k for k in REQUIRED_CANONICAL if k not in mapping]
    if unmapped_required:
        report_rows.append(
            ("column_mapping_required", "WARNING", "Unmapped required fields: " + ", ".join(unmapped_required))
        )
    else:
        report_rows.append(("column_mapping_required", "OK", "All required fields mapped"))

    df_norm, missing_required = build_normalized_table(df, mapping, report_rows)
    unresolved_rows = []
    if "REF" not in df_norm.columns or (df_norm.get("REF", pd.Series(dtype=str)).astype(str).str.strip() == "").any():
        try:
            unresolved_rows = derive_missing_ref(df_norm, args.ref_fasta, report_rows)
        except Exception as exc:
            report_rows.append(("ref_derivation", "ERROR", str(exc)))
            unresolved_rows = [(-1, str(exc))]

    missing_after_derivation = []
    for field in REQUIRED_CANONICAL:
        if field not in df_norm.columns:
            missing_after_derivation.append(field)
        elif df_norm[field].astype(str).map(clean_text).eq("").all():
            missing_after_derivation.append(field)
    if missing_after_derivation:
        report_rows.append(
            ("required_fields_final", "ERROR", "Missing/unresolved required fields: " + ", ".join(missing_after_derivation))
        )
    else:
        report_rows.append(("required_fields_final", "OK", "Required fields resolved"))

    df_vcf, row_errors = normalize_rows_for_vcf(df_norm, report_rows)
    if df_vcf.empty:
        report_rows.append(("output_rows", "ERROR", "No valid rows available for VCF output"))
    else:
        report_rows.append(("output_rows", "OK", str(len(df_vcf))))

    if args.normalized_table:
        out_tab = Path(args.normalized_table)
        out_tab.parent.mkdir(parents=True, exist_ok=True)
        preferred_cols = ["CHROM", "POS", "ID", "REF", "ALT", "P", "BETA", "OR", "EAF", "SE", "QUAL", "FILTER", "INFO", "GT", "SAMPLE"]
        cols = [c for c in preferred_cols if c in df_norm.columns] + [c for c in df_norm.columns if c not in preferred_cols]
        df_norm[cols].to_csv(out_tab, sep="\t", index=False)
        report_rows.append(("normalized_table", "OK", str(out_tab)))

    write_report(args.report, report_rows)
    print_diagnostics(report_rows, unresolved_rows, row_errors)

    if missing_after_derivation or df_vcf.empty:
        raise SystemExit(
            "Cannot convert input table to VCF. See report for missing fields/row errors."
        )

    has_gt = "GT" in df_vcf.columns
    sample_name = "SAMPLE1"
    if "SAMPLE" in df_vcf.columns and df_vcf["SAMPLE"].astype(str).map(clean_text).ne("").any():
        sample_name = str(df_vcf["SAMPLE"].astype(str).map(clean_text).loc[lambda s: s.ne("")].iloc[0])

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=csv_to_vcf.py\n")
        out.write('##INFO=<ID=SRC,Number=1,Type=String,Description="Source format">\n')
        out.write('##INFO=<ID=PVALUE,Number=1,Type=Float,Description="Association p-value from input table">\n')
        out.write('##INFO=<ID=BETA,Number=1,Type=Float,Description="Association beta/effect from input table">\n')
        out.write('##INFO=<ID=OR,Number=1,Type=Float,Description="Association odds ratio from input table">\n')
        out.write('##INFO=<ID=EAF,Number=1,Type=Float,Description="Effect allele frequency from input table">\n')
        out.write('##INFO=<ID=SE,Number=1,Type=Float,Description="Standard error from input table">\n')
        if has_gt:
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        else:
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for _, row in df_vcf.iterrows():
            chrom = str(row["CHROM"])
            pos = int(row["POS"])
            ref = str(row["REF"])
            alt = str(row["ALT"])
            var_id = clean_text(row.get("ID", "")) or "."
            qual = clean_text(row.get("QUAL", "")) or "."
            flt = clean_text(row.get("FILTER", "")) or "PASS"
            info_tokens = []
            raw_info = clean_text(row.get("INFO", ""))
            if raw_info:
                info_tokens.append(raw_info)
            info_tokens.append("SRC=TABLE")

            pval = to_info_token(row, "P")
            if pval:
                info_tokens.append(f"PVALUE={pval}")
            beta = to_info_token(row, "BETA")
            if beta:
                info_tokens.append(f"BETA={beta}")
            odds_ratio = to_info_token(row, "OR")
            if odds_ratio:
                info_tokens.append(f"OR={odds_ratio}")
            eaf = to_info_token(row, "EAF")
            if eaf:
                info_tokens.append(f"EAF={eaf}")
            se = to_info_token(row, "SE")
            if se:
                info_tokens.append(f"SE={se}")

            info = ";".join([tok for tok in info_tokens if tok]) or "."
            if has_gt:
                gt = clean_text(row.get("GT", "")) or "0/1"
                out.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\tGT\t{gt}\n")
            else:
                out.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\n")


if __name__ == "__main__":
    main()
