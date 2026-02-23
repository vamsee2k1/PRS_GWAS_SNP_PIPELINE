#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Preflight validation for pipeline resources and user inputs."
    )
    p.add_argument("--mode", required=True, choices=["full", "variant_only"])
    p.add_argument("--assay", required=True, choices=["dna", "rna"])
    p.add_argument("--read-type", required=True, choices=["short", "long"])
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--reference-fasta", required=True)
    p.add_argument("--reference-gtf", required=True)
    p.add_argument("--reference-star-index", default="")
    p.add_argument("--expected-build", default="")
    p.add_argument("--variants-input", default="")
    p.add_argument("--prs-weights", default="")
    p.add_argument("--prs-weights-format", default="auto")
    p.add_argument("--gene-sets", default="")
    p.add_argument("--annotation-method", default="overlap")
    p.add_argument("--snpeff-database", default="")
    p.add_argument("--enforce-build-match", action="store_true")
    p.add_argument("--fail-on-warning", action="store_true")
    p.add_argument("--out-report", required=True)
    return p.parse_args()


class Reporter:
    def __init__(self):
        self.rows = []
        self.errors = 0
        self.warnings = 0

    def ok(self, check, detail):
        self.rows.append((check, "OK", detail))

    def warning(self, check, detail):
        self.rows.append((check, "WARNING", detail))
        self.warnings += 1

    def error(self, check, detail):
        self.rows.append((check, "ERROR", detail))
        self.errors += 1

    def write(self, out_path):
        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["check", "status", "detail"])
            w.writerows(self.rows)


def looks_like_gzip(path):
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            _ = fh.readline()
        return True
    except OSError:
        return False


def first_non_comment_line(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\n")
    return ""


def first_non_comment_line_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\n")
    return ""


def detect_vcf_build(path, is_gz):
    opener = gzip.open if is_gz else open
    hints = []
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i > 5000:
                    break
                if not line.startswith("#"):
                    break
                low = line.lower()
                if "assembly=b37" in low or "grch37" in low or "hg19" in low:
                    hints.append("GRCh37")
                if "assembly=grch38" in low or "grch38" in low or "hg38" in low:
                    hints.append("GRCh38")
    except OSError:
        return "unknown"

    if "GRCh38" in hints and "GRCh37" in hints:
        return "mixed_or_ambiguous"
    if "GRCh38" in hints:
        return "GRCh38"
    if "GRCh37" in hints:
        return "GRCh37"
    return "unknown"


def validate_tsv_columns(path, required_cols):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in required_cols if c not in header]
        return missing, header


VARIANT_TABLE_ALIASES = {
    "CHROM": ["CHROM", "#CHROM", "chrom", "chr", "chromosome", "#dbSNP_hg38_chr"],
    "POS": ["POS", "position", "bp", "dbSNP_hg38_position"],
    "REF": ["REF", "ref", "reference", "ref_allele", "other_allele", "allele2"],
    "ALT": ["ALT", "alt", "alternate", "alt_allele", "nonref_allele", "effect_allele", "allele1"],
}


def norm_col_name(name):
    return re.sub(r"[^a-z0-9]+", "", str(name or "").strip().lower())


def suggest_variant_table_mapping(header):
    mapping = {}
    norm_to_orig = {norm_col_name(h): h for h in header}
    for canonical, aliases in VARIANT_TABLE_ALIASES.items():
        hit = None
        for a in aliases:
            if a in header:
                hit = a
                break
        if hit is None:
            for a in aliases:
                cand = norm_to_orig.get(norm_col_name(a))
                if cand:
                    hit = cand
                    break
        if hit:
            mapping[canonical] = hit
    return mapping


def normalize_build_name(name):
    n = str(name or "").strip().lower()
    if n in {"grch38", "hg38"}:
        return "GRCh38"
    if n in {"grch37", "hg19"}:
        return "GRCh37"
    return ""


def main():
    args = parse_args()
    rpt = Reporter()

    # Core references
    fasta = Path(args.reference_fasta)
    if fasta.exists() and fasta.is_file():
        rpt.ok("reference_fasta_exists", str(fasta))
    else:
        rpt.error("reference_fasta_exists", f"Missing FASTA: {fasta}")

    fai = Path(f"{args.reference_fasta}.fai")
    if fai.exists() and fai.is_file():
        rpt.ok("reference_fasta_index_exists", str(fai))
    else:
        rpt.error(
            "reference_fasta_index_exists",
            f"Missing FASTA index (.fai): {fai}",
        )

    gtf = Path(args.reference_gtf)
    if gtf.exists() and gtf.is_file():
        rpt.ok("reference_gtf_exists", str(gtf))
        line = first_non_comment_line(str(gtf))
        if line and len(line.split("\t")) >= 9:
            rpt.ok("reference_gtf_format", "GTF has >=9 tab-separated columns")
        else:
            rpt.error(
                "reference_gtf_format",
                "GTF does not appear to contain valid 9-column feature rows",
            )
    else:
        rpt.error("reference_gtf_exists", f"Missing GTF: {gtf}")

    if args.assay == "rna" and args.read_type == "short":
        if not args.reference_star_index:
            rpt.error(
                "star_index_path",
                "RNA short-read mode requires reference.star_index",
            )
        else:
            star_dir = Path(args.reference_star_index)
            if star_dir.exists() and star_dir.is_dir():
                marker_files = [star_dir / "Genome", star_dir / "SA", star_dir / "chrName.txt"]
                if all(f.exists() for f in marker_files):
                    rpt.ok("star_index_ready", f"STAR index appears ready: {star_dir}")
                else:
                    rpt.warning(
                        "star_index_ready",
                        (
                            "STAR index directory exists but standard index files are missing; "
                            "pipeline will attempt to build index"
                        ),
                    )
            else:
                rpt.warning(
                    "star_index_ready",
                    "STAR index directory not found; pipeline will attempt to build index",
                )

    # Sample/metadata checks for full mode
    if args.mode != "variant_only":
        sample_path = Path(args.samplesheet)
        if sample_path.exists() and sample_path.is_file():
            rpt.ok("samplesheet_exists", str(sample_path))
            missing, _ = validate_tsv_columns(str(sample_path), ["sample", "fastq_1", "fastq_2"])
            if missing:
                rpt.error(
                    "samplesheet_columns",
                    f"samplesheet missing required columns: {', '.join(missing)}",
                )
            else:
                rpt.ok("samplesheet_columns", "Required columns present")

            try:
                with open(sample_path, "r", encoding="utf-8", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                if not rows:
                    rpt.error("samplesheet_rows", "samplesheet contains no data rows")
                else:
                    missing_fastq = 0
                    for row in rows:
                        r1 = str(row.get("fastq_1", "")).strip()
                        if not r1 or not Path(r1).exists():
                            missing_fastq += 1
                    if missing_fastq > 0:
                        rpt.error(
                            "samplesheet_fastq_paths",
                            f"{missing_fastq} sample rows have missing fastq_1 paths",
                        )
                    else:
                        rpt.ok("samplesheet_fastq_paths", "All fastq_1 paths exist")
            except OSError as exc:
                rpt.error("samplesheet_read", f"Could not parse samplesheet: {exc}")
        else:
            rpt.error("samplesheet_exists", f"Missing samplesheet: {sample_path}")

        metadata_path = Path(args.metadata)
        if metadata_path.exists() and metadata_path.is_file():
            rpt.ok("metadata_exists", str(metadata_path))
            req_meta_cols = ["sample"] + (["condition"] if args.assay == "rna" else [])
            missing, _ = validate_tsv_columns(str(metadata_path), req_meta_cols)
            if missing:
                rpt.error(
                    "metadata_columns",
                    f"metadata missing required columns: {', '.join(missing)}",
                )
            else:
                rpt.ok("metadata_columns", "Required columns present")
        else:
            rpt.error("metadata_exists", f"Missing metadata: {metadata_path}")

    # Variant input checks for variant-only mode or external override
    if args.mode == "variant_only" or args.variants_input:
        if not args.variants_input:
            rpt.error("variants_input_present", "variants_input is required in variant_only mode")
        else:
            vpath = Path(args.variants_input)
            if vpath.exists() and vpath.is_file():
                rpt.ok("variants_input_exists", str(vpath))
            else:
                rpt.error("variants_input_exists", f"Missing variants_input: {vpath}")

            low = str(vpath).lower()
            is_vcf = low.endswith(".vcf")
            is_vcfgz = low.endswith(".vcf.gz")
            is_tsv = low.endswith(".tsv")
            is_csv = low.endswith(".csv")
            if not (is_vcf or is_vcfgz or is_tsv or is_csv):
                rpt.error(
                    "variants_input_format",
                    "variants_input must end with .vcf, .vcf.gz, .csv, or .tsv",
                )
            else:
                rpt.ok("variants_input_format", f"Accepted extension detected for {vpath.name}")

            if vpath.exists() and is_vcfgz:
                if looks_like_gzip(str(vpath)):
                    rpt.ok("variants_vcfgz_integrity", "VCF.GZ appears readable")
                else:
                    rpt.error("variants_vcfgz_integrity", "VCF.GZ failed gzip integrity read")

                tbi = Path(f"{vpath}.tbi")
                csi = Path(f"{vpath}.csi")
                if tbi.exists() or csi.exists():
                    rpt.ok("variants_vcfgz_index", "Found .tbi or .csi index")
                else:
                    rpt.warning(
                        "variants_vcfgz_index",
                        "No .tbi/.csi found (pipeline can still reindex, but slower)",
                    )

            if vpath.exists() and is_vcf:
                first_data = first_non_comment_line(str(vpath))
                if first_data:
                    rpt.ok("variants_vcf_has_records", "VCF has at least one non-header line")
                else:
                    rpt.warning("variants_vcf_has_records", "VCF has no non-header records")

            if vpath.exists() and is_vcfgz:
                first_data = first_non_comment_line_gz(str(vpath))
                if first_data:
                    rpt.ok("variants_vcfgz_has_records", "VCF.GZ has at least one non-header line")
                else:
                    rpt.warning("variants_vcfgz_has_records", "VCF.GZ has no non-header records")

            if vpath.exists() and (is_csv or is_tsv):
                delim = "," if is_csv else "\t"
                try:
                    with open(vpath, "r", encoding="utf-8", errors="replace") as fh:
                        reader = csv.DictReader(fh, delimiter=delim)
                        header = reader.fieldnames or []
                    mapping = suggest_variant_table_mapping(header)
                    if mapping:
                        mapped_desc = ", ".join(f"{k}<-{v}" for k, v in sorted(mapping.items()))
                        rpt.ok("variants_table_column_mapping", mapped_desc)
                    else:
                        rpt.warning(
                            "variants_table_column_mapping",
                            "No known column aliases detected in raw table header",
                        )

                    need = {"CHROM", "POS", "REF", "ALT"}
                    missing = sorted(need - set(mapping))
                    if not missing:
                        rpt.ok("variants_table_columns", "Required columns present/mappable")
                    elif missing == ["REF"]:
                        if fasta.exists() and fai.exists():
                            rpt.warning(
                                "variants_table_columns",
                                "REF column missing, but importer can derive REF from reference FASTA for SNP rows",
                            )
                        else:
                            rpt.error(
                                "variants_table_columns",
                                "REF column missing and reference FASTA/FAI not available for derivation",
                            )
                    else:
                        rpt.error(
                            "variants_table_columns",
                            (
                                "Missing required columns after alias mapping: "
                                + ", ".join(missing)
                                + ". Required canonical fields: CHROM, POS, REF, ALT"
                            ),
                        )
                except OSError as exc:
                    rpt.error("variants_table_read", f"Could not read table: {exc}")

            if vpath.exists() and (is_vcf or is_vcfgz):
                detected = detect_vcf_build(str(vpath), is_vcfgz)
                expected = normalize_build_name(args.expected_build)
                rpt.ok("variants_detected_build", detected)
                if expected and detected not in {"unknown", "mixed_or_ambiguous"}:
                    if expected == detected:
                        rpt.ok(
                            "variants_reference_build_match",
                            f"Detected {detected}; expected {expected}",
                        )
                    else:
                        msg = f"Detected {detected}; expected {expected}"
                        if args.enforce_build_match:
                            rpt.error("variants_reference_build_match", msg)
                        else:
                            rpt.warning("variants_reference_build_match", msg)
                elif expected:
                    rpt.warning(
                        "variants_reference_build_match",
                        (
                            f"Expected build {expected}, but VCF build could not be confidently "
                            f"detected ({detected})"
                        ),
                    )

    # PRS weights
    if args.prs_weights:
        pw = Path(args.prs_weights)
        if pw.exists() and pw.is_file():
            rpt.ok("prs_weights_exists", str(pw))
            low = str(pw).lower()
            try:
                if low.endswith(".gz"):
                    fh = gzip.open(pw, "rt", encoding="utf-8", errors="replace")
                else:
                    fh = open(pw, "r", encoding="utf-8", errors="replace")
                with fh:
                    header = None
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        header = line.rstrip("\n").split("\t")
                        break
                if header:
                    hset = set(header)
                    if {"ID", "A1", "BETA"}.issubset(hset):
                        rpt.ok("prs_weights_columns", "Simple PRS columns detected (ID/A1/BETA)")
                    elif {"effect_allele", "effect_weight"}.issubset(hset):
                        rpt.ok("prs_weights_columns", "PGS Catalog columns detected")
                    else:
                        rpt.warning(
                            "prs_weights_columns",
                            "Weights header does not match common simple/PGS formats",
                        )
                else:
                    rpt.error("prs_weights_columns", "Could not read weights header")
            except OSError as exc:
                rpt.error("prs_weights_read", f"Could not read PRS weights: {exc}")
        else:
            rpt.error("prs_weights_exists", f"Missing PRS weights: {pw}")
    else:
        rpt.warning("prs_weights_exists", "No PRS weights configured; PRS outputs will be skipped")

    # Gene sets
    if args.gene_sets:
        gs = Path(args.gene_sets)
        if gs.exists() and gs.is_file():
            rpt.ok("gene_sets_exists", str(gs))
            try:
                valid_line = False
                with open(gs, "r", encoding="utf-8", errors="replace") as fh:
                    for i, line in enumerate(fh):
                        if i > 5000:
                            break
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 3:
                            valid_line = True
                            break
                if valid_line:
                    rpt.ok("gene_sets_format", "GMT-style line with >=3 tab-separated fields detected")
                else:
                    rpt.error("gene_sets_format", "No valid GMT-style lines detected")
            except OSError as exc:
                rpt.error("gene_sets_read", f"Could not read gene sets file: {exc}")
        else:
            rpt.error("gene_sets_exists", f"Missing gene sets file: {gs}")
    else:
        rpt.warning("gene_sets_exists", "No gene set file configured; enrichment steps will fail")

    # Annotation settings
    method = str(args.annotation_method or "").strip().lower()
    if method not in {"overlap", "snpeff"}:
        rpt.error("annotation_method", f"Unsupported annotation.method: {args.annotation_method}")
    elif method == "snpeff":
        if args.snpeff_database.strip():
            rpt.ok("snpeff_database", f"snpEff database configured: {args.snpeff_database.strip()}")
        else:
            rpt.error(
                "snpeff_database",
                "annotation.method=snpeff requires annotation.snpeff_database",
            )
    else:
        rpt.ok("annotation_method", "Using overlap-based annotation")

    rpt.write(args.out_report)

    if rpt.errors > 0:
        raise SystemExit(1)
    if args.fail_on_warning and rpt.warnings > 0:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
