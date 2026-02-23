#!/usr/bin/env python3
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def parse_args():
    parser = argparse.ArgumentParser(description="Compute PRS from VCF and PRS model weights")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ")
    parser.add_argument("--weights", required=True, help="PRS model file")
    parser.add_argument(
        "--weights-format",
        default="auto",
        choices=["auto", "pgs_catalog", "simple_tsv"],
        help="Weight file format",
    )
    parser.add_argument("--out-sscore", required=True, help="Output PRS score table")
    parser.add_argument("--out-qc", required=True, help="Output PRS QC metrics table")
    parser.add_argument(
        "--prefer-ds",
        action="store_true",
        help="Prefer dosage field (DS) over genotype (GT) when available",
    )
    parser.add_argument(
        "--expected-build",
        default="",
        help="Expected genome build for the input VCF (e.g., GRCh38, GRCh37)",
    )
    parser.add_argument(
        "--allow-ambiguous-strand",
        action="store_true",
        help="Allow A/T and C/G SNPs during strand harmonization (default: skip ambiguous SNPs)",
    )
    return parser.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def normalize_chr(chrom):
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def complement_base(base):
    b = base.upper()
    if b in COMPLEMENT:
        return COMPLEMENT[b]
    return b


def detect_weights_format(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            hset = set(header)
            if {"chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"}.issubset(hset):
                return "pgs_catalog"
            if {"ID", "A1", "BETA"}.issubset(hset):
                return "simple_tsv"
            break
    raise ValueError("Could not detect weights format. Supported: pgs_catalog, simple_tsv")


def load_simple_weights(path):
    by_id = {}
    by_pos = defaultdict(list)

    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"ID", "A1", "BETA"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError("Simple weights must contain columns: ID, A1, BETA")

        for row in reader:
            var_id = row["ID"].strip()
            if not var_id:
                continue
            try:
                beta = float(row["BETA"])
            except ValueError:
                continue

            effect = row["A1"].strip().upper()
            record = {
                "model_id": var_id,
                "effect_allele": effect,
                "other_allele": None,
                "beta": beta,
                "rsid": None,
                "chr": None,
                "pos": None,
            }

            parts = var_id.split(":")
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
                chrom = normalize_chr(chrom)
                try:
                    ipos = int(pos)
                except ValueError:
                    ipos = None
                if ipos is not None:
                    record["chr"] = chrom
                    record["pos"] = ipos
                    record["other_allele"] = ref.upper() if effect == alt.upper() else alt.upper()
                    by_pos[(chrom, ipos)].append(record)
            elif var_id.lower().startswith("rs"):
                record["rsid"] = var_id

            by_id[var_id] = record

    return by_id, by_pos, "simple_tsv"


def load_pgs_catalog_weights(path):
    by_id = {}
    by_pos = defaultdict(list)
    pgs_meta = {
        "pgs_id": "",
        "pgs_name": "",
        "trait_reported": "",
        "trait_mapped": "",
        "genome_build": "",
        "variants_number": "",
        "citation": "",
    }

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                clean = line.strip().lstrip("#")
                if "=" in clean:
                    k, v = clean.split("=", 1)
                    if k in pgs_meta:
                        pgs_meta[k] = v.strip()
                continue
            header = line.rstrip("\n").split("\t")
            reader = csv.DictReader(fh, delimiter="\t", fieldnames=header)
            break
        else:
            raise ValueError("PGS file appears empty")

        required = {"chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"}
        if not required.issubset(set(header)):
            raise ValueError("PGS Catalog file missing required columns")

        for row in reader:
            try:
                beta = float(row["effect_weight"])
                pos = int(row["chr_position"])
            except (ValueError, TypeError):
                continue

            chrom = normalize_chr(row["chr_name"])
            effect = row["effect_allele"].upper()
            other = row["other_allele"].upper()
            rsid = (row.get("rsID") or "").strip()

            model_id = rsid if rsid else f"{chrom}:{pos}:{other}:{effect}"
            rec = {
                "model_id": model_id,
                "effect_allele": effect,
                "other_allele": other,
                "beta": beta,
                "rsid": rsid if rsid else None,
                "chr": chrom,
                "pos": pos,
            }

            by_id[model_id] = rec
            by_pos[(chrom, pos)].append(rec)

    return by_id, by_pos, "pgs_catalog", pgs_meta


def gt_alt_dosage(gt):
    if gt in {".", "./.", ".|."}:
        return None
    alleles = gt.replace("|", "/").split("/")
    dose = 0.0
    for a in alleles:
        if a == ".":
            return None
        if a == "1":
            dose += 1.0
        elif a != "0":
            return None
    return dose


def parse_sample_dose(sample_field, fmt_keys, prefer_ds=True):
    values = sample_field.split(":")
    fdict = {k: values[i] if i < len(values) else "" for i, k in enumerate(fmt_keys)}

    if prefer_ds and "DS" in fdict and fdict["DS"] not in {"", "."}:
        ds = fdict["DS"].split(",")[0]
        try:
            return float(ds), "DS"
        except ValueError:
            pass

    if "GT" in fdict:
        gt = fdict["GT"]
        dose = gt_alt_dosage(gt)
        if dose is not None:
            return dose, "GT"

    if "DS" in fdict and fdict["DS"] not in {"", "."}:
        ds = fdict["DS"].split(",")[0]
        try:
            return float(ds), "DS"
        except ValueError:
            return None, "NONE"

    return None, "NONE"


def alleles_match_with_optional_strand(model_eff, model_other, ref, alt):
    candidates = [
        (model_eff, model_other, False),
        (complement_base(model_eff), complement_base(model_other), True),
    ]

    for eff, other, strand_flip in candidates:
        if {eff, other} == {ref, alt}:
            return eff, other, strand_flip
    return None, None, None


def is_ambiguous_snp(ref, alt):
    if len(ref) != 1 or len(alt) != 1:
        return False
    pair = {ref.upper(), alt.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def choose_model_record(model_candidates, ref, alt, vcf_ids):
    ref = ref.upper()
    alt = alt.upper()

    for rid in vcf_ids:
        for rec in model_candidates:
            if rec.get("rsid") and rec["rsid"] == rid:
                eff, other, strand_flip = alleles_match_with_optional_strand(
                    rec["effect_allele"], rec["other_allele"] or ref, ref, alt
                )
                if eff is not None:
                    return rec, eff, strand_flip, "rsid"

    for rec in model_candidates:
        other = rec["other_allele"]
        if other is None:
            if rec["effect_allele"] in {ref, alt}:
                return rec, rec["effect_allele"], False, "position"
            continue

        eff, _other, strand_flip = alleles_match_with_optional_strand(
            rec["effect_allele"], other, ref, alt
        )
        if eff is not None:
            return rec, eff, strand_flip, "position"

    return None, None, None, ""


def effect_dose_from_alt_dose(effect_allele, ref, alt, alt_dose):
    if effect_allele == alt:
        return alt_dose
    if effect_allele == ref:
        return 2.0 - alt_dose
    return None


def main():
    args = parse_args()

    weights_format = args.weights_format
    if weights_format == "auto":
        weights_format = detect_weights_format(args.weights)

    if weights_format == "pgs_catalog":
        by_id, by_pos, fmt_name, pgs_meta = load_pgs_catalog_weights(args.weights)
    else:
        by_id, by_pos, fmt_name = load_simple_weights(args.weights)
        pgs_meta = {
            "pgs_id": "",
            "pgs_name": "",
            "trait_reported": "",
            "trait_mapped": "",
            "genome_build": "",
            "variants_number": str(len(by_id)),
            "citation": "",
        }

    model_size = len(by_id)
    matched_model_ids = set()

    samples = []
    scores = {}
    sample_variant_hits = defaultdict(int)
    ds_calls = 0
    gt_calls = 0
    skipped_multiallelic = 0
    skipped_nomatch = 0
    skipped_unparsable = 0
    skipped_ambiguous = 0
    strand_flipped_matches = 0
    matched_by_rsid = 0
    matched_by_position = 0

    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                scores = {s: 0.0 for s in samples}
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos_str, rid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            chrom = normalize_chr(chrom)
            ref = ref.upper()
            alt = alt.upper()

            if "," in alt:
                skipped_multiallelic += 1
                continue

            if not args.allow_ambiguous_strand and is_ambiguous_snp(ref, alt):
                skipped_ambiguous += 1
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                skipped_unparsable += 1
                continue

            vcf_ids = [x for x in rid.split(";") if x and x != "."]
            model_candidates = by_pos.get((chrom, pos), [])

            for vid in vcf_ids:
                if vid in by_id and by_id[vid] not in model_candidates:
                    model_candidates.append(by_id[vid])

            rec, effect_allele, strand_flip, match_mode = choose_model_record(model_candidates, ref, alt, vcf_ids)
            if rec is None:
                skipped_nomatch += 1
                continue

            matched_model_ids.add(rec["model_id"])
            if match_mode == "rsid":
                matched_by_rsid += 1
            else:
                matched_by_position += 1
            if strand_flip:
                strand_flipped_matches += 1

            fmt_keys = parts[8].split(":")
            for sample, sample_field in zip(samples, parts[9:]):
                alt_dose, source = parse_sample_dose(sample_field, fmt_keys, prefer_ds=args.prefer_ds)
                if alt_dose is None:
                    continue
                if source == "DS":
                    ds_calls += 1
                elif source == "GT":
                    gt_calls += 1

                eff_dose = effect_dose_from_alt_dose(effect_allele, ref, alt, alt_dose)
                if eff_dose is None:
                    continue

                scores[sample] += eff_dose * rec["beta"]
                sample_variant_hits[sample] += 1

    Path(args.out_sscore).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_sscore, "w", encoding="utf-8") as out:
        out.write("IID\tSCORE1_SUM\tN_MATCHED_VARIANTS\tMODEL_VARIANTS\tMATCH_RATE\n")
        for sample in samples:
            n_hit = sample_variant_hits.get(sample, 0)
            rate = (n_hit / model_size * 100.0) if model_size else 0.0
            out.write(f"{sample}\t{scores.get(sample, 0.0):.8f}\t{n_hit}\t{model_size}\t{rate:.3f}\n")

    with open(args.out_qc, "w", encoding="utf-8") as qc:
        expected_build = str(args.expected_build).strip()
        model_build = str(pgs_meta.get("genome_build", "")).strip()
        if expected_build and model_build:
            build_match = str(expected_build.lower() == model_build.lower())
        else:
            build_match = "unknown"

        qc.write("metric\tvalue\n")
        qc.write(f"weights_format\t{fmt_name}\n")
        qc.write(f"model_variants\t{model_size}\n")
        qc.write(f"matched_model_variants\t{len(matched_model_ids)}\n")
        qc.write(f"matched_by_rsid\t{matched_by_rsid}\n")
        qc.write(f"matched_by_position\t{matched_by_position}\n")
        match_rate = (len(matched_model_ids) / model_size * 100.0) if model_size else 0.0
        qc.write(f"matched_model_rate_pct\t{match_rate:.3f}\n")
        qc.write(f"parsed_ds_calls\t{ds_calls}\n")
        qc.write(f"parsed_gt_calls\t{gt_calls}\n")
        qc.write(f"strand_flipped_matches\t{strand_flipped_matches}\n")
        qc.write(f"skipped_multiallelic\t{skipped_multiallelic}\n")
        qc.write(f"skipped_ambiguous_strand\t{skipped_ambiguous}\n")
        qc.write(f"skipped_nomatch\t{skipped_nomatch}\n")
        qc.write(f"skipped_unparsable\t{skipped_unparsable}\n")
        qc.write(f"expected_genome_build\t{expected_build}\n")
        qc.write(f"build_match\t{build_match}\n")
        qc.write(f"pgs_id\t{pgs_meta.get('pgs_id', '')}\n")
        qc.write(f"pgs_name\t{pgs_meta.get('pgs_name', '')}\n")
        qc.write(f"trait_reported\t{pgs_meta.get('trait_reported', '')}\n")
        qc.write(f"trait_mapped\t{pgs_meta.get('trait_mapped', '')}\n")
        qc.write(f"model_genome_build\t{pgs_meta.get('genome_build', '')}\n")
        qc.write(f"citation\t{pgs_meta.get('citation', '')}\n")


if __name__ == "__main__":
    main()
