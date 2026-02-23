#!/usr/bin/env python3
import argparse
import re


KEEP_FEATURES = {"gene", "exon", "CDS", "UTR", "three_prime_utr", "five_prime_utr"}
NC_CHROM_RE = re.compile(r"^NC_0*([0-9]+)\\.")


def parse_args():
    parser = argparse.ArgumentParser(description="Convert GTF to BED features used for variant annotation")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out-bed", required=True)
    return parser.parse_args()


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


def parse_attributes(attr_text):
    attrs = {}
    for token in attr_text.strip().split(";"):
        token = token.strip()
        if not token:
            continue
        if " " not in token:
            attrs[token] = ""
            continue
        key, value = token.split(" ", 1)
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def main():
    args = parse_args()

    with open(args.gtf, "r", encoding="utf-8") as in_fh, open(args.out_bed, "w", encoding="utf-8") as out_fh:
        for line in in_fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            chrom_raw, _src, feature, start_raw, end_raw, _score, _strand, _frame, attrs_raw = parts
            if feature not in KEEP_FEATURES:
                continue

            try:
                start = int(start_raw)
                end = int(end_raw)
            except ValueError:
                continue

            attrs = parse_attributes(attrs_raw)
            gene_id = attrs.get("gene_id", attrs.get("gene", ""))
            gene_name = attrs.get("gene_name", attrs.get("Name", gene_id))
            chrom = normalize_chrom(chrom_raw)
            start0 = max(start - 1, 0)
            out_fh.write(
                f"{chrom}\t{start0}\t{end}\t{feature}\t{gene_id}\t{gene_name}\n"
            )


if __name__ == "__main__":
    main()
