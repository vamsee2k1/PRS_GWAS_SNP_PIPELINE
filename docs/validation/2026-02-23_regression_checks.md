# Validation Summary (2026-02-23)

## Scope

Validated two recent changes:

1. Preflight alignment-index checks for DNA short-read mode (`bwa-mem2` / `bwa` / `auto` fallback)
2. Regression tests and real-file validation for `workflow/scripts/csv_to_vcf.py` (GWAS/TSV -> VCF conversion)

## Commands Run

### Unit tests

```bash
python -m unittest -v tests.test_csv_to_vcf
```

### Real example-file conversion test (ADVP GWAS TSV)

```bash
PATH=.snakemake/conda/719d46ce0a478f43e34991fd6377429c_/bin:$PATH \
  .snakemake/conda/719d46ce0a478f43e34991fd6377429c_/bin/python \
  workflow/scripts/csv_to_vcf.py \
  --input INPUT_TEST_FILES/advp.variant.records.hg38.tsv \
  --output /tmp/advp_example_test.vcf \
  --ref-fasta resources/reference/GRCh38.fa \
  --report /tmp/advp_example_test.report.tsv
```

### Preflight index check (GIAB-style full DNA config paths)

```bash
python3 workflow/scripts/preflight_validate_resources.py \
  --mode full --assay dna --read-type short \
  --samplesheet config/samples.giab_hg002.tsv \
  --metadata config/metadata.giab_hg002.tsv \
  --reference-fasta resources/reference/GRCh38.fa \
  --reference-gtf resources/reference/Homo_sapiens.GRCh38.110.gtf \
  --reference-star-index resources/reference/star_grch38_index \
  --dna-short-aligner auto \
  --expected-build GRCh38 \
  --prs-weights resources/prs/models/PGS002280.txt.gz \
  --prs-weights-format pgs_catalog \
  --gene-sets resources/enrichment/real/Reactome_2022.gmt \
  --annotation-method overlap \
  --enforce-build-match \
  --out-report /tmp/preflight_giab_indexcheck.tsv
```

## Results

### A) `csv_to_vcf.py` unit tests

- Status: **PASS** (`5/5`)
- Coverage includes:
  - ADVP alias mapping (`#dbSNP_hg38_chr`, `dbSNP_hg38_position`, `nonref_allele`, `Top SNP`, `P-value`)
  - REF derivation + contig normalization (`19` -> `chr19`)
  - invalid/non-SNP ALT handling (`NR`)
  - row validation drop behavior
  - missing-REF warning path

### B) Real ADVP example file conversion

Input file:
- `INPUT_TEST_FILES/advp.variant.records.hg38.tsv`

Observed conversion summary (from generated report):
- Input rows: `6346`
- `REF` derived from reference FASTA: `2910` rows
- `REF` unresolved (non-SNP/ambiguous/etc.): `3436` rows
- Invalid rows dropped during VCF normalization: `4365`
- Valid rows output to VCF: `1981`

Interpretation:
- This confirms the TSV/GWAS alias mapping and REF derivation flow works on the real ADVP-style example file.
- Unresolved rows are expected because ADVP includes values like `NR` that are not valid nucleotide ALT alleles.

### C) Preflight aligner-index detection (new change)

Observed preflight checks (`alignment.dna_short_aligner=auto`):
- `dna_short_aligner`: `OK (auto)`
- `bwa_mem2_index_ready`: `WARNING` (missing `.bwt.2bit.64`)
- `bwa_index_ready`: `OK` (classic BWA index found)
- `dna_short_alignment_index`: `WARNING` (`auto` will fall back to classic BWA)

Interpretation:
- Preflight now catches the previously confusing full-mode alignment failure condition *before* alignment starts and explains that `auto` will use the BWA fallback when available.

## Attached Validation Artifacts

- `advp_example_csv_to_vcf_report.tsv`
- `preflight_giab_indexcheck.tsv`
- `unittest_test_csv_to_vcf.txt`
