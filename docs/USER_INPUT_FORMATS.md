# User Upload Formats

This file defines what users should upload for each pipeline mode.

## 1. Full DNA/RNA Pipeline (Recommended)

### Required

1. Raw reads in FASTQ format:
- `.fastq.gz` (preferred)
- paired-end: `R1` and `R2`
- long-read single-end is supported

2. Sample sheet (`TSV`): `config/samples.tsv`
Required columns:
- `sample`
- `fastq_1`
- `fastq_2` (set to blank/NA for single-end)
- `condition` (for RNA differential expression)

3. Metadata (`TSV`): `config/metadata.tsv`
Required columns:
- `sample`
- `condition`
Optional columns:
- `batch`, `sex`, `age`, covariates

4. Reference files:
- Genome FASTA (`.fa` / `.fasta`)
- Annotation GTF (`.gtf`) for RNA
- STAR index directory (RNA short-read mode)

### Optional

- External variant input (`VCF/VCF.GZ/CSV/TSV`) for PRS override
- PRS model file:
  - Preferred: official PGS Catalog scoring file (`.txt` or `.txt.gz`)
  - Also supported: simple `TSV` with columns below
- Simple PRS TSV columns:
  - `ID` (variant ID, e.g. `1:55550:A:G`)
  - `A1` (effect allele)
  - `BETA` (effect size)
- Gene set file (`GMT`) for enrichment

## 2. Variant-Only Mode

Set `run.mode: variant_only` in `config/config.yaml`.
Set `run.variant_data_mode` to either:
- `vcf_interpretation` for genotype/callset VCF
- `gwas_summary` for association-style inputs
Shorthand is also supported by setting `run.mode` directly to `vcf_interpretation` or `gwas_summary`.

### Required

1. Variants input (`paths.variants_input`) in one of:
- `.vcf`
- `.vcf.gz`
- `.csv`
- `.tsv`

2. PRS model file for scoring:
- PGS Catalog file (`.txt`/`.txt.gz`) or
- Simple `TSV` (`ID`, `A1`, `BETA`)

### CSV Variant Format

For CSV input, required columns are:
- `CHROM`
- `POS`
- `REF`
- `ALT`

Optional columns:
- `ID`, `QUAL`, `FILTER`, `INFO`, `GT`, `SAMPLE`
- Optional GWAS columns supported during CSV/TSV to VCF conversion: `P`, `PVALUE`, `PVAL`, `BETA`, `EFFECT`, `OR`, `EAF`, `AF`, `SE`
- For association-style variant plots (`volcano`/`manhattan`), include p-value/effect in INFO when available
  using common keys such as `P`/`PVALUE` and `BETA` or `OR` (OR is log-transformed internally).
- If p-value/effect are not present, the pipeline now uses QUAL/AF-derived proxy ranking so plots remain informative.
- For ClinVar-style clinical VCFs, `CLNSIG` is also used as a proxy signal for prioritization plots.

Template example: `data/variants/input_template.csv`

## 3. Reference Genome Choice

- Preferred: **GRCh38** (latest widely used human reference)
- Legacy compatibility: **hg19/GRCh37**

Use the same genome build across:
- alignment reference
- annotations
- variant resources
- PRS weights

## 4. Preflight Validation (Recommended)

Before running a full workflow, run:

```bash
snakemake --use-conda --cores 1 --until preflight_resources
```

This validates key resources and inputs (FASTA/GTF, variant file format, PRS model, GMT file, and build consistency).

Config controls:
- `validation.enforce_build_match: true` (default) fails fast on detected VCF vs reference build mismatch.
- `validation.fail_on_warning: true` turns warnings into hard failures for stricter production gates.
