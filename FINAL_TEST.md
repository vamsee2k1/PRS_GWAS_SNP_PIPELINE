# FINAL TEST

Date: 2026-02-27
Project: Reproducible Bioinformatics Pipeline (3-mode publish readiness rerun)
Status in docs: canonical real-dataset examples replacing older toy/example run references.

## Scope
Re-tested the pipeline end-to-end in all 3 target modes using fresh output directories:
1. AD PRS (variant-only, vcf_interpretation) with a different high-PRS-focused dataset.
2. GWAS mode (variant-only, gwas_summary) with a different GWAS input.
3. Full mode (DNA short-read FASTQ) with a different full-mode input.

All runs completed successfully (100%).

## Configs And Inputs Used

### 1) AD PRS final test
- Config: `config/final_tests/final_test_ad_high_prs.yaml`
- Output dir: `results_final_test_ad_high_prs_20260227`
- Input variants: `data/variants/alzheimers/final_tests/1kg_chr19_high_prs_top40.vcf.gz`
- PRS weights: `resources/prs/models/PGS002280.txt.gz`
- Dataset preparation: top 40 PRS samples were selected and subset from `data/variants/alzheimers/1kg_chr19.GRCh38.vcf.gz` using:
  - sample list: `config/final_tests/high_prs_top40_samples.txt`

### 2) GWAS final test
- Config: `config/final_tests/final_test_gwas_advp.yaml`
- Output dir: `results_final_test_gwas_advp_20260227`
- Input GWAS file: `INPUT_TEST_FILES/advp.variant.records.hg38.tsv`
- PRS intentionally disabled for GWAS summary mode in this test (`prs_weights: ""`).

### 3) Full mode final test
- Config: `config/final_tests/final_test_full_giab.yaml`
- Output dir: `results_final_test_full_giab_20260227`
- Samplesheet: `config/samples.giab_hg002.tsv`
- Metadata: `config/metadata.giab_hg002.tsv`
- Input FASTQ: `data/giab_hg002/HG002_sub_R1.fastq.gz`, `data/giab_hg002/HG002_sub_R2.fastq.gz`
- PRS enabled with `resources/prs/models/PGS002280.txt.gz`

## Execution Status

### AD PRS run
- Preflight: `0 ERROR`, `0 WARNING`
- Run manifest: present
- Status: PASS

Key metrics:
- `variant_qc_metrics.tsv`:
  - total_variants: `1,625,698`
  - sample_count: `40`
  - ts/tv: `2.35086`
- `prs_qc.tsv`:
  - model_variants: `83`
  - matched_model_variants: `3`
  - matched_model_rate_pct: `3.614`
  - build_match: `True`
- `prs_scores.tsv` (high end observed):
  - top PRS in this set: `0.3782`

### GWAS run
- Preflight: `0 ERROR`, `2 WARNING`
- Run manifest: present
- Status: PASS

Preflight warnings:
- REF column missing in raw table (handled by reference-based derivation).
- PRS weights not configured (expected for this GWAS test).

Key metrics:
- `external.input.field_report.tsv`:
  - input_rows: `6346`
  - output_rows: `1981`
  - dropped invalid rows: `4365`
- `variant_qc_metrics.tsv`:
  - total_variants: `1981`
  - ts/tv: `2.22114`
- `variant_association_hits.tsv`:
  - reported-signal hits present (APOE/TOMM40-region signals observed in top rows).

### Full mode run
- Preflight: `0 ERROR`, `2 WARNING`
- Run manifest: present
- Status: PASS

Preflight warnings:
- Missing `bwa-mem2` index file for current reference path.
- Auto-aligner fallback to classic `bwa` (expected behavior, run completed successfully).

Key metrics:
- `qc/qc_before_after.tsv`:
  - reads_retention_pct: `94.167`
  - bases_retention_pct: `93.92684527027026`
- `depth/depth_summary.tsv`:
  - mean_depth: `0.08902475635873401`
  - max_depth: `2224`
- `variants/variant_qc_metrics.tsv`:
  - total_variants: `13,520`
  - snp_count: `12,556`
  - indel_count: `905`
  - ts/tv: `1.53657`
- `prs/prs_qc.tsv`:
  - matched_model_variants: `0`
  - matched_model_rate_pct: `0.000`
  - build_match: `True`
- `prs/prs_scores.tsv` note correctly reports non-interpretable PRS due no callable AD loci in this sample's called set.

## Final Verdict
The pipeline passed full rerun in all three requested modes with fresh outputs and manifests:
- `results_final_test_ad_high_prs_20260227`
- `results_final_test_gwas_advp_20260227`
- `results_final_test_full_giab_20260227`

No hard failures were observed. Preflight warnings were expected/data-context-dependent and did not block execution.
