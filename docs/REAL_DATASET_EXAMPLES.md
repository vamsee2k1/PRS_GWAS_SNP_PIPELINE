# Real Dataset Examples (Canonical)

This document replaces older toy/example-run references with real dataset reruns completed on 2026-02-27.

## Canonical Example Runs

1. AD PRS (`variant_only`, `vcf_interpretation`)
- Config: `config/final_tests/final_test_ad_high_prs.yaml`
- Output: `results_final_test_ad_high_prs_20260227`
- Input: `data/variants/alzheimers/final_tests/1kg_chr19_high_prs_top40.vcf.gz`
- Note: this VCF is a derived subset of a real 1000G chr19 VCF, selecting top-PRS samples from prior AD PRS output.

2. GWAS (`variant_only`, `gwas_summary`)
- Config: `config/final_tests/final_test_gwas_advp.yaml`
- Output: `results_final_test_gwas_advp_20260227`
- Input: `INPUT_TEST_FILES/advp.variant.records.hg38.tsv` (real GWAS summary-style table)

3. Full mode (DNA short-read FASTQ)
- Config: `config/final_tests/final_test_full_giab.yaml`
- Output: `results_final_test_full_giab_20260227`
- Inputs: `data/giab_hg002/HG002_sub_R1.fastq.gz`, `data/giab_hg002/HG002_sub_R2.fastq.gz`

Run-level summary and pass status are documented in `FINAL_TEST.md`.

## How The Plots Are Generated

## QC and Depth Plots

- `qc_before_after.png`
  - Source: `workflow/scripts/qc_compare.py`
  - Inputs: `trimmed/*.fastp.json`
  - Meaning: visual read/base retention and Q20/Q30 improvement after trimming.

- `multiqc_report.html`
  - Source: `multiqc` over raw/trimmed FastQC and fastp outputs.
  - Meaning: consolidated quality dashboard across preprocessing artifacts.

- `depth_distribution.png`
  - Source: `workflow/scripts/depth_summary.py`
  - Inputs: mosdepth summaries/distributions + optional `samtools depth`.
  - Meaning: coverage profile; for WGS-like data, low median depth indicates sparse coverage in subset/test runs.

## Variant QC Plots

- `variant_type_counts.png`, `variant_qual_distribution.png`, `variant_dp_distribution.png`,
  `variant_gq_distribution.png`, `variant_het_allele_balance.png`, `variant_missingness_by_sample.png`
  - Source: `workflow/scripts/variant_qc_stats.py`
  - Meaning:
    - variant class composition (SNP/indel/etc.),
    - confidence/coverage distribution (QUAL/DP/GQ),
    - genotype-balance quality for heterozygous calls,
    - sample-level missingness burden.

## Variant Association/Interpretation Plots

- `variant_volcano.png`, `variant_manhattan.png`, `variant_qq.png`, `variant_heatmap_top.png`,
  `variant_functional_classes.png`, `variant_chromosome_counts.png`
  - Source: `workflow/scripts/variant_visualization.py`
  - Meaning:
    - functional and chromosome plots summarize annotation distribution.
    - volcano/manhattan/qq reflect association signal quality and ranking.
    - heatmap shows top-ranked variant genotype dosage when genotype columns exist,
      otherwise falls back to numeric feature heatmap.

Association signal mode:
- If p-value/effect columns are present and valid: `reported` mode.
- If not present: fallback proxy mode from QUAL/AF (and clinical significance if available).
- Thresholding uses `thresholds.variant_assoc_pvalue` and `thresholds.variant_assoc_abs_effect`.

## Enrichment Plots

- `variant_enrichment_dotplot.png` and `variant_enrichment.tsv`
  - Source: `workflow/scripts/variant_enrichment.R` (clusterProfiler enricher)
  - Input genes: `variant_genes.tsv`
  - Meaning: pathway-level overrepresentation from genes mapped by variant-feature overlap.

## PRS Plots

- `prs_distribution.png`, `prs_scores.tsv`, `prs_qc.tsv`
  - Sources: `workflow/scripts/prs_from_vcf.py`, `workflow/scripts/prs_report.py`
  - Meaning:
    - distribution: cohort-relative PRS spread,
    - scores table: per-sample PRS, percentile, z-score,
    - QC table: model coverage, harmonization counts, build check, skipped categories.

## Biological Significance Notes By Example

## AD PRS Example (`results_final_test_ad_high_prs_20260227`)

- Biological interpretation:
  - This run demonstrates relative AD genetic risk stratification within a selected high-PRS-enriched subset.
  - Higher PRS values indicate higher burden under the model weighting, not diagnosis.
- Critical caveat:
  - `matched_model_variants` is low (3 of 83), so absolute interpretation is limited; use as demo of scoring workflow and harmonization behavior.

## GWAS Example (`results_final_test_gwas_advp_20260227`)

- Biological interpretation:
  - Top associations include known AD loci (for example APOE/TOMM40-region hits in top rows of `variant_association_hits.tsv`).
  - Enrichment output reflects biology from mapped genes among retained variants.
- Critical caveat:
  - Row-level cleaning is expected in heterogeneous public summary tables; review `external.input.field_report.tsv` for dropped rows and REF derivation outcomes.

## Full Mode Example (`results_final_test_full_giab_20260227`)

- Biological interpretation:
  - Primarily an end-to-end technical validation on benchmark-style HG002 subset data.
  - Variant QC, depth, and enrichment outputs validate pipeline mechanics and reproducibility.
- Critical caveat:
  - This is not a disease cohort PRS study; AD PRS output here reports zero matched loci and is correctly marked non-interpretable.

## Reproducibility

For each canonical run, preserve:
- `docs/preflight_checks.tsv`
- `docs/run_manifest.txt`
- key result tables and plots in `variants/`, `prs/`, `qc/`, and `depth/` as applicable.

These three run folders are the recommended examples to reference in publication/readme materials:
- `results_final_test_ad_high_prs_20260227`
- `results_final_test_gwas_advp_20260227`
- `results_final_test_full_giab_20260227`
