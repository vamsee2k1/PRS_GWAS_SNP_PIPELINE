# Alzheimer Disease PRS (Production Configuration)

This pipeline is configured to use a real Alzheimer disease polygenic score model:

- Model: `PGS002280` (`GRS83_AD`)
- Trait: Alzheimer’s disease
- Genome build: GRCh38
- Source model file: `resources/prs/models/PGS002280.txt.gz`

## PRS Formula

For each sample:

`PRS = Σ (effect_allele_dosage_i × effect_weight_i)`

Where for each variant `i`:
- `effect_weight_i` comes from the PGS model file.
- `effect_allele_dosage_i` is computed from VCF `DS` (preferred) or `GT`.

## What the pipeline does

1. Parses the official PGS Catalog file.
2. Harmonizes model alleles with VCF REF/ALT alleles.
3. Computes per-sample PRS from matched variants.
4. Writes score + coverage metrics:
- `results/prs/plink.sscore`
- `results/prs/prs_scores.tsv`
- `results/prs/prs_qc.tsv`
- `results/prs/prs_distribution.png`

`prs_scores.tsv` now includes an `alzheimers_prs_note` column. If `N_MATCHED_VARIANTS = 0`, the note explicitly states that PRS is shown as `0` due to no callable Alzheimer loci coverage.

`prs_qc.tsv` also reports harmonization diagnostics including:
- `matched_by_rsid`
- `matched_by_position`
- `strand_flipped_matches`
- `skipped_ambiguous_strand`
- `expected_genome_build`
- `build_match`

## Accuracy-critical requirements

1. Use genome build-matched data (GRCh38 for this model).
2. Prefer high-quality imputed genotype VCF with dosage (`DS`) fields.
3. Ensure ancestry match to model development population when possible.
4. Apply strict sample/variant QC before PRS interpretation.
5. Treat PRS as risk stratification, not diagnosis.

## Run command (variant-only)

```bash
snakemake --use-conda --cores 8 --configfile config/config.alzheimers_prs.yaml
```

Set `paths.variants_input` in `config/config.alzheimers_prs.yaml` to your cohort VCF/VCF.GZ.
