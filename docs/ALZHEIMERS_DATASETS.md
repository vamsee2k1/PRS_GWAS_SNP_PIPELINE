# Alzheimer Datasets Added

## 1) Real disease-annotated VCF source (downloaded)

- Source: ClinVar GRCh38 VCF
- URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
- Local file: `data/variants/alzheimers/clinvar.vcf.gz`

## 2) Alzheimer-related subset used for pipeline run

- Local file: `data/variants/alzheimers/clinvar_alzheimers_sites.vcf`
- Built by filtering ClinVar entries where annotation contains "Alzheimer".
- This is real disease-associated variant data (not synthetic), but it is site-level ClinVar data.

## Important note for PRS

- ClinVar site VCF does not contain per-sample genotype columns (`GT`/`DS`).
- Therefore, this dataset is suitable for disease-variant processing/filtering, not individual-level PRS calculation.
- For real Alzheimer PRS per sample, use `config/config.alzheimers_prs.yaml` with a real cohort genotype/imputation VCF in GRCh38.
