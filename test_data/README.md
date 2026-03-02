# Test Data Pack (Real Public Dataset-Derived)

This folder contains small, runnable test files derived from public datasets so users can quickly validate the 3 pipeline modes.

## 1) GWAS Summary Test Data

- Dataset label: `NIAGADS ADVP Alzheimer GWAS (hg38-style table)`
- Folder: `test_data/niagads_advp_ad_gwas_hg38/`
- File: `niagads_advp_ad_variant_records_hg38.tsv`
- Config: `config/test_data/test_niagads_advp_ad_gwas_hg38.yaml`

```yaml
output_dir: results_test_data_niagads_advp_ad_gwas_hg38

run:
  mode: variant_only
  variant_data_mode: gwas_summary
  assay: dna
  read_type: short
  threads: 4

reference:
  genome_build: GRCh38

paths:
  variants_input: test_data/niagads_advp_ad_gwas_hg38/niagads_advp_ad_variant_records_hg38.tsv
  prs_weights: ""
  gene_sets: resources/enrichment/real/Reactome_2022.gmt

validation:
  enforce_build_match: true
  fail_on_warning: false
```

Run:

```bash
./run_pipeline.sh --use-conda --cores 8 --configfile config/test_data/test_niagads_advp_ad_gwas_hg38.yaml
```

Reference:
- https://advp.niagads.org/

## 2) Variant VCF Test Data

- Dataset label: `1000 Genomes chr19 AD high-PRS top-40 subset (GRCh38)`
- Folder: `test_data/onekg_chr19_ad_high_prs_top40_grch38/`
- Files: `onekg_chr19_ad_high_prs_top40_grch38.vcf.gz` and `.tbi`
- Config: `config/test_data/test_onekg_chr19_ad_high_prs_top40_grch38.yaml`

```yaml
output_dir: results_test_data_onekg_chr19_ad_high_prs_top40_grch38

run:
  mode: variant_only
  variant_data_mode: vcf_interpretation
  assay: dna
  read_type: short
  threads: 4

reference:
  genome_build: GRCh38

paths:
  variants_input: test_data/onekg_chr19_ad_high_prs_top40_grch38/onekg_chr19_ad_high_prs_top40_grch38.vcf.gz
  prs_weights: resources/prs/models/PGS002280.txt.gz
  gene_sets: resources/enrichment/real/Reactome_2022.gmt

validation:
  enforce_build_match: true
  fail_on_warning: false
```

Run:

```bash
./run_pipeline.sh --use-conda --cores 8 --configfile config/test_data/test_onekg_chr19_ad_high_prs_top40_grch38.yaml
```

References:
- https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- https://www.pgscatalog.org/score/PGS002280/

## 3) Full FASTQ Test Data

- Dataset label: `GIAB/NIST7035 FASTQ subset (GRCh38-aligned workflow test)`
- Folder: `test_data/giab_nist7035_fastq_subset_grch38/`
- Files:
  - `giab_nist7035_subset_R1.fastq.gz`
  - `giab_nist7035_subset_R2.fastq.gz`
  - `samples.tsv`
  - `metadata.tsv`
- Config: `config/test_data/test_giab_nist7035_fastq_subset_grch38.yaml`

```yaml
output_dir: results_test_data_giab_nist7035_fastq_subset_grch38
samplesheet: test_data/giab_nist7035_fastq_subset_grch38/samples.tsv
metadata: test_data/giab_nist7035_fastq_subset_grch38/metadata.tsv

run:
  mode: full
  assay: dna
  read_type: short
  threads: 4
  use_external_variants_for_prs: false
  enable_depth: false

reference:
  genome_build: GRCh38

alignment:
  dna_short_aligner: auto

paths:
  variants_input: ""
  prs_weights: ""
  gene_sets: resources/enrichment/real/Reactome_2022.gmt

validation:
  enforce_build_match: true
  fail_on_warning: false
```

Run:

```bash
./run_pipeline.sh --use-conda --cores 8 --configfile config/test_data/test_giab_nist7035_fastq_subset_grch38.yaml
```

Reference:
- https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/

## Notes

- These are lightweight subsets for reproducible testing, not full-cohort production analyses.
- Main mode templates for custom user data remain in `README.md` Section 4.3.
