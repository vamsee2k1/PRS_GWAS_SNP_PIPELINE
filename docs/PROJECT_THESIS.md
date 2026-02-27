# Project Thesis: Reproducible Multi-Mode Bioinformatics Pipeline

## 1. Problem Statement

Modern bioinformatics workflows are often fragmented across scripts, ad hoc notebooks, and dataset-specific logic. This makes it difficult to:

- reproduce results across machines and users
- validate new input formats safely
- benchmark methods with public reference datasets
- document exact pipeline behavior for publication

This project addresses that gap with a Snakemake-based pipeline that supports three operational modes:

- `full`: raw sequencing reads (`FASTQ`) to alignment, variant calling, annotation, QC, and optional PRS
- `variant_only`: external `VCF/VCF.GZ` to normalization, annotation, QC, enrichment, and optional PRS
- `gwas_summary` (via `variant_only` + `run.variant_data_mode: gwas_summary`): GWAS-style summary TSV/CSV ingestion, harmonization, plotting, and gene/pathway interpretation

## 2. Core Design Goals

### Reproducibility

- Workflow orchestration with Snakemake DAGs
- Rule-scoped conda environments in `envs/`
- Deterministic output tree (`results/...`)
- Run manifest generation (`results/.../docs/run_manifest.txt`)
- Preflight validation before expensive compute

### Practical Input Tolerance

- Accept `VCF`, `VCF.GZ`, `CSV`, `TSV`
- Alias-aware header normalization for raw uploads (including GWAS-style tables)
- Missing-field reporting instead of silent failure
- `REF` allele derivation from reference FASTA where possible

### Publication Readiness

- Explicit citation metadata (`CITATION.cff`)
- Open-source license and reuse notice
- Contribution guidance and code review expectations
- Machine- and human-readable run metadata

## 3. Architecture Overview

The pipeline is intentionally modular. A single `Snakefile` resolves mode and branches into the appropriate subgraph.

### Shared stages

- Config load + mode resolution
- Preflight resource/input validation
- Final variant interpretation outputs (annotation, QC, enrichment)
- Optional PRS branch (if PRS weights configured)
- Run manifest

### Mode-specific branches

#### A. `full` mode

Inputs:
- sample sheet + metadata
- raw FASTQ files
- reference resources

Stages:
- FastQC (raw)
- fastp trimming/filtering
- alignment (`bwa-mem2`, `bwa`, or `minimap2` depending config/mode)
- BAM indexing / duplicate marking (DNA short-read)
- alignment QC (`samtools flagstat`, `samtools stats`)
- depth summaries (`mosdepth`, `samtools depth`)
- variant calling (`bcftools mpileup/call/filter`)
- downstream variant interpretation

#### B. `variant_only` mode

Inputs:
- external variants (`VCF/VCF.GZ/CSV/TSV`)

Stages:
- normalize/import external variants
- optional CSV/TSV -> VCF conversion
- filtering / multiallelic handling
- final VCF canonicalization
- downstream variant interpretation

#### C. `gwas_summary` mode (implemented as `variant_only`)

Inputs:
- GWAS summary-style `TSV/CSV` (e.g., ADVP/NIAGADS-like)

Stages:
- alias mapping (`CHROM`, `POS`, `ALT`, etc.)
- `REF` derivation (reference FASTA)
- diagnostics (`external.input.field_report.tsv`)
- conversion of valid rows to VCF-like representation
- GWAS-style visualizations (Manhattan/QQ/volcano)
- variant-gene mapping + enrichment

## 4. Engineering Decisions and Refinements

### 4.1 Input normalization for real-world TSV/CSV uploads

Raw public datasets often do not match strict VCF-like column schemas. The pipeline now supports:

- common aliases for coordinates/alleles
- row-level rejection with diagnostics for invalid alleles
- partial conversion when some rows are usable

This avoids failing an entire run due to mixed-quality source tables.

### 4.2 DNA short-read aligner fallback for local systems

Some systems can stall on `bwa-mem2 index` for large references (e.g., GRCh38). To improve reliability:

- `alignment.dna_short_aligner: auto` now prefers `bwa-mem2` when its index exists
- automatically falls back to classic `bwa` when classic BWA index files exist
- supports `minimap2_sr` as an explicit override for certain test workflows

This preserves correctness while improving local developer ergonomics.

### 4.3 Depth computation controls for test vs production runs

Whole-genome `samtools depth -a` outputs can be extremely large and slow to summarize. The pipeline now exposes:

- `run.enable_depth` to disable depth outputs for quick validation runs
- `depth.emit_all_positions` to avoid `-a` during tests

This separates production depth reporting from rapid functional testing.

## 5. Validation Strategy (Three-Mode Test Plan)

### Mode 1: GWAS summary

Recommended test:
- ADVP Alzheimer GWAS summary TSV (`advp.variant.records.hg38.tsv`)

What it validates:
- raw TSV ingestion
- alias mapping
- missing-field diagnostics
- `REF` derivation
- GWAS-style variant visualization
- gene mapping + enrichment

Expected PRS outcome:
- no individual PRS scores from summary statistics alone (no sample genotypes)

### Mode 2: Variant-only VCF

Recommended test:
- 1000 Genomes Phase 3 chromosome VCF (chr19 for faster local tests)

What it validates:
- external VCF normalization/filtering
- genotype-aware QC metrics
- feature intersection and annotation
- enrichment
- PRS harmonization logic

Reference-build caveat:
- 1000G Phase 3 is typically GRCh37; disable strict build matching for mechanics tests or use GRCh37 resources

### Mode 3: Full FASTQ pipeline

Recommended tests:
- GIAB HG002 subset (benchmark-style)
- smaller local subset FASTQ for rapid path validation

What it validates:
- preflight
- QC + trimming
- alignment + indexing + duplicate marking
- calling + filtering
- annotation/QC/enrichment
- manifest generation

### Completed real-data rerun (2026-02-27)

These three canonical runs were completed successfully and are the current publish-ready examples:

- AD PRS high-PRS-focused variant-only rerun:
  - `results_final_test_ad_high_prs_20260227`
  - Config: `config/final_tests/final_test_ad_high_prs.yaml`
- GWAS summary rerun:
  - `results_final_test_gwas_advp_20260227`
  - Config: `config/final_tests/final_test_gwas_advp.yaml`
- Full DNA FASTQ rerun:
  - `results_final_test_full_giab_20260227`
  - Config: `config/final_tests/final_test_full_giab.yaml`

Detailed metrics and interpretation:
- `FINAL_TEST.md`
- `docs/REAL_DATASET_EXAMPLES.md`

## 6. Environments and Portability

This project uses layered environments:

- `envs/workflow.yaml`: Snakemake launcher environment
- `envs/alignment.yaml`: aligners + `samtools`
- `envs/qc.yaml`: `fastp`, `fastqc`, `multiqc`
- `envs/variant.yaml`: `bcftools`, `bedtools`, `mosdepth`, `snpEff`, etc.
- `envs/prs.yaml`: Python PRS scoring/reporting stack
- `envs/python_plot.yaml`: plotting/statistics scripts
- `envs/rna.yaml`, `envs/r_stats.yaml`: RNA assembly/counting + R analyses

The pipeline is designed to run in:

- local macOS/Linux development environments
- workstation or server environments
- HPC nodes (with Snakemake profiles, if added)
- CI-style dry-runs for config validation

## 7. Limitations and Known Constraints

- PRS requires per-sample genotype/dosage fields; GWAS summary tables cannot produce individual scores by themselves
- Biological interpretation is only valid when reference build and annotation build match
- Whole-genome depth reports can be heavy in storage/compute unless test-mode depth controls are used
- Public datasets vary in schema quality and may require per-source alias expansion over time

## 8. Future Work

- Add explicit preflight checks for alignment index readiness (`bwa-mem2` and `bwa`)
- Add optional profile support for HPC/SLURM execution
- Add small synthetic fixture datasets for CI
- Add stricter production mode to fail if any raw TSV rows are dropped during conversion
- Add containerized deployment (Docker/Singularity) examples

## 9. Reproducibility and Citation

This repository is open source and intended for reuse, benchmarking, and extension. If used in a publication, benchmark, or report:

- cite the repository release (see `CITATION.cff`)
- report exact config file(s) used
- report reference resources/build versions
- report pipeline version/commit or release tag

That combination is necessary for meaningful scientific reproducibility.
